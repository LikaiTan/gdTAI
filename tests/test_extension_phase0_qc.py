import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

pd.options.future.infer_string = False


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "workflows/intake/qc_extension_h5ads.py"
SPEC = importlib.util.spec_from_file_location("qc_extension_h5ads", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def valid_obs() -> pd.DataFrame:
    barcodes = ["AAACCCGG-1", "AAACCGGT-1", "AAACCTTG-1", "AAAGGGTT-1"]
    cores = [value.split("-")[0] for value in barcodes]
    obs = pd.DataFrame(
        {
            "sample_id": ["sample_1", "sample_1", "sample_2", "sample_2"],
            "library_id": ["library_1", "library_1", "library_2", "library_2"],
            "donor_id": ["donor_1", "donor_1", "donor_2", "donor_2"],
            "tissue": ["blood", "blood", "tumor", "tumor"],
            "diagnosis": ["cancer"] * 4,
            "technology_simple": ["10x 5'"] * 4,
            "source_accession": ["GSETEST"] * 4,
            "barcode": barcodes,
            "barcode_core": cores,
            "tcr_schema_provenance": ["embedded_airr_tra_trb"] * 4,
        },
        index=pd.Index(
            [f"library_{1 if index < 2 else 2}:{core}" for index, core in enumerate(cores)],
            name="cell_id",
        ),
    )
    for chain in MODULE.CHAINS:
        for suffix in MODULE.STRING_SUFFIXES:
            obs[f"{chain}_{suffix}"] = ""
        for suffix in MODULE.NUMERIC_SUFFIXES:
            obs[f"{chain}_{suffix}"] = np.zeros(len(obs), dtype=np.int64)

    for chain, row, cdr3 in (
        ("TRA", 0, "CAVRDSNYQLIW"),
        ("TRB", 0, "CASSLGQETQYF"),
        ("TRG", 1, "CALWEVQELGKKIKVF"),
        ("TRD", 1, "CACDTGGY"),
        ("TRB", 2, "CASSQETQYF"),
    ):
        obs.iloc[row, obs.columns.get_loc(f"{chain}_cdr3")] = cdr3
        obs.iloc[row, obs.columns.get_loc(f"{chain}_v")] = f"{chain}V1"
        obs.iloc[row, obs.columns.get_loc(f"{chain}_j")] = f"{chain}J1"
        obs.iloc[row, obs.columns.get_loc(f"{chain}_umis")] = 2
        obs.iloc[row, obs.columns.get_loc(f"{chain}_reads")] = 10

    for chain in MODULE.CHAINS:
        obs[f"has_{chain}"] = obs[f"{chain}_cdr3"].ne("")
    obs["has_TRA_TRB_paired"] = obs["has_TRA"] & obs["has_TRB"]
    obs["has_TRG_TRD_paired"] = obs["has_TRG"] & obs["has_TRD"]
    obs["has_any_ab_tcr"] = obs["has_TRA"] | obs["has_TRB"]
    obs["has_any_gd_tcr"] = obs["has_TRG"] | obs["has_TRD"]
    obs["TCRseq"] = np.where(obs["has_any_ab_tcr"] | obs["has_any_gd_tcr"], "yes", "no")
    return obs


def write_h5ad(
    path: Path,
    *,
    obs: pd.DataFrame | None = None,
    matrix: np.ndarray | sparse.spmatrix | None = None,
    external_join: bool = False,
) -> None:
    obs = valid_obs() if obs is None else obs
    obs = obs.copy()
    for column in obs.columns:
        if pd.api.types.is_string_dtype(obs[column].dtype):
            obs[column] = obs[column].astype(object)
    obs.index = pd.Index(obs.index.astype(str).to_numpy(dtype=object), name=obs.index.name)
    if matrix is None:
        matrix = sparse.csr_matrix(
            np.array(
                [
                    [5, 1, 0, 2],
                    [3, 0, 4, 1],
                    [2, 1, 0, 3],
                    [1, 0, 0, 1],
                ],
                dtype=np.int32,
            )
        )
    var = pd.DataFrame(
        index=pd.Index(np.asarray(["CD3D", "MT-CO1", "TRDC", "TRBC1"], dtype=object), name="gene_symbol")
    )
    adata = ad.AnnData(X=matrix, obs=obs, var=var)
    if external_join:
        adata.obs["tcr_schema_provenance"] = "productive_contigs"
        adata.uns["tcr_join_summaries"] = json.dumps(
            [
                {
                    "sample_id": "sample_1",
                    "expression_cells": 2,
                    "matched_cells": 2,
                    "unmatched_tcr_cells": 1,
                },
                {
                    "sample_id": "sample_2",
                    "expression_cells": 2,
                    "matched_cells": 1,
                    "unmatched_tcr_cells": 0,
                },
            ]
        )
    adata.write_h5ad(path)


class ExtensionPhase0QCTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.h5ad = self.root / "test_cohort.h5ad"

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def entry(self) -> object:
        return MODULE.ManifestEntry("TEST_COHORT", self.h5ad, 4, 2)

    def write_manifest(self) -> Path:
        path = self.root / "built_h5ads_manifest.csv"
        pd.DataFrame(
            [
                {
                    "cohort_id": "TEST_COHORT",
                    "h5ad_path": str(self.h5ad),
                    "build_status": "built",
                    "n_cells": 4,
                }
            ]
        ).to_csv(path, index=False)
        return path

    def test_valid_gate_writes_canonical_review_artifacts_without_mutation(self) -> None:
        write_h5ad(self.h5ad)
        manifest = self.write_manifest()
        output_root = self.root / "Integrated_dataset"
        before = self.h5ad.stat()

        return_code = MODULE.main(
            ["--manifest", str(manifest), "--output-root", str(output_root)]
        )

        self.assertEqual(return_code, 0)
        after = self.h5ad.stat()
        self.assertEqual((before.st_size, before.st_mtime_ns), (after.st_size, after.st_mtime_ns))
        summary_path = (
            output_root / "logs/extension_intake/extension_phase0_qc_summary.json"
        )
        summary = json.loads(summary_path.read_text())
        self.assertEqual(summary["gate_status"], "PASS_REVIEW_REQUIRED")
        self.assertTrue(summary["review_required"])
        self.assertFalse(summary["merge_approved"])
        self.assertEqual(summary["cell_count"], 4)
        cohort = pd.read_csv(
            output_root / "tables/extension_intake/extension_phase0_cohort_summary.csv"
        ).iloc[0]
        self.assertEqual(int(cohort["tcr_positive_cells"]), 3)
        self.assertEqual(int(cohort["paired_TRA_TRB_cells"]), 1)
        self.assertEqual(int(cohort["paired_TRG_TRD_cells"]), 1)
        self.assertTrue(
            (output_root / "figures/extension_intake/extension_phase0_cohort_overview.png").is_file()
        )

    def test_dry_run_writes_nothing(self) -> None:
        write_h5ad(self.h5ad)
        manifest = self.write_manifest()
        output_root = self.root / "no_outputs"
        return_code = MODULE.main(
            [
                "--manifest",
                str(manifest),
                "--output-root",
                str(output_root),
                "--dry-run",
            ]
        )
        self.assertEqual(return_code, 0)
        self.assertFalse(output_root.exists())

    def test_duplicate_sample_barcode_key_fails_closed(self) -> None:
        obs = valid_obs()
        obs.loc[obs.index[1], ["barcode", "barcode_core"]] = obs.loc[
            obs.index[0], ["barcode", "barcode_core"]
        ].to_numpy()
        write_h5ad(self.h5ad, obs=obs)
        with self.assertRaisesRegex(MODULE.Phase0QCError, "not unique"):
            MODULE.audit_h5ad(self.entry())

    def test_inconsistent_tcr_boolean_fails_closed(self) -> None:
        obs = valid_obs()
        obs.loc[obs.index[3], "has_TRD"] = True
        write_h5ad(self.h5ad, obs=obs)
        with self.assertRaisesRegex(MODULE.Phase0QCError, "has_TRD disagrees"):
            MODULE.audit_h5ad(self.entry())

    def test_command_failure_writes_fail_record_and_returns_nonzero(self) -> None:
        obs = valid_obs()
        obs.loc[obs.index[3], "has_TRD"] = True
        write_h5ad(self.h5ad, obs=obs)
        manifest = self.write_manifest()
        output_root = self.root / "failed_review"

        return_code = MODULE.main(
            ["--manifest", str(manifest), "--output-root", str(output_root)]
        )

        self.assertEqual(return_code, 2)
        summary = json.loads(
            (
                output_root
                / "logs/extension_intake/extension_phase0_qc_summary.json"
            ).read_text()
        )
        self.assertEqual(summary["gate_status"], "FAIL")
        self.assertEqual(summary["error_count"], 1)
        self.assertFalse(summary["merge_approved"])

    def test_dense_or_noninteger_x_fails_closed(self) -> None:
        dense_path = self.root / "dense.h5ad"
        write_h5ad(dense_path, matrix=np.ones((4, 4), dtype=np.int32))
        dense_entry = MODULE.ManifestEntry("TEST_COHORT", dense_path, 4, 2)
        with self.assertRaisesRegex(MODULE.Phase0QCError, "dense X"):
            MODULE.audit_h5ad(dense_entry)

        matrix = sparse.csr_matrix(np.ones((4, 4), dtype=np.float32))
        matrix.data[0] = 0.5
        write_h5ad(self.h5ad, matrix=matrix)
        with self.assertRaisesRegex(MODULE.Phase0QCError, "non-integer"):
            MODULE.audit_h5ad(self.entry())

    def test_external_join_coverage_is_checked_and_reported(self) -> None:
        write_h5ad(self.h5ad, external_join=True)
        summary, libraries, issues = MODULE.audit_h5ad(self.entry())
        self.assertEqual(summary["join_mode"], "sample_id+barcode_core")
        self.assertAlmostEqual(summary["join_coverage"], 0.75)
        self.assertEqual(int(summary["unmatched_tcr_cells"]), 1)
        self.assertEqual(len(libraries), 2)
        self.assertTrue(any(item["check"] == "unmatched_tcr_cells" for item in issues))


if __name__ == "__main__":
    unittest.main()
