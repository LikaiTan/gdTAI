import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite

pd.options.future.infer_string = False


ROOT = Path(__file__).resolve().parents[1]
INTAKE_DIR = ROOT / "workflows" / "intake"
if str(INTAKE_DIR) not in sys.path:
    sys.path.insert(0, str(INTAKE_DIR))

import build_extension_h5ads as builder  # noqa: E402
import stage_extension_cohorts as stager  # noqa: E402
import validate_extension_cohorts as validator  # noqa: E402


class ExtensionCohortIntakeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        report, cohorts, libraries = validator.validate_extension_manifests()
        if not report.passed:
            raise AssertionError(report.to_dict())
        cls.cohorts = cohorts
        cls.libraries = libraries
        cls.cohort_by_id = {row["cohort_id"]: row for row in cohorts}

    def test_manifests_encode_required_cohorts_and_policies(self) -> None:
        report, _, _ = validator.validate_extension_manifests()
        self.assertTrue(report.passed, report.to_dict())
        self.assertEqual(report.included_cohort_count, 8)

        gse121 = self.cohort_by_id["GSE121636_GSE121637"]
        self.assertEqual(gse121["gex_accessions"], "GSE121636")
        self.assertEqual(gse121["vdj_accessions"], "GSE121637")
        self.assertEqual(gse121["intake_role"], "sealed_mixed_t_nk_holdout")
        self.assertEqual(gse121["sealed"], "true")

        gse315 = self.cohort_by_id["GSE315928"]
        self.assertEqual(gse315["intake_role"], "sealed_abt_negative_holdout")
        self.assertEqual(gse315["tcr_schema"], "embedded_airr_tra_trb")

        gse159 = self.cohort_by_id["GSE159251"]
        self.assertEqual(gse159["tcr_schema"], "partial_embedded_tra_trb_cdr3")

        gse294 = self.cohort_by_id["GSE294273_GSE294274"]
        self.assertEqual(gse294["default_tissue"], "lymph_node_metastasis")
        self.assertEqual(gse294["default_diagnosis"], "melanoma")
        self.assertNotIn("esophageal", " ".join(gse294.values()).lower())

        gse169 = self.cohort_by_id["GSE169246"]
        self.assertEqual(gse169["builder_adapter"], "gse169246_processed_h5ad")
        self.assertEqual(gse169["tcr_schema"], "partial_embedded_paired_tra_trb_cdr3")
        self.assertEqual(gse169["intake_role"], "extension_candidate_pending_phase0")
        self.assertEqual(gse169["stage_enabled"], "false")
        self.assertEqual(gse169["build_enabled"], "false")

        duplicate = self.cohort_by_id["TAN2021_GDT_DUPLICATE"]
        self.assertEqual(duplicate["duplicate_of"], "GDT_2020AUG_woCOV")
        self.assertEqual(duplicate["stage_enabled"], "false")
        self.assertEqual(duplicate["build_enabled"], "false")

    def _write_cohort_rows(self, path: Path, rows: list[dict[str, str]]) -> None:
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=validator.COHORT_COLUMNS, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)

    def test_duplicate_checks_cover_cohort_accession_publication_and_bioproject(self) -> None:
        cases = {
            "cohort": ("cohort_key", self.cohort_by_id["GSE114724"]["cohort_key"]),
            "accession": ("gex_accessions", "GSE114724"),
            "publication": ("publication_ids", "PMID:29961579"),
            "bioproject": ("bioproject_ids", "PRJNA472381"),
        }
        with tempfile.TemporaryDirectory() as temporary:
            for code, (field, duplicate_value) in cases.items():
                with self.subTest(code=code):
                    rows = [dict(row) for row in self.cohorts]
                    target = next(row for row in rows if row["cohort_id"] == "GSE292700")
                    target[field] = duplicate_value
                    path = Path(temporary) / f"{code}.csv"
                    self._write_cohort_rows(path, rows)
                    report, _, _ = validator.validate_extension_manifests(
                        cohort_path=path,
                        library_path=validator.DEFAULT_LIBRARIES,
                        canonical_registry_path=validator.DEFAULT_CANONICAL_REGISTRY,
                        shared_metasheet_path=validator.DEFAULT_SHARED_METASHEET,
                    )
                    self.assertIn(f"duplicate_{code}", {issue.code for issue in report.errors})

    def test_shared_metasheet_is_filtered_by_linked_accessions(self) -> None:
        columns, rows = stager.filter_shared_metasheet(
            validator.DEFAULT_SHARED_METASHEET,
            self.cohort_by_id["GSE294273_GSE294274"],
        )
        self.assertIn("gsm", columns)
        self.assertEqual(len(rows), 24)
        self.assertEqual({row["gse"] for row in rows}, {"GSE294273", "GSE294274"})
        self.assertTrue(all("melanoma" in row["cancer"].lower() for row in rows))

    def test_gse169246_selected_source_remains_stage_and_build_disabled(self) -> None:
        cohort = [self.cohort_by_id["GSE169246"]]
        libraries = [row for row in self.libraries if row["cohort_id"] == "GSE169246"]
        with tempfile.TemporaryDirectory() as temporary:
            stage_summary = stager.stage_extension_cohorts(
                Path(temporary),
                Path(temporary) / "stage",
                cohort,
                libraries,
                validator.DEFAULT_SHARED_METASHEET,
                dry_run=True,
            )
            self.assertEqual(stage_summary["cohorts"], [])
            build_summary = builder.build_extension_cohorts(
                Path(temporary) / "stage",
                Path(temporary) / "output",
                cohort,
                dry_run=True,
            )
            self.assertEqual(build_summary["plans"], [])

    def test_stage_snapshots_filtered_shared_metadata_without_local_metasheet(self) -> None:
        cohort = self.cohort_by_id["GSE294273_GSE294274"]
        libraries = [row for row in self.libraries if row["cohort_id"] == cohort["cohort_id"]]
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            gex = root / "sources/GSE294273/extracted"
            vdj = root / "sources/GSE294274/extracted"
            gex.mkdir(parents=True)
            vdj.mkdir(parents=True)
            for suffix in ("matrix.mtx.gz", "barcodes.tsv.gz", "features.tsv.gz"):
                (gex / f"GSM8901614_HM002_{suffix}").write_bytes(b"source")
            (vdj / "GSM8901626_HM002_filtered_contig_annotations.csv.gz").write_bytes(b"source")
            stage_root = root / "stage"
            summary = stager.stage_extension_cohorts(
                root / "sources",
                stage_root,
                [cohort],
                libraries,
                validator.DEFAULT_SHARED_METASHEET,
            )
            self.assertEqual(len(summary["staged"]), 1)
            snapshot = pd.read_csv(stage_root / cohort["cohort_id"] / "sample_metasheet.csv")
            self.assertEqual(len(snapshot), 24)
            self.assertEqual(set(snapshot["gse"]), {"GSE294273", "GSE294274"})
            manifest = json.loads(
                (stage_root / cohort["cohort_id"] / "stage_manifest.json").read_text()
            )
            self.assertEqual(manifest["metasheet"]["row_count"], 24)
            self.assertEqual(manifest["join_key"], "sample_id+barcode_core")

    def test_biological_pair_keys_are_deterministic(self) -> None:
        pairs = (
            ("HM002_GEX", "HM002_TCR", "HM002"),
            ("GU0700_T", "GU0700_T [VDJ sequencing]", "GU0700_T"),
            ("BRCA2-mutant tumors", "BRCA2-mutant tumors (TCR)", "BRCA2-mutant_tumors"),
            (
                "scRNA-seq of HNSCC-infiltrating T cells; 17B",
                "scVDJ-seq of HNSCC-infiltrating T cells; 17B",
                "17B",
            ),
        )
        for gex, tcr, expected in pairs:
            with self.subTest(gex=gex):
                self.assertEqual(builder.biological_sample_key(gex), expected)
                self.assertEqual(builder.biological_sample_key(tcr), expected)

    def test_productive_contigs_use_support_ranking_and_full_chain_schema(self) -> None:
        frame = pd.DataFrame(
            {
                "barcode": ["AACCGGTT-1"] * 5,
                "chain": ["TRA", "TRA", "TRB", "TRG", "TRD"],
                "productive": [True] * 5,
                "is_cell": [True] * 5,
                "high_confidence": [True] * 5,
                "cdr3": ["LOW", "HIGH", "BETA", "GAMMA", "DELTA"],
                "v_gene": ["TRAV1", "TRAV2", "TRBV1", "TRGV1", "TRDV1"],
                "j_gene": ["TRAJ1", "TRAJ2", "TRBJ1", "TRGJ1", "TRDJ1"],
                "umis": [1, 10, 4, 3, 2],
                "reads": [100, 20, 40, 30, 20],
            }
        )
        result = builder.collapse_productive_contigs(frame, "S1")
        self.assertEqual(result.loc[0, "TRA_cdr3"], "HIGH")
        self.assertEqual(result.loc[0, "TRB_cdr3"], "BETA")
        self.assertEqual(result.loc[0, "TRG_cdr3"], "GAMMA")
        self.assertEqual(result.loc[0, "TRD_cdr3"], "DELTA")
        self.assertTrue(set(builder.TCR_COLUMNS).issubset(result.columns))

        duplicated = pd.concat([frame.iloc[[0]], frame.iloc[[0]]], ignore_index=True)
        with self.assertRaises(builder.IntakeBuildError):
            builder.collapse_productive_contigs(duplicated, "S1")

    def test_tcr_join_never_uses_barcode_only_and_rejects_duplicate_keys(self) -> None:
        adata = ad.AnnData(
            X=sparse.csr_matrix(np.array([[1, 0], [0, 2]], dtype=np.int64)),
            obs=pd.DataFrame(
                {
                    "sample_id": ["S1", "S2"],
                    "barcode_core": ["AACCGGTT", "AACCGGTT"],
                },
                index=["L1:AACCGGTT", "L2:AACCGGTT"],
            ),
        )
        tcr = pd.DataFrame(
            {
                "sample_id": ["S1", "S2"],
                "barcode_core": ["AACCGGTT", "AACCGGTT"],
                "TRA_cdr3": ["TRA1", "TRA2"],
            }
        )
        for column in builder.TCR_COLUMNS:
            if column not in tcr:
                tcr[column] = 0 if column.endswith(("_umis", "_reads")) else ""
        summary = builder.join_tcr_by_sample_barcode(adata, tcr)
        self.assertEqual(summary["matched_cells"], 2)
        self.assertEqual(adata.obs["TRA_cdr3"].tolist(), ["TRA1", "TRA2"])

        duplicated = pd.concat([tcr.iloc[[0]], tcr.iloc[[0]]], ignore_index=True)
        with self.assertRaises(builder.IntakeBuildError):
            builder.join_tcr_by_sample_barcode(adata, duplicated)
        with self.assertRaises(builder.IntakeBuildError):
            builder.join_tcr_by_sample_barcode(adata, tcr.drop(columns="sample_id"))

    def test_gse159251_partial_schema_and_notcr_sentinel(self) -> None:
        parsed = builder.parse_embedded_tra_trb(
            pd.Series(["TRA_A|TRB_A", "notcr", "TRA_ONLY|", "|TRB_ONLY"])
        )
        self.assertEqual(parsed["TRA_cdr3"].tolist(), ["TRA_A", "", "TRA_ONLY", ""])
        self.assertEqual(parsed["TRB_cdr3"].tolist(), ["TRB_A", "", "", "TRB_ONLY"])
        self.assertTrue(parsed["TRG_cdr3"].eq("").all())
        self.assertTrue(parsed["TRD_cdr3"].eq("").all())
        with self.assertRaises(builder.IntakeBuildError):
            builder.parse_embedded_tra_trb(pd.Series(["TRA|TRB|EXTRA"]))

    def test_gse315928_airr_mapping_selects_productive_supported_chains(self) -> None:
        obs = pd.DataFrame(
            {
                "IR_VJ_1_locus": ["TRA"],
                "IR_VJ_1_productive": [True],
                "IR_VJ_1_junction_aa": ["TRA_LOW"],
                "IR_VJ_1_duplicate_count": [1],
                "IR_VJ_1_consensus_count": [100],
                "IR_VJ_2_locus": ["TRA"],
                "IR_VJ_2_productive": [True],
                "IR_VJ_2_junction_aa": ["TRA_HIGH"],
                "IR_VJ_2_duplicate_count": [5],
                "IR_VJ_2_consensus_count": [10],
                "IR_VDJ_1_locus": ["TRB"],
                "IR_VDJ_1_productive": [True],
                "IR_VDJ_1_junction_aa": ["TRB_ONE"],
                "IR_VDJ_1_duplicate_count": [3],
                "IR_VDJ_1_consensus_count": [30],
                "clonotype": ["clone1"],
            },
            index=["cell1"],
        )
        mapped = builder.map_embedded_airr_tra_trb(obs)
        self.assertEqual(mapped.loc["cell1", "TRA_cdr3"], "TRA_HIGH")
        self.assertEqual(mapped.loc["cell1", "TRB_cdr3"], "TRB_ONE")
        self.assertEqual(mapped.loc["cell1", "TRA_clone_id"], "clone1")
        self.assertEqual(mapped.loc["cell1", "TRG_cdr3"], "")
        self.assertEqual(mapped.loc["cell1", "TRD_cdr3"], "")

    def _make_synthetic_gse294_stage(self, root: Path) -> Path:
        cohort_id = "GSE294273_GSE294274"
        stage = root / cohort_id
        gex_dir = stage / "sources/GSE294273_GEX"
        tcr_dir = stage / "sources/GSE294274_VDJ"
        gex_dir.mkdir(parents=True)
        tcr_dir.mkdir(parents=True)
        prefix = "GSM8901614_HM002_"
        matrix_path = gex_dir / f"{prefix}matrix.mtx"
        mmwrite(matrix_path, sparse.coo_matrix(np.array([[1, 0], [0, 2]], dtype=np.int64)))
        (gex_dir / f"{prefix}barcodes.tsv").write_text("AACCGGTT-1\nTTGGCCAA-1\n")
        (gex_dir / f"{prefix}features.tsv").write_text("g1\tGENE1\ng2\tGENE2\n")
        tcr_path = tcr_dir / "GSM8901626_HM002_filtered_contig_annotations.csv"
        pd.DataFrame(
            {
                "barcode": ["AACCGGTT-1", "AACCGGTT-1", "TTGGCCAA-1", "TTGGCCAA-1"],
                "chain": ["TRA", "TRB", "TRG", "TRD"],
                "productive": [True, True, True, True],
                "is_cell": [True, True, True, True],
                "high_confidence": [True, True, True, True],
                "cdr3": ["TRA_A", "TRB_A", "TRG_B", "TRD_B"],
                "umis": [3, 4, 5, 6],
                "reads": [30, 40, 50, 60],
            }
        ).to_csv(tcr_path, index=False)
        metasheet = pd.DataFrame(
            [
                {
                    "gse": "GSE294273",
                    "gsm": "GSM8901614",
                    "sample_title": "HM002_GEX",
                    "modality": "GEX",
                    "cancer": "melanoma (LN metastasis)",
                    "source": "lymph node",
                },
                {
                    "gse": "GSE294274",
                    "gsm": "GSM8901626",
                    "sample_title": "HM002_TCR",
                    "modality": "TCR",
                    "cancer": "melanoma LN metastasis",
                    "source": "lymph node",
                },
            ]
        )
        metasheet.to_csv(stage / "sample_metasheet.csv", index=False)
        manifest = {
            "schema_version": 1,
            "cohort_id": cohort_id,
            "join_key": "sample_id+barcode_core",
            "metasheet": {},
            "sources": [
                {
                    "source_role": "gex",
                    "staged_path": str(matrix_path.relative_to(stage)),
                },
                {
                    "source_role": "vdj",
                    "staged_path": str(tcr_path.relative_to(stage)),
                },
            ],
        }
        (stage / "stage_manifest.json").write_text(json.dumps(manifest))
        return stage

    def test_synthetic_staged_h5ad_passes_qc_and_serialized_validation(self) -> None:
        cohort = self.cohort_by_id["GSE294273_GSE294274"]
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            stage = self._make_synthetic_gse294_stage(root)
            adata, summary = builder.build_one_cohort(
                stage, cohort, root, builder.DEFAULT_R_EXPORTER
            )
            self.assertTrue(sparse.isspmatrix_csr(adata.X))
            self.assertEqual(summary["samples"], 1)
            self.assertTrue(adata.obs_names.is_unique)
            self.assertTrue(adata.obs["tissue"].eq("lymph_node_metastasis").all())
            self.assertTrue(adata.obs["diagnosis"].eq("melanoma").all())
            self.assertTrue(set(builder.TCR_COLUMNS).issubset(adata.obs.columns))
            self.assertTrue(set(builder.TCR_FLAGS).issubset(adata.obs.columns))
            self.assertEqual(int(adata.obs["has_TRA_TRB_paired"].sum()), 1)
            self.assertEqual(int(adata.obs["has_TRG_TRD_paired"].sum()), 1)

            output = root / f"{cohort['cohort_id']}.h5ad"
            adata.write_h5ad(output)
            issues = validator.validate_built_h5ad(output, cohort)
            self.assertEqual(issues, [])
            qc_path = root / f"{cohort['cohort_id']}.standalone_qc.json"
            qc_path.write_text(json.dumps(summary), encoding="utf-8")
            resumed = builder.validated_existing_build(output, qc_path, cohort)
            self.assertTrue(resumed["resumed_existing"])
            self.assertEqual(resumed["cells"], adata.n_obs)

    def test_standalone_qc_rejects_dense_or_fractional_x(self) -> None:
        obs = pd.DataFrame(
            {
                "sample_id": ["S1"],
                "library_id": ["L1"],
                "donor_id": ["D1"],
                "tissue": ["blood"],
                "diagnosis": ["healthy"],
                "barcode_core": ["AACCGGTT"],
            },
            index=["L1:AACCGGTT"],
        )
        dense = ad.AnnData(X=np.array([[1]], dtype=np.int64), obs=obs.copy())
        with self.assertRaises(builder.IntakeBuildError):
            builder.standalone_sample_qc(dense, "dense")
        fractional = ad.AnnData(X=sparse.csr_matrix([[1.5]]), obs=obs.copy())
        with self.assertRaises(builder.IntakeBuildError):
            builder.standalone_sample_qc(fractional, "fractional")

    def test_built_manifest_is_atomic_and_checksum_pinned(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output_root = root / "built_h5ads"
            output_root.mkdir()
            h5ad = output_root / "TEST.h5ad"
            qc = output_root / "TEST.standalone_qc.json"
            h5ad.write_bytes(b"immutable-h5ad-fixture")
            qc.write_text("{}\n", encoding="utf-8")
            path = builder.write_built_manifest(
                output_root,
                [
                    {
                        "cohort_id": "TEST",
                        "h5ad": str(h5ad),
                        "standalone_qc": str(qc),
                        "cells": 17,
                        "genes": 23,
                    }
                ],
            )
            frame = pd.read_csv(path)
            self.assertEqual(frame.loc[0, "build_status"], "built")
            self.assertEqual(int(frame.loc[0, "n_cells"]), 17)
            self.assertEqual(frame.loc[0, "sha256"], builder._sha256(h5ad))
            with self.assertRaises(FileExistsError):
                builder.write_built_manifest(output_root, [])

    def test_r_exporter_encodes_nested_gzip_and_partial_provenance(self) -> None:
        source = builder.DEFAULT_R_EXPORTER.read_text(encoding="utf-8")
        self.assertIn("layers < 2L", source)
        self.assertIn("raw_sparse_integer_counts", source)
        self.assertIn("partial_embedded_tra_trb_cdr3", source)
        self.assertIn("notcr", source)


if __name__ == "__main__":
    unittest.main()
