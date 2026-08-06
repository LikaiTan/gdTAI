from pathlib import Path
import sys
import unittest

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]
INTAKE_DIR = ROOT / "workflows" / "intake"
if str(INTAKE_DIR) not in sys.path:
    sys.path.insert(0, str(INTAKE_DIR))

import filter_extension_tnk_cells as tnk_filter  # noqa: E402


class ExtensionTnkFilterTests(unittest.TestCase):
    def test_selection_keeps_t_nk_and_productive_tcr_but_removes_non_tnk(self) -> None:
        genes = ["CD3D", "NKG7", "GNLY", "KLRD1", "MS4A1", "LST1", "TRAC", "IL7R"]
        matrix = sparse.csr_matrix(
            np.asarray(
                [
                    [3, 0, 0, 0, 0, 0, 2, 0],
                    [0, 4, 3, 2, 0, 0, 0, 0],
                    [0, 0, 0, 0, 5, 0, 0, 0],
                    [0, 0, 0, 0, 0, 4, 0, 0],
                    [0, 0, 0, 0, 5, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 5],
                ],
                dtype=np.int64,
            )
        )
        obs = pd.DataFrame(
            {
                "sample_id": ["S1"] * 6,
                "library_id": ["L1"] * 6,
                "donor_id": ["D1"] * 6,
                "tissue": ["tumor"] * 6,
                "diagnosis": ["cancer"] * 6,
                "major_cluster": ["Tcell", "NKcell", "Bcell", "Myeloid", "Bcell", ""],
                "TRA_cdr3": ["", "", "", "", "CAVR", ""],
                "TRB_cdr3": ["", "", "", "", "CASS", ""],
            },
            index=[f"cell{i}" for i in range(6)],
        )
        adata = ad.AnnData(X=matrix, obs=obs, var=pd.DataFrame(index=genes))
        decisions = tnk_filter.compute_decisions(adata, "TEST")
        self.assertEqual(decisions["keep"].tolist(), [True, True, False, False, True, False])
        self.assertEqual(int(decisions["productive_tcr_evidence"].sum()), 1)
        self.assertEqual(decisions.loc["cell4", "selection_reason"], "productive_tcr")
        self.assertEqual(int(decisions.loc["cell5", "t_hits"]), 1)
        self.assertEqual(int(decisions.loc["cell5", "core_t_hits"]), 0)

    def test_gse169246_metadata_propagates_within_sample_without_guessing(self) -> None:
        obs = pd.DataFrame(
            {
                "sample_id": ["P001.Pre", "P001.Pre", "P002.Post"],
                "patient_id": ["P001", "", ""],
                "tissue": ["breast", "", ""],
                "timepoint_group": ["Pre-treatment", "", ""],
                "TRA_cdr3": ["", "", ""],
                "TRB_cdr3": ["", "", ""],
                "raw_barcode": ["AAAA-1", "CCCC-1", "GGGG-1"],
            },
            index=["a", "b", "c"],
        )
        result = tnk_filter.harmonize_gse169246_obs(obs)
        self.assertEqual(result.loc["b", "tissue_harmonized"], "breast")
        self.assertEqual(result.loc["c", "tissue_harmonized"], "unresolved")
        self.assertEqual(result.loc["c", "donor_id"], "P002")
        self.assertEqual(result.loc["a", "specimen_context"], "primary_tumor")
        self.assertEqual(result.loc["a", "library_id"], "P001.Pre")

    def test_gse169246_distinguishes_libraries_with_reused_barcodes(self) -> None:
        obs = pd.DataFrame(
            {
                "sample_id": ["P002.Post", "P002.Post"],
                "patient_id": ["P002", "P002"],
                "tissue": ["breast", "breast"],
                "raw_barcode": [
                    "AAAGATGAGCCTATGT.Post_P002_b",
                    "AAAGATGAGCCTATGT.Post_P002_t",
                ],
            },
            index=["b_cell", "t_cell"],
        )
        result = tnk_filter.harmonize_gse169246_obs(obs)
        self.assertEqual(result["sample_id"].nunique(), 1)
        self.assertEqual(result["library_id"].nunique(), 2)
        self.assertEqual(
            result["library_id"].tolist(), ["Post_P002_b", "Post_P002_t"]
        )

    def test_common_harmonization_uses_exact_source_accession(self) -> None:
        obs = pd.DataFrame(
            {
                "sample_id": ["S1"],
                "library_id": ["L1"],
                "donor_id": ["D1"],
                "tissue": ["tumor"],
                "tissue_harmonized": ["breast_tumor"],
                "diagnosis": ["breast_cancer"],
                "specimen_context": ["primary_tumor"],
                "source_accession": ["GSE121636"],
            },
            index=["cell1"],
        )
        adata = ad.AnnData(
            X=sparse.csr_matrix(np.ones((1, 1), dtype=np.int64)),
            obs=obs,
            var=pd.DataFrame(index=["CD3D"]),
        )
        tnk_filter.harmonize_common_obs(adata, "GSE121636_GSE121637")
        self.assertEqual(adata.obs.loc["cell1", "source_gse_id"], "GSE121636")
        self.assertEqual(
            adata.obs.loc["cell1", "extension_cohort_id"], "GSE121636_GSE121637"
        )

    def test_source_contract_includes_all_new_cohorts(self) -> None:
        self.assertEqual(len(tnk_filter.ALL_COHORTS), 8)
        self.assertIn("GSE169246", tnk_filter.ALL_COHORTS)
        self.assertIn("GSE292700", tnk_filter.ALL_COHORTS)


if __name__ == "__main__":
    unittest.main()
