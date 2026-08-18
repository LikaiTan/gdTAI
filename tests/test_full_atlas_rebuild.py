from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

from workflows.integration.rebuild_full_atlas import (
    align_matrix,
    hvg_excluded,
    read_config,
    standardize_obs,
)


ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "configs/datasets/full_atlas_rebuild.json"


def test_frozen_contract_totals() -> None:
    config = read_config(CONFIG)
    assert len(config["inputs"]) == 16
    assert sum(int(row["expected_cells_effective"]) for row in config["inputs"]) == 5_933_312
    assert config["expected_historical_cells"] + config["expected_extension_cells"] + config["expected_copd_cells"] == 5_933_312
    assert config["source_policy"]["propagate_repaired_tcr_sidecar"] is False
    assert all(row["expected_sha256"] != "PENDING_FULL_SHA256" for row in config["inputs"])


def test_sparse_alignment_preserves_values_and_target_order() -> None:
    matrix = sparse.csr_matrix(np.asarray([[1, 2, 3], [4, 5, 6]], dtype=np.float32))
    aligned, missing = align_matrix(matrix, pd.Index(["B", "A", "unused"]), pd.Index(["A", "B", "C"]))
    np.testing.assert_array_equal(aligned.toarray(), np.asarray([[2, 1, 0], [5, 4, 0]], dtype=np.float32))
    assert missing == 1


def test_standardized_metadata_derives_tcr_flags_and_unique_ids() -> None:
    obs = pd.DataFrame(
        {
            "source_gse_id": ["GSEX", "GSEX"],
            "sample_id": ["S1", "S1"],
            "library_id": ["L1", "L1"],
            "TRA_cdr3": ["CAVA", ""],
            "TRB_cdr3": ["CASB", ""],
            "TRD_cdr3_aa": ["", "CALD"],
            "productive_TRG": [False, True],
        },
        index=["AAAC-1", "AAAG-1"],
    )
    row = {
        "cohort_id": "GSEX",
        "fallback_source_gse_id": "GSEX",
        "atlas_input_role": "new_extension",
        "model_evaluation_role_frozen": "locked_validation_cohort",
    }
    result = standardize_obs(obs, row)
    assert result.index.tolist() == ["GSEX::AAAC-1", "GSEX::AAAG-1"]
    assert result.index.is_unique
    assert result["has_TRA_TRB_paired"].tolist() == [True, False]
    assert result["has_TRG_TRD_paired"].tolist() == [False, True]
    assert result["has_any_ab_tcr"].tolist() == [True, False]
    assert result["has_any_gd_tcr"].tolist() == [False, True]
    assert result["integration_batch"].astype(str).tolist() == ["GSEX::library_id::L1", "GSEX::library_id::L1"]


def test_hvg_exclusion_is_specific_to_tcr_segments() -> None:
    assert hvg_excluded("TRDV2")
    assert hvg_excluded("TRBJ1-1")
    assert hvg_excluded("MT-CO1")
    assert hvg_excluded("RPL13")
    assert not hvg_excluded("TRDC")
    assert not hvg_excluded("TRGC1")
    assert not hvg_excluded("CD3D")
