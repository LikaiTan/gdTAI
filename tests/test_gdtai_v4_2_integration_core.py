from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from workflows.gdtai.gdtai_v4_2_integration_core import (
    apply_current_atlas_recovery,
    attach_recovery_metadata,
    cluster_consensus_votes,
    common_eligible_genes,
    ensure_no_locked_cohorts,
    is_hvg_excluded,
    make_integration_batch,
    minimum_free_gib_for_stage,
    mixed_boundary_mask,
    read_sparse_rows_genes,
    recovery_contract_sha256,
    resource_contract_sha256,
    source_balanced_sample_indices,
    validate_preflight_approval,
    validate_project_data_approval,
    validate_pseudolabel_contract,
    validate_recovery_preflight,
    validate_resource_preflight,
    validate_role_separation,
)


ROOT = Path(__file__).resolve().parents[1]
EXECUTION_CONFIG = ROOT / "configs/models/gdtai/v4_2_integration_execution.json"
if hasattr(pd.options, "future") and hasattr(pd.options.future, "infer_string"):
    pd.options.future.infer_string = False


def pseudolabel_contract() -> dict:
    return {
        "minimum_run_agreement": 0.8,
        "minimum_anchor_nk_purity": 0.95,
        "maximum_productive_tcr_contamination": 0.02,
        "minimum_independent_anchors_per_cluster": 20,
        "minimum_primary_nk_anchors_per_cluster": 10,
        "minimum_development_sources": 3,
        "maximum_single_source_fraction": 0.7,
        "marker_thresholds_may_define_truth": False,
        "may_control_validation_guardrails": False,
    }


def test_actual_preflight_approval_and_role_contract_are_bound() -> None:
    config = json.loads(EXECUTION_CONFIG.read_text())
    result = validate_preflight_approval(config)
    assert result["approval"]["approved"] is True
    assert result["summary"]["result"] == "PASS_REVIEW_REQUIRED"
    assert result["roles"]["include_in_integration_fit"].sum() == 6
    assert result["roles"]["allow_locked_evaluation"].sum() == 3


def test_locked_roles_fail_closed() -> None:
    roles = pd.DataFrame(
        {
            "cohort_id": ["development", "locked"],
            "include_in_integration_fit": [True, True],
            "include_in_cluster_label_design": [True, False],
            "allow_pseudolabel_training": [True, False],
            "allow_model_tuning": [True, False],
            "allow_locked_evaluation": [False, True],
        }
    )
    assert "a locked cohort has a development permission" in validate_role_separation(roles)
    frame = pd.DataFrame({"input_cohort_id": ["development", "locked"]})
    with pytest.raises(RuntimeError, match="Locked cohorts leaked"):
        ensure_no_locked_cohorts(frame, roles)


def test_project_data_approval_is_absent_and_fails_closed(tmp_path: Path) -> None:
    config_path = tmp_path / "execution.json"
    runner = tmp_path / "runner.py"
    core = tmp_path / "core.py"
    runner.write_text("runner\n")
    core.write_text("core\n")
    config = {"protocol_version": "test", "execution_approval": str(tmp_path / "missing.json")}
    config_path.write_text(json.dumps(config))
    with pytest.raises(PermissionError, match="approval is absent"):
        validate_project_data_approval(config_path, config, runner, core)


def test_recovery_replaces_only_a_missing_current_atlas(tmp_path: Path) -> None:
    source = tmp_path / "cleaned.h5ad"
    source.touch()
    original = tmp_path / "integrated.h5ad"
    roles = pd.DataFrame(
        {
            "cohort_id": ["current_atlas", "extension"],
            "path": [str(original), str(tmp_path / "extension.h5ad")],
            "expected_sha256": ["old", "extension"],
            "expected_cells": [7, 3],
        }
    )
    config = {
        "current_atlas_recovery": {
            "active": True,
            "require_original_missing": True,
            "source_h5ad": str(source),
            "source_sha256": "replacement",
            "expected_effective_cells": 7,
        }
    }
    recovered = apply_current_atlas_recovery(config, roles)
    current = recovered.loc[recovered["cohort_id"].eq("current_atlas")].iloc[0]
    assert current["path"] == str(source)
    assert current["expected_sha256"] == "replacement"
    assert bool(current["recovery_active"])
    assert not bool(recovered.loc[recovered["cohort_id"].eq("extension"), "recovery_active"].iloc[0])

    original.touch()
    with pytest.raises(RuntimeError, match="forbidden while the original"):
        apply_current_atlas_recovery(config, roles)


def test_recovery_metadata_join_is_source_and_cell_specific(tmp_path: Path) -> None:
    metadata = tmp_path / "metadata.csv"
    pd.DataFrame(
        {
            "source_gse_id": ["A", "B"],
            "original_cell_id": ["cell-1", "cell-1"],
            "donor_id": ["donor-a", "donor-b"],
        }
    ).to_csv(metadata, index=False)
    config = {
        "current_atlas_recovery": {
            "metadata_sources": [{"path": str(metadata), "sha256": "not-used-here"}],
            "metadata_columns": ["donor_id"],
        }
    }
    frame = pd.DataFrame(
        {
            "source_gse_id": ["B", "A"],
            "source_original_cell_id": ["cell-1", "cell-1"],
        },
        index=["row-b", "row-a"],
    )
    joined = attach_recovery_metadata(frame, config)
    assert joined["donor_id"].tolist() == ["donor-b", "donor-a"]


def test_recovery_preflight_is_bound_to_the_exact_contract(tmp_path: Path) -> None:
    summary = tmp_path / "summary.json"
    recovery = {
        "active": True,
        "source_h5ad": str(tmp_path / "cleaned.h5ad"),
        "source_sha256": "source",
        "expected_raw_cells": 10,
        "expected_effective_cells": 9,
        "expected_genes": 4,
        "row_exclusion_manifest": str(tmp_path / "exclude.csv"),
        "row_exclusion_manifest_sha256": "exclude",
        "expected_intersecting_exclusions": 1,
        "metadata_sources": [{"path": str(tmp_path / "metadata.csv"), "sha256": "metadata"}],
        "metadata_columns": ["donor_id"],
        "expected_metadata_audit": {"donor_id": {"n_missing": 0, "n_unique": 2}},
        "recovery_preflight_summary": str(summary),
    }
    config = {"current_atlas_recovery": recovery}
    summary.write_text(
        json.dumps(
            {
                "result": "PASS_REVIEW_REQUIRED",
                "recovery_contract_sha256": recovery_contract_sha256(config),
            }
        )
    )
    assert validate_recovery_preflight(config)["result"] == "PASS_REVIEW_REQUIRED"

    config["current_atlas_recovery"]["expected_genes"] = 5
    with pytest.raises(RuntimeError, match="does not match"):
        validate_recovery_preflight(config)


def test_stage_specific_resource_floor_is_bound_and_fail_closed(tmp_path: Path) -> None:
    summary = tmp_path / "resource_summary.json"
    config = {
        "resources": {
            "ssd_root": "/ssd/test",
            "minimum_ssd_free_gib": 300,
            "minimum_ssd_free_gib_by_stage": {
                "prepare": 300,
                "fit": 300,
                "cluster": 150,
                "consensus": 150,
            },
            "maximum_ram_gib": 800,
            "maximum_cpu_cores": 80,
            "resource_preflight_summary": str(summary),
        }
    }
    summary.write_text(
        json.dumps(
            {
                "result": "PASS_REVIEW_REQUIRED",
                "resource_contract_sha256": resource_contract_sha256(config),
            }
        )
    )
    assert minimum_free_gib_for_stage(config, "prepare") == 300
    assert minimum_free_gib_for_stage(config, "cluster") == 150
    assert validate_resource_preflight(config)["result"] == "PASS_REVIEW_REQUIRED"

    config["resources"]["minimum_ssd_free_gib_by_stage"]["cluster"] = 149
    with pytest.raises(RuntimeError, match="does not match"):
        validate_resource_preflight(config)

    config["resources"]["minimum_ssd_free_gib_by_stage"]["cluster"] = 301
    with pytest.raises(ValueError, match="Invalid SSD"):
        minimum_free_gib_for_stage(config, "cluster")


def test_hvg_exclusions_are_specific_to_forbidden_families() -> None:
    excluded = ["TRDC", "TRDV1", "TRGV9", "TRBJ2-7", "TRAC", "MT-CO1", "RPL3", "RPS18", "IGHM"]
    retained = ["CD3D", "CD4", "NKG7", "GNLY", "KLRD1", "TYROBP", "FCER1G", "MALAT1"]
    assert all(is_hvg_excluded(gene) for gene in excluded)
    assert not any(is_hvg_excluded(gene) for gene in retained)
    common = common_eligible_genes([set(excluded + retained), set(excluded + retained + ["X"])])
    assert common == sorted(retained)


def test_batch_hierarchy_uses_library_then_sample_then_source() -> None:
    obs = pd.DataFrame(
        {
            "source_gse_id": ["A", "A", "B"],
            "library_id": ["L1", "", None],
            "sample_id": ["S1", "S2", ""],
        }
    )
    result = make_integration_batch(obs, ["library_id", "sample_id", "source_gse_id"])
    assert result["integration_batch"].tolist() == [
        "A::library_id::L1",
        "A::sample_id::S2",
        "B::source_gse_id::B",
    ]
    assert result["integration_batch_level"].tolist() == ["library_id", "sample_id", "source_gse_id"]


def test_source_balanced_sampling_is_deterministic_and_capped() -> None:
    sources = np.asarray(["A"] * 100 + ["B"] * 9 + ["C"] * 31)
    first = source_balanced_sample_indices(sources, cap_per_source=10, seed=11)
    second = source_balanced_sample_indices(sources, cap_per_source=10, seed=11)
    assert np.array_equal(first, second)
    counts = pd.Series(sources[first]).value_counts().to_dict()
    assert counts == {"A": 10, "B": 9, "C": 10}


@pytest.mark.parametrize("matrix_format", ["csr", "csc"])
def test_sparse_h5ad_gene_and_row_extraction(tmp_path: Path, matrix_format: str) -> None:
    dense = np.arange(60, dtype=np.float32).reshape(10, 6)
    matrix = sparse.csr_matrix(dense) if matrix_format == "csr" else sparse.csc_matrix(dense)
    adata = ad.AnnData(X=matrix, var=pd.DataFrame(index=["A", "B", "C", "D", "E", "F"]))
    path = tmp_path / f"{matrix_format}.h5ad"
    adata.write_h5ad(path)
    rows = np.asarray([0, 3, 7, 9], dtype=np.int64)
    observed = read_sparse_rows_genes(path, "X", ["F", "B", "D"], rows=rows, row_chunk_size=3)
    expected = dense[rows][:, [5, 1, 3]]
    assert sparse.isspmatrix_csr(observed)
    assert np.array_equal(observed.toarray(), expected)


def test_boundary_uses_only_independent_anchor_mixing() -> None:
    partitions = np.asarray(
        [
            [0, 0], [0, 0], [0, 0], [1, 2], [1, 2], [2, 1], [2, 1],
        ],
        dtype=np.int32,
    )
    nk = np.asarray([True, False, False, True, False, False, False])
    tcr = np.asarray([False, True, False, False, False, True, False])
    mask, counts = mixed_boundary_mask(partitions, nk, tcr, minimum_mixed_runs=1)
    assert mask[:3].all()
    assert not mask[3:].any()
    assert counts[:3].tolist() == [2, 2, 2]


def test_consensus_requires_anchor_purity_tcr_control_and_source_spread() -> None:
    n_cells = 90
    partitions = np.tile(np.repeat([0, 1, 2], 30)[:, None], (1, 5)).astype(np.int32)
    nk = np.zeros(n_cells, dtype=bool)
    tcr = np.zeros(n_cells, dtype=bool)
    candidate = np.zeros(n_cells, dtype=bool)
    nk[:20] = True
    candidate[20:30] = True
    tcr[30:50] = True
    candidate[50:60] = True
    nk[60:80] = True
    candidate[80:90] = True
    sources = np.asarray(["anchor"] * n_cells, dtype=object)
    sources[20:30] = np.asarray(["A", "B", "C", "A", "B", "C", "A", "B", "C", "A"])
    sources[50:60] = np.asarray(["A", "B", "C", "A", "B", "C", "A", "B", "C", "A"])
    sources[80:90] = "A"
    cell, clusters = cluster_consensus_votes(
        partitions, nk, tcr, candidate, sources, pseudolabel_contract()
    )
    assert cell.loc[20:29, "selected_pseudo_nk"].all()
    assert not cell.loc[50:59, "selected_pseudo_nk"].any()
    assert not cell.loc[80:89, "selected_pseudo_nk"].any()
    run0 = clusters[clusters["run"].eq("run_0")].set_index("cluster")
    assert bool(run0.loc[0, "qualifies_as_nk"])
    assert not bool(run0.loc[1, "qualifies_as_nk"])
    assert not bool(run0.loc[2, "qualifies_as_nk"])


def test_marker_threshold_or_guardrail_use_is_rejected() -> None:
    contract = pseudolabel_contract()
    contract["marker_thresholds_may_define_truth"] = True
    with pytest.raises(ValueError, match="marker thresholds"):
        validate_pseudolabel_contract(contract)
    contract = pseudolabel_contract()
    contract["may_control_validation_guardrails"] = True
    with pytest.raises(ValueError, match="validation guardrails"):
        validate_pseudolabel_contract(contract)
