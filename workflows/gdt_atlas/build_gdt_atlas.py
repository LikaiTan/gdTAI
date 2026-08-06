#!/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python
"""Build a de novo gdT atlas from an approved rebuilt integrated object.

The workflow deliberately fails before reading expression data unless the
rebuilt H5AD, approval marker, four-profile gdTAI decision, passed sealed
holdout, and all pinned digests agree. It never mutates the rebuilt input.
"""

# ruff: noqa: E402

from __future__ import annotations

# TNK_WORKFLOW_BOOTSTRAP
import sys as _tnk_sys
from pathlib import Path as _TnkPath

_TNK_PROJECT_ROOT = _TnkPath(__file__).resolve().parents[2]
for _tnk_path in (
    _TNK_PROJECT_ROOT,
    _TNK_PROJECT_ROOT / "src",
    _TNK_PROJECT_ROOT / "workflows" / "gdt_atlas",
):
    _tnk_value = str(_tnk_path)
    if _tnk_value not in _tnk_sys.path:
        _tnk_sys.path.insert(0, _tnk_value)

import argparse
import hashlib
import importlib.util
import json
import logging
import os
import re
import shutil
import sys
import uuid
from dataclasses import dataclass
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path
from types import ModuleType
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd
from scipy import sparse

from tnk_atlas.provenance import atomic_write_json, sha256_file


PROJECT_ROOT = _TNK_PROJECT_ROOT
DEFAULT_CONFIG = PROJECT_ROOT / "configs/gdt_atlas/atlas_rebuild.json"
TRUE_STRINGS = frozenset({"true", "1", "yes", "y", "t"})
FALSE_STRINGS = frozenset({"false", "0", "no", "n", "f"})
SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")
SUPPORTED_ENGINEERED_FEATURES = frozenset(
    {
        "any_TRDV",
        "any_TRDJ",
        "any_TRG",
        "any_ab_TCR_gene",
        "TRDC_only",
        "TRDC_plus_TRDV",
        "TRDC_plus_TRDJ",
        "CD3_score",
        "NK_score",
        "gdT_TCR_score",
        "abT_TCR_score",
        "NK_minus_CD3_score",
        "TRDC_log1p",
        "TRDV_score",
        "TRDJ_score",
        "TRG_score",
    }
)
FROZEN_NK_CONTROL_GENES = {
    "NKG7",
    "GNLY",
    "PRF1",
    "GZMB",
    "GZMH",
    "KLRD1",
    "KLRF1",
    "FCGR3A",
    "NCAM1",
    "TYROBP",
    "FCER1G",
    "CST7",
}
REQUIRED_FROZEN_PROFILE_SPECS = {
    "v2_high_f1": ("gdtai_v2", "high_f1"),
    "v2_high_purity": ("gdtai_v2", "high_purity"),
    "v3_round14_balanced": ("gdtai_v3", "fixed_threshold"),
    "v3_round12_high_purity": ("gdtai_v3_round12", "fixed_threshold"),
}
REQUIRED_FROZEN_PROFILE_IDS = frozenset(REQUIRED_FROZEN_PROFILE_SPECS)
REQUIRED_SEALED_HOLDOUT_IDS = frozenset({"GSE315928", "GSE121636"})
REQUIRED_SELECTION_GUARDRAIL_CHECKS = frozenset(
    {
        "gold_recall_pass",
        "pooled_abt_fpr_pass",
        "pooled_strict_nk_fpr_pass",
        "large_cohort_fpr_pass",
        "labeled_fp_fraction_pass",
    }
)


class AtlasConfigError(ValueError):
    """Raised when the atlas config weakens or contradicts the contract."""


class PrerequisiteError(RuntimeError):
    """Raised when an input, approval, profile, or completed build is absent."""


@dataclass(frozen=True)
class FrozenProfileSelection:
    comparison_workflow: Path
    comparison_workflow_sha256: str
    evaluation_config: Path
    evaluation_config_sha256: str
    evaluation_input_manifest: Path
    evaluation_input_manifest_sha256: str
    selection_decision: Path
    selection_decision_sha256: str
    selection_decision_file_sha256: str
    holdout_status: Path
    holdout_status_sha256: str
    selected_profile_id: str
    selected_model_id: str
    selected_mode: str
    selected_model_path: Path
    selected_model_sha256: str
    threshold_contract: dict[str, Any]
    threshold_contract_sha256: str
    decision: dict[str, Any]
    holdout: dict[str, Any]
    evaluator: ModuleType
    profiles: tuple[Any, ...]
    selected_profile: Any


@dataclass(frozen=True)
class ApprovedInput:
    input_h5ad: Path
    approval_marker: Path
    input_sha256: str
    input_size_bytes: int
    config_sha256: str
    frozen_selection: FrozenProfileSelection
    approval: dict[str, Any]
    h5ad_contract: dict[str, Any]

    @property
    def model_path(self) -> Path:
        return self.frozen_selection.selected_model_path

    @property
    def model_sha256(self) -> str:
        return self.frozen_selection.selected_model_sha256


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument(
        "--preflight-only",
        action="store_true",
        help="Validate every prerequisite and print the contract without building the atlas.",
    )
    return parser.parse_args(argv)


def _require(mapping: Mapping[str, Any], key: str, context: str) -> Any:
    if key not in mapping:
        raise AtlasConfigError(f"Missing `{context}.{key}` in atlas config.")
    return mapping[key]


def _relative_project_path(value: str, *, context: str) -> Path:
    path = Path(value)
    if path.is_absolute() or ".." in path.parts:
        raise AtlasConfigError(
            f"`{context}` must be a project-relative path without `..`: {value!r}"
        )
    return path


def resolve_path(value: str | Path, root: Path = PROJECT_ROOT) -> Path:
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (root / path).resolve()


def load_config(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    config_path = path.expanduser().resolve()
    if not config_path.is_file():
        raise FileNotFoundError(config_path)
    try:
        config = json.loads(config_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise AtlasConfigError(f"Invalid JSON in {config_path}: {exc}") from exc
    if not isinstance(config, dict):
        raise AtlasConfigError("Atlas config must be a JSON object.")
    validate_config(config)
    return config


def validate_config(config: Mapping[str, Any]) -> None:
    if config.get("schema_version") != 1:
        raise AtlasConfigError("Only atlas config schema_version 1 is supported.")
    workflow_id = str(config.get("workflow_id", "")).strip()
    if not workflow_id:
        raise AtlasConfigError("`workflow_id` must be nonempty.")

    execution = _require(config, "execution", "config")
    prerequisites = _require(config, "prerequisites", "config")
    profile_selection = _require(config, "frozen_profile_selection", "config")
    selection = _require(config, "selection", "config")
    feature_selection = _require(config, "feature_selection", "config")
    scvi_config = _require(config, "scvi", "config")
    clustering = _require(config, "embedding_clustering", "config")
    analysis = _require(config, "analysis", "config")
    for name, section in [
        ("execution", execution),
        ("prerequisites", prerequisites),
        ("frozen_profile_selection", profile_selection),
        ("selection", selection),
        ("feature_selection", feature_selection),
        ("scvi", scvi_config),
        ("embedding_clustering", clustering),
        ("analysis", analysis),
    ]:
        if not isinstance(section, Mapping):
            raise AtlasConfigError(f"`{name}` must be a JSON object.")

    path_keys = [
        "rebuilt_integrated_h5ad",
        "approval_marker",
        "atlas_h5ad",
        "build_manifest",
        "build_table_dir",
        "build_log_dir",
        "scvi_model_dir",
        "analysis_manifest",
        "analysis_table_dir",
    ]
    for key in path_keys:
        _relative_project_path(
            str(_require(execution, key, "execution")), context=f"execution.{key}"
        )
    if bool(execution.get("allow_overwrite", False)):
        raise AtlasConfigError(
            "`execution.allow_overwrite` must remain false for the approved atlas workflow."
        )
    if execution["rebuilt_integrated_h5ad"] == execution["atlas_h5ad"]:
        raise AtlasConfigError("The rebuilt input and atlas output paths must differ.")

    if prerequisites.get("digest_algorithm") != "sha256":
        raise AtlasConfigError("The prerequisite digest algorithm must be full SHA256.")
    for key in [
        "require_config_sha256",
        "require_input_sha256",
        "require_input_size_bytes",
        "require_csr_counts",
        "require_nonnegative_integer_counts",
        "require_unique_obs_names",
        "require_unique_var_names",
    ]:
        if prerequisites.get(key) is not True:
            raise AtlasConfigError(f"`prerequisites.{key}` must remain true.")

    profile_path_keys = [
        "comparison_workflow",
        "evaluation_config",
        "evaluation_input_manifest",
        "selection_decision",
        "holdout_status",
    ]
    for key in profile_path_keys:
        _relative_project_path(
            str(_require(profile_selection, key, "frozen_profile_selection")),
            context=f"frozen_profile_selection.{key}",
        )
    forbidden_selected_values = {
        "profile_id",
        "model_id",
        "model_path",
        "model_sha256",
        "threshold",
    }
    present_forbidden = sorted(forbidden_selected_values & set(profile_selection))
    if present_forbidden:
        raise AtlasConfigError(
            "The atlas config must not preselect a gdTAI profile or threshold; "
            f"remove: {present_forbidden}."
        )
    required_profile_rows = profile_selection.get("required_profiles")
    if not isinstance(required_profile_rows, list):
        raise AtlasConfigError(
            "`frozen_profile_selection.required_profiles` must be a list."
        )
    configured_profile_specs = {
        str(row.get("profile_id", "")): (
            str(row.get("model_id", "")),
            str(row.get("mode", "")),
        )
        for row in required_profile_rows
        if isinstance(row, Mapping)
    }
    if (
        len(required_profile_rows) != len(REQUIRED_FROZEN_PROFILE_SPECS)
        or configured_profile_specs != REQUIRED_FROZEN_PROFILE_SPECS
    ):
        raise AtlasConfigError(
            "The atlas must wait for the exact four-profile V2/V3 model-mode comparison."
        )
    required_holdouts = {
        str(value) for value in profile_selection.get("required_sealed_holdout_ids", [])
    }
    if required_holdouts != REQUIRED_SEALED_HOLDOUT_IDS:
        raise AtlasConfigError(
            "The atlas gate must retain both prescribed sealed holdouts."
        )
    if profile_selection.get("required_selection_status") != "selected":
        raise AtlasConfigError(
            "`frozen_profile_selection.required_selection_status` must be `selected`."
        )
    if profile_selection.get("required_holdout_promotion_status") != "holdout_passed":
        raise AtlasConfigError("The frozen-profile gate must require `holdout_passed`.")
    if profile_selection.get("required_normalization") != (
        "log1p(raw_counts_per_10000)"
    ) or not np.isclose(
        float(profile_selection.get("normalization_target_sum", np.nan)), 10_000.0
    ):
        raise AtlasConfigError(
            "Selected-profile inference must retain log1p raw counts per 10,000."
        )
    if profile_selection.get("counts_layer") != prerequisites.get("counts_layer"):
        raise AtlasConfigError(
            "Selected-profile inference and scVI prerequisites must use the same raw-count layer."
        )
    if int(profile_selection.get("chunk_size", 0)) <= 0:
        raise AtlasConfigError(
            "Selected-profile inference chunk size must be positive."
        )
    if int(profile_selection.get("maximum_missing_model_genes", -1)) < 0:
        raise AtlasConfigError(
            "Maximum missing selected-model genes cannot be negative."
        )
    output_columns = [
        str(profile_selection.get(key, "")).strip()
        for key in [
            "score_obs_column",
            "threshold_obs_column",
            "prediction_obs_column",
            "annotation_obs_column",
            "annotation_source_obs_column",
        ]
    ]
    if any(not value for value in output_columns) or len(set(output_columns)) != len(
        output_columns
    ):
        raise AtlasConfigError(
            "Selected-profile score, threshold, prediction, and annotation output columns "
            "must be nonempty and distinct."
        )

    if selection.get("allow_silver_only_fn_addback") is not False:
        raise AtlasConfigError("Silver-only FN add-back is prohibited.")
    if selection.get("hard_trdv_exclusion") is not False:
        raise AtlasConfigError("Hard TRDV exclusion is prohibited.")
    if selection.get("nk_gene_only_exclusion") is not False:
        raise AtlasConfigError("NK-gene-alone exclusion is prohibited.")
    if selection.get("require_both_productive_cdr3_for_gold_addback") is not True:
        raise AtlasConfigError(
            "Gold add-back must require both productive TRG and TRD CDR3 fields."
        )
    if selection.get("hard_exclusions_apply_to_gold_addback") is not True:
        raise AtlasConfigError(
            "Hard exclusions must also apply to primary-gold add-back candidates."
        )
    if int(selection.get("minimum_selected_cells", 0)) <= 0:
        raise AtlasConfigError("`selection.minimum_selected_cells` must be positive.")
    required_selector_names = {
        "paired_alpha_beta",
        "any_alpha_beta",
        "paired_gamma_delta",
        "trg_cdr3",
        "trd_cdr3",
        "truth_class",
        "high_confidence_nk",
        "doublet",
        "non_t_contaminant",
        "severe_qc_failure",
    }
    columns = selection.get("columns", {})
    if not isinstance(columns, Mapping) or not required_selector_names.issubset(
        columns
    ):
        missing = sorted(required_selector_names - set(columns))
        raise AtlasConfigError(f"Missing required selection column mappings: {missing}")
    hard_flag_columns = {
        str(columns[name])
        for name in [
            "high_confidence_nk",
            "doublet",
            "non_t_contaminant",
            "severe_qc_failure",
        ]
    }
    if any(
        re.match(r"^(NKG7|GNLY|KLRD1|NCR1|FCER1G|TYROBP)$", column, re.I)
        for column in hard_flag_columns
    ):
        raise AtlasConfigError(
            "Selection hard flags must be upstream decisions, not individual NK genes."
        )
    if (
        not selection.get("primary_gold_values")
        or not selection.get("alpha_beta_gold_values")
        or not selection.get("silver_only_values")
    ):
        raise AtlasConfigError(
            "Primary-gold, alpha-beta-gold, and silver-only truth values must be configured."
        )
    if selection.get("known_alpha_beta_only_requires_paired_trab_or_gold") is not True:
        raise AtlasConfigError(
            "Known alpha-beta-only exclusion must require paired TRA/TRB or reviewed alpha-beta gold."
        )
    if (
        selection.get("known_alpha_beta_only_requires_no_productive_gdt_evidence")
        is not True
    ):
        raise AtlasConfigError(
            "Known alpha-beta-only exclusion must require absence of productive gamma-delta evidence."
        )

    expected_tcr_prefixes = {"TRA", "TRB", "TRG", "TRD"}
    hvg_prefixes = {
        str(value).upper()
        for value in feature_selection.get("tcr_gene_prefixes_excluded", [])
    }
    de_config = analysis.get("donor_pseudobulk_de", {})
    de_prefixes = {
        str(value).upper() for value in de_config.get("tcr_gene_prefixes_excluded", [])
    }
    if not expected_tcr_prefixes.issubset(hvg_prefixes):
        raise AtlasConfigError(
            "All TRA/TRB/TRG/TRD genes must be excluded from HVG selection."
        )
    if not expected_tcr_prefixes.issubset(de_prefixes):
        raise AtlasConfigError(
            "All TRA/TRB/TRG/TRD genes must be excluded from primary DE."
        )
    report_prefixes = {
        str(value).upper()
        for value in analysis.get("separate_tcr_expression", {}).get(
            "gene_prefixes", []
        )
    }
    if not {"TRDV", "TRGV"}.issubset(report_prefixes):
        raise AtlasConfigError("TRDV and TRGV expression must be reported separately.")
    targeted_rule = analysis.get("targeted_cohort_exclusion", {})
    if targeted_rule.get("exclude_from_abundance_inference") is not True:
        raise AtlasConfigError(
            "Targeted cohorts must remain excluded from abundance inference."
        )
    if not str(targeted_rule.get("required_boolean_column", "")).strip():
        raise AtlasConfigError(
            "A required upstream targeted-cohort decision column must be configured."
        )
    if (
        de_config.get("enabled") is not True
        or de_config.get("comparison") != "cluster_vs_rest_within_donor"
    ):
        raise AtlasConfigError(
            "Primary DE must remain enabled as donor-paired cluster-vs-rest pseudobulk DE."
        )
    if (
        analysis.get("random_effects_meta_analysis", {}).get("method")
        != "dersimonian_laird"
    ):
        raise AtlasConfigError(
            "The configured cross-study abundance model must use DerSimonian-Laird random effects."
        )

    required_obs = set(selection.get("required_obs_columns", []))
    required_from_contract = set(str(value) for value in columns.values())
    required_from_contract.update(
        str(analysis[key])
        for key in [
            "study_column",
            "sample_column",
            "library_column",
            "donor_column",
            "condition_column",
            "tissue_column",
        ]
    )
    required_from_contract.add(str(targeted_rule["required_boolean_column"]))
    required_from_contract.add(str(feature_selection["batch_key"]))
    required_from_contract.add(str(scvi_config["batch_key"]))
    required_from_contract.update(
        str(value) for value in scvi_config.get("categorical_covariate_keys", [])
    )
    if not required_from_contract.issubset(required_obs):
        raise AtlasConfigError(
            "`selection.required_obs_columns` is missing configured build/analysis metadata: "
            f"{sorted(required_from_contract - required_obs)}"
        )

    seeds = [int(value) for value in clustering.get("seeds", [])]
    resolutions = [float(value) for value in clustering.get("leiden_resolutions", [])]
    if len(set(seeds)) < 2 or len(set(resolutions)) < 2:
        raise AtlasConfigError(
            "At least two UMAP/Leiden seeds and two Leiden resolutions are required."
        )
    if int(clustering.get("primary_seed", -1)) not in seeds:
        raise AtlasConfigError("`primary_seed` must be included in `seeds`.")
    primary_resolution = float(clustering.get("primary_resolution", -1.0))
    if not any(np.isclose(primary_resolution, value) for value in resolutions):
        raise AtlasConfigError(
            "`primary_resolution` must be included in `leiden_resolutions`."
        )
    if int(feature_selection.get("n_top_genes", 0)) <= 0:
        raise AtlasConfigError("`feature_selection.n_top_genes` must be positive.")
    if int(feature_selection.get("candidate_pool_size", 0)) < int(
        feature_selection["n_top_genes"]
    ):
        raise AtlasConfigError(
            "The HVG candidate pool cannot be smaller than the final HVG count."
        )


def validate_runtime_environment(config: Mapping[str, Any]) -> None:
    expected = str(config["execution"]["required_environment"])
    conda_env = os.environ.get("CONDA_DEFAULT_ENV", "")
    executable_parts = set(Path(sys.executable).resolve().parts)
    if conda_env != expected and expected not in executable_parts:
        raise PrerequisiteError(
            f"This workflow requires `{expected}`; current CONDA_DEFAULT_ENV={conda_env!r} "
            f"and executable={sys.executable!r}."
        )


def _parse_approved_at(value: Any) -> datetime:
    text = str(value).strip().replace("Z", "+00:00")
    try:
        parsed = datetime.fromisoformat(text)
    except ValueError as exc:
        raise PrerequisiteError(
            "Approval marker `approved_at_utc` must be an ISO-8601 timestamp."
        ) from exc
    if parsed.tzinfo is None:
        raise PrerequisiteError(
            "Approval marker `approved_at_utc` must include a timezone."
        )
    return parsed.astimezone(timezone.utc)


def _stable_json_sha256(value: Any) -> str:
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _read_json_object(path: Path, label: str) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise PrerequisiteError(f"Cannot read {label} {path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise PrerequisiteError(f"{label.capitalize()} must be a JSON object.")
    return payload


def verify_selection_decision_digest(
    decision: Mapping[str, Any], required_status: str = "selected"
) -> str:
    digest = str(decision.get("decision_sha256", ""))
    if not SHA256_PATTERN.fullmatch(digest):
        raise PrerequisiteError(
            "The gdTAI selection decision lacks a valid canonical SHA256 digest."
        )
    unsigned = dict(decision)
    unsigned.pop("decision_sha256", None)
    observed = _stable_json_sha256(unsigned)
    if observed != digest:
        raise PrerequisiteError(
            "The gdTAI selection decision digest is invalid; atlas preflight is closed."
        )
    if decision.get("status") != required_status:
        raise PrerequisiteError(
            f"The gdTAI selection decision did not pass: {decision.get('status')!r}."
        )
    selected = str(decision.get("selected_profile", "")).strip()
    identities = decision.get("profile_identity")
    if (
        not selected
        or not isinstance(identities, Mapping)
        or selected not in identities
    ):
        raise PrerequisiteError(
            "The gdTAI selection decision does not identify a selected frozen profile."
        )
    ranking = decision.get("selection_ranking")
    if not isinstance(ranking, list):
        raise PrerequisiteError("The gdTAI selection decision lacks its ranking audit.")
    selected_rows = [
        row
        for row in ranking
        if isinstance(row, Mapping) and row.get("selected") is True
    ]
    if (
        len(selected_rows) != 1
        or str(selected_rows[0].get("profile_id", "")) != selected
        or selected_rows[0].get("eligible") is not True
        or any(
            selected_rows[0].get(key) is not True
            for key in REQUIRED_SELECTION_GUARDRAIL_CHECKS
        )
    ):
        raise PrerequisiteError(
            "The checksum-pinned ranking does not contain exactly one selected profile "
            "that passed every T/NK guardrail."
        )
    ranking_profile_ids = {
        str(row.get("profile_id", "")) for row in ranking if isinstance(row, Mapping)
    }
    if ranking_profile_ids != set(identities):
        raise PrerequisiteError(
            "The checksum-pinned ranking does not cover every frozen profile identity."
        )
    return digest


def validate_passed_holdout_status(
    decision: Mapping[str, Any],
    holdout: Mapping[str, Any],
    required_promotion_status: str = "holdout_passed",
) -> None:
    selected = str(decision.get("selected_profile", ""))
    expected_digest = str(decision.get("decision_sha256", ""))
    if holdout.get("selected_profile") != selected:
        raise PrerequisiteError(
            "The sealed-holdout status names a different selected gdTAI profile."
        )
    if holdout.get("selection_decision_sha256") != expected_digest:
        raise PrerequisiteError(
            "The sealed-holdout status is not linked to the current checksum-pinned selection decision."
        )
    if holdout.get("holdout_schema_failure") is not False:
        raise PrerequisiteError("The sealed holdout has a schema failure.")
    if holdout.get("holdout_veto") is not False:
        raise PrerequisiteError("The sealed holdout vetoed the selected gdTAI profile.")
    if holdout.get("promotion_status") != required_promotion_status:
        raise PrerequisiteError(
            "The sealed holdout has not passed; atlas preflight remains closed."
        )


def _load_gdtai_evaluator(path: Path, workflow_sha256: str) -> ModuleType:
    module_name = f"_atlas_frozen_profile_evaluator_{workflow_sha256}"
    cached = sys.modules.get(module_name)
    if isinstance(cached, ModuleType):
        return cached
    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise PrerequisiteError(f"Cannot import gdTAI comparison workflow {path}.")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    try:
        spec.loader.exec_module(module)
    except Exception as exc:
        sys.modules.pop(module_name, None)
        raise PrerequisiteError(
            f"Cannot load gdTAI comparison workflow {path}: {exc}"
        ) from exc
    if sha256_file(path) != workflow_sha256:
        sys.modules.pop(module_name, None)
        raise PrerequisiteError(
            "The gdTAI comparison workflow changed while it was being loaded."
        )
    required = {
        "CANONICAL_OBS_COLUMNS",
        "FeatureSpec",
        "feature_schema",
        "load_profiles",
        "normalize_annotation",
        "predict_profiles",
        "profile_identity",
        "v2_thresholds",
        "verify_decision",
    }
    missing = sorted(name for name in required if not hasattr(module, name))
    if missing:
        raise PrerequisiteError(
            f"gdTAI comparison workflow lacks required interfaces: {missing}"
        )
    return module


def _probability_threshold(value: Any, context: str) -> float:
    try:
        threshold = float(value)
    except (TypeError, ValueError) as exc:
        raise PrerequisiteError(f"{context} is not a numeric threshold.") from exc
    if not np.isfinite(threshold) or not 0.0 < threshold <= 1.0:
        raise PrerequisiteError(f"{context} must be in (0, 1].")
    return threshold


def _threshold_contract(profile: Any) -> dict[str, Any]:
    if profile.model_id == "gdtai_v2":
        operating_modes = profile.payload.get("operating_modes", {})
        if profile.mode not in operating_modes:
            raise PrerequisiteError(
                f"Selected V2 mode {profile.mode!r} is absent from its frozen artifact."
            )
        mode = operating_modes[profile.mode]
        if profile.mode == "high_f1":
            return {
                "kind": "fixed",
                "threshold": _probability_threshold(
                    mode.get("threshold"), "Selected V2 high-F1 threshold"
                ),
            }
        if profile.mode != "high_purity":
            raise PrerequisiteError(f"Unsupported selected V2 mode: {profile.mode!r}.")
        rules = mode.get("annotation_thresholds")
        required = {
            "other_threshold",
            "gdt_threshold",
            "cd8_threshold",
            "cd4_threshold",
            "treg_threshold",
            "nk_threshold",
        }
        if not isinstance(rules, Mapping) or not required.issubset(rules):
            raise PrerequisiteError(
                "Selected V2 high-purity artifact lacks annotation-specific thresholds."
            )
        normalized: dict[str, float | str] = {}
        for key in sorted(required):
            value = rules[key]
            normalized[key] = (
                "disabled"
                if str(value) == "disabled"
                else _probability_threshold(value, f"Selected V2 threshold {key}")
            )
        return {"kind": "annotation_specific", "thresholds": normalized}
    payload = profile.payload
    return {
        "kind": "fixed",
        "threshold": _probability_threshold(
            payload.get("threshold"), "Selected V3 fixed threshold"
        ),
    }


def selected_profile_thresholds(
    evaluator: ModuleType, profile: Any, annotations: np.ndarray
) -> np.ndarray:
    annotations = np.asarray(annotations, dtype=object)
    if profile.model_id == "gdtai_v2":
        thresholds = np.asarray(
            evaluator.v2_thresholds(profile, annotations), dtype=np.float32
        )
    else:
        thresholds = np.full(
            annotations.size,
            _probability_threshold(
                profile.payload.get("threshold"), "Selected V3 fixed threshold"
            ),
            dtype=np.float32,
        )
    valid = (
        np.isfinite(thresholds) & (thresholds > 0.0) & (thresholds <= 1.0)
    ) | np.isposinf(thresholds)
    if thresholds.shape != annotations.shape or not valid.all():
        raise PrerequisiteError(
            "The selected frozen profile produced invalid or misaligned thresholds."
        )
    return thresholds


def validate_frozen_profile_selection(
    config: Mapping[str, Any], root: Path = PROJECT_ROOT
) -> FrozenProfileSelection:
    selection = config["frozen_profile_selection"]
    paths = {
        key: resolve_path(selection[key], root)
        for key in [
            "comparison_workflow",
            "evaluation_config",
            "evaluation_input_manifest",
            "selection_decision",
            "holdout_status",
        ]
    }
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise PrerequisiteError(
            "Four-profile gdTAI decision prerequisites are absent; atlas preflight is closed: "
            f"{missing}"
        )
    for key in ("selection_decision", "holdout_status"):
        if paths[key].is_symlink():
            raise PrerequisiteError(
                f"The {key.replace('_', ' ')} must be a regular file, not a symlink."
            )

    evaluation_config = _read_json_object(
        paths["evaluation_config"], "gdTAI evaluation config"
    )
    profile_rows = evaluation_config.get("profiles")
    if not isinstance(profile_rows, list):
        raise PrerequisiteError("gdTAI evaluation config lacks a profile list.")
    profile_specs = {
        str(row.get("profile_id", "")): (
            str(row.get("model_id", "")),
            str(row.get("mode", "")),
        )
        for row in profile_rows
        if isinstance(row, Mapping)
    }
    required_specs = {
        str(row["profile_id"]): (str(row["model_id"]), str(row["mode"]))
        for row in selection["required_profiles"]
    }
    if (
        len(profile_rows) != len(REQUIRED_FROZEN_PROFILE_IDS)
        or profile_specs != required_specs
    ):
        raise PrerequisiteError(
            "The upstream evaluation is not the required four-profile V2/V3 comparison."
        )
    configured_holdouts = set(evaluation_config.get("sealed_holdouts", {}))
    if configured_holdouts != set(selection["required_sealed_holdout_ids"]):
        raise PrerequisiteError(
            "The upstream evaluation config does not contain the required sealed holdouts."
        )
    if evaluation_config.get("normalization") != selection["required_normalization"]:
        raise PrerequisiteError(
            "The upstream evaluation normalization differs from the atlas inference contract."
        )
    try:
        manifest = pd.read_csv(paths["evaluation_input_manifest"], dtype=str).fillna("")
    except Exception as exc:
        raise PrerequisiteError(
            f"Cannot read gdTAI evaluation input manifest: {exc}"
        ) from exc
    required_manifest_columns = {"dataset_id", "h5ad_path", "sealed_holdout"}
    if not required_manifest_columns.issubset(manifest.columns):
        raise PrerequisiteError(
            "gdTAI evaluation input manifest lacks required decision/holdout columns."
        )
    sealed = manifest["sealed_holdout"].str.strip().str.lower().isin(TRUE_STRINGS)
    manifest_holdouts = set(manifest.loc[sealed, "dataset_id"].astype(str))
    if manifest_holdouts != set(selection["required_sealed_holdout_ids"]):
        raise PrerequisiteError(
            "gdTAI evaluation input manifest does not seal exactly the required holdout cohorts."
        )

    decision = _read_json_object(
        paths["selection_decision"], "gdTAI selection decision"
    )
    decision_sha = verify_selection_decision_digest(
        decision, str(selection["required_selection_status"])
    )
    if (
        resolve_path(str(decision.get("config_path", "")), root)
        != paths["evaluation_config"]
    ):
        raise PrerequisiteError(
            "The checksum-pinned decision points to a different gdTAI evaluation config."
        )
    if (
        resolve_path(str(decision.get("input_manifest_path", "")), root)
        != paths["evaluation_input_manifest"]
    ):
        raise PrerequisiteError(
            "The checksum-pinned decision points to a different evaluation input manifest."
        )

    workflow_sha = sha256_file(paths["comparison_workflow"])
    evaluator = _load_gdtai_evaluator(paths["comparison_workflow"], workflow_sha)
    try:
        profiles = tuple(evaluator.load_profiles(evaluation_config))
        verified_decision = evaluator.verify_decision(
            paths["selection_decision"],
            list(profiles),
            paths["evaluation_config"],
            paths["evaluation_input_manifest"],
        )
    except Exception as exc:
        raise PrerequisiteError(
            f"The frozen-profile decision failed evaluator verification: {exc}"
        ) from exc
    for profile in profiles:
        observed_model_sha = sha256_file(Path(profile.artifact).resolve())
        if observed_model_sha != str(profile.sha256):
            raise PrerequisiteError(
                f"Frozen profile {profile.profile_id!r} changed while being loaded."
            )
    if verified_decision.get("decision_sha256") != decision_sha:
        raise PrerequisiteError(
            "Evaluator verification returned a different decision digest."
        )
    if set(evaluator.profile_identity(list(profiles))) != set(profile_specs):
        raise PrerequisiteError(
            "Loaded frozen profile identities differ from the four-profile decision."
        )

    holdout = _read_json_object(paths["holdout_status"], "sealed-holdout status")
    validate_passed_holdout_status(
        decision,
        holdout,
        str(selection["required_holdout_promotion_status"]),
    )
    selected_profile_id = str(decision["selected_profile"])
    matches = [p for p in profiles if p.profile_id == selected_profile_id]
    if len(matches) != 1:
        raise PrerequisiteError(
            "The selected frozen profile is not unique in the loaded evaluation profiles."
        )
    selected_profile = matches[0]
    threshold_contract = _threshold_contract(selected_profile)
    return FrozenProfileSelection(
        comparison_workflow=paths["comparison_workflow"],
        comparison_workflow_sha256=workflow_sha,
        evaluation_config=paths["evaluation_config"],
        evaluation_config_sha256=sha256_file(paths["evaluation_config"]),
        evaluation_input_manifest=paths["evaluation_input_manifest"],
        evaluation_input_manifest_sha256=sha256_file(
            paths["evaluation_input_manifest"]
        ),
        selection_decision=paths["selection_decision"],
        selection_decision_sha256=decision_sha,
        selection_decision_file_sha256=sha256_file(paths["selection_decision"]),
        holdout_status=paths["holdout_status"],
        holdout_status_sha256=sha256_file(paths["holdout_status"]),
        selected_profile_id=selected_profile_id,
        selected_model_id=str(selected_profile.model_id),
        selected_mode=str(selected_profile.mode),
        selected_model_path=Path(selected_profile.artifact).resolve(),
        selected_model_sha256=str(selected_profile.sha256),
        threshold_contract=threshold_contract,
        threshold_contract_sha256=_stable_json_sha256(threshold_contract),
        decision=decision,
        holdout=holdout,
        evaluator=evaluator,
        profiles=profiles,
        selected_profile=selected_profile,
    )


def validate_approval_marker(
    config: Mapping[str, Any],
    config_path: Path,
    frozen_selection: FrozenProfileSelection,
    root: Path = PROJECT_ROOT,
) -> tuple[Path, Path, str, int, str, dict[str, Any]]:
    execution = config["execution"]
    prerequisites = config["prerequisites"]
    input_h5ad = resolve_path(execution["rebuilt_integrated_h5ad"], root)
    marker_path = resolve_path(execution["approval_marker"], root)
    missing = []
    if not input_h5ad.is_file():
        missing.append(f"rebuilt integrated H5AD: {input_h5ad}")
    if not marker_path.is_file():
        missing.append(f"approval marker: {marker_path}")
    if missing:
        raise PrerequisiteError(
            "Required approved rebuild prerequisites are absent: " + "; ".join(missing)
        )
    if marker_path.is_symlink():
        raise PrerequisiteError(
            "The approval marker must be a regular file, not a symlink."
        )
    try:
        marker = json.loads(marker_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise PrerequisiteError(
            f"Cannot read approval marker {marker_path}: {exc}"
        ) from exc
    if not isinstance(marker, dict):
        raise PrerequisiteError("Approval marker must be a JSON object.")
    required_keys = {
        "schema_version",
        "decision",
        "purpose",
        "workflow_id",
        "approved_by",
        "approved_at_utc",
        "input_h5ad",
        "input_sha256",
        "input_size_bytes",
        "config_sha256",
        "gdtai_comparison_workflow_sha256",
        "gdtai_selection_decision_sha256",
        "gdtai_selection_decision_file_sha256",
        "gdtai_holdout_status_sha256",
        "gdtai_selected_profile_id",
        "gdtai_selected_model_id",
        "gdtai_selected_mode",
        "gdtai_selected_model_sha256",
        "gdtai_threshold_contract_sha256",
    }
    absent = sorted(required_keys - set(marker))
    if absent:
        raise PrerequisiteError(
            f"Approval marker is incomplete; missing keys: {absent}"
        )
    expected_values = {
        "schema_version": prerequisites["approval_schema_version"],
        "decision": prerequisites["approval_decision"],
        "purpose": prerequisites["approval_purpose"],
        "workflow_id": config["workflow_id"],
        "gdtai_comparison_workflow_sha256": frozen_selection.comparison_workflow_sha256,
        "gdtai_selection_decision_sha256": frozen_selection.selection_decision_sha256,
        "gdtai_selection_decision_file_sha256": frozen_selection.selection_decision_file_sha256,
        "gdtai_holdout_status_sha256": frozen_selection.holdout_status_sha256,
        "gdtai_selected_profile_id": frozen_selection.selected_profile_id,
        "gdtai_selected_model_id": frozen_selection.selected_model_id,
        "gdtai_selected_mode": frozen_selection.selected_mode,
        "gdtai_selected_model_sha256": frozen_selection.selected_model_sha256,
        "gdtai_threshold_contract_sha256": frozen_selection.threshold_contract_sha256,
    }
    for key, expected in expected_values.items():
        if marker.get(key) != expected:
            raise PrerequisiteError(
                f"Approval marker `{key}` mismatch: expected {expected!r}, observed {marker.get(key)!r}."
            )
    if not str(marker["approved_by"]).strip():
        raise PrerequisiteError("Approval marker `approved_by` must be nonempty.")
    _parse_approved_at(marker["approved_at_utc"])
    marker_input = resolve_path(str(marker["input_h5ad"]), root)
    if marker_input != input_h5ad:
        raise PrerequisiteError(
            f"Approval marker names a different H5AD: {marker_input}; configured input is {input_h5ad}."
        )
    config_sha = sha256_file(config_path.resolve())
    if marker["config_sha256"] != config_sha:
        raise PrerequisiteError(
            f"Approval marker config SHA256 mismatch: expected current config {config_sha}, "
            f"observed {marker['config_sha256']}."
        )
    stat = input_h5ad.stat()
    if int(marker["input_size_bytes"]) != int(stat.st_size):
        raise PrerequisiteError(
            "Approved input size no longer matches the rebuilt H5AD."
        )
    approved_sha = str(marker["input_sha256"])
    if not SHA256_PATTERN.fullmatch(approved_sha):
        raise PrerequisiteError(
            "Approval marker `input_sha256` is not a lowercase SHA256 digest."
        )
    observed_sha = sha256_file(input_h5ad)
    if observed_sha != approved_sha:
        raise PrerequisiteError(
            f"Approved input SHA256 mismatch: expected {approved_sha}, observed {observed_sha}."
        )
    return input_h5ad, marker_path, observed_sha, int(stat.st_size), config_sha, marker


def _decode_strings(values: np.ndarray) -> np.ndarray:
    flat = np.asarray(values)
    if flat.dtype.kind in {"S", "O"}:
        return (
            pd.Series(flat, copy=False)
            .map(
                lambda value: value.decode("utf-8")
                if isinstance(value, (bytes, np.bytes_))
                else str(value)
            )
            .to_numpy(dtype=object)
        )
    return flat.astype(str)


def _matrix_shape(matrix: Any) -> tuple[int, int]:
    shape = (
        matrix.shape
        if hasattr(matrix, "shape") and len(matrix.shape) == 2
        else matrix.attrs.get("shape")
    )
    if shape is None:
        raise PrerequisiteError(
            "Cannot determine count-matrix shape from the rebuilt H5AD."
        )
    return int(shape[0]), int(shape[1])


def _matrix_encoding(matrix: Any) -> str:
    import h5py

    if isinstance(matrix, h5py.Dataset):
        return "dense"
    value = matrix.attrs.get("encoding-type", "")
    if isinstance(value, bytes):
        value = value.decode("utf-8")
    if not value and {"data", "indices", "indptr"}.issubset(matrix.keys()):
        value = "csr_matrix"
    return str(value)


def _sample_count_values(data: Any, maximum: int = 100_000) -> np.ndarray:
    n_values = int(data.shape[0])
    if n_values == 0:
        return np.asarray([], dtype=np.float64)
    if n_values <= maximum:
        return np.asarray(data[:], dtype=np.float64)
    block = maximum // 3
    pieces = [data[:block], data[n_values // 2 : n_values // 2 + block], data[-block:]]
    return np.concatenate([np.asarray(piece, dtype=np.float64) for piece in pieces])


def inspect_h5ad_contract(path: Path, config: Mapping[str, Any]) -> dict[str, Any]:
    import h5py

    prerequisites = config["prerequisites"]
    required_obs = list(config["selection"]["required_obs_columns"])
    counts_key = str(prerequisites["counts_layer"])
    with h5py.File(path, "r") as handle:
        if "obs" not in handle or "var" not in handle or "layers" not in handle:
            raise PrerequisiteError(
                "Rebuilt input is not a valid H5AD with obs, var, and layers groups."
            )
        if counts_key not in handle["layers"]:
            raise PrerequisiteError(
                f"Rebuilt input must contain layers[{counts_key!r}] with raw counts."
            )
        missing_obs = sorted(set(required_obs) - set(handle["obs"].keys()))
        if missing_obs:
            raise PrerequisiteError(
                f"Rebuilt input is missing required obs columns: {missing_obs}"
            )
        if "_index" not in handle["obs"] or "_index" not in handle["var"]:
            raise PrerequisiteError(
                "Rebuilt input must contain explicit obs and var indexes."
            )
        obs_names = _decode_strings(handle["obs"]["_index"][:])
        var_names = _decode_strings(handle["var"]["_index"][:])
        counts = handle["layers"][counts_key]
        n_obs, n_vars = _matrix_shape(counts)
        if n_obs != len(obs_names) or n_vars != len(var_names):
            raise PrerequisiteError(
                "Counts-layer dimensions do not match obs/var indexes."
            )
        if n_obs < int(prerequisites["minimum_cells"]) or n_vars < int(
            prerequisites["minimum_genes"]
        ):
            raise PrerequisiteError(
                f"Rebuilt input is unexpectedly small ({n_obs:,} cells x {n_vars:,} genes)."
            )
        if (
            prerequisites["require_unique_obs_names"]
            and pd.Index(obs_names).has_duplicates
        ):
            raise PrerequisiteError("Rebuilt input has duplicated obs names.")
        if (
            prerequisites["require_unique_var_names"]
            and pd.Index(var_names).has_duplicates
        ):
            raise PrerequisiteError("Rebuilt input has duplicated var names.")
        encoding = _matrix_encoding(counts)
        if prerequisites["require_csr_counts"] and encoding != "csr_matrix":
            raise PrerequisiteError(f"Raw counts must be CSR; found {encoding!r}.")
        data = counts["data"] if encoding == "csr_matrix" else counts
        sampled = _sample_count_values(data)
        if sampled.size:
            if not np.isfinite(sampled).all() or np.any(sampled < 0):
                raise PrerequisiteError(
                    "Raw counts contain non-finite or negative sampled values."
                )
            if prerequisites["require_nonnegative_integer_counts"] and not np.allclose(
                sampled, np.rint(sampled), rtol=0.0, atol=1e-6
            ):
                raise PrerequisiteError(
                    "Raw counts are not integer-like in the deterministic validation sample."
                )
    return {
        "n_obs": n_obs,
        "n_vars": n_vars,
        "counts_layer": counts_key,
        "counts_encoding": encoding,
        "sampled_count_values": int(sampled.size),
        "required_obs_columns": required_obs,
    }


def validate_prerequisites(
    config: Mapping[str, Any], config_path: Path, root: Path = PROJECT_ROOT
) -> ApprovedInput:
    frozen_selection = validate_frozen_profile_selection(config, root)
    input_h5ad, marker_path, input_sha, input_size, config_sha, marker = (
        validate_approval_marker(config, config_path, frozen_selection, root)
    )
    contract = inspect_h5ad_contract(input_h5ad, config)
    return ApprovedInput(
        input_h5ad=input_h5ad,
        approval_marker=marker_path,
        input_sha256=input_sha,
        input_size_bytes=input_size,
        config_sha256=config_sha,
        frozen_selection=frozen_selection,
        approval=marker,
        h5ad_contract=contract,
    )


def strict_bool_array(
    values: pd.Series | np.ndarray | Sequence[Any], column: str
) -> np.ndarray:
    series = pd.Series(values, copy=False)
    if pd.api.types.is_bool_dtype(series.dtype):
        if series.isna().any():
            raise PrerequisiteError(
                f"Boolean decision column `{column}` contains missing values."
            )
        return series.to_numpy(dtype=bool)
    if pd.api.types.is_numeric_dtype(series.dtype):
        numeric = pd.to_numeric(series, errors="coerce")
        invalid = numeric.isna() | ~numeric.isin([0, 1])
        if invalid.any():
            examples = series.loc[invalid].astype(str).head(5).tolist()
            raise PrerequisiteError(
                f"Boolean decision column `{column}` has invalid values: {examples}"
            )
        return numeric.to_numpy(dtype=np.int8).astype(bool)
    normalized = series.astype("string").str.strip().str.lower()
    valid = normalized.isin(TRUE_STRINGS | FALSE_STRINGS)
    if not bool(valid.all()):
        examples = series.loc[~valid].astype(str).head(5).tolist()
        raise PrerequisiteError(
            f"Boolean decision column `{column}` has invalid or missing values: {examples}"
        )
    return normalized.isin(TRUE_STRINGS).to_numpy(dtype=bool)


def _normalized_truth(values: pd.Series | Sequence[Any]) -> np.ndarray:
    return (
        pd.Series(values, copy=False)
        .astype("string")
        .fillna("")
        .str.strip()
        .str.lower()
        .to_numpy(dtype=object)
    )


def _valid_cdr3(
    values: pd.Series | Sequence[Any], invalid_tokens: Sequence[str]
) -> np.ndarray:
    invalid = {str(value).strip().lower() for value in invalid_tokens}
    cleaned = (
        pd.Series(values, copy=False)
        .astype("string")
        .fillna("")
        .str.strip()
        .str.lower()
    )
    return (~cleaned.isin(invalid)).to_numpy(dtype=bool)


def build_selection_table(
    obs: pd.DataFrame,
    scores: np.ndarray,
    thresholds: np.ndarray,
    config: Mapping[str, Any],
) -> pd.DataFrame:
    selection = config["selection"]
    columns = selection["columns"]
    missing = sorted(set(selection["required_obs_columns"]) - set(obs.columns))
    if missing:
        raise PrerequisiteError(
            f"Selection metadata is missing required columns: {missing}"
        )
    scores = np.asarray(scores, dtype=np.float32)
    if (
        scores.ndim != 1
        or scores.shape[0] != obs.shape[0]
        or not np.isfinite(scores).all()
    ):
        raise PrerequisiteError(
            "Frozen gdTAI scores must be finite and aligned one-to-one with obs rows."
        )
    thresholds = np.asarray(thresholds, dtype=np.float32)
    valid_thresholds = (
        np.isfinite(thresholds) & (thresholds > 0.0) & (thresholds <= 1.0)
    ) | np.isposinf(thresholds)
    if (
        thresholds.ndim != 1
        or thresholds.shape[0] != obs.shape[0]
        or not valid_thresholds.all()
    ):
        raise PrerequisiteError(
            "Selected-profile thresholds must be in (0, 1] or positive infinity "
            "and aligned one-to-one with obs rows."
        )

    paired_ab = strict_bool_array(
        obs[columns["paired_alpha_beta"]], columns["paired_alpha_beta"]
    )
    any_ab = strict_bool_array(
        obs[columns["any_alpha_beta"]], columns["any_alpha_beta"]
    )
    paired_gd = strict_bool_array(
        obs[columns["paired_gamma_delta"]], columns["paired_gamma_delta"]
    )
    high_confidence_nk = strict_bool_array(
        obs[columns["high_confidence_nk"]], columns["high_confidence_nk"]
    )
    doublet = strict_bool_array(obs[columns["doublet"]], columns["doublet"])
    non_t = strict_bool_array(
        obs[columns["non_t_contaminant"]], columns["non_t_contaminant"]
    )
    severe_qc = strict_bool_array(
        obs[columns["severe_qc_failure"]], columns["severe_qc_failure"]
    )
    truth = _normalized_truth(obs[columns["truth_class"]])
    primary_values = {
        str(value).strip().lower() for value in selection["primary_gold_values"]
    }
    alpha_beta_values = {
        str(value).strip().lower() for value in selection["alpha_beta_gold_values"]
    }
    silver_values = {
        str(value).strip().lower() for value in selection["silver_only_values"]
    }
    primary_gold = np.isin(truth, list(primary_values))
    alpha_beta_gold = np.isin(truth, list(alpha_beta_values))
    silver_only = np.isin(truth, list(silver_values))
    if np.any(primary_gold & (alpha_beta_gold | silver_only)) or np.any(
        alpha_beta_gold & silver_only
    ):
        raise PrerequisiteError(
            "Truth labels cannot overlap primary-gold, alpha-beta-gold, and silver-only classes."
        )
    valid_trg = _valid_cdr3(obs[columns["trg_cdr3"]], selection["invalid_cdr3_tokens"])
    valid_trd = _valid_cdr3(obs[columns["trd_cdr3"]], selection["invalid_cdr3_tokens"])

    predicted = scores >= thresholds
    productive_gd_evidence = paired_gd | valid_trg | valid_trd
    strict_alpha_beta_evidence = (paired_ab & any_ab) | alpha_beta_gold
    known_alpha_beta_only = strict_alpha_beta_evidence & ~productive_gd_evidence
    primary_gold_paired = primary_gold & paired_gd & valid_trg & valid_trd
    eligible_gold_fn = (~predicted) & primary_gold_paired
    silver_only_fn = (~predicted) & silver_only
    hard_exclusion = (
        known_alpha_beta_only | high_confidence_nk | doublet | non_t | severe_qc
    )
    selected = (predicted | eligible_gold_fn) & ~hard_exclusion
    gold_fn_added = eligible_gold_fn & ~hard_exclusion

    if np.any(gold_fn_added & predicted):
        raise AssertionError(
            "Gold FN add-back unexpectedly overlaps frozen-profile predictions."
        )
    if np.any(gold_fn_added & ~primary_gold_paired):
        raise AssertionError(
            "Gold FN add-back includes a cell without paired primary-gold TRG/TRD evidence."
        )
    if np.any(selected & silver_only_fn):
        raise AssertionError("A silver-only false negative entered the final atlas.")
    if np.any(selected & hard_exclusion):
        raise AssertionError("A hard-excluded cell entered the final atlas.")

    reason = np.full(obs.shape[0], "profile_negative_not_added", dtype=object)
    reason[predicted] = "frozen_profile_prediction"
    reason[silver_only_fn] = "silver_only_fn_not_added"
    reason[eligible_gold_fn] = "primary_gold_paired_trg_trd_fn_added"
    reason[known_alpha_beta_only] = "excluded_known_alpha_beta_only"
    reason[high_confidence_nk] = "excluded_high_confidence_nk"
    reason[non_t] = "excluded_non_t_contaminant"
    reason[severe_qc] = "excluded_severe_qc_failure"
    reason[doublet] = "excluded_doublet"

    return pd.DataFrame(
        {
            "gdtai_frozen_score": scores,
            "gdtai_frozen_threshold": thresholds,
            "gdtai_frozen_prediction": predicted,
            "selection_primary_gold": primary_gold,
            "selection_alpha_beta_gold": alpha_beta_gold,
            "selection_silver_only": silver_only,
            "selection_any_alpha_beta": any_ab,
            "selection_paired_gamma_delta": paired_gd,
            "selection_valid_trg_cdr3": valid_trg,
            "selection_valid_trd_cdr3": valid_trd,
            "selection_primary_gold_paired_trg_trd": primary_gold_paired,
            "selection_known_alpha_beta_only": known_alpha_beta_only,
            "selection_high_confidence_nk": high_confidence_nk,
            "selection_doublet": doublet,
            "selection_non_t_contaminant": non_t,
            "selection_severe_qc_failure": severe_qc,
            "selection_hard_exclusion": hard_exclusion,
            "selection_eligible_gold_fn": eligible_gold_fn,
            "selection_silver_only_fn_not_added": silver_only_fn,
            "selection_gold_fn_added": gold_fn_added,
            "gdt_atlas_selected": selected,
            "gdt_atlas_selection_reason": reason,
        },
        index=obs.index,
    )


def summarize_selection(selection: pd.DataFrame) -> pd.DataFrame:
    masks = [
        ("input_cells", np.ones(selection.shape[0], dtype=bool)),
        ("frozen_profile_predictions", selection["gdtai_frozen_prediction"]),
        ("known_alpha_beta_only_removed", selection["selection_known_alpha_beta_only"]),
        ("high_confidence_nk_removed", selection["selection_high_confidence_nk"]),
        ("doublets_removed", selection["selection_doublet"]),
        ("non_t_contaminants_removed", selection["selection_non_t_contaminant"]),
        ("severe_qc_failures_removed", selection["selection_severe_qc_failure"]),
        ("silver_only_fn_not_added", selection["selection_silver_only_fn_not_added"]),
        ("primary_gold_paired_trg_trd_fn_added", selection["selection_gold_fn_added"]),
        ("final_atlas_cells", selection["gdt_atlas_selected"]),
    ]
    n_input = max(selection.shape[0], 1)
    return pd.DataFrame(
        {
            "category": [name for name, _ in masks],
            "n_cells": [int(np.asarray(mask, dtype=bool).sum()) for _, mask in masks],
            "fraction_of_input": [
                float(np.asarray(mask, dtype=bool).sum() / n_input) for _, mask in masks
            ],
        }
    )


def _sum_gene_prefixes(
    x_gene: np.ndarray, gene_names: Sequence[str], prefixes: tuple[str, ...]
) -> np.ndarray:
    columns = [
        index
        for index, gene in enumerate(gene_names)
        if str(gene).upper().startswith(prefixes)
    ]
    if not columns:
        return np.zeros(x_gene.shape[0], dtype=np.float32)
    return x_gene[:, columns].sum(axis=1, dtype=np.float32)


def _mean_named_genes(
    x_gene: np.ndarray, gene_names: Sequence[str], targets: set[str]
) -> np.ndarray:
    columns = [
        index for index, gene in enumerate(gene_names) if str(gene).upper() in targets
    ]
    if not columns:
        return np.zeros(x_gene.shape[0], dtype=np.float32)
    return x_gene[:, columns].mean(axis=1, dtype=np.float32)


def append_frozen_engineered_features(
    x_gene: np.ndarray, gene_names: Sequence[str], engineered_names: Sequence[str]
) -> np.ndarray:
    names = [str(name) for name in engineered_names]
    unsupported = sorted(set(names) - SUPPORTED_ENGINEERED_FEATURES)
    if unsupported:
        raise PrerequisiteError(
            f"Frozen gdTAI payload contains unsupported engineered features: {unsupported}"
        )
    lookup = {str(gene).upper(): index for index, gene in enumerate(gene_names)}
    trdc = (
        x_gene[:, lookup["TRDC"]]
        if "TRDC" in lookup
        else np.zeros(x_gene.shape[0], dtype=np.float32)
    )
    trdv = _sum_gene_prefixes(x_gene, gene_names, ("TRDV",))
    trdj = _sum_gene_prefixes(x_gene, gene_names, ("TRDJ", "TRDD"))
    trg = _sum_gene_prefixes(x_gene, gene_names, ("TRGV", "TRGJ", "TRGC"))
    ab = _sum_gene_prefixes(
        x_gene, gene_names, ("TRAV", "TRAJ", "TRAC", "TRBV", "TRBJ", "TRBC")
    )
    cd3 = _mean_named_genes(x_gene, gene_names, {"CD3D", "CD3E", "CD3G"})
    nk = _mean_named_genes(
        x_gene,
        gene_names,
        FROZEN_NK_CONTROL_GENES,
    )
    values: dict[str, np.ndarray] = {
        "any_TRDV": trdv > 0,
        "any_TRDJ": trdj > 0,
        "any_TRG": trg > 0,
        "any_ab_TCR_gene": ab > 0,
        "TRDC_only": (trdc > 0) & (trdv <= 0) & (trdj <= 0),
        "TRDC_plus_TRDV": (trdc > 0) & (trdv > 0),
        "TRDC_plus_TRDJ": (trdc > 0) & (trdj > 0),
        "CD3_score": cd3,
        "NK_score": nk,
        "gdT_TCR_score": trdc + trdv + trdj + trg,
        "abT_TCR_score": ab,
        "NK_minus_CD3_score": nk - cd3,
        "TRDC_log1p": trdc,
        "TRDV_score": trdv,
        "TRDJ_score": trdj,
        "TRG_score": trg,
    }
    engineered = np.column_stack([values[name] for name in names]).astype(
        np.float32, copy=False
    )
    return np.column_stack([x_gene, engineered]).astype(np.float32, copy=False)


def _read_csr_chunk(
    matrix: Any, start: int, end: int, n_vars: int
) -> sparse.csr_matrix:
    absolute_indptr = np.asarray(matrix["indptr"][start : end + 1], dtype=np.int64)
    data_start = int(absolute_indptr[0])
    data_end = int(absolute_indptr[-1])
    local_indptr = absolute_indptr - data_start
    data = np.asarray(matrix["data"][data_start:data_end])
    indices = np.asarray(matrix["indices"][data_start:data_end])
    return sparse.csr_matrix((data, indices, local_indptr), shape=(end - start, n_vars))


def frozen_profile_inference(
    approved: ApprovedInput,
    config: Mapping[str, Any],
    obs: pd.DataFrame,
) -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    pd.DataFrame,
    dict[str, Any],
]:
    import h5py

    frozen = approved.frozen_selection
    evaluator = frozen.evaluator
    selected_profile = frozen.selected_profile
    profile_config = config["frozen_profile_selection"]
    union_genes, _ = evaluator.feature_schema(list(frozen.profiles))
    selected_payload = (
        selected_profile.payload["base_model"]
        if selected_profile.model_id == "gdtai_v2"
        else selected_profile.payload
    )
    selected_genes = [str(value) for value in selected_payload["gene_names"]]
    counts_layer = str(profile_config["counts_layer"])
    with h5py.File(approved.input_h5ad, "r") as handle:
        var_names = _decode_strings(handle["var"]["_index"][:]).astype(str)
        gene_lookup = {gene: index for index, gene in enumerate(var_names)}
        missing_selected_genes = [
            gene for gene in selected_genes if gene not in gene_lookup
        ]
        if len(missing_selected_genes) > int(
            profile_config["maximum_missing_model_genes"]
        ):
            raise PrerequisiteError(
                "Rebuilt input is missing "
                f"{len(missing_selected_genes)} selected frozen-model genes: "
                f"{missing_selected_genes[:20]}"
            )
        available_union_genes = [gene for gene in union_genes if gene in gene_lookup]
        union_gene_indices = np.asarray(
            [gene_lookup[gene] for gene in available_union_genes], dtype=np.int64
        )
        union_spec = evaluator.FeatureSpec(
            gene_names=available_union_genes,
            gene_indices=union_gene_indices.astype(np.int32),
            gene_feature_names=[
                f"{gene}_log1p_cp10k" for gene in available_union_genes
            ],
            engineered_feature_names=[],
            model_feature_names=[
                f"{gene}_log1p_cp10k" for gene in available_union_genes
            ],
            gene_to_col={
                gene: index for index, gene in enumerate(available_union_genes)
            },
            engineered_to_col={},
        )
        matrix = handle["layers"][counts_layer]
        n_obs, n_vars = _matrix_shape(matrix)
        if obs.shape[0] != n_obs:
            raise PrerequisiteError(
                "Selected-profile inference metadata is not aligned with the count matrix."
            )
        scores = np.empty(n_obs, dtype=np.float32)
        thresholds = np.empty(n_obs, dtype=np.float32)
        annotations = np.empty(n_obs, dtype=object)
        annotation_sources = np.empty(n_obs, dtype=object)
        canonical_obs = {
            column: (
                obs[column].to_numpy(copy=False)
                if column in obs.columns
                else np.full(n_obs, "", dtype=object)
            )
            for column in evaluator.CANONICAL_OBS_COLUMNS
        }
        chunk_size = int(profile_config["chunk_size"])
        target_sum = float(profile_config["normalization_target_sum"])
        for start in range(0, n_obs, chunk_size):
            end = min(n_obs, start + chunk_size)
            counts = _read_csr_chunk(matrix, start, end, n_vars)
            totals = np.asarray(counts.sum(axis=1), dtype=np.float64).ravel()
            if np.any(totals <= 0):
                raise PrerequisiteError(
                    f"Raw counts contain zero-library cells in rows {start}:{end}."
                )
            x_union = (
                counts[:, union_gene_indices].toarray().astype(np.float32, copy=False)
            )
            x_union *= (target_sum / totals).astype(np.float32)[:, None]
            np.log1p(x_union, out=x_union)
            chunk_obs = {
                column: values[start:end] for column, values in canonical_obs.items()
            }
            chunk_annotation, chunk_annotation_source = evaluator.normalize_annotation(
                chunk_obs, x_union, union_spec
            )
            prediction = evaluator.predict_profiles(
                [selected_profile], x_union, union_spec, chunk_annotation
            )
            chunk_score, chunk_prediction = prediction[selected_profile.profile_id]
            chunk_score = np.asarray(chunk_score, dtype=np.float32)
            chunk_threshold = selected_profile_thresholds(
                evaluator, selected_profile, chunk_annotation
            )
            if not np.isfinite(chunk_score).all():
                raise PrerequisiteError(
                    f"Frozen gdTAI produced non-finite scores in rows {start}:{end}."
                )
            if not np.array_equal(
                np.asarray(chunk_prediction, dtype=bool),
                chunk_score >= chunk_threshold,
            ):
                raise PrerequisiteError(
                    "Selected-profile predictions disagree with their frozen thresholds."
                )
            scores[start:end] = chunk_score
            thresholds[start:end] = chunk_threshold
            annotations[start:end] = chunk_annotation
            annotation_sources[start:end] = chunk_annotation_source
            logging.info(
                "Selected frozen gdTAI inference: %s / %s cells",
                f"{end:,}",
                f"{n_obs:,}",
            )
    selected_gene_set = set(selected_genes)
    availability = pd.DataFrame(
        {
            "feature_index": np.arange(len(union_genes), dtype=int),
            "gene": union_genes,
            "required_by_selected_profile": [
                gene in selected_gene_set for gene in union_genes
            ],
            "available_in_rebuilt_h5ad": [gene in gene_lookup for gene in union_genes],
            "h5ad_gene_index": [gene_lookup.get(gene, "") for gene in union_genes],
        }
    )
    finite_thresholds = thresholds[np.isfinite(thresholds)]
    details = {
        "profile_id": frozen.selected_profile_id,
        "model_id": frozen.selected_model_id,
        "mode": frozen.selected_mode,
        "model_path": str(frozen.selected_model_path),
        "model_sha256": approved.model_sha256,
        "threshold_contract": frozen.threshold_contract,
        "threshold_contract_sha256": frozen.threshold_contract_sha256,
        "minimum_finite_threshold": (
            float(finite_thresholds.min()) if finite_thresholds.size else None
        ),
        "maximum_finite_threshold": (
            float(finite_thresholds.max()) if finite_thresholds.size else None
        ),
        "n_disabled_threshold_cells": int(np.isposinf(thresholds).sum()),
        "n_model_genes": len(selected_genes),
        "n_profile_union_genes": len(union_genes),
        "missing_model_genes": missing_selected_genes,
        "selection_decision": str(frozen.selection_decision),
        "selection_decision_sha256": frozen.selection_decision_sha256,
        "selection_decision_file_sha256": frozen.selection_decision_file_sha256,
        "holdout_status": str(frozen.holdout_status),
        "holdout_status_sha256": frozen.holdout_status_sha256,
        "comparison_workflow_sha256": frozen.comparison_workflow_sha256,
    }
    return (
        scores,
        thresholds,
        annotations,
        annotation_sources,
        availability,
        details,
    )


def read_source_metadata(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    import anndata as ad

    source = ad.read_h5ad(path, backed="r")
    try:
        obs = source.obs.copy()
        var = source.var.copy()
    finally:
        source.file.close()
    obs.index = obs.index.astype(str)
    var.index = var.index.astype(str)
    return obs, var


def extract_selected_counts(
    path: Path, layer: str, selected_rows: np.ndarray, block_rows: int = 10_000
) -> sparse.csr_matrix:
    import h5py

    rows = np.asarray(selected_rows, dtype=np.int64)
    if rows.size == 0 or np.any(np.diff(rows) <= 0):
        raise PrerequisiteError(
            "Selected source rows must be nonempty, unique, and strictly increasing."
        )
    with h5py.File(path, "r") as handle:
        matrix = handle["layers"][layer]
        n_obs, n_vars = _matrix_shape(matrix)
        if rows[0] < 0 or rows[-1] >= n_obs:
            raise PrerequisiteError("Selected source row is outside the rebuilt H5AD.")
        source_indptr = np.asarray(matrix["indptr"][:], dtype=np.int64)
        row_nnz = source_indptr[rows + 1] - source_indptr[rows]
        output_indptr = np.empty(rows.size + 1, dtype=np.int64)
        output_indptr[0] = 0
        np.cumsum(row_nnz, out=output_indptr[1:])
        total_nnz = int(output_indptr[-1])
        output_data = np.empty(total_nnz, dtype=matrix["data"].dtype)
        output_indices = np.empty(total_nnz, dtype=matrix["indices"].dtype)
        copied_rows = 0
        for source_block_start in range(0, n_obs, block_rows):
            source_block_end = min(n_obs, source_block_start + block_rows)
            selected_start = int(np.searchsorted(rows, source_block_start, side="left"))
            selected_end = int(np.searchsorted(rows, source_block_end, side="left"))
            if selected_end <= selected_start:
                continue
            block_data_start = int(source_indptr[source_block_start])
            block_data_end = int(source_indptr[source_block_end])
            block_data = np.asarray(matrix["data"][block_data_start:block_data_end])
            block_indices = np.asarray(
                matrix["indices"][block_data_start:block_data_end]
            )
            for output_row in range(selected_start, selected_end):
                source_row = int(rows[output_row])
                source_start = int(source_indptr[source_row] - block_data_start)
                source_end = int(source_indptr[source_row + 1] - block_data_start)
                output_start = int(output_indptr[output_row])
                output_end = int(output_indptr[output_row + 1])
                output_data[output_start:output_end] = block_data[
                    source_start:source_end
                ]
                output_indices[output_start:output_end] = block_indices[
                    source_start:source_end
                ]
            copied_rows = selected_end
            if source_block_start % (block_rows * 25) == 0:
                logging.info(
                    "Count extraction: %s / %s selected cells",
                    f"{copied_rows:,}",
                    f"{rows.size:,}",
                )
        if copied_rows != rows.size:
            raise RuntimeError(
                f"Count extraction copied {copied_rows:,} of {rows.size:,} selected rows."
            )
    output = sparse.csr_matrix(
        (output_data, output_indices, output_indptr), shape=(rows.size, n_vars)
    )
    output.sum_duplicates()
    output.sort_indices()
    return output


def is_tcr_gene(gene: str, prefixes: Sequence[str]) -> bool:
    upper = str(gene).upper()
    return upper.startswith(tuple(str(prefix).upper() for prefix in prefixes))


def select_hvgs(adata: Any, config: Mapping[str, Any], table_dir: Path) -> list[str]:
    import scanpy as sc

    hvg_config = config["feature_selection"]
    candidate_pool = min(int(hvg_config["candidate_pool_size"]), adata.n_vars)
    target = int(hvg_config["n_top_genes"])
    prefixes = list(hvg_config["tcr_gene_prefixes_excluded"])
    if sum(not is_tcr_gene(gene, prefixes) for gene in adata.var_names) < target:
        raise PrerequisiteError(
            "The rebuilt object does not contain enough non-TCR genes for HVG selection."
        )
    sc.pp.highly_variable_genes(
        adata,
        layer=config["prerequisites"]["counts_layer"],
        flavor=str(hvg_config["flavor"]),
        n_top_genes=candidate_pool,
        batch_key=str(hvg_config["batch_key"]),
        subset=False,
    )
    tcr_mask = np.asarray(
        [is_tcr_gene(gene, prefixes) for gene in adata.var_names], dtype=bool
    )
    ranks = pd.to_numeric(adata.var.get("highly_variable_rank"), errors="coerce")
    ranked = pd.DataFrame(
        {
            "gene": adata.var_names.astype(str),
            "rank": ranks.to_numpy(dtype=float),
            "candidate": adata.var["highly_variable"]
            .fillna(False)
            .to_numpy(dtype=bool),
            "excluded_tcr_gene": tcr_mask,
        }
    )
    selected = (
        ranked.loc[ranked["candidate"] & ~ranked["excluded_tcr_gene"]]
        .sort_values(["rank", "gene"], kind="mergesort")["gene"]
        .head(target)
        .tolist()
    )
    if len(selected) != target:
        raise PrerequisiteError(
            f"HVG selection produced only {len(selected)} non-TCR genes; required exactly {target}."
        )
    adata.var["highly_variable_gdt_atlas"] = adata.var_names.isin(selected)
    adata.var["excluded_from_hvg_tcr_gene"] = tcr_mask
    ranked["selected_for_scvi"] = ranked["gene"].isin(selected)
    ranked.to_csv(table_dir / "hvg_selection_audit.csv.gz", index=False)
    pd.DataFrame({"gene": selected}).to_csv(
        table_dir / "scvi_hvg_genes.csv", index=False
    )
    ranked.loc[ranked["excluded_tcr_gene"], ["gene"]].to_csv(
        table_dir / "tcr_genes_excluded_from_hvg.csv", index=False
    )
    return selected


def train_scvi(
    adata: Any,
    hvg_genes: Sequence[str],
    config: Mapping[str, Any],
    model_dir: Path,
    table_dir: Path,
) -> dict[str, Any]:
    import scvi
    import torch

    settings = config["scvi"]
    seed = int(settings["training_seed"])
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    scvi.settings.seed = seed
    hvg = adata[:, list(hvg_genes)].copy()
    categorical_covariates = [
        str(key) for key in settings.get("categorical_covariate_keys", [])
    ]
    for key in [settings["batch_key"], *categorical_covariates]:
        if key not in hvg.obs:
            raise PrerequisiteError(
                f"scVI metadata key `{key}` is missing from the gdT subset."
            )
        if (
            hvg.obs[key].isna().any()
            or hvg.obs[key].astype("string").str.strip().eq("").any()
        ):
            raise PrerequisiteError(
                f"scVI metadata key `{key}` contains missing or blank values."
            )
        hvg.obs[key] = hvg.obs[key].astype(str).astype("category")
    scvi.model.SCVI.setup_anndata(
        hvg,
        layer=config["prerequisites"]["counts_layer"],
        batch_key=str(settings["batch_key"]),
        categorical_covariate_keys=categorical_covariates or None,
    )
    model = scvi.model.SCVI(
        hvg,
        n_hidden=int(settings["n_hidden"]),
        n_latent=int(settings["n_latent"]),
        n_layers=int(settings["n_layers"]),
        dropout_rate=float(settings["dropout_rate"]),
        dispersion=str(settings["dispersion"]),
        gene_likelihood=str(settings["gene_likelihood"]),
    )
    model.train(
        max_epochs=int(settings["max_epochs"]),
        batch_size=int(settings["batch_size"]),
        early_stopping=bool(settings["early_stopping"]),
        early_stopping_patience=int(settings["early_stopping_patience"]),
        accelerator=settings["accelerator"],
        devices=settings["devices"],
    )
    latent = np.asarray(model.get_latent_representation(), dtype=np.float32)
    if latent.shape[0] != adata.n_obs or not np.isfinite(latent).all():
        raise RuntimeError("scVI returned an invalid latent representation.")
    latent_key = str(settings["latent_key"])
    adata.obsm[latent_key] = latent
    temporary_model_dir = model_dir.parent / f".{model_dir.name}.tmp-{uuid.uuid4().hex}"
    try:
        model.save(temporary_model_dir, overwrite=False, save_anndata=False)
        temporary_model_dir.replace(model_dir)
    except Exception:
        shutil.rmtree(temporary_model_dir, ignore_errors=True)
        raise
    history = getattr(model, "history", None)
    if history is not None:
        try:
            pd.DataFrame(history).to_csv(
                table_dir / "scvi_training_history.csv", index=False
            )
        except Exception as exc:
            logging.warning("Could not serialize scVI training history: %s", exc)
    return {
        "latent_key": latent_key,
        "latent_dimensions": int(latent.shape[1]),
        "training_seed": seed,
        "cuda_available": bool(torch.cuda.is_available()),
        "scvi_version": str(scvi.__version__),
        "model_dir": str(model_dir),
    }


def _resolution_token(value: float) -> str:
    text = f"{float(value):g}".replace("-", "m").replace(".", "p")
    return text


def run_embeddings_and_clustering(
    adata: Any, config: Mapping[str, Any], table_dir: Path
) -> tuple[pd.DataFrame, pd.DataFrame, list[str], list[str]]:
    import scanpy as sc
    from sklearn.metrics import adjusted_rand_score

    settings = config["embedding_clustering"]
    latent_key = config["scvi"]["latent_key"]
    seeds = [int(value) for value in settings["seeds"]]
    resolutions = [float(value) for value in settings["leiden_resolutions"]]
    sc.pp.neighbors(
        adata,
        n_neighbors=int(settings["n_neighbors"]),
        metric=str(settings["metric"]),
        use_rep=str(latent_key),
        random_state=int(settings["primary_seed"]),
    )
    umap_keys: list[str] = []
    cluster_keys: list[str] = []
    audit_rows: list[dict[str, Any]] = []
    for seed in seeds:
        sc.tl.umap(
            adata,
            min_dist=float(settings["umap_min_dist"]),
            spread=float(settings["umap_spread"]),
            random_state=seed,
        )
        umap_key = f"X_umap_s{seed}"
        adata.obsm[umap_key] = np.asarray(adata.obsm["X_umap"], dtype=np.float32).copy()
        umap_keys.append(umap_key)
        for resolution in resolutions:
            key = f"leiden_s{seed}_r{_resolution_token(resolution)}"
            sc.tl.leiden(
                adata,
                resolution=resolution,
                key_added=key,
                random_state=seed,
                flavor="igraph",
                directed=False,
            )
            adata.obs[key] = adata.obs[key].astype(str).astype("category")
            counts = adata.obs[key].value_counts()
            cluster_keys.append(key)
            audit_rows.append(
                {
                    "seed": seed,
                    "resolution": resolution,
                    "cluster_key": key,
                    "n_clusters": int(counts.size),
                    "minimum_cluster_cells": int(counts.min()),
                    "median_cluster_cells": float(counts.median()),
                    "maximum_cluster_cells": int(counts.max()),
                }
            )
    primary_seed = int(settings["primary_seed"])
    primary_resolution = float(settings["primary_resolution"])
    primary_key = f"leiden_s{primary_seed}_r{_resolution_token(primary_resolution)}"
    output_key = str(settings["primary_cluster_key"])
    adata.obs[output_key] = adata.obs[primary_key].astype(str).astype("category")
    adata.obsm["X_umap"] = adata.obsm[f"X_umap_s{primary_seed}"].copy()

    stability_rows: list[dict[str, Any]] = []
    for resolution in resolutions:
        keys = [f"leiden_s{seed}_r{_resolution_token(resolution)}" for seed in seeds]
        for left, right in combinations(keys, 2):
            stability_rows.append(
                {
                    "resolution": resolution,
                    "cluster_key_a": left,
                    "cluster_key_b": right,
                    "adjusted_rand_index": float(
                        adjusted_rand_score(
                            adata.obs[left].astype(str), adata.obs[right].astype(str)
                        )
                    ),
                }
            )
    audit = pd.DataFrame(audit_rows)
    stability = pd.DataFrame(stability_rows)
    audit.to_csv(table_dir / "leiden_seed_resolution_audit.csv", index=False)
    stability.to_csv(table_dir / "leiden_seed_stability.csv", index=False)
    return audit, stability, umap_keys, cluster_keys


def _prepare_output_paths(
    config: Mapping[str, Any], root: Path = PROJECT_ROOT
) -> dict[str, Path]:
    execution = config["execution"]
    paths = {
        key: resolve_path(execution[key], root)
        for key in [
            "atlas_h5ad",
            "build_manifest",
            "build_table_dir",
            "build_log_dir",
            "scvi_model_dir",
        ]
    }
    collisions = [str(path) for path in paths.values() if path.exists()]
    if collisions:
        raise PrerequisiteError(
            "Atlas outputs already exist and overwrite is disabled; move or approve their disposition first: "
            + "; ".join(collisions)
        )
    return paths


def _setup_logging(log_dir: Path) -> None:
    log_dir.mkdir(parents=True, exist_ok=False)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[
            logging.FileHandler(
                log_dir / "build_gdt_atlas.log", mode="x", encoding="utf-8"
            ),
            logging.StreamHandler(),
        ],
        force=True,
    )


def _write_h5ad_atomic(adata: Any, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.parent / f".{output.name}.tmp-{uuid.uuid4().hex}"
    try:
        adata.write_h5ad(temporary, compression="gzip")
        temporary.replace(output)
    finally:
        temporary.unlink(missing_ok=True)


def build_atlas(
    config: Mapping[str, Any], config_path: Path, approved: ApprovedInput
) -> dict[str, Any]:
    import anndata as ad
    import scanpy as sc

    paths = _prepare_output_paths(config)
    paths["build_table_dir"].mkdir(parents=True, exist_ok=False)
    paths["scvi_model_dir"].parent.mkdir(parents=True, exist_ok=True)
    _setup_logging(paths["build_log_dir"])
    source_stat_before = approved.input_h5ad.stat()
    logging.info("Approved rebuilt input: %s", approved.input_h5ad)
    logging.info(
        "Selected frozen gdTAI profile: %s (%s)",
        approved.frozen_selection.selected_profile_id,
        approved.model_path,
    )

    obs, var = read_source_metadata(approved.input_h5ad)
    (
        scores,
        thresholds,
        annotations,
        annotation_sources,
        availability,
        inference_info,
    ) = frozen_profile_inference(approved, config, obs)
    availability.to_csv(
        paths["build_table_dir"] / "frozen_model_gene_availability.csv", index=False
    )
    selection = build_selection_table(obs, scores, thresholds, config)
    summary = summarize_selection(selection)
    summary.to_csv(paths["build_table_dir"] / "selection_summary.csv", index=False)
    audit = pd.concat(
        [
            pd.DataFrame({"cell_id": obs.index.astype(str)}, index=obs.index),
            obs[
                [
                    config["analysis"]["study_column"],
                    config["analysis"]["donor_column"],
                    config["analysis"]["sample_column"],
                ]
            ],
            selection,
        ],
        axis=1,
    )
    audit.to_csv(paths["build_table_dir"] / "cell_selection_audit.csv.gz", index=False)
    selected_rows = np.flatnonzero(
        selection["gdt_atlas_selected"].to_numpy(dtype=bool)
    ).astype(np.int64)
    if selected_rows.size < int(config["selection"]["minimum_selected_cells"]):
        raise PrerequisiteError(
            f"Curation retained only {selected_rows.size:,} cells, below the configured atlas minimum."
        )
    counts = extract_selected_counts(
        approved.input_h5ad,
        str(config["prerequisites"]["counts_layer"]),
        selected_rows,
    )
    selected_obs = obs.iloc[selected_rows].copy()
    selected_selection = selection.iloc[selected_rows]
    profile_config = config["frozen_profile_selection"]
    score_column = str(profile_config["score_obs_column"])
    threshold_column = str(profile_config["threshold_obs_column"])
    prediction_column = str(profile_config["prediction_obs_column"])
    annotation_column = str(profile_config["annotation_obs_column"])
    annotation_source_column = str(profile_config["annotation_source_obs_column"])
    selected_obs[score_column] = selected_selection["gdtai_frozen_score"].to_numpy(
        dtype=np.float32
    )
    selected_obs[threshold_column] = selected_selection[
        "gdtai_frozen_threshold"
    ].to_numpy(dtype=np.float32)
    selected_obs[prediction_column] = selected_selection[
        "gdtai_frozen_prediction"
    ].to_numpy(dtype=bool)
    selected_obs[annotation_column] = annotations[selected_rows]
    selected_obs[annotation_source_column] = annotation_sources[selected_rows]
    for column in [
        "selection_primary_gold_paired_trg_trd",
        "selection_gold_fn_added",
        "gdt_atlas_selection_reason",
    ]:
        selected_obs[column] = selected_selection[column].to_numpy()
    adata = ad.AnnData(X=counts.copy(), obs=selected_obs, var=var.copy())
    adata.layers[str(config["prerequisites"]["counts_layer"])] = counts
    hvg_genes = select_hvgs(adata, config, paths["build_table_dir"])
    scvi_info = train_scvi(
        adata,
        hvg_genes,
        config,
        paths["scvi_model_dir"],
        paths["build_table_dir"],
    )
    cluster_audit, stability, umap_keys, cluster_keys = run_embeddings_and_clustering(
        adata, config, paths["build_table_dir"]
    )
    sc.pp.normalize_total(adata, target_sum=10_000.0)
    sc.pp.log1p(adata)
    adata.uns["gdt_atlas_build"] = {
        "workflow_id": str(config["workflow_id"]),
        "config_sha256": approved.config_sha256,
        "source_h5ad": str(approved.input_h5ad),
        "source_sha256": approved.input_sha256,
        "approval_marker": str(approved.approval_marker),
        "gdtai_profile_id": approved.frozen_selection.selected_profile_id,
        "gdtai_model_id": approved.frozen_selection.selected_model_id,
        "gdtai_mode": approved.frozen_selection.selected_mode,
        "gdtai_model_sha256": approved.model_sha256,
        "gdtai_threshold_contract_sha256": approved.frozen_selection.threshold_contract_sha256,
        "gdtai_selection_decision_sha256": approved.frozen_selection.selection_decision_sha256,
        "gdtai_holdout_status_sha256": approved.frozen_selection.holdout_status_sha256,
        "selection_rule": (
            "post-holdout selected frozen gdTAI positive or primary-gold paired productive TRG/TRD false negative; "
            "minus alpha-beta-only, strict high-confidence NK, doublet, non-T contaminant, and severe QC flags"
        ),
        "silver_only_fn_addback": False,
        "hard_trdv_exclusion": False,
        "nk_gene_only_exclusion": False,
        "tcr_genes_excluded_from_hvg": True,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    _write_h5ad_atomic(adata, paths["atlas_h5ad"])
    output_sha = sha256_file(paths["atlas_h5ad"])
    source_stat_after = approved.input_h5ad.stat()
    if (
        source_stat_before.st_size != source_stat_after.st_size
        or source_stat_before.st_mtime_ns != source_stat_after.st_mtime_ns
    ):
        raise RuntimeError(
            "The rebuilt source H5AD changed while the atlas was being built."
        )
    manifest = {
        "schema_version": 1,
        "status": "complete",
        "workflow_id": config["workflow_id"],
        "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        "config_path": str(config_path.resolve()),
        "config_sha256": approved.config_sha256,
        "input_h5ad": str(approved.input_h5ad),
        "input_sha256": approved.input_sha256,
        "input_size_bytes": approved.input_size_bytes,
        "approval_marker": str(approved.approval_marker),
        "approved_by": approved.approval["approved_by"],
        "gdtai": inference_info,
        "selection": {
            row.category: int(row.n_cells) for row in summary.itertuples(index=False)
        },
        "atlas_h5ad": str(paths["atlas_h5ad"]),
        "atlas_sha256": output_sha,
        "atlas_size_bytes": int(paths["atlas_h5ad"].stat().st_size),
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "n_hvgs": len(hvg_genes),
        "scvi": scvi_info,
        "umap_keys": umap_keys,
        "cluster_keys": cluster_keys,
        "primary_cluster_key": config["embedding_clustering"]["primary_cluster_key"],
        "n_clustering_runs": int(cluster_audit.shape[0]),
        "n_seed_stability_comparisons": int(stability.shape[0]),
        "source_unchanged_by_size_and_mtime": True,
    }
    atomic_write_json(paths["build_manifest"], manifest)
    logging.info("Completed de novo gdT atlas: %s", paths["atlas_h5ad"])
    return manifest


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    config_path = args.config.expanduser().resolve()
    config = load_config(config_path)
    validate_runtime_environment(config)
    approved = validate_prerequisites(config, config_path)
    if args.preflight_only:
        print(
            json.dumps(
                {
                    "status": "approved_preflight_passed",
                    "workflow_id": config["workflow_id"],
                    "config_sha256": approved.config_sha256,
                    "input_h5ad": str(approved.input_h5ad),
                    "input_sha256": approved.input_sha256,
                    "approval_marker": str(approved.approval_marker),
                    "selected_profile_id": approved.frozen_selection.selected_profile_id,
                    "selected_model_id": approved.frozen_selection.selected_model_id,
                    "selected_mode": approved.frozen_selection.selected_mode,
                    "model_sha256": approved.model_sha256,
                    "selection_decision_sha256": approved.frozen_selection.selection_decision_sha256,
                    "holdout_status_sha256": approved.frozen_selection.holdout_status_sha256,
                    "threshold_contract": approved.frozen_selection.threshold_contract,
                    "h5ad_contract": approved.h5ad_contract,
                },
                indent=2,
                sort_keys=True,
            )
        )
        return 0
    build_atlas(config, config_path, approved)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
