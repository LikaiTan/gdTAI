#!/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python
"""Run donor-aware biology analyses on a completed de novo gdT atlas build."""

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
import json
import logging
import math
import shutil
import uuid
from datetime import datetime, timezone
from pathlib import Path
from statistics import NormalDist
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd
from scipy import sparse, stats

from tnk_atlas.provenance import atomic_write_json, sha256_file

try:
    from .build_gdt_atlas import (
        DEFAULT_CONFIG,
        PROJECT_ROOT,
        ApprovedInput,
        AtlasConfigError,
        PrerequisiteError,
        is_tcr_gene,
        load_config,
        resolve_path,
        strict_bool_array,
        validate_prerequisites,
        validate_runtime_environment,
    )
except ImportError:
    from build_gdt_atlas import (  # type: ignore[no-redef]
        DEFAULT_CONFIG,
        PROJECT_ROOT,
        ApprovedInput,
        AtlasConfigError,
        PrerequisiteError,
        is_tcr_gene,
        load_config,
        resolve_path,
        strict_bool_array,
        validate_prerequisites,
        validate_runtime_environment,
    )


INVALID_UNIT_VALUES = frozenset({"", "na", "nan", "none", "null", "unknown", "<na>"})


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument(
        "--preflight-only",
        action="store_true",
        help="Validate the approved rebuild and completed atlas build without writing analysis outputs.",
    )
    return parser.parse_args(argv)


def validate_completed_build(
    config: Mapping[str, Any], approved: ApprovedInput, root: Path = PROJECT_ROOT
) -> tuple[Path, Path, dict[str, Any]]:
    execution = config["execution"]
    manifest_path = resolve_path(execution["build_manifest"], root)
    atlas_path = resolve_path(execution["atlas_h5ad"], root)
    missing = [str(path) for path in [manifest_path, atlas_path] if not path.is_file()]
    if missing:
        raise PrerequisiteError(
            f"A completed de novo gdT atlas build is absent: {missing}"
        )
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise PrerequisiteError(
            f"Cannot read build manifest {manifest_path}: {exc}"
        ) from exc
    expected = {
        "status": "complete",
        "workflow_id": config["workflow_id"],
        "config_sha256": approved.config_sha256,
        "input_sha256": approved.input_sha256,
    }
    for key, value in expected.items():
        if manifest.get(key) != value:
            raise PrerequisiteError(
                f"Build manifest `{key}` mismatch: expected {value!r}, observed {manifest.get(key)!r}."
            )
    if resolve_path(str(manifest.get("atlas_h5ad", "")), root) != atlas_path:
        raise PrerequisiteError(
            "Build manifest atlas path does not match the configured atlas path."
        )
    observed_sha = sha256_file(atlas_path)
    if manifest.get("atlas_sha256") != observed_sha:
        raise PrerequisiteError(
            f"Atlas SHA256 no longer matches the completed build manifest: observed {observed_sha}."
        )
    gdtai = manifest.get("gdtai", {})
    frozen = approved.frozen_selection
    expected_gdtai = {
        "profile_id": frozen.selected_profile_id,
        "model_id": frozen.selected_model_id,
        "mode": frozen.selected_mode,
        "model_sha256": frozen.selected_model_sha256,
        "selection_decision_sha256": frozen.selection_decision_sha256,
        "selection_decision_file_sha256": frozen.selection_decision_file_sha256,
        "holdout_status_sha256": frozen.holdout_status_sha256,
        "threshold_contract_sha256": frozen.threshold_contract_sha256,
    }
    for key, expected_value in expected_gdtai.items():
        if gdtai.get(key) != expected_value:
            raise PrerequisiteError(
                f"Build manifest gdTAI `{key}` mismatch: expected {expected_value!r}, "
                f"observed {gdtai.get(key)!r}."
            )
    return atlas_path, manifest_path, manifest


def benjamini_hochberg(p_values: Sequence[float]) -> np.ndarray:
    values = np.asarray(p_values, dtype=float)
    adjusted = np.full(values.shape, np.nan, dtype=float)
    finite = np.isfinite(values)
    if not finite.any():
        return adjusted
    finite_indices = np.flatnonzero(finite)
    finite_values = values[finite]
    order = np.argsort(finite_values, kind="mergesort")
    ranked = finite_values[order]
    n = ranked.size
    corrected = ranked * n / np.arange(1, n + 1, dtype=float)
    corrected = np.minimum.accumulate(corrected[::-1])[::-1]
    corrected = np.clip(corrected, 0.0, 1.0)
    adjusted_indices = finite_indices[order]
    adjusted[adjusted_indices] = corrected
    return adjusted


def random_effects_dersimonian_laird(
    effects: Sequence[float], variances: Sequence[float], confidence_level: float = 0.95
) -> dict[str, float | int]:
    y = np.asarray(effects, dtype=float)
    v = np.asarray(variances, dtype=float)
    valid = np.isfinite(y) & np.isfinite(v) & (v > 0)
    y = y[valid]
    v = v[valid]
    k = int(y.size)
    if k == 0:
        return {
            "k_studies": 0,
            "pooled_effect": np.nan,
            "standard_error": np.nan,
            "ci_low": np.nan,
            "ci_high": np.nan,
            "z_score": np.nan,
            "p_value": np.nan,
            "tau_squared": np.nan,
            "i_squared": np.nan,
            "q_statistic": np.nan,
        }
    fixed_weights = 1.0 / v
    fixed_mean = float(np.sum(fixed_weights * y) / np.sum(fixed_weights))
    q = float(np.sum(fixed_weights * np.square(y - fixed_mean)))
    df = k - 1
    c_value = float(
        np.sum(fixed_weights) - np.sum(np.square(fixed_weights)) / np.sum(fixed_weights)
    )
    tau_squared = max(0.0, (q - df) / c_value) if k > 1 and c_value > 0 else 0.0
    random_weights = 1.0 / (v + tau_squared)
    pooled = float(np.sum(random_weights * y) / np.sum(random_weights))
    standard_error = float(math.sqrt(1.0 / np.sum(random_weights)))
    z_score = pooled / standard_error if standard_error > 0 else np.nan
    p_value = (
        float(2.0 * stats.norm.sf(abs(z_score))) if np.isfinite(z_score) else np.nan
    )
    alpha = 1.0 - float(confidence_level)
    critical = NormalDist().inv_cdf(1.0 - alpha / 2.0)
    i_squared = max(0.0, (q - df) / q) * 100.0 if k > 1 and q > 0 else 0.0
    return {
        "k_studies": k,
        "pooled_effect": pooled,
        "standard_error": standard_error,
        "ci_low": pooled - critical * standard_error,
        "ci_high": pooled + critical * standard_error,
        "z_score": z_score,
        "p_value": p_value,
        "tau_squared": tau_squared,
        "i_squared": i_squared,
        "q_statistic": q,
    }


def _clean_identifier(values: pd.Series | Sequence[Any]) -> pd.Series:
    cleaned = pd.Series(values, copy=False).astype("string").fillna("").str.strip()
    invalid = cleaned.str.lower().isin(INVALID_UNIT_VALUES)
    return cleaned.mask(invalid, "")


def targeted_cohort_mask(obs: pd.DataFrame, config: Mapping[str, Any]) -> np.ndarray:
    analysis = config["analysis"]
    rule = analysis["targeted_cohort_exclusion"]
    study_column = str(analysis["study_column"])
    flag_column = str(rule["required_boolean_column"])
    missing = [column for column in [study_column, flag_column] if column not in obs]
    if missing:
        raise PrerequisiteError(
            f"Targeted-cohort exclusion metadata is missing: {missing}"
        )
    explicit = strict_bool_array(obs[flag_column], flag_column)
    sources = {str(value) for value in rule.get("source_ids", [])}
    source_match = (
        obs[study_column].astype("string").fillna("").isin(sources).to_numpy(dtype=bool)
    )
    return explicit | source_match


def _stable_group_value(values: pd.Series) -> str:
    cleaned = _clean_identifier(values)
    unique = sorted(set(cleaned[cleaned.ne("")].astype(str)))
    if not unique:
        return ""
    if len(unique) > 1:
        return "__MULTIPLE__"
    return unique[0]


def _composition_for_units(
    frame: pd.DataFrame,
    unit_columns: list[str],
    cluster_column: str,
    minimum_cells: int,
    metadata_columns: list[str],
) -> pd.DataFrame:
    valid = np.ones(frame.shape[0], dtype=bool)
    for column in unit_columns:
        valid &= (
            frame[column]
            .astype("string")
            .fillna("")
            .str.strip()
            .ne("")
            .to_numpy(dtype=bool)
        )
    work = frame.loc[valid].copy()
    if work.empty:
        columns = [
            *unit_columns,
            cluster_column,
            "cluster_cells",
            "unit_cells",
            "cluster_fraction",
        ]
        return pd.DataFrame(columns=columns)
    work[cluster_column] = work[cluster_column].astype(str)
    counts = (
        work.groupby([*unit_columns, cluster_column], observed=True, dropna=False)
        .size()
        .rename("cluster_cells")
        .reset_index()
    )
    totals = (
        work.groupby(unit_columns, observed=True, dropna=False)
        .size()
        .rename("unit_cells")
        .reset_index()
    )
    cluster_levels = pd.DataFrame(
        {cluster_column: sorted(work[cluster_column].astype(str).unique())}
    )
    skeleton = totals[unit_columns].merge(cluster_levels, how="cross")
    counts = skeleton.merge(
        counts,
        on=[*unit_columns, cluster_column],
        how="left",
        validate="one_to_one",
    )
    counts["cluster_cells"] = counts["cluster_cells"].fillna(0).astype(np.int64)
    metadata_agg = {
        column: (column, _stable_group_value)
        for column in metadata_columns
        if column not in unit_columns
    }
    metadata = (
        work.groupby(unit_columns, observed=True, dropna=False)
        .agg(**metadata_agg)
        .reset_index()
        if metadata_agg
        else totals[unit_columns].copy()
    )
    targeted = (
        work.groupby(unit_columns, observed=True, dropna=False)["_targeted_cohort"]
        .any()
        .rename("targeted_cohort")
        .reset_index()
    )
    out = counts.merge(totals, on=unit_columns, validate="many_to_one")
    out = out.merge(metadata, on=unit_columns, validate="many_to_one")
    out = out.merge(targeted, on=unit_columns, validate="many_to_one")
    out["cluster_fraction"] = out["cluster_cells"] / out["unit_cells"]
    out["minimum_cell_requirement_met"] = out["unit_cells"] >= int(minimum_cells)
    out["eligible_for_abundance_inference"] = (
        out["minimum_cell_requirement_met"] & ~out["targeted_cohort"]
    )
    return out.sort_values(
        [*unit_columns, cluster_column], kind="mergesort"
    ).reset_index(drop=True)


def build_composition_tables(
    obs: pd.DataFrame, config: Mapping[str, Any]
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    analysis = config["analysis"]
    composition = analysis["composition"]
    study = str(analysis["study_column"])
    sample = str(analysis["sample_column"])
    library = str(analysis["library_column"])
    donor = str(analysis["donor_column"])
    condition = str(analysis["condition_column"])
    tissue = str(analysis["tissue_column"])
    cluster = str(analysis["cluster_column"])
    required = [study, sample, library, donor, condition, tissue, cluster]
    missing = sorted(set(required) - set(obs.columns))
    if missing:
        raise PrerequisiteError(f"Atlas composition metadata is missing: {missing}")
    frame = obs[required].copy()
    frame[study] = _clean_identifier(frame[study])
    sample_values = _clean_identifier(frame[sample])
    library_values = _clean_identifier(frame[library])
    sample_values = sample_values.mask(sample_values.eq(""), library_values)
    frame["_sample_unit"] = frame[study] + "||" + sample_values
    donor_values = _clean_identifier(frame[donor])
    frame["_donor_unit"] = frame[study] + "||" + donor_values
    frame.loc[sample_values.eq(""), "_sample_unit"] = ""
    frame.loc[donor_values.eq(""), "_donor_unit"] = ""
    frame["_targeted_cohort"] = targeted_cohort_mask(obs, config)
    metadata = [study, donor, condition, tissue]
    sample_table = _composition_for_units(
        frame,
        [study, "_sample_unit"],
        cluster,
        int(composition["minimum_cells_per_sample"]),
        metadata,
    )
    donor_table = _composition_for_units(
        frame,
        [study, "_donor_unit"],
        cluster,
        int(composition["minimum_cells_per_donor"]),
        [condition, tissue],
    )
    donor_condition_table = _composition_for_units(
        frame,
        [study, "_donor_unit", condition],
        cluster,
        int(composition["minimum_cells_per_donor"]),
        [tissue],
    )
    audit = pd.DataFrame(
        [
            {
                "unit_type": "sample",
                "n_units": int(sample_table["_sample_unit"].nunique())
                if not sample_table.empty
                else 0,
                "n_targeted_units": int(
                    sample_table.loc[
                        sample_table["targeted_cohort"], "_sample_unit"
                    ].nunique()
                )
                if not sample_table.empty
                else 0,
                "n_inference_eligible_units": int(
                    sample_table.loc[
                        sample_table["eligible_for_abundance_inference"], "_sample_unit"
                    ].nunique()
                )
                if not sample_table.empty
                else 0,
            },
            {
                "unit_type": "donor",
                "n_units": int(donor_table["_donor_unit"].nunique())
                if not donor_table.empty
                else 0,
                "n_targeted_units": int(
                    donor_table.loc[
                        donor_table["targeted_cohort"], "_donor_unit"
                    ].nunique()
                )
                if not donor_table.empty
                else 0,
                "n_inference_eligible_units": int(
                    donor_table.loc[
                        donor_table["eligible_for_abundance_inference"], "_donor_unit"
                    ].nunique()
                )
                if not donor_table.empty
                else 0,
            },
        ]
    )
    return sample_table, donor_table, donor_condition_table, audit


def _logit_fraction(
    cluster_cells: pd.Series, total_cells: pd.Series, offset: float
) -> np.ndarray:
    fraction = (cluster_cells.to_numpy(dtype=float) + offset) / (
        total_cells.to_numpy(dtype=float) + 2.0 * offset
    )
    return np.log(fraction / (1.0 - fraction))


def run_within_study_contrasts(
    donor_condition: pd.DataFrame, config: Mapping[str, Any]
) -> tuple[pd.DataFrame, pd.DataFrame]:
    analysis = config["analysis"]
    study_column = str(analysis["study_column"])
    donor_column = "_donor_unit"
    cluster_column = str(analysis["cluster_column"])
    minimum = int(analysis["composition"]["minimum_donors_per_group"])
    offset = float(analysis["composition"]["continuity_offset"])
    rows: list[dict[str, Any]] = []
    plan_rows: list[dict[str, Any]] = []
    for contrast in analysis.get("within_study_contrasts", []):
        contrast_id = str(contrast.get("contrast_id", "")).strip()
        enabled = bool(contrast.get("enabled", False))
        group_column = str(contrast.get("group_column", ""))
        case_values = {
            str(value).strip().lower() for value in contrast.get("case_values", [])
        }
        reference_values = {
            str(value).strip().lower() for value in contrast.get("reference_values", [])
        }
        study_ids = {str(value) for value in contrast.get("study_ids", [])}
        plan_rows.append(
            {
                "contrast_id": contrast_id,
                "enabled": enabled,
                "group_column": group_column,
                "case_values": "|".join(sorted(case_values)),
                "reference_values": "|".join(sorted(reference_values)),
                "study_ids": "|".join(sorted(study_ids)),
                "paired_by_donor": bool(contrast.get("paired_by_donor", False)),
                "note": str(contrast.get("note", "")),
            }
        )
        if not enabled:
            continue
        if (
            not contrast_id
            or not group_column
            or not case_values
            or not reference_values
        ):
            raise AtlasConfigError(
                f"Enabled contrast {contrast_id!r} has an incomplete definition."
            )
        if case_values & reference_values:
            raise AtlasConfigError(
                f"Enabled contrast {contrast_id!r} has overlapping case/reference values."
            )
        if group_column not in donor_condition:
            raise PrerequisiteError(
                f"Contrast group column `{group_column}` is absent from donor composition."
            )
        eligible = donor_condition[
            donor_condition["eligible_for_abundance_inference"]
        ].copy()
        if study_ids:
            eligible = eligible[eligible[study_column].astype(str).isin(study_ids)]
        normalized_group = (
            eligible[group_column].astype("string").fillna("").str.strip().str.lower()
        )
        eligible["_contrast_group"] = np.where(
            normalized_group.isin(case_values),
            "case",
            np.where(normalized_group.isin(reference_values), "reference", ""),
        )
        eligible = eligible[eligible["_contrast_group"].ne("")]
        for (study_id, cluster_id), group in eligible.groupby(
            [study_column, cluster_column], observed=True, sort=True
        ):
            case = group[group["_contrast_group"].eq("case")].copy()
            reference = group[group["_contrast_group"].eq("reference")].copy()
            paired = bool(contrast.get("paired_by_donor", False))
            if paired:
                case = case.set_index(donor_column)
                reference = reference.set_index(donor_column)
                shared = case.index.intersection(reference.index)
                if shared.size < minimum:
                    continue
                case_logit = _logit_fraction(
                    case.loc[shared, "cluster_cells"],
                    case.loc[shared, "unit_cells"],
                    offset,
                )
                reference_logit = _logit_fraction(
                    reference.loc[shared, "cluster_cells"],
                    reference.loc[shared, "unit_cells"],
                    offset,
                )
                differences = case_logit - reference_logit
                effect = float(np.mean(differences))
                standard_error = float(
                    np.std(differences, ddof=1) / math.sqrt(shared.size)
                )
                n_case = n_reference = int(shared.size)
                mean_case = float(case.loc[shared, "cluster_fraction"].mean())
                mean_reference = float(reference.loc[shared, "cluster_fraction"].mean())
            else:
                n_case = int(case[donor_column].nunique())
                n_reference = int(reference[donor_column].nunique())
                shared_donors = set(case[donor_column].astype(str)) & set(
                    reference[donor_column].astype(str)
                )
                if shared_donors:
                    raise PrerequisiteError(
                        f"Contrast {contrast_id!r} treats {len(shared_donors)} repeated donors as independent; "
                        "set paired_by_donor=true or correct the study mapping."
                    )
                if n_case < minimum or n_reference < minimum:
                    continue
                case_logit = _logit_fraction(
                    case["cluster_cells"], case["unit_cells"], offset
                )
                reference_logit = _logit_fraction(
                    reference["cluster_cells"], reference["unit_cells"], offset
                )
                effect = float(np.mean(case_logit) - np.mean(reference_logit))
                standard_error = float(
                    math.sqrt(
                        np.var(case_logit, ddof=1) / n_case
                        + np.var(reference_logit, ddof=1) / n_reference
                    )
                )
                mean_case = float(case["cluster_fraction"].mean())
                mean_reference = float(reference["cluster_fraction"].mean())
            if not np.isfinite(standard_error) or standard_error <= 0:
                continue
            z_score = effect / standard_error
            rows.append(
                {
                    "contrast_id": contrast_id,
                    study_column: str(study_id),
                    cluster_column: str(cluster_id),
                    "paired_by_donor": paired,
                    "n_case_donors": n_case,
                    "n_reference_donors": n_reference,
                    "mean_case_fraction": mean_case,
                    "mean_reference_fraction": mean_reference,
                    "effect_log_odds": effect,
                    "standard_error": standard_error,
                    "variance": standard_error**2,
                    "z_score": z_score,
                    "p_value": float(2.0 * stats.norm.sf(abs(z_score))),
                }
            )
    results = pd.DataFrame(rows)
    if not results.empty:
        results["fdr_bh"] = benjamini_hochberg(results["p_value"])
    plan = pd.DataFrame(plan_rows)
    return results, plan


def run_random_effects_meta_analysis(
    study_results: pd.DataFrame, config: Mapping[str, Any]
) -> pd.DataFrame:
    if study_results.empty:
        return pd.DataFrame(
            columns=[
                "contrast_id",
                config["analysis"]["cluster_column"],
                "k_studies",
                "pooled_effect",
                "standard_error",
                "ci_low",
                "ci_high",
                "z_score",
                "p_value",
                "tau_squared",
                "i_squared",
                "q_statistic",
                "fdr_bh",
            ]
        )
    settings = config["analysis"]["random_effects_meta_analysis"]
    if settings["method"] != "dersimonian_laird":
        raise AtlasConfigError(
            "Only DerSimonian-Laird random-effects meta-analysis is implemented."
        )
    minimum = int(settings["minimum_studies"])
    cluster_column = str(config["analysis"]["cluster_column"])
    rows: list[dict[str, Any]] = []
    for (contrast_id, cluster_id), group in study_results.groupby(
        ["contrast_id", cluster_column], observed=True, sort=True
    ):
        if group.shape[0] < minimum:
            continue
        result = random_effects_dersimonian_laird(
            group["effect_log_odds"],
            group["variance"],
            float(settings["confidence_level"]),
        )
        rows.append(
            {
                "contrast_id": str(contrast_id),
                cluster_column: str(cluster_id),
                **result,
            }
        )
    output = pd.DataFrame(rows)
    if not output.empty:
        output["fdr_bh"] = benjamini_hochberg(output["p_value"])
    return output


def _donor_codes(
    obs: pd.DataFrame, config: Mapping[str, Any]
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    analysis = config["analysis"]
    study = _clean_identifier(obs[str(analysis["study_column"])])
    donor = _clean_identifier(obs[str(analysis["donor_column"])])
    valid = study.ne("") & donor.ne("")
    donor_key = (study + "||" + donor).to_numpy(dtype=object)
    codes, uniques = pd.factorize(donor_key[valid.to_numpy(dtype=bool)], sort=True)
    return (
        valid.to_numpy(dtype=bool),
        codes.astype(np.int64),
        np.asarray(uniques, dtype=object),
    )


def donor_pseudobulk_cluster_vs_rest(
    adata: Any, config: Mapping[str, Any]
) -> tuple[pd.DataFrame, pd.DataFrame]:
    settings = config["analysis"]["donor_pseudobulk_de"]
    if not bool(settings.get("enabled", False)):
        return pd.DataFrame(), pd.DataFrame([{"status": "disabled"}])
    if settings.get("comparison") != "cluster_vs_rest_within_donor":
        raise AtlasConfigError(
            "Only donor-paired cluster-vs-rest pseudobulk DE is implemented."
        )
    counts_key = str(config["prerequisites"]["counts_layer"])
    if counts_key not in adata.layers:
        raise PrerequisiteError(
            f"Atlas is missing layers[{counts_key!r}] for pseudobulk DE."
        )
    counts = sparse.csr_matrix(adata.layers[counts_key])
    obs = adata.obs
    valid_donor, donor_codes, donor_names = _donor_codes(obs, config)
    if bool(
        config["analysis"]["targeted_cohort_exclusion"].get(
            "exclude_from_primary_de", False
        )
    ):
        valid_donor &= ~targeted_cohort_mask(obs, config)
        donor_codes, donor_names = pd.factorize(
            (
                _clean_identifier(
                    obs.loc[valid_donor, str(config["analysis"]["study_column"])]
                )
                + "||"
                + _clean_identifier(
                    obs.loc[valid_donor, str(config["analysis"]["donor_column"])]
                )
            ).to_numpy(dtype=object),
            sort=True,
        )
        donor_codes = donor_codes.astype(np.int64)
        donor_names = np.asarray(donor_names, dtype=object)
    selected_cells = np.flatnonzero(valid_donor)
    if selected_cells.size == 0 or donor_names.size == 0:
        raise PrerequisiteError(
            "No cells with valid donor identifiers are available for pseudobulk DE."
        )
    cell_counts = counts[selected_cells]
    n_donors = donor_names.size
    donor_map = sparse.csr_matrix(
        (
            np.ones(selected_cells.size, dtype=np.float32),
            (donor_codes, np.arange(selected_cells.size, dtype=np.int64)),
        ),
        shape=(n_donors, selected_cells.size),
    )
    prefixes = list(settings["tcr_gene_prefixes_excluded"])
    non_tcr = np.asarray(
        [not is_tcr_gene(gene, prefixes) for gene in adata.var_names], dtype=bool
    )
    candidate_indices = np.flatnonzero(non_tcr)
    total_by_gene = np.asarray(cell_counts[:, candidate_indices].sum(axis=0)).ravel()
    keep_gene = total_by_gene >= int(settings["minimum_total_counts_per_gene"])
    gene_indices = candidate_indices[keep_gene]
    gene_names = adata.var_names[gene_indices].astype(str).to_numpy(dtype=object)
    if gene_indices.size == 0:
        raise PrerequisiteError("No non-TCR genes pass the pseudobulk count filter.")
    expression = cell_counts[:, gene_indices].tocsr()
    total_expression_by_donor = (donor_map @ expression).tocsr()
    cell_library = np.asarray(cell_counts.sum(axis=1), dtype=np.float64).ravel()
    total_library_by_donor = np.asarray(
        donor_map @ cell_library[:, None], dtype=np.float64
    ).ravel()
    cluster_column = str(config["analysis"]["cluster_column"])
    clusters = (
        obs.iloc[selected_cells][cluster_column].astype(str).to_numpy(dtype=object)
    )
    minimum_cells = int(settings["minimum_cells_per_donor_cluster"])
    minimum_donors = int(settings["minimum_paired_donors"])
    target = float(settings["normalization_target"])
    pseudocount = float(settings["pseudocount"])
    gene_chunk_size = int(settings["gene_chunk_size"])
    rows: list[pd.DataFrame] = []
    audit_rows: list[dict[str, Any]] = []
    donor_total_cells = np.bincount(donor_codes, minlength=n_donors)
    for cluster_id in sorted(set(clusters)):
        cluster_mask = clusters == cluster_id
        cluster_donor_map = donor_map[:, cluster_mask]
        cluster_expression = (cluster_donor_map @ expression[cluster_mask]).tocsr()
        cluster_library = np.asarray(
            cluster_donor_map @ cell_library[cluster_mask, None], dtype=np.float64
        ).ravel()
        cluster_cells = np.asarray(cluster_donor_map.sum(axis=1)).ravel().astype(int)
        rest_cells = donor_total_cells - cluster_cells
        valid = (cluster_cells >= minimum_cells) & (rest_cells >= minimum_cells)
        donor_count = int(valid.sum())
        audit_rows.append(
            {
                "cluster": str(cluster_id),
                "paired_donors": donor_count,
                "minimum_cells_per_side": minimum_cells,
                "tested": donor_count >= minimum_donors,
            }
        )
        if donor_count < minimum_donors:
            continue
        rest_library = total_library_by_donor - cluster_library
        for start in range(0, gene_indices.size, gene_chunk_size):
            end = min(gene_indices.size, start + gene_chunk_size)
            cluster_dense = (
                cluster_expression[valid, start:end]
                .toarray()
                .astype(np.float64, copy=False)
            )
            total_dense = (
                total_expression_by_donor[valid, start:end]
                .toarray()
                .astype(np.float64, copy=False)
            )
            rest_dense = np.maximum(total_dense - cluster_dense, 0.0)
            cluster_log_cpm = np.log2(
                (cluster_dense + pseudocount)
                / (cluster_library[valid, None] + 1.0)
                * target
            )
            rest_log_cpm = np.log2(
                (rest_dense + pseudocount) / (rest_library[valid, None] + 1.0) * target
            )
            differences = cluster_log_cpm - rest_log_cpm
            effects = np.mean(differences, axis=0)
            standard_errors = np.std(differences, axis=0, ddof=1) / math.sqrt(
                donor_count
            )
            t_statistics = np.divide(
                effects,
                standard_errors,
                out=np.full(effects.shape, np.nan, dtype=float),
                where=standard_errors > 0,
            )
            p_values = 2.0 * stats.t.sf(np.abs(t_statistics), df=donor_count - 1)
            rows.append(
                pd.DataFrame(
                    {
                        "cluster": str(cluster_id),
                        "gene": gene_names[start:end],
                        "n_paired_donors": donor_count,
                        "mean_log2_cpm_difference": effects,
                        "standard_error": standard_errors,
                        "t_statistic": t_statistics,
                        "p_value": p_values,
                        "donor_detection_fraction_cluster": np.mean(
                            cluster_dense > 0, axis=0
                        ),
                        "donor_detection_fraction_rest": np.mean(
                            rest_dense > 0, axis=0
                        ),
                    }
                )
            )
    output = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    if not output.empty:
        output["fdr_bh"] = output.groupby("cluster", observed=True)[
            "p_value"
        ].transform(benjamini_hochberg)
        output = output.sort_values(
            ["cluster", "fdr_bh", "p_value", "gene"], kind="mergesort"
        )
    return output, pd.DataFrame(audit_rows)


def report_trdv_trgv_expression(
    adata: Any, config: Mapping[str, Any]
) -> tuple[pd.DataFrame, pd.DataFrame]:
    settings = config["analysis"]["separate_tcr_expression"]
    prefixes = tuple(str(value).upper() for value in settings["gene_prefixes"])
    gene_indices = np.asarray(
        [
            index
            for index, gene in enumerate(adata.var_names)
            if str(gene).upper().startswith(prefixes)
        ],
        dtype=np.int64,
    )
    if gene_indices.size == 0:
        raise PrerequisiteError(
            "No TRDV/TRGV genes are present for the required separate expression report."
        )
    counts_key = str(config["prerequisites"]["counts_layer"])
    counts = sparse.csr_matrix(adata.layers[counts_key])
    valid_donor, donor_codes, donor_names = _donor_codes(adata.obs, config)
    selected_cells = np.flatnonzero(valid_donor)
    cluster_column = str(config["analysis"]["cluster_column"])
    clusters = (
        adata.obs.iloc[selected_cells][cluster_column]
        .astype(str)
        .to_numpy(dtype=object)
    )
    group_keys = np.asarray(
        [
            f"{donor_names[donor_codes[index]]}||{clusters[index]}"
            for index in range(selected_cells.size)
        ],
        dtype=object,
    )
    group_codes, group_names = pd.factorize(group_keys, sort=True)
    group_map = sparse.csr_matrix(
        (
            np.ones(selected_cells.size, dtype=np.float32),
            (group_codes, np.arange(selected_cells.size, dtype=np.int64)),
        ),
        shape=(len(group_names), selected_cells.size),
    )
    selected_counts = counts[selected_cells]
    tcr_counts = selected_counts[:, gene_indices].tocsr()
    group_counts = (group_map @ tcr_counts).toarray().astype(np.float64, copy=False)
    group_detection = (group_map @ (tcr_counts > 0).astype(np.float32)).toarray()
    group_cells = np.asarray(group_map.sum(axis=1)).ravel()
    cell_library = np.asarray(selected_counts.sum(axis=1), dtype=np.float64).ravel()
    group_library = np.asarray(
        group_map @ cell_library[:, None], dtype=np.float64
    ).ravel()
    target = float(settings["normalization_target"])
    cpm = group_counts / np.maximum(group_library[:, None], 1.0) * target
    detection = group_detection / np.maximum(group_cells[:, None], 1.0)
    genes = adata.var_names[gene_indices].astype(str).to_numpy(dtype=object)
    records: list[dict[str, Any]] = []
    minimum_cells = int(settings["minimum_cells_per_donor_cluster"])
    for group_index, group_name in enumerate(group_names):
        donor_key, cluster_id = str(group_name).rsplit("||", 1)
        study_id, donor_id = donor_key.split("||", 1)
        for gene_index, gene in enumerate(genes):
            records.append(
                {
                    "source_gse_id": study_id,
                    "donor_id": donor_id,
                    "donor_unit": donor_key,
                    "cluster": cluster_id,
                    "gene": gene,
                    "gene_family": "TRDV"
                    if str(gene).upper().startswith("TRDV")
                    else "TRGV",
                    "n_cells": int(group_cells[group_index]),
                    "eligible_group_size": bool(
                        group_cells[group_index] >= minimum_cells
                    ),
                    "counts": float(group_counts[group_index, gene_index]),
                    "cpm": float(cpm[group_index, gene_index]),
                    "cell_detection_fraction": float(
                        detection[group_index, gene_index]
                    ),
                }
            )
    donor_long = pd.DataFrame(records)
    eligible = donor_long[donor_long["eligible_group_size"]]
    summary = (
        eligible.groupby(
            ["cluster", "gene", "gene_family"], observed=True, as_index=False
        )
        .agg(
            n_donors=("donor_unit", "nunique"),
            median_donor_cpm=("cpm", "median"),
            mean_donor_cpm=("cpm", "mean"),
            median_cell_detection_fraction=("cell_detection_fraction", "median"),
        )
        .sort_values(
            ["cluster", "gene_family", "median_donor_cpm", "gene"],
            ascending=[True, True, False, True],
        )
    )
    return donor_long, summary


def _analysis_plan(config: Mapping[str, Any]) -> dict[str, Any]:
    analysis = config["analysis"]
    return {
        "schema_version": 1,
        "abundance_unit": "donor and sample; cells are never inferential replicates",
        "abundance_estimand": "within-gdT cluster fraction",
        "targeted_cohorts": {
            "descriptive_reporting": "included",
            "abundance_inference": "excluded",
            "rule": analysis["targeted_cohort_exclusion"],
        },
        "within_study_contrasts": analysis["within_study_contrasts"],
        "meta_analysis": {
            **analysis["random_effects_meta_analysis"],
            "input": "study-specific donor-level log-odds effects",
        },
        "primary_de": {
            **analysis["donor_pseudobulk_de"],
            "replicate": "donor",
            "tcr_genes_in_primary_results": False,
        },
        "tcr_expression": {
            **analysis["separate_tcr_expression"],
            "separate_from_primary_de": True,
            "selection_use": False,
        },
    }


def _prepare_analysis_outputs(config: Mapping[str, Any]) -> tuple[Path, Path, Path]:
    execution = config["execution"]
    output_dir = resolve_path(execution["analysis_table_dir"])
    manifest = resolve_path(execution["analysis_manifest"])
    log_dir = resolve_path(execution["build_log_dir"])
    if output_dir.exists() or manifest.exists():
        raise PrerequisiteError(
            "Analysis outputs already exist and overwrite is disabled: "
            f"{output_dir if output_dir.exists() else manifest}"
        )
    if not log_dir.is_dir():
        raise PrerequisiteError(
            "Build log directory is absent despite a completed build manifest."
        )
    temporary = output_dir.parent / f".{output_dir.name}.tmp-{uuid.uuid4().hex}"
    temporary.mkdir(parents=True, exist_ok=False)
    return output_dir, temporary, manifest


def _setup_logging(config: Mapping[str, Any]) -> Path:
    log_dir = resolve_path(config["execution"]["build_log_dir"])
    path = (
        log_dir
        / f"analyze_gdt_atlas_{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}_{uuid.uuid4().hex[:8]}.log"
    )
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[
            logging.FileHandler(path, mode="x", encoding="utf-8"),
            logging.StreamHandler(),
        ],
        force=True,
    )
    return path


def run_analysis(
    config: Mapping[str, Any],
    config_path: Path,
    approved: ApprovedInput,
    atlas_path: Path,
    build_manifest: dict[str, Any],
) -> dict[str, Any]:
    import anndata as ad

    output_dir, temporary_dir, manifest_path = _prepare_analysis_outputs(config)
    log_path = _setup_logging(config)
    try:
        logging.info("Reading completed gdT atlas: %s", atlas_path)
        adata = ad.read_h5ad(atlas_path)
        build_uns = adata.uns.get("gdt_atlas_build", {})
        if build_uns.get("config_sha256") != approved.config_sha256:
            raise PrerequisiteError(
                "Atlas-embedded config digest does not match the approved config."
            )
        if build_uns.get("source_sha256") != approved.input_sha256:
            raise PrerequisiteError(
                "Atlas-embedded source digest does not match the approved rebuilt object."
            )
        sample, donor, donor_condition, composition_audit = build_composition_tables(
            adata.obs, config
        )
        sample.to_csv(temporary_dir / "sample_cluster_composition.csv.gz", index=False)
        donor.to_csv(temporary_dir / "donor_cluster_composition.csv.gz", index=False)
        donor_condition.to_csv(
            temporary_dir / "donor_condition_cluster_composition.csv.gz", index=False
        )
        composition_audit.to_csv(
            temporary_dir / "composition_inference_eligibility.csv", index=False
        )
        sample.loc[sample["eligible_for_abundance_inference"]].to_csv(
            temporary_dir / "sample_cluster_composition_inference_only.csv.gz",
            index=False,
        )
        donor.loc[donor["eligible_for_abundance_inference"]].to_csv(
            temporary_dir / "donor_cluster_composition_inference_only.csv.gz",
            index=False,
        )

        study_contrasts, contrast_plan = run_within_study_contrasts(
            donor_condition, config
        )
        study_contrasts.to_csv(
            temporary_dir / "within_study_abundance_contrasts.csv", index=False
        )
        contrast_plan.to_csv(
            temporary_dir / "within_study_contrast_plan.csv", index=False
        )
        meta = run_random_effects_meta_analysis(study_contrasts, config)
        meta.to_csv(temporary_dir / "random_effects_meta_analysis.csv", index=False)

        de, de_audit = donor_pseudobulk_cluster_vs_rest(adata, config)
        de.to_csv(temporary_dir / "donor_pseudobulk_primary_de.csv.gz", index=False)
        de_audit.to_csv(temporary_dir / "donor_pseudobulk_de_audit.csv", index=False)
        excluded_tcr = pd.DataFrame(
            {
                "gene": [
                    str(gene)
                    for gene in adata.var_names
                    if is_tcr_gene(
                        str(gene),
                        config["analysis"]["donor_pseudobulk_de"][
                            "tcr_gene_prefixes_excluded"
                        ],
                    )
                ]
            }
        )
        excluded_tcr.to_csv(
            temporary_dir / "tcr_genes_excluded_from_primary_de.csv", index=False
        )
        if not de.empty and de["gene"].isin(set(excluded_tcr["gene"])).any():
            raise AssertionError(
                "Primary pseudobulk DE contains a prohibited TCR gene."
            )

        tcr_donor, tcr_summary = report_trdv_trgv_expression(adata, config)
        tcr_donor.to_csv(
            temporary_dir / "trdv_trgv_donor_cluster_expression.csv.gz", index=False
        )
        tcr_summary.to_csv(
            temporary_dir / "trdv_trgv_cluster_expression_summary.csv", index=False
        )
        atomic_write_json(temporary_dir / "analysis_plan.json", _analysis_plan(config))
        temporary_dir.replace(output_dir)
    except Exception:
        shutil.rmtree(temporary_dir, ignore_errors=True)
        raise

    manifest = {
        "schema_version": 1,
        "status": "complete",
        "workflow_id": config["workflow_id"],
        "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        "config_path": str(config_path.resolve()),
        "config_sha256": approved.config_sha256,
        "input_h5ad": str(approved.input_h5ad),
        "input_sha256": approved.input_sha256,
        "atlas_h5ad": str(atlas_path),
        "atlas_sha256": build_manifest["atlas_sha256"],
        "analysis_table_dir": str(output_dir),
        "analysis_log": str(log_path),
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "sample_units": int(sample["_sample_unit"].nunique())
        if not sample.empty
        else 0,
        "donor_units": int(donor["_donor_unit"].nunique()) if not donor.empty else 0,
        "abundance_inference_excludes_targeted_cohorts": True,
        "within_study_effects": int(study_contrasts.shape[0]),
        "random_effects_results": int(meta.shape[0]),
        "primary_de_rows": int(de.shape[0]),
        "primary_de_excludes_tcr_genes": True,
        "trdv_trgv_expression_rows": int(tcr_donor.shape[0]),
    }
    atomic_write_json(manifest_path, manifest)
    logging.info("Completed donor-aware gdT atlas analysis: %s", output_dir)
    return manifest


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    config_path = args.config.expanduser().resolve()
    config = load_config(config_path)
    validate_runtime_environment(config)
    approved = validate_prerequisites(config, config_path)
    atlas_path, build_manifest_path, build_manifest = validate_completed_build(
        config, approved
    )
    if args.preflight_only:
        print(
            json.dumps(
                {
                    "status": "completed_build_preflight_passed",
                    "workflow_id": config["workflow_id"],
                    "input_sha256": approved.input_sha256,
                    "atlas_h5ad": str(atlas_path),
                    "atlas_sha256": build_manifest["atlas_sha256"],
                    "build_manifest": str(build_manifest_path),
                },
                indent=2,
                sort_keys=True,
            )
        )
        return 0
    run_analysis(config, config_path, approved, atlas_path, build_manifest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
