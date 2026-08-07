"""Nested project-data orchestration for the approved gdTAI V4.1-GPU Gate C."""

from __future__ import annotations

import hashlib
import fcntl
import json
import logging
import math
import os
import subprocess
import time
from dataclasses import asdict
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd

import run_gdtai_v4_nested_evaluation as legacy
from gdtai_v4_gpu_core import (
    AtomicCandidateCheckpoint,
    GpuFit,
    decision_dict,
    fit_gpu_binary_estimator,
    json_safe,
    sha256_file,
    stage1_threshold_frontier,
    stage2_threshold_frontier,
)
from gdtai_v4_nested_core import (
    DERIVED_FEATURE_ORDER,
    PlattCalibrator,
    balanced_training_weights,
    candidate_id,
    fit_platt_calibrator,
    stable_seed,
)


PROJECT_ROOT = Path(__file__).resolve().parents[2]


class ProjectGpuLock:
    def __init__(self, config: Mapping[str, Any]):
        self.config = config
        self.handle: Any = None

    def __enter__(self) -> "ProjectGpuLock":
        path = Path("/tmp/gdtai-v41-gpu.lock")
        self.handle = path.open("w")
        try:
            fcntl.flock(self.handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise RuntimeError("Another gdTAI V4.1-GPU process holds the project GPU lock") from error
        self.handle.write(f"pid={os.getpid()}\n")
        self.handle.flush()
        self._validate_competing_clients()
        return self

    def _validate_competing_clients(self) -> None:
        query = subprocess.run(
            ["nvidia-smi", "--query-compute-apps=pid,used_memory", "--format=csv,noheader,nounits"],
            capture_output=True, text=True, check=True, timeout=30,
        )
        competitors = []
        for line in query.stdout.splitlines():
            fields = [field.strip() for field in line.split(",")]
            if len(fields) != 2:
                continue
            pid, memory = int(fields[0]), int(fields[1])
            if pid != os.getpid() and memory > int(self.config["environment"]["maximum_competing_gpu_memory_mib"]):
                competitors.append({"pid": pid, "used_memory_mib": memory})
        if competitors:
            raise RuntimeError(f"Competing GPU clients exceed the frozen memory limit: {competitors}")
        utilization = subprocess.run(
            ["nvidia-smi", "--query-gpu=utilization.gpu", "--format=csv,noheader,nounits"],
            capture_output=True, text=True, check=True, timeout=30,
        )
        initial = int(utilization.stdout.strip().splitlines()[0])
        if initial > int(self.config["environment"]["maximum_competing_gpu_utilization_percent"]):
            time.sleep(60)
            repeated = subprocess.run(
                ["nvidia-smi", "--query-gpu=utilization.gpu", "--format=csv,noheader,nounits"],
                capture_output=True, text=True, check=True, timeout=30,
            )
            final = int(repeated.stdout.strip().splitlines()[0])
            if final > int(self.config["environment"]["maximum_competing_gpu_utilization_percent"]):
                raise RuntimeError(f"Sustained competing GPU utilization is {final}%")

    def __exit__(self, *_: Any) -> None:
        if self.handle is not None:
            fcntl.flock(self.handle.fileno(), fcntl.LOCK_UN)
            self.handle.close()


def resolve(value: str | Path) -> Path:
    path = Path(value)
    return path if path.is_absolute() else PROJECT_ROOT / path


def array_sha256(values: np.ndarray) -> str:
    array = np.ascontiguousarray(values)
    return hashlib.sha256(memoryview(array)).hexdigest()


def safe_name(value: str) -> str:
    return hashlib.sha256(value.encode()).hexdigest()[:16] + "_" + "".join(
        character if character.isalnum() or character in "-_" else "_" for character in value
    )[:80]


def stage_candidates(config: Mapping[str, Any], stage: str) -> list[tuple[str, dict[str, Any]]]:
    spec = config["models"][stage]
    candidates = [("torch_ridge", {"l2_strength": float(value)}) for value in spec["torch_ridge_l2"]]
    candidates.extend(
        ("xgboost", dict(spec["xgboost_fixed"]) | dict(parameters))
        for parameters in spec["xgboost"]
    )
    return candidates


def comparator_candidates(config: Mapping[str, Any]) -> list[tuple[str, dict[str, Any]]]:
    return [
        ("torch_ridge", {"l2_strength": float(value)})
        for value in config["models"]["stage2"]["torch_ridge_l2"]
    ]


def fit_gpu(
    matrix: np.ndarray,
    labels: np.ndarray,
    weights: np.ndarray,
    names: Sequence[str],
    gene_count: int,
    family: str,
    parameters: Mapping[str, Any],
    config: Mapping[str, Any],
    seed: int,
) -> GpuFit:
    return fit_gpu_binary_estimator(
        matrix, labels, weights, names, gene_count, family, parameters,
        float(config["feature_policy"]["minimum_detection_fraction"]),
        int(config["feature_policy"]["maximum_retained_genes"]), seed,
        config["models"]["torch_fixed"], config["models"]["xgboost_fixed"],
    )


def checkpoint_contract(
    config: Mapping[str, Any], outer_fold_id: str, stage: str, model_id: str,
    train_rows: np.ndarray, assignments: np.ndarray,
) -> dict[str, Any]:
    cache = json.loads((resolve(config["outputs"]["cache_dir"]) / "feature_cache_manifest.json").read_text())
    return {
        "protocol_version": config["protocol_version"],
        "cell_manifest_sha256": config["preflight"]["cell_manifest_sha256"],
        "split_manifest_sha256": config["preflight"]["split_manifest_sha256"],
        "feature_manifest_sha256": config["preflight"]["feature_manifest_sha256"],
        "gene_feature_sha256": cache["gene_feature_sha256"],
        "derived_feature_sha256": cache["derived_feature_sha256"],
        "outer_fold_id": outer_fold_id,
        "stage": stage,
        "candidate_id": model_id,
        "train_rows_sha256": array_sha256(np.asarray(train_rows, dtype=np.int64)),
        "assignments_sha256": array_sha256(np.asarray(assignments, dtype=np.int8)),
        "config_sha256": config["_gate_c_code_sha256"]["config"],
        "runner_sha256": config["_gate_c_code_sha256"]["runner"],
        "core_sha256": config["_gate_c_code_sha256"]["core"],
        "project_runner_sha256": config["_gate_c_code_sha256"]["project_runner"],
        "environment": config["_gpu_environment"],
    }


def crossfit_gpu(
    *, matrix: np.ndarray, feature_names: Sequence[str], gene_feature_count: int,
    train_rows: np.ndarray, assignments: np.ndarray, labels: np.ndarray,
    sources: np.ndarray, reliability: np.ndarray, family: str,
    parameters: Mapping[str, Any], config: Mapping[str, Any], outer_fold_id: str,
    stage: str,
) -> legacy.CrossfitResult:
    model_id = candidate_id(family, parameters)
    checkpoint = AtomicCandidateCheckpoint(
        resolve(config["outputs"]["checkpoint_dir"]) / outer_fold_id / stage,
        safe_name(model_id),
        checkpoint_contract(config, outer_fold_id, stage, model_id, train_rows, assignments),
    )
    cached = checkpoint.load()
    if cached is not None:
        metadata, arrays = cached["metadata"], cached["arrays"]
        calibrator = PlattCalibrator(float(metadata["calibrator_intercept"]), float(metadata["calibrator_slope"]))
        logging.info("Resumed %s/%s from atomic checkpoint", stage, model_id)
        return legacy.CrossfitResult(
            candidate_id=model_id, family=family, parameters=dict(parameters), train_rows=train_rows,
            calibrated_oof=arrays["calibrated_oof"], calibrated_control=np.asarray([], dtype=np.float64),
            calibrator=calibrator, fold_models=[],
            retained_feature_counts=[int(value) for value in arrays["retained_feature_counts"]],
            retained_gene_counts=[int(value) for value in arrays["retained_gene_counts"]],
            fold_iterations=[int(value) for value in arrays["fold_iterations"]],
            fold_converged=[bool(value) for value in arrays["fold_converged"]],
        )

    raw_oof = np.full(train_rows.size, np.nan, dtype=np.float64)
    retained_features, retained_genes, iterations, converged = [], [], [], []
    model_root = checkpoint.path / "fold_models"
    for fold in (0, 1, 2):
        fit_mask = assignments != fold
        validation_mask = assignments == fold
        weights = balanced_training_weights(labels[fit_mask], sources[fit_mask], reliability[fit_mask])
        model = fit_gpu(
            np.asarray(matrix[train_rows[fit_mask]], dtype=np.float32), labels[fit_mask], weights,
            feature_names, gene_feature_count, family, parameters, config,
            stable_seed(int(config["random_seed"]), outer_fold_id, stage, model_id, fold),
        )
        if not model.converged:
            logging.warning("Ineligible fit %s/%s/fold%d: %s", stage, model_id, fold, model.convergence_reason)
        raw_oof[validation_mask] = model.predict_probability(np.asarray(matrix[train_rows[validation_mask]], dtype=np.float32))
        model.export(model_root / f"fold_{fold}")
        retained_features.append(len(model.selected_feature_names))
        retained_genes.append(int((model.selected_columns < gene_feature_count).sum()))
        iterations.append(model.n_iterations)
        converged.append(model.converged)
        del model
    if not np.isfinite(raw_oof).all():
        raise RuntimeError(f"Incomplete OOF probabilities for {stage}/{model_id}")
    weights = balanced_training_weights(labels, sources, reliability)
    calibrator = fit_platt_calibrator(
        raw_oof, labels, weights,
        stable_seed(int(config["random_seed"]), outer_fold_id, stage, model_id, "calibration"),
    )
    calibrated = calibrator.predict(raw_oof)
    checkpoint.save(
        {
            "calibrated_oof": calibrated,
            "retained_feature_counts": np.asarray(retained_features, dtype=np.int16),
            "retained_gene_counts": np.asarray(retained_genes, dtype=np.int16),
            "fold_iterations": np.asarray(iterations, dtype=np.int16),
            "fold_converged": np.asarray(converged, dtype=bool),
        },
        {
            "candidate_id": model_id, "family": family, "parameters": dict(parameters),
            "calibrator_intercept": calibrator.intercept, "calibrator_slope": calibrator.slope,
            "all_folds_converged": bool(all(converged)),
        },
    )
    return legacy.CrossfitResult(
        candidate_id=model_id, family=family, parameters=dict(parameters), train_rows=train_rows,
        calibrated_oof=calibrated, calibrated_control=np.asarray([], dtype=np.float64), calibrator=calibrator,
        fold_models=[], retained_feature_counts=retained_features, retained_gene_counts=retained_genes,
        fold_iterations=iterations, fold_converged=converged,
    )


def save_frontier(frontier: pd.DataFrame, config: Mapping[str, Any], outer: str, stage: str, model_id: str, mode: str) -> Path:
    path = resolve(config["outputs"]["table_dir"]) / "threshold_frontiers" / outer / stage
    path.mkdir(parents=True, exist_ok=True)
    output = path / f"{safe_name(model_id)}__{mode}.parquet"
    temporary = output.with_suffix(".tmp.parquet")
    frontier.to_parquet(temporary, index=False, compression="zstd")
    os.replace(temporary, output)
    return output


def candidate_table_row(result: legacy.CrossfitResult, decisions: Mapping[str, Any], heldout: str, stage: str) -> dict[str, Any]:
    row = {
        "heldout_source": heldout, "stage": stage, "candidate_id": result.candidate_id,
        "family": result.family, "parameters_json": json.dumps(result.parameters, sort_keys=True),
        "mean_retained_features": float(np.mean(result.retained_feature_counts)),
        "min_retained_features": int(min(result.retained_feature_counts)),
        "max_retained_features": int(max(result.retained_feature_counts)),
        "mean_retained_genes": float(np.mean(result.retained_gene_counts)),
        "all_folds_converged": bool(all(result.fold_converged)),
        "converged_fold_count": int(sum(result.fold_converged)),
        "maximum_fold_iterations": int(max(result.fold_iterations)),
        "calibrator_intercept": result.calibrator.intercept,
        "calibrator_slope": result.calibrator.slope,
    }
    for mode, decision in decisions.items():
        row.update({f"{mode}_{key}": value for key, value in decision_dict(decision).items()})
    return row


def inner_assignments_with_hard_nk(
    cells: pd.DataFrame, splits: pd.DataFrame, outer_fold_id: str,
    base_rows: np.ndarray, hard_nk_rows: np.ndarray,
) -> np.ndarray:
    base_assign = legacy.inner_assignments(cells, splits, outer_fold_id, "stage2", base_rows)
    frozen = splits[(splits.outer_fold_id == outer_fold_id) & (splits.stage == "stage1")]
    mapping = frozen.set_index("group_key")["inner_fold"]
    hard = cells.loc[hard_nk_rows, "group_key"].map(mapping)
    if hard.isna().any():
        raise RuntimeError("Hard NK group is absent from frozen Stage-1 split")
    return np.concatenate([base_assign, hard.to_numpy(dtype=np.int8)])


def select_stage1(
    cells: pd.DataFrame, splits: pd.DataFrame, gene_matrix: np.ndarray,
    feature_names: list[str], heldout: str, outer_fold_id: str,
    config: Mapping[str, Any],
) -> dict[str, Any]:
    indices = [feature_names.index(gene) for gene in config["feature_policy"]["stage1_genes"]]
    matrix = np.asarray(gene_matrix[:, indices], dtype=np.float32)
    names = [feature_names[index] for index in indices]
    eligible = cells.stage1_role.isin(["t_positive", "nk_negative"]) & cells.source_gse_id.ne(heldout)
    train_rows = cells.index[eligible].to_numpy(dtype=np.int64)
    assignments = legacy.inner_assignments(cells, splits, outer_fold_id, "stage1", train_rows)
    local = cells.loc[train_rows]
    labels = local.stage1_role.eq("t_positive").to_numpy(dtype=np.int8)
    sources = local.source_gse_id.to_numpy(dtype=object)
    reliability = local.stage1_weight.to_numpy(dtype=np.float32)
    primary = [value for value in legacy.config_primary_sources(config) if value != heldout]
    records, results = [], {}
    for family, parameters in stage_candidates(config, "stage1"):
        result = crossfit_gpu(
            matrix=matrix, feature_names=names, gene_feature_count=len(names), train_rows=train_rows,
            assignments=assignments, labels=labels, sources=sources, reliability=reliability,
            family=family, parameters=parameters, config=config, outer_fold_id=outer_fold_id, stage="stage1",
        )
        frontier, decision = stage1_threshold_frontier(
            result.calibrated_oof, sources,
            local.truth_class.eq("gdT_primary").to_numpy(bool),
            (local.has_any_ab_tcr.astype(bool) & local.stage1_role.eq("t_positive")).to_numpy(bool),
            local.stage1_role.eq("nk_negative").to_numpy(bool), primary, config["stage1_guardrails"],
        )
        save_frontier(frontier, config, outer_fold_id, "stage1", result.candidate_id, "stage1")
        records.append({"result": result, "decision": decision})
        results[result.candidate_id] = result
    table = pd.DataFrame([candidate_table_row(row["result"], {"stage1": row["decision"]}, heldout, "stage1") for row in records])
    eligible_records = [row for row in records if row["decision"].passed and all(row["result"].fold_converged)]
    if not eligible_records:
        return {"selected": None, "table": table, "matrix": matrix, "names": names, "train_rows": train_rows}
    selected_row = sorted(
        eligible_records,
        key=lambda row: (-row["decision"].objective_value, -row["decision"].macro_recall,
                         np.mean(row["result"].retained_gene_counts), row["result"].candidate_id),
    )[0]
    selected = selected_row["result"]
    weights = balanced_training_weights(labels, sources, reliability)
    final = fit_gpu(matrix[train_rows], labels, weights, names, len(names), selected.family,
                    selected.parameters, config, stable_seed(config["random_seed"], outer_fold_id, "stage1", "outer"))
    if not final.converged:
        raise RuntimeError(f"Selected Stage-1 outer fit did not converge: {final.convergence_reason}")
    final.export(resolve(config["outputs"]["checkpoint_dir"]) / outer_fold_id / "selected_stage1")
    probability = selected.calibrator.predict(final.predict_probability(matrix))
    probability[train_rows] = selected.calibrated_oof
    return {
        "selected": selected, "decision": selected_row["decision"], "final": final,
        "calibrator": selected.calibrator, "probability": probability, "table": table,
        "matrix": matrix, "names": names, "train_rows": train_rows,
    }


def evaluate_stage2_grid(
    *, cells: pd.DataFrame, matrix: np.ndarray, names: list[str], gene_count: int,
    base_rows: np.ndarray, hard_nk_rows: np.ndarray, assignments: np.ndarray,
    control_rows: np.ndarray, p1: np.ndarray, p1_threshold: float,
    heldout: str, outer_fold_id: str, config: Mapping[str, Any],
    stage_name: str, candidates: Sequence[tuple[str, dict[str, Any]]],
) -> tuple[list[dict[str, Any]], dict[str, legacy.CrossfitResult]]:
    train_rows = np.concatenate([base_rows, hard_nk_rows])
    local = cells.loc[train_rows]
    labels = np.concatenate([
        cells.loc[base_rows, "stage2_role"].eq("positive").to_numpy(dtype=np.int8),
        np.zeros(hard_nk_rows.size, dtype=np.int8),
    ])
    sources = local.source_gse_id.to_numpy(dtype=object)
    reliability = np.concatenate([
        cells.loc[base_rows, "truth_reliability"].to_numpy(dtype=np.float32),
        cells.loc[hard_nk_rows, "stage1_weight"].to_numpy(dtype=np.float32),
    ])
    primary_mask = local.truth_class.isin(["gdT_primary", "abT_primary"]).to_numpy(bool)
    primary_rows = train_rows[primary_mask]
    primary = cells.loc[primary_rows]
    hard_position = {int(row): base_rows.size + index for index, row in enumerate(hard_nk_rows)}
    records, results = [], {}
    for family, parameters in candidates:
        result = crossfit_gpu(
            matrix=matrix, feature_names=names, gene_feature_count=gene_count, train_rows=train_rows,
            assignments=assignments, labels=labels, sources=sources, reliability=reliability,
            family=family, parameters=parameters, config=config, outer_fold_id=outer_fold_id, stage=stage_name,
        )
        control_probability = np.zeros(control_rows.size, dtype=np.float64)
        for index, row in enumerate(control_rows):
            position = hard_position.get(int(row))
            if position is not None:
                control_probability[index] = result.calibrated_oof[position]
        decisions = {}
        for mode, mode_spec in config["operating_modes"].items():
            frontier, decision = stage2_threshold_frontier(
                result.calibrated_oof[primary_mask],
                primary.truth_class.eq("gdT_primary").to_numpy(dtype=np.int8),
                primary.source_gse_id.to_numpy(dtype=object),
                primary.truth_class.eq("abT_primary").to_numpy(bool),
                control_probability, p1[primary_rows], p1[control_rows], p1_threshold,
                primary.exclusion_union.to_numpy(bool), cells.loc[control_rows, "exclusion_union"].to_numpy(bool),
                mode, mode_spec,
            )
            save_frontier(frontier, config, outer_fold_id, stage_name, result.candidate_id, mode)
            decisions[mode] = decision
        records.append({"result": result, "decisions": decisions})
        results[result.candidate_id] = result
    return records, results


def select_stage2_record(records: Sequence[Mapping[str, Any]]) -> Mapping[str, Any] | None:
    eligible = [row for row in records if row["decisions"]["balanced"].passed and all(row["result"].fold_converged)]
    if not eligible:
        return None
    return sorted(eligible, key=lambda row: (
        -row["decisions"]["balanced"].objective_value,
        row["decisions"]["balanced"].strict_nk_fpr,
        row["decisions"]["balanced"].paired_abt_fpr,
        np.mean(row["result"].retained_gene_counts), row["result"].candidate_id,
    ))[0]


def model_importance(model: GpuFit, heldout: str) -> pd.DataFrame:
    if model.family == "torch_ridge":
        coefficient = model.model["coefficient"].detach().cpu().numpy().astype(float)
        importance = np.abs(coefficient)
        direction = np.sign(coefficient).astype(np.int8)
        method = "standardized_gpu_ridge_coefficient"
    else:
        booster = model.model.get_booster()
        gain = booster.get_score(importance_type="gain")
        importance = np.asarray([gain.get(f"f{index}", 0.0) for index in range(len(model.selected_feature_names))])
        direction = np.zeros(len(importance), dtype=np.int8)
        method = "gpu_xgboost_gain"
    normalized = importance / importance.sum() if importance.sum() else np.zeros_like(importance)
    table = pd.DataFrame({
        "heldout_source": heldout, "family": model.family, "method": method,
        "feature": model.selected_feature_names, "importance": importance,
        "normalized_importance": normalized, "direction": direction,
    })
    table["rank"] = table.importance.rank(method="first", ascending=False).astype(int)
    return table.sort_values("rank").reset_index(drop=True)


def fit_selected_outer(
    matrix: np.ndarray, names: list[str], gene_count: int, cells: pd.DataFrame,
    base_rows: np.ndarray, hard_rows: np.ndarray, selected: legacy.CrossfitResult,
    config: Mapping[str, Any], outer: str, stage: str,
) -> GpuFit:
    rows = np.concatenate([base_rows, hard_rows])
    labels = np.concatenate([
        cells.loc[base_rows, "stage2_role"].eq("positive").to_numpy(np.int8),
        np.zeros(hard_rows.size, dtype=np.int8),
    ])
    sources = cells.loc[rows, "source_gse_id"].to_numpy(object)
    reliability = np.concatenate([
        cells.loc[base_rows, "truth_reliability"].to_numpy(np.float32),
        cells.loc[hard_rows, "stage1_weight"].to_numpy(np.float32),
    ])
    weights = balanced_training_weights(labels, sources, reliability)
    model = fit_gpu(matrix[rows], labels, weights, names, gene_count, selected.family,
                    selected.parameters, config, stable_seed(config["random_seed"], outer, stage, "outer"))
    if not model.converged:
        raise RuntimeError(f"Selected {stage} outer fit did not converge: {model.convergence_reason}")
    model.export(resolve(config["outputs"]["checkpoint_dir"]) / outer / f"selected_{stage}")
    return model


def evaluate_comparator(
    *, model_id: str, matrix: np.ndarray, names: list[str], cells: pd.DataFrame,
    splits: pd.DataFrame, heldout: str, outer: str, config: Mapping[str, Any],
) -> dict[str, Any]:
    base_mask = cells.stage2_role.isin(["positive", "negative"]) & cells.source_gse_id.ne(heldout)
    base_rows = cells.index[base_mask].to_numpy(np.int64)
    control_rows = cells.index[cells.stage1_role.eq("nk_negative") & cells.source_gse_id.ne(heldout)].to_numpy(np.int64)
    hard_rows = control_rows.copy()
    assignments = inner_assignments_with_hard_nk(cells, splits, outer, base_rows, hard_rows)
    ones = np.ones(cells.shape[0], dtype=np.float64)
    records, _ = evaluate_stage2_grid(
        cells=cells, matrix=matrix, names=names, gene_count=len(names), base_rows=base_rows,
        hard_nk_rows=hard_rows, assignments=assignments, control_rows=control_rows,
        p1=ones, p1_threshold=0.0, heldout=heldout, outer_fold_id=outer, config=config,
        stage_name=model_id, candidates=comparator_candidates(config),
    )
    table = pd.DataFrame([candidate_table_row(row["result"], row["decisions"], heldout, model_id) for row in records])
    chosen = select_stage2_record(records)
    output = {"candidate_table": table, "selected": None, "metrics": [], "predictions": [], "control_metrics": []}
    if chosen is None:
        return output
    selected = chosen["result"]
    model = fit_selected_outer(matrix, names, len(names), cells, base_rows, hard_rows, selected, config, outer, model_id)
    probability = selected.calibrator.predict(model.predict_probability(matrix))
    output["selected"] = {
        "candidate_id": selected.candidate_id, "family": selected.family,
        "parameters": selected.parameters, "model": model, "calibrator": selected.calibrator,
        **{mode: decision_dict(decision) for mode, decision in chosen["decisions"].items()},
    }
    for mode, decision in chosen["decisions"].items():
        if not decision.passed:
            continue
        metrics, predictions = legacy.evaluate_outer_predictions(
            cells=cells, heldout_source=heldout, model_id=model_id, mode=mode,
            stage1_probability=ones, stage2_probability=probability,
            stage1_threshold=0.0, stage2_threshold=decision.threshold,
        )
        output["metrics"].append(metrics); output["predictions"].append(predictions)
        output["control_metrics"].append(legacy.heldout_nk_control_metrics(
            cells=cells, heldout_source=heldout, model_id=model_id, mode=mode,
            stage1_probability=ones, stage2_probability=probability,
            stage1_threshold=0.0, stage2_threshold=decision.threshold,
        ))
    return output


def run_outer_fold(
    *, cells: pd.DataFrame, splits: pd.DataFrame, gene_matrix: np.ndarray,
    derived_matrix: np.ndarray, legacy_score: np.ndarray, feature_names: list[str],
    heldout_source: str, config: Mapping[str, Any],
) -> dict[str, Any]:
    outer_index = legacy.config_primary_sources(config).index(heldout_source)
    outer = f"outer_{outer_index}_{heldout_source}"
    logging.info("Starting V4.1-GPU outer fold %s", outer)
    stage1 = select_stage1(cells, splits, gene_matrix, feature_names, heldout_source, outer, config)
    result: dict[str, Any] = {
        "outer_fold_id": outer, "heldout_source": heldout_source,
        "stage1_candidate_table": stage1["table"], "stage2_candidate_table": pd.DataFrame(),
        "stage1_selected": None, "stage2_selected": None,
        "comparator_candidate_tables": {}, "comparator_selected": {},
        "metrics": [], "predictions": [], "control_metrics": [],
        "feature_importance": pd.DataFrame(),
    }
    if stage1["selected"] is None:
        result["stage1_heldout_guardrails"] = {
            "heldout_source": heldout_source, "pass": False, "reason": "No eligible Stage-1 candidate",
            "gdt_recall": math.nan, "abt_recall": math.nan, "strict_nk_call_rate": math.nan,
        }
    else:
        p1 = stage1["probability"]
        p1_threshold = float(stage1["decision"].threshold)
        result["stage1_selected"] = {
            "candidate_id": stage1["selected"].candidate_id,
            "family": stage1["selected"].family,
            "parameters": stage1["selected"].parameters,
            "threshold": p1_threshold,
        }
        result["stage1_final_model"] = stage1["final"]
        result["stage1_calibrator"] = stage1["calibrator"]
        result["stage1_heldout_guardrails"] = legacy.heldout_stage1_guardrails(cells, heldout_source, p1, p1_threshold, config)
        base_rows = cells.index[cells.stage2_role.isin(["positive", "negative"]) & cells.source_gse_id.ne(heldout_source)].to_numpy(np.int64)
        control_rows = cells.index[cells.stage1_role.eq("nk_negative") & cells.source_gse_id.ne(heldout_source)].to_numpy(np.int64)
        hard_rows = control_rows[p1[control_rows] >= p1_threshold]
        assignments = inner_assignments_with_hard_nk(cells, splits, outer, base_rows, hard_rows)
        stage2_matrix = np.column_stack([gene_matrix, derived_matrix, p1]).astype(np.float32, copy=False)
        stage2_names = [*feature_names, *DERIVED_FEATURE_ORDER, "stage1_probability"]
        records, _ = evaluate_stage2_grid(
            cells=cells, matrix=stage2_matrix, names=stage2_names, gene_count=len(feature_names),
            base_rows=base_rows, hard_nk_rows=hard_rows, assignments=assignments,
            control_rows=control_rows, p1=p1, p1_threshold=p1_threshold, heldout=heldout_source,
            outer_fold_id=outer, config=config, stage_name="stage2_v4_gpu",
            candidates=stage_candidates(config, "stage2"),
        )
        result["stage2_candidate_table"] = pd.DataFrame([
            candidate_table_row(row["result"], row["decisions"], heldout_source, "stage2_v4_gpu") for row in records
        ])
        chosen = select_stage2_record(records)
        if chosen is not None:
            selected = chosen["result"]
            model = fit_selected_outer(stage2_matrix, stage2_names, len(feature_names), cells, base_rows, hard_rows, selected, config, outer, "stage2_v4_gpu")
            p2 = selected.calibrator.predict(model.predict_probability(stage2_matrix))
            result["stage2_selected"] = {
                "candidate_id": selected.candidate_id, "family": selected.family,
                "parameters": selected.parameters, "hard_nk_count": int(hard_rows.size),
                **{mode: decision_dict(decision) for mode, decision in chosen["decisions"].items()},
            }
            result["stage2_final_model"] = model
            result["stage2_calibrator"] = selected.calibrator
            result["feature_importance"] = model_importance(model, heldout_source)
            for mode, decision in chosen["decisions"].items():
                if not decision.passed:
                    continue
                metrics, predictions = legacy.evaluate_outer_predictions(
                    cells=cells, heldout_source=heldout_source, model_id="v4_nested_selected", mode=mode,
                    stage1_probability=p1, stage2_probability=p2, stage1_threshold=p1_threshold,
                    stage2_threshold=decision.threshold,
                )
                result["metrics"].append(metrics); result["predictions"].append(predictions)
                result["control_metrics"].append(legacy.heldout_nk_control_metrics(
                    cells=cells, heldout_source=heldout_source, model_id="v4_nested_selected", mode=mode,
                    stage1_probability=p1, stage2_probability=p2, stage1_threshold=p1_threshold,
                    stage2_threshold=decision.threshold,
                ))
            del stage2_matrix

    compact_genes = ["TRDC", "TRDV1", "TRDV2", "TRDV3", "TRAC", "TRBC1", "TRBC2"]
    compact_indices = [feature_names.index(gene) for gene in compact_genes]
    tcr_indices = [index for index, gene in enumerate(feature_names) if gene != "TRAT1" and gene.startswith(("TRA", "TRB", "TRG", "TRD"))]
    for model_id, indices in (("compact_7gene_logistic", compact_indices), ("v2_like_tcr_logistic", tcr_indices)):
        comparator = evaluate_comparator(
            model_id=model_id, matrix=np.asarray(gene_matrix[:, indices], dtype=np.float32),
            names=[feature_names[index] for index in indices], cells=cells, splits=splits,
            heldout=heldout_source, outer=outer, config=config,
        )
        result["comparator_candidate_tables"][model_id] = comparator["candidate_table"]
        result["comparator_selected"][model_id] = comparator["selected"]
        result["metrics"].extend(comparator["metrics"]); result["predictions"].extend(comparator["predictions"])
        result["control_metrics"].extend(comparator["control_metrics"])
    legacy_result = legacy.evaluate_legacy_comparator(cells=cells, legacy_score=legacy_score, heldout_source=heldout_source, config=config)
    result["comparator_selected"]["legacy_trd_minus_trab"] = legacy_result["selected"]
    result["metrics"].extend(legacy_result["metrics"]); result["predictions"].extend(legacy_result["predictions"])
    result["control_metrics"].extend(legacy_result["control_metrics"])
    logging.info("Completed V4.1-GPU outer fold %s", outer)
    return result


def save_outer_result(result: Mapping[str, Any], config: Mapping[str, Any]) -> None:
    table_dir = resolve(config["outputs"]["table_dir"])
    table_dir.mkdir(parents=True, exist_ok=True)
    fold = str(result["outer_fold_id"])
    result["stage1_candidate_table"].to_csv(table_dir / f"{fold}_stage1_candidates.csv", index=False)
    if not result["stage2_candidate_table"].empty:
        result["stage2_candidate_table"].to_csv(table_dir / f"{fold}_stage2_candidates.csv", index=False)
    for model_id, table in result["comparator_candidate_tables"].items():
        table.to_csv(table_dir / f"{fold}_{model_id}_candidates.csv", index=False)
    if result["metrics"]:
        pd.DataFrame(result["metrics"]).to_csv(table_dir / f"{fold}_metrics.csv", index=False)
    if result["control_metrics"]:
        pd.DataFrame(result["control_metrics"]).to_csv(table_dir / f"{fold}_negative_controls.csv", index=False)
    if result["predictions"]:
        pd.concat(result["predictions"], ignore_index=True).to_csv(
            table_dir / f"{fold}_predictions.csv.gz", index=False,
            compression={"method": "gzip", "mtime": 0},
        )
    if not result["feature_importance"].empty:
        result["feature_importance"].to_csv(table_dir / f"{fold}_feature_importance.csv", index=False)
    pd.DataFrame([result["stage1_heldout_guardrails"]]).to_csv(table_dir / f"{fold}_stage1_heldout_guardrails.csv", index=False)
    summary = {
        "outer_fold_id": fold, "heldout_source": result["heldout_source"],
        "stage1_selected": result["stage1_selected"], "stage2_selected": result["stage2_selected"],
        "comparator_selected": {
            key: ({name: value for name, value in selected.items() if name not in {"model", "calibrator"}} if isinstance(selected, dict) else selected)
            for key, selected in result["comparator_selected"].items()
        },
        "metrics": result["metrics"], "control_metrics": result["control_metrics"],
        "stage1_heldout_guardrails": result["stage1_heldout_guardrails"],
    }
    legacy.write_json(json_safe(summary), table_dir / f"{fold}_summary.json")


def run_project_evaluation(config: Mapping[str, Any]) -> None:
    started = time.monotonic()
    log_dir = resolve(config["outputs"]["log_dir"])
    log_dir.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.FileHandler(log_dir / "gpu_nested_evaluation.log", mode="a"), logging.StreamHandler()],
        force=True,
    )
    legacy.run_outer_fold = run_outer_fold
    legacy.save_outer_result = save_outer_result
    with ProjectGpuLock(config):
        output = legacy.run_nested_evaluation(config, [])
    elapsed = time.monotonic() - started
    output["runtime_seconds"] = elapsed
    output["gpu_environment"] = config["_gpu_environment"]
    legacy.write_json(json_safe(output), log_dir / "nested_evaluation_summary.json")
    if elapsed > float(config["environment"]["maximum_wall_seconds"]):
        raise RuntimeError(f"Gate C exceeded the frozen wall-time limit: {elapsed:.1f}s")
    report = legacy.render_nested_report(config, no_pdf=True)
    html_path = Path(report["html"])
    document = html_path.read_text()
    document = document.replace("gdTAI V4 Nested Cross-Dataset Evaluation", "gdTAI V4.1-GPU Nested Cross-Dataset Evaluation")
    document = document.replace(
        "Stage 1 is a soft T-lineage elastic-net model. Stage 2 compares elastic-net and histogram gradient boosting",
        "Stage 1 is a multigene T-lineage gate comparing deterministic GPU ridge and GPU XGBoost. Stage 2 compares the same GPU families with fold-local Stage-1-passing strict-NK hard negatives",
    )
    html_path.write_text(document)
    pdf_path = html_path.parent / "gdtai_v4_gpu_nested_evaluation_report.pdf"
    completed = subprocess.run(
        ["google-chrome", "--headless", "--no-sandbox", "--disable-gpu",
         f"--print-to-pdf={pdf_path}", "--no-pdf-header-footer", html_path.resolve().as_uri()],
        capture_output=True, text=True, timeout=180, check=False,
    )
    if completed.returncode != 0 or not pdf_path.exists():
        raise RuntimeError(f"Gate C Chrome PDF rendering failed: {completed.stderr}")
