"""Deterministic GPU fitting and audit utilities for gdTAI V4.1-GPU.

This module owns no project-data selection. It exposes GPU-only model fitting,
portable exports, complete threshold frontiers, and hash-bound checkpoints for
the guarded launcher.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import tempfile
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd

from gdtai_v4_nested_core import (
    ThresholdDecision,
    count_greater_equal,
    f_beta_from_counts,
    fold_feature_mask,
    threshold_candidates,
)


class GpuEnvironmentError(RuntimeError):
    """Raised when the frozen direct-CUDA execution contract is not met."""


def sha256_file(path: Path, chunk_size: int = 16 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_sha256(value: Any) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def json_safe(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, np.ndarray):
        return json_safe(value.tolist())
    if isinstance(value, (np.integer, np.bool_)):
        return value.item()
    if isinstance(value, (float, np.floating)):
        return float(value) if math.isfinite(float(value)) else None
    return value


def atomic_write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=path.parent, delete=False) as handle:
        json.dump(json_safe(payload), handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
        temporary = Path(handle.name)
    os.replace(temporary, path)


def atomic_save_npz(path: Path, **arrays: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(suffix=".npz", dir=path.parent, delete=False) as handle:
        temporary = Path(handle.name)
    try:
        np.savez_compressed(temporary, **arrays)
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def validate_gpu_environment(spec: Mapping[str, Any]) -> dict[str, Any]:
    required = {
        "CUDA_VISIBLE_DEVICES": str(spec["cuda_visible_devices"]),
        "CUBLAS_WORKSPACE_CONFIG": str(spec["cublas_workspace_config"]),
        "PYTHONHASHSEED": str(spec["pythonhashseed"]),
    }
    failures = [f"{key} must equal {value!r}" for key, value in required.items() if os.environ.get(key) != value]
    mps = os.environ.get("CUDA_MPS_PIPE_DIRECTORY", "")
    prefix = str(spec["mps_path_prefix"])
    forbidden = str(spec["forbidden_mps_path"])
    if not mps.startswith(prefix) or mps == prefix:
        failures.append(f"CUDA_MPS_PIPE_DIRECTORY must be a unique path beginning {prefix!r}")
    if Path(mps or "/nonexistent").resolve() == Path(forbidden).resolve():
        failures.append("The shared NVIDIA MPS directory is forbidden")
    if mps and (Path(mps) / "control").exists():
        failures.append("The isolated CUDA directory already contains an MPS control socket")
    if failures:
        raise GpuEnvironmentError("; ".join(failures))

    import torch
    import xgboost

    if not torch.cuda.is_available():
        raise GpuEnvironmentError("PyTorch cannot initialize CUDA; CPU fallback is forbidden")
    torch.use_deterministic_algorithms(True)
    torch.backends.cudnn.benchmark = False
    return {
        "cuda_device_name": torch.cuda.get_device_name(0),
        "cuda_device_count": torch.cuda.device_count(),
        "torch_version": torch.__version__,
        "xgboost_version": xgboost.__version__,
        "cuda_visible_devices": os.environ["CUDA_VISIBLE_DEVICES"],
        "cuda_mps_pipe_directory": mps,
        "cublas_workspace_config": os.environ["CUBLAS_WORKSPACE_CONFIG"],
        "pythonhashseed": os.environ["PYTHONHASHSEED"],
    }


def weighted_standardize(matrix: np.ndarray, weight: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x = np.asarray(matrix, dtype=np.float32)
    w = np.asarray(weight, dtype=np.float64)
    if x.ndim != 2 or w.shape != (x.shape[0],):
        raise ValueError("Weighted standardization shapes disagree")
    if not np.isfinite(x).all() or not np.isfinite(w).all() or (w < 0).any() or w.sum() <= 0:
        raise ValueError("Weighted standardization requires finite data and positive total weight")
    normalized = w / w.sum()
    mean = np.sum(x.astype(np.float64) * normalized[:, None], axis=0)
    variance = np.sum(np.square(x.astype(np.float64) - mean) * normalized[:, None], axis=0)
    scale = np.sqrt(np.maximum(variance, 0.0))
    scale[scale == 0] = 1.0
    transformed = ((x - mean.astype(np.float32)) / scale.astype(np.float32)).astype(np.float32)
    return transformed, mean.astype(np.float32), scale.astype(np.float32)


@dataclass
class GpuFit:
    family: str
    parameters: dict[str, Any]
    selected_columns: np.ndarray
    selected_feature_names: list[str]
    mean: np.ndarray
    scale: np.ndarray
    model: Any
    n_iterations: int
    converged: bool
    convergence_reason: str
    gradient_inf: float
    relative_terminal_loss_change: float
    runtime_seconds: float
    peak_gpu_memory_bytes: int

    def transform(self, matrix: np.ndarray) -> np.ndarray:
        x = np.asarray(matrix, dtype=np.float32)[:, self.selected_columns]
        return ((x - self.mean) / self.scale).astype(np.float32, copy=False)

    def predict_probability(self, matrix: np.ndarray, batch_size: int = 262144) -> np.ndarray:
        x = self.transform(matrix)
        if self.family == "torch_ridge":
            import torch

            coefficient = self.model["coefficient"]
            intercept = self.model["intercept"]
            output: list[np.ndarray] = []
            with torch.no_grad():
                for start in range(0, x.shape[0], batch_size):
                    tensor = torch.as_tensor(x[start : start + batch_size], device="cuda")
                    output.append(torch.sigmoid(tensor @ coefficient + intercept).cpu().numpy())
            return np.concatenate(output).astype(np.float64) if output else np.asarray([], dtype=np.float64)
        if self.family == "xgboost":
            import cupy as cp

            output = []
            for start in range(0, x.shape[0], batch_size):
                prediction = self.model.predict_proba(cp.asarray(x[start : start + batch_size]))[:, 1]
                output.append(cp.asnumpy(prediction) if isinstance(prediction, cp.ndarray) else np.asarray(prediction))
            return np.concatenate(output).astype(np.float64) if output else np.asarray([], dtype=np.float64)
        raise ValueError(f"Unknown GPU family: {self.family}")

    def metadata(self) -> dict[str, Any]:
        gradient_inf = self.gradient_inf if math.isfinite(self.gradient_inf) else None
        terminal_change = (
            self.relative_terminal_loss_change
            if math.isfinite(self.relative_terminal_loss_change)
            else None
        )
        return {
            "family": self.family,
            "parameters": self.parameters,
            "selected_columns": self.selected_columns.tolist(),
            "selected_feature_names": self.selected_feature_names,
            "n_iterations": self.n_iterations,
            "converged": self.converged,
            "convergence_reason": self.convergence_reason,
            "gradient_inf": gradient_inf,
            "relative_terminal_loss_change": terminal_change,
            "runtime_seconds": self.runtime_seconds,
            "peak_gpu_memory_bytes": self.peak_gpu_memory_bytes,
        }

    def export(self, directory: Path) -> dict[str, str]:
        directory.mkdir(parents=True, exist_ok=True)
        metadata_path = directory / "model.json"
        if self.family == "torch_ridge":
            model_path = directory / "model.npz"
            atomic_save_npz(
                model_path,
                coefficient=self.model["coefficient"].detach().cpu().numpy(),
                intercept=np.asarray([self.model["intercept"].detach().cpu().item()], dtype=np.float32),
                mean=self.mean,
                scale=self.scale,
                selected_columns=self.selected_columns,
            )
        else:
            model_path = directory / "model.ubj"
            temporary = directory / "model.tmp.ubj"
            self.model.save_model(temporary)
            os.replace(temporary, model_path)
            atomic_save_npz(
                directory / "preprocessing.npz",
                mean=self.mean,
                scale=self.scale,
                selected_columns=self.selected_columns,
            )
        metadata = self.metadata() | {"model_file": model_path.name}
        atomic_write_json(metadata_path, metadata)
        return {"model": sha256_file(model_path), "metadata": sha256_file(metadata_path)}


def fit_gpu_binary_estimator(
    matrix: np.ndarray,
    labels: np.ndarray,
    sample_weight: np.ndarray,
    feature_names: Sequence[str],
    gene_feature_count: int,
    family: str,
    parameters: Mapping[str, Any],
    minimum_detection_fraction: float,
    maximum_retained_genes: int,
    seed: int,
    torch_spec: Mapping[str, Any],
    xgboost_spec: Mapping[str, Any],
) -> GpuFit:
    import torch

    x = np.asarray(matrix, dtype=np.float32)
    y = np.asarray(labels, dtype=np.float32)
    weight = np.asarray(sample_weight, dtype=np.float32)
    if np.unique(y).size != 2 or x.shape[0] != y.size or y.shape != weight.shape:
        raise ValueError("GPU fitting requires aligned binary labels and weights")
    mask = fold_feature_mask(x, gene_feature_count, minimum_detection_fraction, maximum_retained_genes)
    selected = np.flatnonzero(mask)
    scaled, mean, scale = weighted_standardize(x[:, selected], weight)
    names = [str(feature_names[index]) for index in selected]
    torch.cuda.empty_cache()
    torch.cuda.reset_peak_memory_stats()
    started = time.monotonic()
    gradient_inf = math.nan
    terminal_change = math.nan

    if family == "torch_ridge":
        torch.manual_seed(int(seed))
        torch.cuda.manual_seed_all(int(seed))
        tx = torch.as_tensor(scaled, device="cuda")
        ty = torch.as_tensor(y, device="cuda")
        tw = torch.as_tensor(weight / weight.mean(), device="cuda")
        coefficient = torch.zeros(tx.shape[1], device="cuda", requires_grad=True)
        intercept = torch.zeros((), device="cuda", requires_grad=True)
        maximum = int(torch_spec["maximum_iterations"])
        optimizer = torch.optim.LBFGS(
            [coefficient, intercept], max_iter=maximum, tolerance_grad=0.0,
            tolerance_change=0.0, line_search_fn="strong_wolfe",
        )
        losses: list[float] = []

        def closure() -> Any:
            optimizer.zero_grad(set_to_none=True)
            logits = tx @ coefficient + intercept
            per_cell = torch.nn.functional.binary_cross_entropy_with_logits(logits, ty, reduction="none")
            loss = (per_cell * tw).mean() + float(parameters["l2_strength"]) * coefficient.square().sum()
            loss.backward()
            losses.append(float(loss.detach().cpu()))
            return loss

        optimizer.step(closure)
        state = optimizer.state[coefficient]
        n_iterations = int(state.get("n_iter", maximum))
        with torch.no_grad():
            logits = tx @ coefficient + intercept
        loss = closure()
        del loss, logits, tx, ty, tw
        gradient_inf = max(
            float(coefficient.grad.detach().abs().max().cpu()),
            float(intercept.grad.detach().abs().cpu()),
        )
        if len(losses) >= 2:
            terminal_change = abs(losses[-1] - losses[-2]) / max(abs(losses[-2]), 1.0)
        finite = bool(torch.isfinite(coefficient).all() and torch.isfinite(intercept))
        converged = bool(
            n_iterations < maximum
            and finite
            and gradient_inf <= float(torch_spec["gradient_inf_maximum"])
            and terminal_change <= float(torch_spec["relative_terminal_loss_change_maximum"])
        )
        reason = "PASS" if converged else (
            f"n_iter={n_iterations}; finite={finite}; gradient_inf={gradient_inf:.3g}; "
            f"relative_terminal_loss_change={terminal_change:.3g}"
        )
        model = {"coefficient": coefficient.detach(), "intercept": intercept.detach()}
    elif family == "xgboost":
        import cupy as cp
        import xgboost as xgb

        fixed = dict(xgboost_spec)
        fixed.update(dict(parameters))
        fixed["random_state"] = int(seed)
        fixed.setdefault("device", "cuda")
        fixed.setdefault("tree_method", "hist")
        fixed.setdefault("objective", "binary:logistic")
        fixed.setdefault("eval_metric", "logloss")
        model = xgb.XGBClassifier(**fixed)
        model.fit(cp.asarray(scaled), cp.asarray(y.astype(np.int8)), sample_weight=cp.asarray(weight))
        n_iterations = int(fixed["n_estimators"])
        prediction = model.predict_proba(cp.asarray(scaled[: min(4096, scaled.shape[0])]))[:, 1]
        probe = cp.asnumpy(prediction) if isinstance(prediction, cp.ndarray) else np.asarray(prediction)
        converged = bool(np.isfinite(probe).all())
        reason = "PASS" if converged else "Non-finite XGBoost probabilities"
    else:
        raise ValueError(f"Unknown GPU family: {family}")

    torch.cuda.synchronize()
    return GpuFit(
        family=family,
        parameters=dict(parameters),
        selected_columns=selected,
        selected_feature_names=names,
        mean=mean,
        scale=scale,
        model=model,
        n_iterations=n_iterations,
        converged=converged,
        convergence_reason=reason,
        gradient_inf=gradient_inf,
        relative_terminal_loss_change=terminal_change,
        runtime_seconds=float(time.monotonic() - started),
        peak_gpu_memory_bytes=int(torch.cuda.max_memory_allocated()),
    )


def predict_portable_export(directory: Path, matrix: np.ndarray, batch_size: int = 262144) -> np.ndarray:
    """Load a portable NPZ/UBJ export and score without pickle."""
    metadata = json.loads((directory / "model.json").read_text())
    family = str(metadata["family"])
    x = np.asarray(matrix, dtype=np.float32)
    if family == "torch_ridge":
        import torch

        with np.load(directory / "model.npz") as archive:
            selected = archive["selected_columns"]
            scaled = (x[:, selected] - archive["mean"]) / archive["scale"]
            coefficient = torch.as_tensor(archive["coefficient"], device="cuda")
            intercept = torch.as_tensor(archive["intercept"][0], device="cuda")
        output = []
        with torch.no_grad():
            for start in range(0, scaled.shape[0], batch_size):
                values = torch.as_tensor(scaled[start : start + batch_size], device="cuda")
                output.append(torch.sigmoid(values @ coefficient + intercept).cpu().numpy())
        return np.concatenate(output).astype(np.float64) if output else np.asarray([], dtype=np.float64)
    if family == "xgboost":
        import cupy as cp
        import xgboost as xgb

        with np.load(directory / "preprocessing.npz") as archive:
            selected = archive["selected_columns"]
            scaled = (x[:, selected] - archive["mean"]) / archive["scale"]
        model = xgb.Booster(params={"device": "cuda"})
        model.load_model(directory / "model.ubj")
        model.set_param({"device": "cuda"})
        output = []
        for start in range(0, scaled.shape[0], batch_size):
            values = model.inplace_predict(cp.asarray(scaled[start : start + batch_size]))
            output.append(cp.asnumpy(values) if isinstance(values, cp.ndarray) else np.asarray(values))
        return np.concatenate(output).astype(np.float64) if output else np.asarray([], dtype=np.float64)
    raise ValueError(f"Unknown portable model family: {family}")


def stage1_threshold_frontier(
    probability: np.ndarray,
    sources: np.ndarray,
    gdt_mask: np.ndarray,
    abt_mask: np.ndarray,
    strict_nk_mask: np.ndarray,
    primary_sources: Sequence[str],
    guardrails: Mapping[str, Any],
) -> tuple[pd.DataFrame, ThresholdDecision]:
    score = np.asarray(probability, dtype=np.float64)
    source = np.asarray(sources, dtype=object)
    gdt = np.asarray(gdt_mask, dtype=bool)
    abt = np.asarray(abt_mask, dtype=bool)
    nk = np.asarray(strict_nk_mask, dtype=bool)
    thresholds = threshold_candidates(score)
    data: dict[str, Any] = {"threshold": thresholds}
    gdt_recalls, abt_recalls, nk_passages = [], [], []
    missing: list[str] = []
    for value in primary_sources:
        positive = score[gdt & (source == value)]
        alpha_beta = score[abt & (source == value)]
        if positive.size == 0 or alpha_beta.size == 0:
            missing.append(str(value))
            continue
        gdt_recall = count_greater_equal(positive, thresholds) / positive.size
        abt_recall = count_greater_equal(alpha_beta, thresholds) / alpha_beta.size
        data[f"gdt_recall__{value}"] = gdt_recall
        data[f"abt_recall__{value}"] = abt_recall
        gdt_recalls.append(gdt_recall)
        abt_recalls.append(abt_recall)
    for value in sorted(pd.unique(source[nk]).tolist()):
        values = score[nk & (source == value)]
        passage = count_greater_equal(values, thresholds) / values.size
        data[f"nk_passage__{value}"] = passage
        nk_passages.append(passage)
    if missing or not nk_passages:
        reason = "Missing Stage-1 strata: " + ", ".join(missing or ["strict_NK"])
        empty = pd.DataFrame(data)
        empty["eligible"] = False
        return empty, ThresholdDecision("stage1", False, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, len(thresholds), 0, hashlib.sha256(thresholds.tobytes()).hexdigest(), reason)
    data["minimum_source_gdt_recall"] = np.min(np.vstack(gdt_recalls), axis=0)
    data["minimum_source_abt_recall"] = np.min(np.vstack(abt_recalls), axis=0)
    data["macro_nk_passage"] = np.mean(np.vstack(nk_passages), axis=0)
    data["maximum_source_nk_passage"] = np.max(np.vstack(nk_passages), axis=0)
    frontier = pd.DataFrame(data)
    frontier["eligible"] = (
        frontier["minimum_source_gdt_recall"].ge(float(guardrails["gdt_recall_per_source"]))
        & frontier["minimum_source_abt_recall"].ge(float(guardrails["abt_recall_per_source"]))
        & frontier["maximum_source_nk_passage"].le(float(guardrails["maximum_nk_passage_per_source"]))
    )
    valid = frontier.index[frontier["eligible"]].to_numpy()
    digest = hashlib.sha256(thresholds.tobytes()).hexdigest()
    if valid.size == 0:
        decision = ThresholdDecision("stage1", False, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, len(thresholds), 0, digest, "No threshold satisfies Stage-1 recall and per-source NK passage guardrails")
    else:
        best = int(valid[-1])
        row = frontier.loc[best]
        decision = ThresholdDecision("stage1", True, float(row.threshold), float(1 - row.macro_nk_passage), math.nan, float((row.minimum_source_gdt_recall + row.minimum_source_abt_recall) / 2), math.nan, math.nan, math.nan, float(row.macro_nk_passage), len(thresholds), int(valid.size), digest, "PASS")
    return frontier, decision


def stage2_threshold_frontier(
    primary_probability: np.ndarray,
    labels: np.ndarray,
    sources: np.ndarray,
    paired_abt_mask: np.ndarray,
    nk_probability: np.ndarray,
    stage1_primary_probability: np.ndarray,
    stage1_nk_probability: np.ndarray,
    stage1_threshold: float,
    primary_exclusion: np.ndarray,
    nk_exclusion: np.ndarray,
    mode: str,
    mode_spec: Mapping[str, Any],
) -> tuple[pd.DataFrame, ThresholdDecision]:
    score = np.asarray(primary_probability, dtype=np.float64)
    y = np.asarray(labels, dtype=np.int8)
    source = np.asarray(sources, dtype=object)
    paired = np.asarray(paired_abt_mask, dtype=bool)
    nk_score = np.asarray(nk_probability, dtype=np.float64)
    effective = np.where((np.asarray(stage1_primary_probability) >= stage1_threshold) & (~np.asarray(primary_exclusion, bool)), score, -np.inf)
    nk_effective = np.where((np.asarray(stage1_nk_probability) >= stage1_threshold) & (~np.asarray(nk_exclusion, bool)), nk_score, -np.inf)
    thresholds = threshold_candidates(score, nk_score)
    data: dict[str, Any] = {"threshold": thresholds}
    precisions, recalls, f1s, f05s = [], [], [], []
    for value in sorted(pd.unique(source).tolist()):
        pos = effective[(source == value) & (y == 1)]
        neg = effective[(source == value) & (y == 0)]
        if not pos.size or not neg.size:
            raise ValueError(f"Source {value} lacks a Stage-2 class")
        tp = count_greater_equal(pos, thresholds).astype(float)
        fp = count_greater_equal(neg, thresholds).astype(float)
        fn = pos.size - tp
        precision = np.divide(tp, tp + fp, out=np.zeros_like(tp), where=(tp + fp) > 0)
        recall = tp / pos.size
        f1 = f_beta_from_counts(tp, fp, fn, 1.0)
        f05 = f_beta_from_counts(tp, fp, fn, 0.5)
        data[f"precision__{value}"], data[f"recall__{value}"] = precision, recall
        precisions.append(precision); recalls.append(recall); f1s.append(f1); f05s.append(f05)
    data["macro_precision"] = np.mean(np.vstack(precisions), axis=0)
    data["macro_recall"] = np.mean(np.vstack(recalls), axis=0)
    data["macro_f1"] = np.mean(np.vstack(f1s), axis=0)
    data["macro_f0_5"] = np.mean(np.vstack(f05s), axis=0)
    data["paired_abt_fpr"] = count_greater_equal(effective[paired], thresholds) / int(paired.sum()) if paired.any() else np.nan
    data["strict_nk_fpr"] = count_greater_equal(nk_effective, thresholds) / nk_effective.size if nk_effective.size else np.nan
    frontier = pd.DataFrame(data)
    frontier["eligible"] = (
        frontier.macro_recall.ge(float(mode_spec["minimum_macro_recall"]))
        & frontier.paired_abt_fpr.le(float(mode_spec["maximum_paired_abt_fpr"]))
        & frontier.strict_nk_fpr.le(float(mode_spec["maximum_strict_nk_fpr"]))
    )
    objective_column = "macro_f1" if mode_spec["objective"] == "f1" else "macro_f0_5"
    valid = frontier.index[frontier.eligible].to_numpy()
    digest = hashlib.sha256(thresholds.tobytes()).hexdigest()
    if not valid.size:
        decision = ThresholdDecision(mode, False, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, len(thresholds), 0, digest, "No threshold satisfies the frozen operating-mode guardrails")
    else:
        objective = frontier.loc[valid, objective_column]
        tied = objective.index[np.isclose(objective, objective.max(), atol=1e-12, rtol=0)].tolist()
        best = sorted(tied, key=lambda index: (frontier.loc[index, "strict_nk_fpr"], frontier.loc[index, "paired_abt_fpr"], -frontier.loc[index, "threshold"]))[0]
        row = frontier.loc[best]
        decision = ThresholdDecision(mode, True, float(row.threshold), float(row[objective_column]), float(row.macro_precision), float(row.macro_recall), float(row.macro_f1), float(row.macro_f0_5), float(row.paired_abt_fpr), float(row.strict_nk_fpr), len(thresholds), int(valid.size), digest, "PASS")
    return frontier, decision


def frontier_diagnostics(frontier: pd.DataFrame, objective: str) -> dict[str, Any]:
    objective_column = "macro_f1" if objective == "f1" else "macro_f0_5"
    unconstrained = frontier.loc[int(frontier[objective_column].idxmax())].to_dict()
    constrained = frontier[frontier.eligible]
    if not constrained.empty:
        near_miss = constrained.loc[int(constrained[objective_column].idxmax())].to_dict()
    else:
        violation_columns = [column for column in ("paired_abt_fpr", "strict_nk_fpr") if column in frontier]
        rank = np.zeros(frontier.shape[0], dtype=np.float64)
        for column in violation_columns:
            rank += frontier[column].rank(method="min", pct=True).to_numpy()
        rank += (1.0 - frontier.get("macro_recall", pd.Series(np.zeros(frontier.shape[0])))).to_numpy()
        near_miss = frontier.loc[int(np.argmin(rank))].to_dict()
    return {"best_unconstrained": unconstrained, "best_or_closest_guardrail": near_miss}


class AtomicCandidateCheckpoint:
    def __init__(self, root: Path, candidate_key: str, contract: Mapping[str, Any]):
        self.path = root / candidate_key
        self.contract = dict(contract)
        self.contract_sha256 = canonical_sha256(self.contract)

    def load(self) -> dict[str, Any] | None:
        metadata_path = self.path / "checkpoint.json"
        probability_path = self.path / "probabilities.npz"
        if not metadata_path.exists() and not probability_path.exists():
            return None
        if not metadata_path.exists() or not probability_path.exists():
            raise RuntimeError(f"Incomplete atomic checkpoint: {self.path}")
        metadata = json.loads(metadata_path.read_text())
        if metadata.get("contract_sha256") != self.contract_sha256 or metadata.get("contract") != self.contract:
            raise RuntimeError(f"Checkpoint contract mismatch: {self.path}")
        if metadata.get("probability_sha256") != sha256_file(probability_path):
            raise RuntimeError(f"Checkpoint probability checksum mismatch: {self.path}")
        with np.load(probability_path) as archive:
            arrays = {name: archive[name] for name in archive.files}
        return {"metadata": metadata, "arrays": arrays}

    def save(self, arrays: Mapping[str, np.ndarray], metadata: Mapping[str, Any]) -> dict[str, Any]:
        self.path.mkdir(parents=True, exist_ok=True)
        probability_path = self.path / "probabilities.npz"
        atomic_save_npz(probability_path, **arrays)
        payload = dict(metadata) | {
            "contract": self.contract,
            "contract_sha256": self.contract_sha256,
            "probability_sha256": sha256_file(probability_path),
            "complete": True,
        }
        atomic_write_json(self.path / "checkpoint.json", payload)
        return payload


def decision_dict(decision: ThresholdDecision) -> dict[str, Any]:
    return asdict(decision)
