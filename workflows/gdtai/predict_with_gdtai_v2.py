#!/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python
"""Apply gdTAI v2.0 with selectable operating modes.

gdTAI v2.0 wraps the same frozen individual TCR-gene classifier with two
operating modes:

- high_f1: use the F1-optimized base threshold.
- high_purity: use the NK-optimized annotation-specific thresholds.

The input H5AD is opened read-only. The script rebuilds the same log1p-CP10K
gene features used by the selected model and writes a compact static report.
"""

from __future__ import annotations

# TNK_WORKFLOW_BOOTSTRAP
import sys as _tnk_sys
from pathlib import Path as _TnkPath

_TNK_PROJECT_ROOT = _TnkPath(__file__).resolve().parents[2]
for _tnk_path in (
    _TNK_PROJECT_ROOT,
    _TNK_PROJECT_ROOT / "src",
    _TNK_PROJECT_ROOT / "workflows" / "integration",
    _TNK_PROJECT_ROOT / "workflows" / "intake",
    _TNK_PROJECT_ROOT / "workflows" / "metadata",
    _TNK_PROJECT_ROOT / "workflows" / "analysis",
    _TNK_PROJECT_ROOT / "workflows" / "gdtai",
    _TNK_PROJECT_ROOT / "workflows" / "gdt_atlas",
    _TNK_PROJECT_ROOT / "workflows" / "reporting",
    _TNK_PROJECT_ROOT / "workflows" / "maintenance",
):
    _tnk_value = str(_tnk_path)
    if _tnk_value not in _tnk_sys.path:
        _tnk_sys.path.insert(0, _tnk_value)


import argparse
import gzip
import html
import json
import logging
import pickle
from pathlib import Path
from typing import Any

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from predict_with_selected_gdt_model import (
    DEFAULT_OBS_COLUMNS,
    apply_death_penalty,
    build_gene_mapping,
    dense_chunk_features,
    extract_csr_features,
    grouped_summary,
    infer_matrix_encoding,
    read_obs_index,
    read_selected_obs,
    sha256_file,
    setup_logging,
)


DEFAULT_V2_MODEL = _TNK_PROJECT_ROOT / "Integrated_dataset/models/gdT_prediction_classifier/gdTAI_v2.0/gdTAI_v2_model.pkl"
STATIC_ROOT = _TNK_PROJECT_ROOT / "gdT_prediction/external_tests"
TABLE_ROOT = _TNK_PROJECT_ROOT / "Integrated_dataset/tables/gdT_prediction/external_tests"
FIGURE_ROOT = _TNK_PROJECT_ROOT / "Integrated_dataset/figures/gdT_prediction/external_tests"
LOG_ROOT = _TNK_PROJECT_ROOT / "Integrated_dataset/logs/gdT_prediction/external_tests"
ANNOTATION_COLUMN_DEFAULT = "simple_annotation_plus6"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Apply gdTAI v2.0 to a local H5AD.")
    parser.add_argument("--input-h5ad", type=Path, required=True)
    parser.add_argument("--model-pkl", type=Path, default=DEFAULT_V2_MODEL)
    parser.add_argument("--mode", choices=["high_f1", "high_purity"], default="high_f1")
    parser.add_argument("--dataset-id", default=None, help="Output dataset ID. Defaults to input H5AD stem.")
    parser.add_argument("--annotation-column", default=ANNOTATION_COLUMN_DEFAULT)
    parser.add_argument("--chunk-size", type=int, default=50_000)
    parser.add_argument("--obs-column", action="append", default=[], help="Extra obs column to copy/summarize.")
    parser.add_argument("--max-umap-points", type=int, default=300_000)
    parser.add_argument("--seed", type=int, default=20260617)
    parser.add_argument("--no-figures", action="store_true")
    return parser.parse_args()


def setup_dirs(output_id: str) -> tuple[Path, Path, Path, Path]:
    table_dir = TABLE_ROOT / output_id
    figure_dir = FIGURE_ROOT / output_id
    log_dir = LOG_ROOT / output_id
    static_dir = STATIC_ROOT / output_id
    for path in [table_dir, figure_dir, log_dir, static_dir]:
        path.mkdir(parents=True, exist_ok=True)
    return table_dir, figure_dir, log_dir, static_dir


def load_v2_model(path: Path) -> dict[str, Any]:
    with path.open("rb") as handle:
        payload = pickle.load(handle)
    required = ["version", "base_model", "operating_modes"]
    missing = [key for key in required if key not in payload]
    if missing:
        raise KeyError(f"gdTAI v2 model pickle missing required keys: {missing}")
    return payload


def normalize_annotation(values: np.ndarray) -> np.ndarray:
    series = pd.Series(values, copy=False).astype("string").fillna("").str.strip().str.upper()
    return series.to_numpy(dtype=object)


def threshold_vector(
    payload: dict[str, Any],
    mode: str,
    n_obs: int,
    annotation_values: np.ndarray | None,
) -> tuple[np.ndarray, np.ndarray]:
    mode_payload = payload["operating_modes"][mode]
    if mode == "high_f1":
        threshold = float(mode_payload["threshold"])
        return np.full(n_obs, threshold, dtype=np.float32), np.full(n_obs, "all_cells", dtype=object)

    thresholds = mode_payload["annotation_thresholds"]
    if annotation_values is None:
        threshold = float(thresholds["other_threshold"])
        return np.full(n_obs, threshold, dtype=np.float32), np.full(n_obs, "missing_annotation_other_threshold", dtype=object)

    annotation = normalize_annotation(annotation_values)
    out = np.full(n_obs, float(thresholds["other_threshold"]), dtype=np.float32)
    label = np.full(n_obs, "other", dtype=object)
    mapping = {
        "GDT_CELL": ("gdt_threshold", "gdT_cell"),
        "CD8_T": ("cd8_threshold", "CD8_T"),
        "CD4_T": ("cd4_threshold", "CD4_T"),
        "TREG": ("treg_threshold", "Treg"),
        "NK_CELL": ("nk_threshold", "NK_cell"),
    }
    for normalized, (key, display) in mapping.items():
        mask = annotation == normalized
        value = thresholds[key]
        out[mask] = np.inf if str(value) == "disabled" else float(value)
        label[mask] = display
    return out, label


def prediction_summary(
    *,
    pred: np.ndarray,
    score: np.ndarray,
    threshold_used: np.ndarray,
    total_cells: int,
    payload: dict[str, Any],
    mode: str,
    annotation_available: bool,
) -> pd.DataFrame:
    mode_payload = payload["operating_modes"][mode]
    finite_thresholds = threshold_used[np.isfinite(threshold_used)]
    return pd.DataFrame(
        [
            {
                "gdtai_version": payload["version"],
                "mode": mode,
                "mode_description": mode_payload["description"],
                "annotation_available": bool(annotation_available),
                "total_cells": int(total_cells),
                "predicted_gdT_cells": int(pred.sum()),
                "predicted_gdT_fraction": float(pred.mean()) if total_cells else np.nan,
                "median_score": float(np.median(score)) if score.size else np.nan,
                "mean_score": float(np.mean(score)) if score.size else np.nan,
                "min_finite_threshold_used": float(np.min(finite_thresholds)) if finite_thresholds.size else np.nan,
                "max_finite_threshold_used": float(np.max(finite_thresholds)) if finite_thresholds.size else np.nan,
            }
        ]
    )


def write_prediction_csv(
    path: Path,
    *,
    obs_index: np.ndarray,
    obs_columns: dict[str, np.ndarray],
    score: np.ndarray,
    raw_score: np.ndarray,
    pred: np.ndarray,
    threshold_used: np.ndarray,
    threshold_annotation: np.ndarray,
    row_sum: np.ndarray,
    n_detected: np.ndarray,
    missing_gene_count: int,
    mode: str,
    chunk_size: int,
) -> None:
    columns = list(obs_columns)
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        for start in range(0, score.size, chunk_size):
            end = min(start + chunk_size, score.size)
            data: dict[str, Any] = {
                "cell_id": obs_index[start:end],
                "gdtai_mode": mode,
                "gdtai_score": score[start:end],
                "gdtai_score_before_penalty": raw_score[start:end],
                "threshold_used": threshold_used[start:end],
                "threshold_annotation": threshold_annotation[start:end],
                "predicted_gdT": pred[start:end],
                "row_sum_x": row_sum[start:end],
                "n_detected_genes_x": n_detected[start:end],
                "missing_model_gene_count": missing_gene_count,
            }
            for column in columns:
                data[column] = obs_columns[column][start:end]
            pd.DataFrame(data).to_csv(handle, index=False, header=(start == 0))


def plot_histogram(score: np.ndarray, pred: np.ndarray, threshold_used: np.ndarray, path: Path, title: str) -> None:
    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    ax.hist(score[~pred], bins=80, alpha=0.75, label="predicted non-gdT", color="#5f6b75")
    ax.hist(score[pred], bins=80, alpha=0.85, label="predicted gdT", color="#b51f2a")
    finite = threshold_used[np.isfinite(threshold_used)]
    if finite.size and np.unique(np.round(finite, 6)).size == 1:
        ax.axvline(float(finite[0]), color="black", lw=1.2, linestyle="--", label=f"threshold {finite[0]:.4f}")
    ax.set_xlabel("gdTAI score")
    ax.set_ylabel("cells")
    ax.set_title(title)
    ax.legend(frameon=False)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def read_umap(handle: h5py.File) -> np.ndarray | None:
    if "obsm" not in handle or "X_umap" not in handle["obsm"]:
        return None
    obj = handle["obsm"]["X_umap"]
    if isinstance(obj, h5py.Dataset):
        arr = obj[:]
    elif isinstance(obj, h5py.Group) and "data" in obj:
        arr = obj["data"][:]
    else:
        return None
    arr = np.asarray(arr)
    if arr.ndim != 2 or arr.shape[1] < 2:
        return None
    return arr[:, :2].astype(np.float32, copy=False)


def plot_umap(handle: h5py.File, pred: np.ndarray, path: Path, max_points: int, seed: int) -> bool:
    umap = read_umap(handle)
    if umap is None:
        return False
    rng = np.random.default_rng(seed)
    idx = np.arange(pred.size)
    if idx.size > max_points:
        positive = idx[pred]
        negative = idx[~pred]
        n_pos = min(positive.size, max_points // 2)
        n_neg = max_points - n_pos
        chosen = []
        if n_pos:
            chosen.append(rng.choice(positive, size=n_pos, replace=False))
        if n_neg and negative.size:
            chosen.append(rng.choice(negative, size=min(n_neg, negative.size), replace=False))
        idx = np.concatenate(chosen) if chosen else np.arange(0)
    idx_neg = idx[~pred[idx]]
    idx_pos = idx[pred[idx]]
    fig, ax = plt.subplots(figsize=(7, 6), constrained_layout=True)
    ax.scatter(umap[idx_neg, 0], umap[idx_neg, 1], s=1, c="#b7bec5", alpha=0.25, linewidths=0)
    ax.scatter(umap[idx_pos, 0], umap[idx_pos, 1], s=2, c="#b51f2a", alpha=0.85, linewidths=0)
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_title("gdTAI v2.0 predicted gdT cells")
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return True


def dataframe_to_markdown(df: pd.DataFrame) -> str:
    view = df.copy()
    for col in view.columns:
        if pd.api.types.is_float_dtype(view[col]):
            view[col] = view[col].map(lambda x: "" if pd.isna(x) else f"{x:.6g}")
    view = view.fillna("").astype(str)
    headers = list(view.columns)
    rows = view.values.tolist()
    out = ["| " + " | ".join(header.replace("|", "\\|") for header in headers) + " |"]
    out.append("| " + " | ".join("---" for _ in headers) + " |")
    for row in rows:
        out.append("| " + " | ".join(str(cell).replace("|", "\\|").replace("\n", " ") for cell in row) + " |")
    return "\n".join(out)


def write_report(
    *,
    path_md: Path,
    path_html: Path,
    dataset_id: str,
    output_id: str,
    input_h5ad: Path,
    model_pkl: Path,
    model_sha: str,
    payload: dict[str, Any],
    mode: str,
    annotation_column: str,
    matrix_encoding: str,
    gene_availability: pd.DataFrame,
    overall: pd.DataFrame,
    summary_paths: list[Path],
    figure_paths: list[Path],
) -> None:
    missing = gene_availability.loc[~gene_availability["available_in_h5ad"], "gene"].astype(str).tolist()
    mode_payload = payload["operating_modes"][mode]
    lines = [
        f"# gdTAI v2.0 external test: {output_id}",
        "",
        f"- input H5AD: `{input_h5ad}`",
        f"- model: `{model_pkl}`",
        f"- model SHA256: `{model_sha}`",
        f"- mode: `{mode}`",
        f"- mode description: {mode_payload['description']}",
        f"- annotation column: `{annotation_column}`",
        f"- matrix encoding: `{matrix_encoding}`",
        f"- model genes available: `{int(gene_availability['available_in_h5ad'].sum())}` / `{gene_availability.shape[0]}`",
        f"- missing model genes: `{len(missing)}`",
        "",
        "## Overall Prediction",
        dataframe_to_markdown(overall),
        "",
        "## Caveats",
        "- This is inference-only; the input H5AD was not modified.",
        "- The model expects count-like `X`; if `X` is normalized/log-transformed, interpret cautiously.",
        "- `high_purity` needs an annotation column. If missing, all cells use the strict `other` threshold.",
        "- Missing model genes were filled with zero.",
        "",
        "## Outputs",
    ]
    for output in summary_paths + figure_paths:
        lines.append(f"- `{output}`")
    if missing:
        lines.extend(["", "## Missing Genes", ", ".join(missing)])
    path_md.write_text("\n".join(lines) + "\n", encoding="utf-8")

    fig_html = "\n".join(
        f"<figure><img src='{html.escape(str(fig.resolve()))}'><figcaption>{html.escape(fig.name)}</figcaption></figure>"
        for fig in figure_paths
    )
    html_text = f"""<!doctype html><html><head><meta charset='utf-8'>
<title>gdTAI v2.0 external test: {html.escape(output_id)}</title>
<style>
body{{font-family:Arial,Helvetica,sans-serif;max-width:1100px;margin:24px auto;color:#1e252b;line-height:1.5}}
code{{background:#eef1f4;padding:1px 4px;border-radius:3px}}
table{{border-collapse:collapse;font-size:13px}}th,td{{border:1px solid #d8dde2;padding:5px 7px;text-align:left}}
img{{max-width:100%;border:1px solid #d8dde2}}figure{{margin:20px 0}}
</style></head><body>
<h1>gdTAI v2.0 external test: {html.escape(output_id)}</h1>
<p><strong>Mode:</strong> <code>{html.escape(mode)}</code> - {html.escape(mode_payload['description'])}</p>
<p><strong>Input:</strong> <code>{html.escape(str(input_h5ad))}</code></p>
<p><strong>Model:</strong> <code>{html.escape(str(model_pkl))}</code></p>
<p><strong>Annotation column:</strong> <code>{html.escape(annotation_column)}</code></p>
<p><strong>Model genes available:</strong> {int(gene_availability['available_in_h5ad'].sum())} / {gene_availability.shape[0]}</p>
<h2>Overall Prediction</h2>{overall.to_html(index=False, escape=True)}
<h2>Caveats</h2><ul>
<li>Inference-only; the input H5AD was not modified.</li>
<li>The model expects count-like X.</li>
<li><code>high_purity</code> needs an annotation column; if missing all cells use the strict other threshold.</li>
<li>Missing model genes were filled with zero.</li>
</ul>
<h2>Figures</h2>{fig_html}
</body></html>"""
    path_html.write_text(html_text, encoding="utf-8")


def main() -> None:
    args = parse_args()
    dataset_id = args.dataset_id or args.input_h5ad.stem
    output_id = f"{dataset_id}_{args.mode}"
    table_dir, figure_dir, log_dir, static_dir = setup_dirs(output_id)
    setup_logging(log_dir)
    payload = load_v2_model(args.model_pkl)
    mode_payload = payload["operating_modes"][args.mode]
    model_sha = sha256_file(args.model_pkl)
    base = payload["base_model"]
    model = base["model_object"]
    model_name = str(base.get("model", ""))
    gene_names = [str(gene) for gene in base["gene_names"]]
    apply_penalty_rule = "penalty" in model_name
    obs_columns_requested = DEFAULT_OBS_COLUMNS + [args.annotation_column] + args.obs_column
    logging.info("Input H5AD: %s", args.input_h5ad)
    logging.info("gdTAI v2.0 model: %s", args.model_pkl)
    logging.info("Mode: %s", args.mode)
    logging.info("Mode description: %s", mode_payload["description"])

    with h5py.File(args.input_h5ad, "r") as handle:
        matrix_encoding = infer_matrix_encoding(handle)
        if matrix_encoding not in {"csr_matrix", "csr_matrix_or_csc_matrix", "dense"}:
            raise TypeError(f"Unsupported X encoding `{matrix_encoding}`. Convert to CSR or provide count-like dense X.")
        n_obs = int(handle["obs"]["_index"].shape[0])
        obs_index = read_obs_index(handle)
        obs_columns = read_selected_obs(handle, obs_columns_requested)
        annotation_values = obs_columns.get(args.annotation_column)
        annotation_available = annotation_values is not None
        threshold_used, threshold_annotation = threshold_vector(payload, args.mode, n_obs, annotation_values)
        gene_idx_to_col, gene_availability = build_gene_mapping(handle, gene_names)
        missing_gene_count = int((~gene_availability["available_in_h5ad"]).sum())
        gene_availability.to_csv(table_dir / "model_gene_availability.csv", index=False)

        score = np.zeros(n_obs, dtype=np.float32)
        raw_score = np.zeros(n_obs, dtype=np.float32)
        pred = np.zeros(n_obs, dtype=bool)
        row_sum = np.zeros(n_obs, dtype=np.float32)
        n_detected = np.zeros(n_obs, dtype=np.int32)
        extract_fn = dense_chunk_features if matrix_encoding == "dense" else extract_csr_features
        for start in range(0, n_obs, args.chunk_size):
            end = min(start + args.chunk_size, n_obs)
            features, chunk_sum, chunk_detected = extract_fn(handle, start, end, gene_idx_to_col, len(gene_names))
            chunk_raw = model.predict_proba(features)[:, 1].astype(np.float32)
            if apply_penalty_rule:
                chunk_score, _ = apply_death_penalty(chunk_raw, features, gene_names)
            else:
                chunk_score = chunk_raw
            raw_score[start:end] = chunk_raw
            score[start:end] = chunk_score
            pred[start:end] = chunk_score >= threshold_used[start:end]
            row_sum[start:end] = chunk_sum
            n_detected[start:end] = chunk_detected
            logging.info("Predicted %s / %s cells", f"{end:,}", f"{n_obs:,}")

        prediction_csv = table_dir / "gdtai_v2_predictions.csv.gz"
        write_prediction_csv(
            prediction_csv,
            obs_index=obs_index,
            obs_columns=obs_columns,
            score=score,
            raw_score=raw_score,
            pred=pred,
            threshold_used=threshold_used,
            threshold_annotation=threshold_annotation,
            row_sum=row_sum,
            n_detected=n_detected,
            missing_gene_count=missing_gene_count,
            mode=args.mode,
            chunk_size=args.chunk_size,
        )
        overall = prediction_summary(
            pred=pred,
            score=score,
            threshold_used=threshold_used,
            total_cells=n_obs,
            payload=payload,
            mode=args.mode,
            annotation_available=annotation_available,
        )
        overall["dataset_id"] = dataset_id
        overall["output_id"] = output_id
        overall["model_sha256"] = model_sha
        overall["missing_model_gene_count"] = missing_gene_count
        overall["matrix_encoding"] = matrix_encoding
        overall.to_csv(table_dir / "prediction_summary_overall.csv", index=False)

        summary_paths = [prediction_csv, table_dir / "prediction_summary_overall.csv", table_dir / "model_gene_availability.csv"]
        for column, values in obs_columns.items():
            group_df = grouped_summary(values, pred, score, np.zeros(n_obs, dtype=bool), column)
            group_path = table_dir / f"prediction_summary_by_{column}.csv"
            group_df.to_csv(group_path, index=False)
            summary_paths.append(group_path)

        figure_paths: list[Path] = []
        if not args.no_figures:
            hist_path = figure_dir / "gdtai_v2_score_histogram.png"
            plot_histogram(score, pred, threshold_used, hist_path, f"gdTAI v2.0 score distribution: {args.mode}")
            figure_paths.append(hist_path)
            umap_path = figure_dir / "gdtai_v2_prediction_umap.png"
            if plot_umap(handle, pred, umap_path, args.max_umap_points, args.seed):
                figure_paths.append(umap_path)

    report_md = log_dir / "gdtai_v2_external_test_report.md"
    report_html = static_dir / "index.html"
    write_report(
        path_md=report_md,
        path_html=report_html,
        dataset_id=dataset_id,
        output_id=output_id,
        input_h5ad=args.input_h5ad,
        model_pkl=args.model_pkl,
        model_sha=model_sha,
        payload=payload,
        mode=args.mode,
        annotation_column=args.annotation_column,
        matrix_encoding=matrix_encoding,
        gene_availability=gene_availability,
        overall=overall,
        summary_paths=summary_paths,
        figure_paths=figure_paths,
    )
    manifest = {
        "dataset_id": dataset_id,
        "output_id": output_id,
        "input_h5ad": str(args.input_h5ad),
        "model_pkl": str(args.model_pkl),
        "model_sha256": model_sha,
        "mode": args.mode,
        "mode_description": mode_payload["description"],
        "annotation_column": args.annotation_column,
        "annotation_available": annotation_available,
        "n_cells": int(n_obs),
        "predicted_gdT_cells": int(pred.sum()),
        "predicted_gdT_fraction": float(pred.mean()) if n_obs else None,
        "missing_model_gene_count": missing_gene_count,
        "table_dir": str(table_dir),
        "figure_dir": str(figure_dir),
        "log_dir": str(log_dir),
        "static_report": str(report_html),
    }
    (log_dir / "prediction_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    logging.info("Wrote gdTAI v2.0 report: %s", report_html)


if __name__ == "__main__":
    main()
