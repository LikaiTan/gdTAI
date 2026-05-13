#!/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python
"""Apply the selected gdT classifier to a local H5AD.

This is an inference-only helper for same-server testing. It never mutates the
input H5AD. It loads the trusted local selected-model pickle, rebuilds the exact
log1p-CP10K gene features used at training time, applies the saved threshold and
FOXP3/CD4/low-CD3 death-penalty rule, then writes prediction tables and a small
static report.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
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

from run_gdt_prediction_package_evaluation import read_obs_column, read_string_dataset


DEFAULT_MODEL = Path(
    "Integrated_dataset/models/gdT_prediction_classifier/gse144469_holdout_tcrgene/selected_model.pkl"
)
TARGET_SUM = 10_000.0
DEFAULT_OBS_COLUMNS = [
    "source_gse_id",
    "sample_id",
    "sampleid",
    "library_id",
    "tissue_corrected",
    "tissue",
    "simple_annotation_plus6",
    "leiden",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Apply selected gdT classifier to a local H5AD.")
    parser.add_argument("--input-h5ad", type=Path, required=True)
    parser.add_argument("--model-pkl", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--dataset-id", default=None, help="Output dataset ID. Defaults to input H5AD stem.")
    parser.add_argument("--chunk-size", type=int, default=50_000)
    parser.add_argument("--obs-column", action="append", default=[], help="Extra obs column to copy/summarize.")
    parser.add_argument("--max-umap-points", type=int, default=300_000)
    parser.add_argument("--seed", type=int, default=20260507)
    parser.add_argument("--no-figures", action="store_true")
    return parser.parse_args()


def setup_dirs(dataset_id: str) -> tuple[Path, Path, Path, Path]:
    table_dir = Path("Integrated_dataset/tables/gdT_prediction/external_tests") / dataset_id
    figure_dir = Path("Integrated_dataset/figures/gdT_prediction/external_tests") / dataset_id
    log_dir = Path("Integrated_dataset/logs/gdT_prediction/external_tests") / dataset_id
    static_dir = Path("gdT_prediction/external_tests") / dataset_id
    for path in [table_dir, figure_dir, log_dir, static_dir]:
        path.mkdir(parents=True, exist_ok=True)
    return table_dir, figure_dir, log_dir, static_dir


def setup_logging(log_dir: Path) -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.FileHandler(log_dir / "predict_with_selected_gdt_model.log", mode="w"), logging.StreamHandler()],
        force=True,
    )


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def load_model(path: Path) -> dict[str, Any]:
    with path.open("rb") as fh:
        payload = pickle.load(fh)
    required = ["model_object", "threshold", "feature_names", "gene_names"]
    missing = [key for key in required if key not in payload]
    if missing:
        raise KeyError(f"Model pickle missing required keys: {missing}")
    return payload


def infer_matrix_encoding(handle: h5py.File) -> str:
    x = handle["X"]
    if isinstance(x, h5py.Dataset):
        return "dense"
    if not isinstance(x, h5py.Group):
        raise TypeError("Unsupported H5AD X object.")
    encoding = x.attrs.get("encoding-type", b"")
    if isinstance(encoding, bytes):
        encoding = encoding.decode("utf-8")
    if encoding:
        return str(encoding)
    keys = set(x.keys())
    if {"data", "indices", "indptr"}.issubset(keys):
        return "csr_matrix_or_csc_matrix"
    raise TypeError(f"Unsupported H5AD X group keys: {sorted(keys)}")


def read_obs_index(handle: h5py.File) -> np.ndarray:
    if "_index" in handle["obs"]:
        return read_string_dataset(handle["obs"]["_index"])
    n_obs = int(handle["X"].attrs.get("shape", [handle["obs"].shape[0]])[0])
    return np.asarray([str(i) for i in range(n_obs)], dtype=object)


def read_selected_obs(handle: h5py.File, columns: list[str]) -> dict[str, np.ndarray]:
    out = {}
    present = set(handle["obs"].keys())
    for column in columns:
        if column in present and column not in out:
            values = read_obs_column(handle, column)
            out[column] = np.asarray(values, dtype=object)
    return out


def feature_column_indices(gene_names: list[str]) -> dict[str, int | None]:
    lookup = {gene: i for i, gene in enumerate(gene_names)}
    return {
        "FOXP3": lookup.get("FOXP3"),
        "CD4": lookup.get("CD4"),
        "CD3D": lookup.get("CD3D"),
        "CD3E": lookup.get("CD3E"),
        "CD3G": lookup.get("CD3G"),
    }


def apply_death_penalty(score: np.ndarray, x_features: np.ndarray, gene_names: list[str]) -> tuple[np.ndarray, np.ndarray]:
    cols = feature_column_indices(gene_names)
    penalty = np.zeros(score.shape[0], dtype=bool)
    if cols["FOXP3"] is not None:
        penalty |= x_features[:, int(cols["FOXP3"])] > 0.25
    if cols["CD4"] is not None:
        penalty |= x_features[:, int(cols["CD4"])] > 0.75
    cd3_cols = [cols[name] for name in ["CD3D", "CD3E", "CD3G"] if cols[name] is not None]
    if cd3_cols:
        penalty |= x_features[:, [int(col) for col in cd3_cols]].sum(axis=1) < 0.25
    out = score.astype(np.float32, copy=True)
    out[penalty] = np.minimum(out[penalty], 0.03)
    return out, penalty


def build_gene_mapping(handle: h5py.File, gene_names: list[str]) -> tuple[dict[int, int], pd.DataFrame]:
    var_names = pd.Index(read_string_dataset(handle["var"]["_index"]), dtype="string").astype(str).tolist()
    gene_to_idx = {gene: idx for idx, gene in enumerate(var_names)}
    rows = []
    mapping: dict[int, int] = {}
    for feature_col, gene in enumerate(gene_names):
        gene_idx = gene_to_idx.get(gene)
        rows.append(
            {
                "feature_index": feature_col,
                "gene": gene,
                "available_in_h5ad": gene_idx is not None,
                "h5ad_gene_index": gene_idx if gene_idx is not None else "",
            }
        )
        if gene_idx is not None:
            mapping[int(gene_idx)] = feature_col
    return mapping, pd.DataFrame(rows)


def extract_csr_features(
    handle: h5py.File,
    start_row: int,
    end_row: int,
    gene_idx_to_col: dict[int, int],
    n_features: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x_group = handle["X"]
    indptr = x_group["indptr"]
    indices = x_group["indices"]
    data = x_group["data"]
    n_rows = end_row - start_row
    features = np.zeros((n_rows, n_features), dtype=np.float32)
    row_sum = np.zeros(n_rows, dtype=np.float32)
    n_detected = np.zeros(n_rows, dtype=np.int32)
    selected = np.asarray(sorted(gene_idx_to_col), dtype=np.int64)
    for out_i, obs_i in enumerate(range(start_row, end_row)):
        row_start = int(indptr[obs_i])
        row_end = int(indptr[obs_i + 1])
        if row_end <= row_start:
            continue
        row_indices = indices[row_start:row_end].astype(np.int64, copy=False)
        row_data = data[row_start:row_end].astype(np.float32, copy=False)
        total = float(np.sum(row_data, dtype=np.float64))
        row_sum[out_i] = total
        n_detected[out_i] = int(np.count_nonzero(row_data > 0))
        if total <= 0 or selected.size == 0:
            continue
        present = np.isin(row_indices, selected, assume_unique=False)
        if not present.any():
            continue
        values = np.log1p(row_data[present] * (TARGET_SUM / total)).astype(np.float32, copy=False)
        for gene_idx, value in zip(row_indices[present], values):
            features[out_i, gene_idx_to_col[int(gene_idx)]] = value
    return features, row_sum, n_detected


def dense_chunk_features(
    handle: h5py.File,
    start_row: int,
    end_row: int,
    gene_idx_to_col: dict[int, int],
    n_features: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x = handle["X"]
    dense = x[start_row:end_row, :]
    row_sum = dense.sum(axis=1).astype(np.float32)
    n_detected = (dense > 0).sum(axis=1).astype(np.int32)
    features = np.zeros((end_row - start_row, n_features), dtype=np.float32)
    for gene_idx, col in gene_idx_to_col.items():
        values = dense[:, gene_idx].astype(np.float32, copy=False)
        with np.errstate(divide="ignore", invalid="ignore"):
            transformed = np.log1p(values * (TARGET_SUM / row_sum))
        transformed[~np.isfinite(transformed)] = 0.0
        features[:, col] = transformed.astype(np.float32, copy=False)
    return features, row_sum, n_detected


def prediction_summary(pred: np.ndarray, score: np.ndarray, penalty: np.ndarray, total_cells: int) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "total_cells": int(total_cells),
                "predicted_gdT_cells": int(pred.sum()),
                "predicted_gdT_fraction": float(pred.mean()) if total_cells else np.nan,
                "median_score": float(np.median(score)) if score.size else np.nan,
                "mean_score": float(np.mean(score)) if score.size else np.nan,
                "death_penalty_cells": int(penalty.sum()),
                "death_penalty_fraction": float(penalty.mean()) if total_cells else np.nan,
            }
        ]
    )


def grouped_summary(values: np.ndarray, pred: np.ndarray, score: np.ndarray, penalty: np.ndarray, key: str) -> pd.DataFrame:
    df = pd.DataFrame(
        {
            key: pd.Series(values, dtype="string").fillna(""),
            "pred": pred.astype(bool),
            "score": score.astype(np.float32),
            "penalty": penalty.astype(bool),
        }
    )
    rows = []
    for value, group in df.groupby(key, sort=True, dropna=False):
        n = int(group.shape[0])
        rows.append(
            {
                key: str(value),
                "total_cells": n,
                "predicted_gdT_cells": int(group["pred"].sum()),
                "predicted_gdT_fraction": float(group["pred"].mean()) if n else np.nan,
                "median_score": float(group["score"].median()) if n else np.nan,
                "death_penalty_cells": int(group["penalty"].sum()),
                "death_penalty_fraction": float(group["penalty"].mean()) if n else np.nan,
            }
        )
    return pd.DataFrame(rows).sort_values("predicted_gdT_cells", ascending=False)


def write_prediction_csv(
    path: Path,
    obs_index: np.ndarray,
    obs_columns: dict[str, np.ndarray],
    score: np.ndarray,
    raw_score: np.ndarray,
    pred: np.ndarray,
    penalty: np.ndarray,
    row_sum: np.ndarray,
    n_detected: np.ndarray,
    missing_gene_count: int,
    chunk_size: int,
) -> None:
    columns = list(obs_columns)
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        for start in range(0, score.size, chunk_size):
            end = min(start + chunk_size, score.size)
            data: dict[str, Any] = {
                "cell_id": obs_index[start:end],
                "gdt_score": score[start:end],
                "gdt_score_before_penalty": raw_score[start:end],
                "predicted_gdT": pred[start:end],
                "death_penalty_applied": penalty[start:end],
                "row_sum_x": row_sum[start:end],
                "n_detected_genes_x": n_detected[start:end],
                "missing_model_gene_count": missing_gene_count,
            }
            for column in columns:
                data[column] = obs_columns[column][start:end]
            df = pd.DataFrame(data)
            df.to_csv(handle, index=False, header=(start == 0))


def plot_histogram(score: np.ndarray, pred: np.ndarray, threshold: float, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    ax.hist(score[~pred], bins=80, alpha=0.75, label="predicted non-gdT", color="#5f6b75")
    ax.hist(score[pred], bins=80, alpha=0.85, label="predicted gdT", color="#b51f2a")
    ax.axvline(threshold, color="black", lw=1.2, linestyle="--", label=f"threshold {threshold:.4f}")
    ax.set_xlabel("selected model score")
    ax.set_ylabel("cells")
    ax.set_title("gdT classifier score distribution")
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
    fig, ax = plt.subplots(figsize=(7, 6), constrained_layout=True)
    idx_neg = idx[~pred[idx]]
    idx_pos = idx[pred[idx]]
    ax.scatter(umap[idx_neg, 0], umap[idx_neg, 1], s=1, c="#b7bec5", alpha=0.25, linewidths=0)
    ax.scatter(umap[idx_pos, 0], umap[idx_pos, 1], s=2, c="#b51f2a", alpha=0.85, linewidths=0)
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_title("Predicted gdT cells")
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return True


def write_report(
    *,
    path_md: Path,
    path_html: Path,
    dataset_id: str,
    input_h5ad: Path,
    model_pkl: Path,
    model_sha: str,
    payload: dict[str, Any],
    matrix_encoding: str,
    gene_availability: pd.DataFrame,
    overall: pd.DataFrame,
    summary_paths: list[Path],
    figure_paths: list[Path],
) -> None:
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

    missing = gene_availability.loc[~gene_availability["available_in_h5ad"], "gene"].astype(str).tolist()
    lines = [
        f"# gdT classifier external test: {dataset_id}",
        "",
        f"- input H5AD: `{input_h5ad}`",
        f"- model: `{model_pkl}`",
        f"- model SHA256: `{model_sha}`",
        f"- model name: `{payload.get('model', '')}`",
        f"- threshold: `{float(payload['threshold']):.16g}`",
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
        "- Missing model genes were filled with zero.",
        "- This is a high-confidence gdT recovery classifier, not a clinical classifier.",
        "",
        "## Outputs",
    ]
    for output in summary_paths + figure_paths:
        lines.append(f"- `{output}`")
    if missing:
        lines.extend(["", "## Missing Genes", ", ".join(missing)])
    path_md.write_text("\n".join(lines) + "\n", encoding="utf-8")
    html_rows = overall.to_html(index=False, escape=True)
    fig_html = "\n".join(
        f"<figure><img src='{html.escape(str(fig.resolve()))}'><figcaption>{html.escape(fig.name)}</figcaption></figure>"
        for fig in figure_paths
    )
    html_text = f"""<!doctype html><html><head><meta charset='utf-8'>
<title>gdT classifier external test: {html.escape(dataset_id)}</title>
<style>
body{{font-family:Arial,Helvetica,sans-serif;max-width:1100px;margin:24px auto;color:#1e252b;line-height:1.5}}
code{{background:#eef1f4;padding:1px 4px;border-radius:3px}}
table{{border-collapse:collapse;font-size:13px}}th,td{{border:1px solid #d8dde2;padding:5px 7px;text-align:left}}
img{{max-width:100%;border:1px solid #d8dde2}}figure{{margin:20px 0}}
</style></head><body>
<h1>gdT classifier external test: {html.escape(dataset_id)}</h1>
<p><strong>Input:</strong> <code>{html.escape(str(input_h5ad))}</code></p>
<p><strong>Model:</strong> <code>{html.escape(str(model_pkl))}</code></p>
<p><strong>Threshold:</strong> <code>{float(payload['threshold']):.16g}</code></p>
<p><strong>Model genes available:</strong> {int(gene_availability['available_in_h5ad'].sum())} / {gene_availability.shape[0]}</p>
<h2>Overall Prediction</h2>{html_rows}
<h2>Caveats</h2><ul>
<li>Inference-only; the input H5AD was not modified.</li>
<li>The model expects count-like X.</li>
<li>Missing model genes were filled with zero.</li>
</ul>
<h2>Figures</h2>{fig_html}
</body></html>"""
    path_html.write_text(html_text, encoding="utf-8")


def main() -> None:
    args = parse_args()
    dataset_id = args.dataset_id or args.input_h5ad.stem
    table_dir, figure_dir, log_dir, static_dir = setup_dirs(dataset_id)
    setup_logging(log_dir)
    logging.info("Input H5AD: %s", args.input_h5ad)
    logging.info("Model pickle: %s", args.model_pkl)
    payload = load_model(args.model_pkl)
    model_sha = sha256_file(args.model_pkl)
    threshold = float(payload["threshold"])
    gene_names = [str(gene) for gene in payload["gene_names"]]
    model = payload["model_object"]
    obs_columns_requested = DEFAULT_OBS_COLUMNS + args.obs_column

    with h5py.File(args.input_h5ad, "r") as handle:
        matrix_encoding = infer_matrix_encoding(handle)
        if matrix_encoding not in {"csr_matrix", "csr_matrix_or_csc_matrix", "dense"}:
            raise TypeError(
                f"Unsupported X encoding `{matrix_encoding}`. Convert to CSR or provide a count-like dense H5AD."
            )
        n_obs = int(handle["obs"]["_index"].shape[0])
        obs_index = read_obs_index(handle)
        obs_columns = read_selected_obs(handle, obs_columns_requested)
        gene_idx_to_col, gene_availability = build_gene_mapping(handle, gene_names)
        missing_gene_count = int((~gene_availability["available_in_h5ad"]).sum())
        gene_availability.to_csv(table_dir / "model_gene_availability.csv", index=False)
        score = np.zeros(n_obs, dtype=np.float32)
        raw_score = np.zeros(n_obs, dtype=np.float32)
        pred = np.zeros(n_obs, dtype=bool)
        penalty = np.zeros(n_obs, dtype=bool)
        row_sum = np.zeros(n_obs, dtype=np.float32)
        n_detected = np.zeros(n_obs, dtype=np.int32)
        extract_fn = dense_chunk_features if matrix_encoding == "dense" else extract_csr_features
        for start in range(0, n_obs, args.chunk_size):
            end = min(start + args.chunk_size, n_obs)
            features, chunk_sum, chunk_detected = extract_fn(
                handle,
                start,
                end,
                gene_idx_to_col,
                len(gene_names),
            )
            chunk_raw = model.predict_proba(features)[:, 1].astype(np.float32)
            chunk_score, chunk_penalty = apply_death_penalty(chunk_raw, features, gene_names)
            raw_score[start:end] = chunk_raw
            score[start:end] = chunk_score
            pred[start:end] = chunk_score >= threshold
            penalty[start:end] = chunk_penalty
            row_sum[start:end] = chunk_sum
            n_detected[start:end] = chunk_detected
            logging.info("Predicted %s / %s cells", f"{end:,}", f"{n_obs:,}")

        prediction_csv = table_dir / "gdt_predictions.csv.gz"
        write_prediction_csv(
            prediction_csv,
            obs_index,
            obs_columns,
            score,
            raw_score,
            pred,
            penalty,
            row_sum,
            n_detected,
            missing_gene_count,
            args.chunk_size,
        )
        overall = prediction_summary(pred, score, penalty, n_obs)
        overall["dataset_id"] = dataset_id
        overall["threshold"] = threshold
        overall["model_sha256"] = model_sha
        overall["missing_model_gene_count"] = missing_gene_count
        overall["matrix_encoding"] = matrix_encoding
        overall.to_csv(table_dir / "prediction_summary_overall.csv", index=False)
        summary_paths = [prediction_csv, table_dir / "prediction_summary_overall.csv", table_dir / "model_gene_availability.csv"]
        for column, values in obs_columns.items():
            group_df = grouped_summary(values, pred, score, penalty, column)
            group_path = table_dir / f"prediction_summary_by_{column}.csv"
            group_df.to_csv(group_path, index=False)
            summary_paths.append(group_path)
        figure_paths: list[Path] = []
        if not args.no_figures:
            hist_path = figure_dir / "prediction_score_histogram.png"
            plot_histogram(score, pred, threshold, hist_path)
            figure_paths.append(hist_path)
            umap_path = figure_dir / "prediction_umap.png"
            if plot_umap(handle, pred, umap_path, args.max_umap_points, args.seed):
                figure_paths.append(umap_path)

    report_md = log_dir / "selected_gdt_model_external_test_report.md"
    report_html = static_dir / "index.html"
    write_report(
        path_md=report_md,
        path_html=report_html,
        dataset_id=dataset_id,
        input_h5ad=args.input_h5ad,
        model_pkl=args.model_pkl,
        model_sha=model_sha,
        payload=payload,
        matrix_encoding=matrix_encoding,
        gene_availability=gene_availability,
        overall=overall,
        summary_paths=summary_paths,
        figure_paths=figure_paths,
    )
    manifest = {
        "dataset_id": dataset_id,
        "input_h5ad": str(args.input_h5ad),
        "model_pkl": str(args.model_pkl),
        "model_sha256": model_sha,
        "threshold": threshold,
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
    logging.info("Wrote prediction report: %s", report_html)


if __name__ == "__main__":
    main()
