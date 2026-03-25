#!/usr/bin/env python3
"""Rebuild count matrices from project selected inputs and optionally repair an H5AD in place.

This script is for Phase 0 rescue of scanpy-processed H5AD files that no longer
retain raw counts internally. It rebuilds the pre-normalization count matrix from
the project's `selected_inputs.csv`, applies the same coarse QC gates used by the
project scaffold, and checks whether the rebuilt count-space matches the current
processed H5AD on obs/var.

Default mode is read-only validation. Use `--write` only after reviewing the QC.
"""

from __future__ import annotations

import argparse
import gzip
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy import io


PRESERVE_UNS_KEYS = {"Channels", "genome", "modality", "uid"}
NONCODING_PATTERNS = ("MALAT1", "NEAT1", "XIST", "LINC", "MIR", "SNOR", "SNORA", "SNORD", "SCARNA", "RNU")


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Rebuild counts from a scanpy project's selected inputs and "
            "validate or repair a processed H5AD."
        )
    )
    parser.add_argument("project_dir", help="Path to analysis_26GSE_V4/scanpy_projects/<GSE>.")
    parser.add_argument(
        "--write",
        action="store_true",
        help="Replace the processed H5AD in place after validation passes.",
    )
    parser.add_argument(
        "--summary-path",
        default="",
        help="Optional markdown path for a repair or dry-run summary.",
    )
    parser.add_argument(
        "--skip-failed-samples",
        action="store_true",
        help=(
            "Skip sample keys listed in outputs/failed_samples.csv. "
            "Use this when a processed H5AD already excludes known-bad source files."
        ),
    )
    return parser.parse_args()


def load_config(config_path: Path) -> dict:
    """Load project config JSON."""
    return json.loads(config_path.read_text(encoding="utf-8"))


def resolve_path(project_dir: Path, raw_path: str) -> Path:
    """Resolve project-relative scaffold paths."""
    path = Path(raw_path)
    if path.is_absolute():
        return path
    return (project_dir.parents[2] / path).resolve()


def resolve_existing_input_path(project_dir: Path, raw_path: str) -> Path:
    """Resolve an input path and fall back to filename-based search under the GSE supplementary folder.

    Several project manifests point at older extracted directories such as
    `extracted_final_*`, while the actual files now live under
    `extracted_filtered_*` or another sibling folder. This fallback keeps the
    original manifest semantics but avoids hard-failing on stale extraction paths.
    """
    resolved = resolve_path(project_dir, raw_path)
    if resolved.exists():
        return resolved

    gse = project_dir.name
    suppl_dir = project_dir.parents[2] / "downloads" / gse / "suppl"
    if suppl_dir.exists():
        matches = sorted({candidate.resolve() for candidate in suppl_dir.rglob(resolved.name) if candidate.is_file()})
        if len(matches) == 1:
            return matches[0]
        if len(matches) > 1:
            preferred = [match for match in matches if "filtered" in match.parent.name.lower()]
            if len(preferred) == 1:
                return preferred[0]
            return matches[0]

    raise FileNotFoundError(f"Could not resolve existing input for manifest path: {raw_path}")


def find_companion(base: Path, stem: str, kind: str) -> Optional[Path]:
    """Locate features/genes/barcodes files adjacent to an MTX input."""
    candidates: List[Path] = []
    variants = [stem, stem.replace("_matrix", f"_{kind}"), stem.replace(".matrix", f".{kind}")]
    for variant in dict.fromkeys(variants):
        candidates.extend(
            [
                base / f"{variant}_{kind}.tsv.gz",
                base / f"{variant}.{kind}.tsv.gz",
                base / f"{variant}_{kind}.tsv",
                base / f"{variant}.{kind}.tsv",
                base / f"{variant}.tsv.gz",
                base / f"{variant}.tsv",
            ]
        )
    candidates.extend([base / f"{kind}.tsv.gz", base / f"{kind}.tsv"])
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def read_tsv_lines(path: Path) -> List[str]:
    """Read lines from a possibly gzipped TSV-like file."""
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as handle:
        return [line.rstrip("\n") for line in handle]


def load_mtx(path: Path) -> ad.AnnData:
    """Load a Matrix Market input and infer orientation from companion metadata."""
    base = path.parent
    stem = path.name.replace(".matrix.mtx.gz", "").replace(".mtx.gz", "").replace(".mtx", "")
    features_path = find_companion(base, stem, "features") or find_companion(base, stem, "genes")
    barcodes_path = find_companion(base, stem, "barcodes")
    if features_path is None or barcodes_path is None:
        raise FileNotFoundError(f"Missing features/barcodes alongside {path}")

    matrix = io.mmread(str(path))
    matrix = matrix.tocsr() if sp.issparse(matrix) else sp.csr_matrix(matrix)

    features = []
    for line in read_tsv_lines(features_path):
        parts = line.split("\t")
        features.append(parts[1] if len(parts) > 1 else parts[0])
    barcodes = [line.split("\t")[0] for line in read_tsv_lines(barcodes_path)]

    n_row, n_col = matrix.shape
    if n_row == len(features) and n_col == len(barcodes):
        matrix = matrix.T
    elif not (n_row == len(barcodes) and n_col == len(features)):
        raise ValueError(
            f"Dimension mismatch for {path}: matrix={matrix.shape}, "
            f"features={len(features)}, barcodes={len(barcodes)}"
        )

    loaded = ad.AnnData(X=matrix)
    loaded.var_names = pd.Index(features).astype(str)
    loaded.obs_names = pd.Index(barcodes).astype(str)
    loaded.var_names_make_unique()
    return loaded


def load_one(path: Path) -> ad.AnnData:
    """Load one selected input source."""
    lower = str(path).lower()
    if lower.endswith(".h5ad"):
        loaded = sc.read_h5ad(path)
    elif lower.endswith(".h5"):
        loaded = sc.read_10x_h5(path)
    elif lower.endswith(".mtx") or lower.endswith(".mtx.gz"):
        loaded = load_mtx(path)
    else:
        raise ValueError(f"Unsupported input format: {path}")
    loaded.var_names_make_unique()
    return loaded


def remove_noncoding(adata: ad.AnnData) -> ad.AnnData:
    """Drop commonly non-coding genes to match the scaffold preprocessing."""
    keep = np.ones(adata.n_vars, dtype=bool)
    for i, gene in enumerate(pd.Index(adata.var_names).astype(str)):
        if any(pattern in gene.upper() for pattern in NONCODING_PATTERNS):
            keep[i] = False
    return adata[:, keep].copy()


def sanitize_dataframe(frame: pd.DataFrame) -> pd.DataFrame:
    """Convert non-numeric columns into HDF5-safe plain strings."""
    clean = frame.copy()
    for column in clean.columns:
        series = clean[column]
        if pd.api.types.is_numeric_dtype(series.dtype) or pd.api.types.is_bool_dtype(series.dtype):
            continue
        as_object = series.astype(object)
        clean[column] = pd.Series(
            ["" if pd.isna(value) else str(value) for value in as_object],
            index=series.index,
            dtype=object,
        )
    return clean


def canonical_barcode(value: object) -> str:
    """Normalize a barcode-like token while removing concat-added numeric tails."""
    text = str(value).strip()
    if text == "" or text.lower() == "nan":
        return ""
    parts = text.split("_")
    if len(parts) > 1:
        text = parts[0]
    if text.count("-") >= 2:
        prefix, tail = text.rsplit("-", 1)
        if tail.isdigit():
            text = prefix
    return text


def build_obs_alignment_key(obs: pd.DataFrame, obs_names: pd.Index) -> pd.Index:
    """Build a stable per-cell alignment key using sample plus barcode-like information.

    This is used when current processed H5AD obs names differ from the rebuilt
    count-space names but the underlying cells are the same.
    """
    sample_series = obs["sample"].astype(str) if "sample" in obs.columns else pd.Series([""] * len(obs), index=obs.index)

    barcode_column = None
    for candidate in ("barcode", "cell_barcode"):
        if candidate in obs.columns:
            barcode_column = candidate
            break

    if barcode_column is not None:
        barcode_series = obs[barcode_column].astype(str).map(canonical_barcode)
    else:
        barcode_series = pd.Series(pd.Index(obs_names).astype(str), index=obs.index).map(canonical_barcode)

    keys = barcode_series + "__" + sample_series
    return pd.Index(keys.astype(str))


def load_failed_sample_keys(project_dir: Path) -> set[str]:
    """Load sample keys that should be skipped during rebuild."""
    failed_path = project_dir / "outputs" / "failed_samples.csv"
    if not failed_path.exists():
        return set()
    failed = pd.read_csv(failed_path)
    if "sample" not in failed.columns:
        return set()
    return set(failed["sample"].astype(str))


def rebuild_counts(project_dir: Path, skip_samples: Optional[set[str]] = None) -> tuple[ad.AnnData, dict]:
    """Rebuild the count-space AnnData from selected inputs and scaffold QC rules."""
    config_path = project_dir / "config" / "config.json"
    config = load_config(config_path)
    manifest_path = resolve_path(project_dir, config["inputs"]["selected_manifest"])
    manifest = pd.read_csv(manifest_path)
    if manifest.empty:
        raise ValueError(f"No selected inputs in {manifest_path}")
    skip_samples = skip_samples or set()

    adatas = []
    per_sample = []
    for row in manifest.itertuples(index=False):
        input_path = resolve_existing_input_path(project_dir, row.path)
        sample = str(getattr(row, "sample_key", input_path.stem))
        if sample in skip_samples:
            continue
        loaded = load_one(input_path)
        loaded.obs["sample"] = sample
        loaded.obs["GSE"] = config["project"]["gse"]
        adatas.append(loaded)
        per_sample.append({"sample": sample, "n_obs": int(loaded.n_obs), "n_vars": int(loaded.n_vars), "path": str(input_path)})

    merged = adatas[0] if len(adatas) == 1 else ad.concat(adatas, join="outer", index_unique="-")
    merged = remove_noncoding(merged)
    merged.var_names_make_unique()

    sc.pp.filter_cells(merged, min_genes=int(config["qc"]["min_genes_per_cell"]))
    sc.pp.filter_genes(merged, min_cells=int(config["qc"]["min_cells_per_gene"]))
    merged.var["mt"] = pd.Index(merged.var_names).str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(merged, qc_vars=["mt"], inplace=True)
    merged = merged[merged.obs["pct_counts_mt"] < float(config["qc"]["max_mito_pct"]), :].copy()

    if not sp.issparse(merged.X):
        merged.X = sp.csr_matrix(merged.X)
    elif not sp.isspmatrix_csr(merged.X):
        merged.X = merged.X.tocsr()

    if merged.X.nnz and np.allclose(merged.X.data, np.round(merged.X.data), atol=1e-6):
        max_value = float(merged.X.data.max())
        min_value = float(merged.X.data.min())
        if min_value >= np.iinfo(np.int32).min and max_value <= np.iinfo(np.int32).max:
            merged.X.data = merged.X.data.astype(np.int32, copy=False)

    metadata = {
        "config_path": str(config_path),
        "manifest_path": str(manifest_path),
        "samples_loaded": len(per_sample),
        "samples_skipped": sorted(skip_samples),
        "per_sample": per_sample,
    }
    return merged, metadata


def compare_axes(current: ad.AnnData, rebuilt: ad.AnnData) -> dict:
    """Return axis comparison details used for dry-run QC."""
    current_obs = pd.Index(current.obs_names.astype(str))
    current_var = pd.Index(current.var_names.astype(str))
    rebuilt_obs = pd.Index(rebuilt.obs_names.astype(str))
    rebuilt_var = pd.Index(rebuilt.var_names.astype(str))
    current_obs_key = build_obs_alignment_key(current.obs, current_obs)
    rebuilt_obs_key = build_obs_alignment_key(rebuilt.obs, rebuilt_obs)

    obs_equal = current_obs.equals(rebuilt_obs)
    var_equal = current_var.equals(rebuilt_var)
    current_vars_subset = current_var.isin(rebuilt_var).all()
    obs_key_equal = current_obs_key.equals(rebuilt_obs_key)
    obs_key_set_equal = current_obs_key.is_unique and rebuilt_obs_key.is_unique and set(current_obs_key) == set(rebuilt_obs_key)
    return {
        "obs_equal": bool(obs_equal),
        "obs_key_equal": bool(obs_key_equal),
        "obs_key_set_equal": bool(obs_key_set_equal),
        "var_equal": bool(var_equal),
        "var_current_subset_of_rebuilt": bool(current_vars_subset),
        "current_n_obs": int(current.n_obs),
        "rebuilt_n_obs": int(rebuilt.n_obs),
        "current_n_vars": int(current.n_vars),
        "rebuilt_n_vars": int(rebuilt.n_vars),
        "obs_only_current": current_obs.difference(rebuilt_obs)[:10].tolist(),
        "obs_only_rebuilt": rebuilt_obs.difference(current_obs)[:10].tolist(),
        "var_only_current": current_var.difference(rebuilt_var)[:10].tolist(),
        "var_only_rebuilt": rebuilt_var.difference(current_var)[:10].tolist(),
    }


def build_repaired(current: ad.AnnData, rebuilt_counts: ad.AnnData) -> ad.AnnData:
    """Create a repaired H5AD using rebuilt counts and current obs/var tables."""
    current_obs = pd.Index(current.obs_names.astype(str))
    current_var = pd.Index(current.var_names.astype(str))
    rebuilt_lookup = rebuilt_counts
    if not current_obs.equals(pd.Index(rebuilt_counts.obs_names.astype(str))):
        current_keys = build_obs_alignment_key(current.obs, current_obs)
        rebuilt_keys = build_obs_alignment_key(rebuilt_counts.obs, pd.Index(rebuilt_counts.obs_names.astype(str)))
        if not current_keys.is_unique or not rebuilt_keys.is_unique:
            raise ValueError("Obs alignment keys are not unique; refusing repair")
        rebuilt_lookup = rebuilt_counts.copy()
        rebuilt_lookup.obs_names = rebuilt_keys
        current_obs = current_keys

    counts = rebuilt_lookup[current_obs, current_var].X
    counts = counts.tocsr() if sp.issparse(counts) else sp.csr_matrix(counts)

    repaired = ad.AnnData(
        X=counts,
        obs=sanitize_dataframe(current.obs),
        var=sanitize_dataframe(current.var),
    )
    kept_uns = {key: current.uns[key] for key in current.uns.keys() if key in PRESERVE_UNS_KEYS}
    removed_uns = [key for key in current.uns.keys() if key not in PRESERVE_UNS_KEYS]
    repaired.uns = kept_uns
    repaired.uns["phase0_raw_rescue"] = {
        "source": "selected_inputs.csv",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "counts_nnz": int(counts.nnz),
        "counts_dtype": str(counts.dtype),
        "removed_obsm_keys": list(current.obsm.keys()),
        "removed_varm_keys": list(current.varm.keys()),
        "removed_uns_keys": removed_uns,
    }
    return repaired


def write_summary(summary_path: Path, project_dir: Path, comparison: dict, mode: str) -> None:
    """Write a concise markdown QC summary."""
    with summary_path.open("w", encoding="utf-8") as handle:
        handle.write("# Selected-Inputs Count Rescue Summary\n\n")
        handle.write(f"- Project: `{project_dir}`\n")
        handle.write(f"- Mode: `{mode}`\n")
        handle.write(f"- Obs axes match exactly: {comparison['obs_equal']}\n")
        handle.write(f"- Obs alignment keys match exactly: {comparison['obs_key_equal']}\n")
        handle.write(f"- Obs alignment key sets match: {comparison['obs_key_set_equal']}\n")
        handle.write(f"- Var axes match exactly: {comparison['var_equal']}\n")
        handle.write(f"- Current vars subset of rebuilt vars: {comparison['var_current_subset_of_rebuilt']}\n")
        handle.write(f"- Current n_obs: {comparison['current_n_obs']}\n")
        handle.write(f"- Rebuilt n_obs: {comparison['rebuilt_n_obs']}\n")
        handle.write(f"- Current n_vars: {comparison['current_n_vars']}\n")
        handle.write(f"- Rebuilt n_vars: {comparison['rebuilt_n_vars']}\n")
        if comparison["obs_only_current"]:
            handle.write(f"- Obs only in current: `{comparison['obs_only_current']}`\n")
        if comparison["obs_only_rebuilt"]:
            handle.write(f"- Obs only in rebuilt: `{comparison['obs_only_rebuilt']}`\n")
        if comparison["var_only_current"]:
            handle.write(f"- Vars only in current: `{comparison['var_only_current']}`\n")
        if comparison["var_only_rebuilt"]:
            handle.write(f"- Vars only in rebuilt: `{comparison['var_only_rebuilt']}`\n")


def main() -> None:
    """Run dry-run validation or in-place repair."""
    args = parse_args()
    project_dir = Path(args.project_dir).resolve()
    gse = project_dir.name
    current_h5ad = project_dir / "outputs" / f"{gse}_scanpy_processed.h5ad"
    if not current_h5ad.exists():
        current_h5ad = project_dir / "outputs" / "scanpy_processed.h5ad"
    if not current_h5ad.exists():
        raise FileNotFoundError(f"Processed H5AD not found under {project_dir / 'outputs'}")

    skip_samples = load_failed_sample_keys(project_dir) if args.skip_failed_samples else set()
    rebuilt, _ = rebuild_counts(project_dir, skip_samples=skip_samples)
    current = ad.read_h5ad(current_h5ad)
    comparison = compare_axes(current, rebuilt)

    if args.summary_path:
        summary_path = Path(args.summary_path).resolve()
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        write_summary(summary_path, project_dir, comparison, "write" if args.write else "dry-run")

    print(json.dumps(comparison, indent=2))

    if not args.write:
        return

    if not (comparison["obs_equal"] or comparison["obs_key_set_equal"]) or not comparison["var_current_subset_of_rebuilt"]:
        raise ValueError("Axis validation failed; refusing to repair in write mode")

    repaired = build_repaired(current, rebuilt)
    temp_path = current_h5ad.with_suffix(current_h5ad.suffix + ".tmp")
    if temp_path.exists():
        temp_path.unlink()
    repaired.write_h5ad(temp_path)
    temp_path.replace(current_h5ad)
    print(f"Repaired {current_h5ad}")


if __name__ == "__main__":
    main()
