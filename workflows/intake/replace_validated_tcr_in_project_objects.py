#!/usr/bin/env python3
"""Propagate approved rebuilt TCR columns into metadata and project milestones.

By default this helper propagates GSEs that already passed the per-GSE
TRD/TRAB validation rule in `tcr_rebuild_phase4_summary.csv`. It can also
include explicitly user-approved repaired GSEs via `--include-gse`, even if
they were unscoreable for that original Phase 4 gate.
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


import logging
import shutil
import argparse
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from anndata._io.specs import write_elem


PROJECT_ROOT = Path(__file__).resolve().parents[2]
REGISTRY_CSV = PROJECT_ROOT / "configs" / "datasets" / "integration_inputs.csv"
SUMMARY_CSV = PROJECT_ROOT / "Integrated_dataset" / "tables" / "tcr_rebuild_phase4" / "tcr_rebuild_phase4_summary.csv"
REPLACEMENT_TABLE = PROJECT_ROOT / "Integrated_dataset" / "tables" / "tcr_rebuild_phase4" / "validated_tcr_replacement_table.csv.gz"
REPLACEMENT_LOG = PROJECT_ROOT / "Integrated_dataset" / "logs" / "validated_tcr_replace.log"
REPLACEMENT_SUMMARY = PROJECT_ROOT / "Integrated_dataset" / "tables" / "tcr_rebuild_phase4" / "validated_tcr_replace_summary.csv"
HARMONIZED_MAIN = PROJECT_ROOT / "analysis_26GSE_V4" / "outputs" / "harmonized_metadata_v4.csv"
HARMONIZED_BACKUP = HARMONIZED_MAIN.with_name(HARMONIZED_MAIN.name + ".bak_before_validated_tcr_replace")

MILESTONE_PATHS = {
    "TNK_candidates": PROJECT_ROOT / "Integrated_dataset" / "TNK_candidates.h5ad",
    "TNK_cleaned": PROJECT_ROOT / "Integrated_dataset" / "TNK_cleaned.h5ad",
    "integrated": PROJECT_ROOT / "high_speed_temp" / "Integrated_dataset" / "integrated.h5ad",
}

STRING_TCR_COLUMNS = [
    "TRA_cdr3",
    "TRA_v",
    "TRA_d",
    "TRA_j",
    "TRA_cdr3_nt",
    "TRA_clone_id",
    "TRB_cdr3",
    "TRB_v",
    "TRB_d",
    "TRB_j",
    "TRB_cdr3_nt",
    "TRB_clone_id",
    "TCRseq",
]
NUMERIC_TCR_COLUMNS = ["TRA_umis", "TRA_reads", "TRB_umis", "TRB_reads"]
TCR_COLUMNS = [*STRING_TCR_COLUMNS[:-1], *NUMERIC_TCR_COLUMNS, "TCRseq"]
METADATA_USECOLS = ["obs_name", "gse_id", *TCR_COLUMNS]


def configure_logging() -> None:
    """Configure console and file logging."""
    handlers = [
        logging.FileHandler(REPLACEMENT_LOG, mode="a", encoding="utf-8"),
        logging.StreamHandler(),
    ]
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=handlers,
        force=True,
    )


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--include-gse",
        nargs="*",
        default=[],
        help="Additional repaired GSEs to propagate even if they are unscoreable for the Phase 4 validation gate.",
    )
    return parser.parse_args()


def clean_text(value: object) -> str:
    """Normalize missing-like values to an empty string."""
    if pd.isna(value):
        return ""
    text = str(value).strip()
    if text.lower() in {"", "nan", "na", "none", "null", "<na>"}:
        return ""
    return text


def load_registry_map() -> dict[str, Path]:
    """Load per-GSE H5AD paths from the registry."""
    registry = pd.read_csv(REGISTRY_CSV, dtype=str)
    return {
        clean_text(row.gse_id): Path(clean_text(row.h5ad_path))
        for row in registry.itertuples(index=False)
        if clean_text(row.gse_id)
    }


def load_validated_gses(extra_gses: list[str] | None = None) -> list[str]:
    """Return the validated GSEs plus any explicitly user-approved extras."""
    summary = pd.read_csv(SUMMARY_CSV)
    validated = summary.loc[summary["validation_pass"].fillna(False), "gse_id"].astype(str).tolist()
    extras = [clean_text(gse) for gse in (extra_gses or []) if clean_text(gse)]
    combined = list(dict.fromkeys([*validated, *extras]))
    if not combined:
        raise RuntimeError("No validated GSEs were available for replacement.")
    return combined


def build_replacement_table(validated_gses: list[str], registry_map: dict[str, Path]) -> pd.DataFrame:
    """Build one validated-only TCR replacement table from per-GSE H5ADs."""
    frames: list[pd.DataFrame] = []
    for gse_id in validated_gses:
        path = registry_map[gse_id]
        logging.info("Collecting validated TCR rows from %s", path)
        adata = ad.read_h5ad(path, backed="r")
        obs = adata.obs.copy()
        for column in STRING_TCR_COLUMNS:
            if column not in obs.columns:
                obs[column] = "no" if column == "TCRseq" else ""
        for column in NUMERIC_TCR_COLUMNS:
            if column not in obs.columns:
                obs[column] = 0
        obs = obs[TCR_COLUMNS].copy()
        obs["obs_name"] = obs.index.astype(str)
        obs["gse_id"] = gse_id
        obs["_index"] = gse_id + "__" + obs["obs_name"].astype(str)
        frames.append(obs[["_index", "obs_name", "gse_id", *TCR_COLUMNS]].reset_index(drop=True))
        adata.file.close()

    replacement = pd.concat(frames, ignore_index=True)
    for column in STRING_TCR_COLUMNS:
        replacement[column] = replacement[column].map(clean_text).astype(object)
    for column in NUMERIC_TCR_COLUMNS:
        replacement[column] = pd.to_numeric(replacement[column], errors="coerce").fillna(0).astype(np.int32)
    replacement.to_csv(REPLACEMENT_TABLE, index=False, compression="gzip")
    return replacement


def update_harmonized_metadata(replacement: pd.DataFrame) -> dict[str, object]:
    """Replace TCR columns in the main harmonized metadata by obs_name."""
    logging.info("Updating harmonized metadata at %s", HARMONIZED_MAIN)
    metadata = pd.read_csv(HARMONIZED_MAIN, dtype=str, low_memory=False)
    if not HARMONIZED_BACKUP.exists():
        shutil.copy2(HARMONIZED_MAIN, HARMONIZED_BACKUP)

    replacement_meta = replacement[["obs_name", "gse_id", *TCR_COLUMNS]].copy()
    replacement_meta = replacement_meta.rename(columns={column: f"{column}__new" for column in TCR_COLUMNS})
    merged = metadata.merge(replacement_meta, on=["obs_name", "gse_id"], how="left")

    updated_rows = 0
    for column in TCR_COLUMNS:
        new_col = f"{column}__new"
        if column in NUMERIC_TCR_COLUMNS:
            base = pd.to_numeric(merged.get(column, 0), errors="coerce").fillna(0).astype(np.int32)
            incoming = pd.to_numeric(merged[new_col], errors="coerce")
            mask = incoming.notna()
            base.loc[mask] = incoming.loc[mask].astype(np.int32)
            merged[column] = base.astype(str)
        else:
            base = merged.get(column, "").map(clean_text)
            incoming = merged[new_col].map(clean_text)
            mask = incoming != ""
            base.loc[mask] = incoming.loc[mask]
            if column == "TCRseq":
                # Explicitly write "no" into validated rows that remain unmatched/blank.
                gse_mask = merged["gse_id"].isin(replacement["gse_id"].unique())
                empty_mask = gse_mask & base.eq("")
                base.loc[empty_mask] = "no"
            merged[column] = base
        updated_rows += int((merged[new_col].map(clean_text) != "").sum()) if column not in NUMERIC_TCR_COLUMNS else int(pd.to_numeric(merged[new_col], errors="coerce").notna().sum())
        merged = merged.drop(columns=[new_col])

    merged.to_csv(HARMONIZED_MAIN, index=False)
    return {
        "target": "harmonized_metadata_v4.csv",
        "rows_total": int(len(metadata)),
        "validated_rows": int(metadata["gse_id"].isin(replacement["gse_id"].unique()).sum()),
        "replacement_rows": int(len(replacement)),
        "column_updates": updated_rows,
    }


def read_string_dataset(obj) -> np.ndarray:
    """Read one H5AD string-like obs element into a Python string array."""
    if isinstance(obj, h5py.Dataset):
        data = obj[:]
        if h5py.check_string_dtype(obj.dtype) is not None:
            return obj.asstr()[:]
        if data.dtype.kind in {"S", "O"}:
            return pd.Series(data).map(lambda x: x.decode("utf-8") if isinstance(x, bytes) else clean_text(x)).to_numpy(dtype=object)
        return data.astype(str)

    if isinstance(obj, h5py.Group) and "codes" in obj and "categories" in obj:
        codes = obj["codes"][:]
        categories = read_string_dataset(obj["categories"])
        out = np.empty(len(codes), dtype=object)
        neg_mask = codes < 0
        out[neg_mask] = ""
        pos_mask = ~neg_mask
        out[pos_mask] = categories[codes[pos_mask]]
        return out

    raise TypeError(f"Unsupported obs element type: {type(obj)}")


def write_obs_column(obs_group: h5py.Group, column_name: str, values: np.ndarray) -> None:
    """Write or replace one obs column in place."""
    column_order = list(obs_group.attrs["column-order"])
    tmp_name = f"__{column_name}_tmp"
    for existing in [tmp_name, column_name]:
        if existing in obs_group:
            del obs_group[existing]

    if values.dtype.kind in {"U", "O"}:
        payload = pd.Series(values).map(clean_text).to_numpy(dtype=object)
    elif values.dtype.kind in {"i", "u"}:
        payload = np.asarray(values, dtype=np.int32)
    else:
        payload = np.asarray(values)

    write_elem(obs_group, tmp_name, payload)
    obs_group.move(tmp_name, column_name)
    if column_name not in column_order:
        column_order.append(column_name)
    obs_group.attrs["column-order"] = np.asarray(column_order, dtype=object)


def update_milestone_h5ad(label: str, path: Path, replacement: pd.DataFrame) -> dict[str, object]:
    """Write validated TCR columns into one milestone H5AD by `_index`."""
    logging.info("Updating %s at %s", label, path)
    replacement_indexed = replacement.set_index("_index", drop=False)
    replacement_gses = replacement["gse_id"].unique()

    with h5py.File(path, "r+") as handle:
        obs_group = handle["obs"]
        cell_index = read_string_dataset(obs_group["_index"])
        source_gse = read_string_dataset(obs_group["source_gse_id"]) if "source_gse_id" in obs_group else np.array([""] * len(cell_index), dtype=object)
        target_mask = np.isin(source_gse, replacement_gses)
        matched_index = pd.Index(cell_index[target_mask])

        stats = {
            "target": label,
            "rows_total": int(len(cell_index)),
            "validated_rows_present": int(target_mask.sum()),
            "replacement_rows_written": 0,
        }

        for column in TCR_COLUMNS:
            if column in STRING_TCR_COLUMNS:
                default_value = "no" if column == "TCRseq" else ""
                full_values = np.full(len(cell_index), default_value, dtype=object)
                mapped = matched_index.map(replacement_indexed[column]).to_numpy(dtype=object)
                mapped = pd.Series(mapped).map(clean_text).to_numpy(dtype=object)
            else:
                full_values = np.zeros(len(cell_index), dtype=np.int32)
                mapped = pd.to_numeric(matched_index.map(replacement_indexed[column]), errors="coerce").fillna(0).to_numpy(dtype=np.int32)

            full_values[target_mask] = mapped
            write_obs_column(obs_group, column, full_values)
            if column == "TCRseq":
                stats["replacement_rows_written"] = int(np.count_nonzero(pd.Series(mapped).map(clean_text).eq("yes")))

    return stats


def write_summary(rows: list[dict[str, object]]) -> None:
    """Write the replacement summary table."""
    pd.DataFrame(rows).to_csv(REPLACEMENT_SUMMARY, index=False)


def main() -> None:
    """Run validated-plus-user-approved TCR propagation."""
    args = parse_args()
    configure_logging()
    registry_map = load_registry_map()
    validated_gses = load_validated_gses(args.include_gse)
    logging.info("Validated/approved GSEs: %s", ", ".join(validated_gses))

    replacement = build_replacement_table(validated_gses, registry_map)
    summary_rows = [update_harmonized_metadata(replacement)]

    for label, path in MILESTONE_PATHS.items():
        summary_rows.append(update_milestone_h5ad(label, path, replacement))

    write_summary(summary_rows)
    logging.info("Wrote replacement table to %s", REPLACEMENT_TABLE)
    logging.info("Wrote replacement summary to %s", REPLACEMENT_SUMMARY)


if __name__ == "__main__":
    main()
