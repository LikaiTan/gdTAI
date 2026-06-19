#!/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python
"""Promote gdTAI v3 and build the final predicted gdT cell-atlas H5AD.

This is a handoff-preparation step only. It does not run downstream gdT-only
integration or phenotype analysis.
"""

from __future__ import annotations

import argparse
import json
import logging
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import anndata as ad
import h5py
import numpy as np
import pandas as pd

from run_gdt_gse144469_holdout_tcrgene_classifier import build_obs_metadata


PROJECT_ROOT = Path(__file__).resolve().parent
DEFAULT_INPUT_H5AD = PROJECT_ROOT / "high_speed_temp/Integrated_dataset/integrated_plus6.h5ad"
V3_CANDIDATE_DIR = PROJECT_ROOT / "Integrated_dataset/models/gdT_prediction_classifier/gdtai_v3_trdc_nk_guard"
PROMOTED_DIR = PROJECT_ROOT / "Integrated_dataset/models/gdT_prediction_classifier/gdTAI_v3.0"
PREDICTION_TABLE = (
    PROJECT_ROOT
    / "Integrated_dataset/tables/gdT_prediction/gdtai_v3_trdc_nk_guard/full_atlas_selected_predicted_cells.csv.gz"
)
TARGET_SELECTION_TABLE = (
    PROJECT_ROOT / "Integrated_dataset/tables/gdT_prediction/gdtai_v3_trdc_nk_guard/full_atlas_target_selection.csv"
)
EXTERNAL_METRICS_TABLE = (
    PROJECT_ROOT / "Integrated_dataset/tables/gdT_prediction/gdtai_v3_trdc_nk_guard/external_primary_metrics.csv"
)
EXTERNAL_FP_TABLE = (
    PROJECT_ROOT / "Integrated_dataset/tables/gdT_prediction/gdtai_v3_trdc_nk_guard/external_false_positive_groups.csv"
)

ATLAS_DIR = PROJECT_ROOT / "Integrated_dataset/gdT_atlas"
TABLE_DIR = PROJECT_ROOT / "Integrated_dataset/tables/gdT_atlas"
LOG_DIR = PROJECT_ROOT / "Integrated_dataset/logs/gdT_atlas"
OUTPUT_H5AD = ATLAS_DIR / "predicted_gdt_cell_atlas.h5ad"
HANDOFF_MD = ATLAS_DIR / "predicted_gdt_cell_atlas.md"
SUMMARY_CSV = TABLE_DIR / "predicted_gdt_cell_atlas_summary.csv"
METADATA_CSV = TABLE_DIR / "predicted_gdt_cell_atlas_metadata.csv.gz"
REMOVED_FP_CSV = TABLE_DIR / "predicted_gdt_cell_atlas_removed_fp_cells.csv.gz"
GOLD_FN_CSV = TABLE_DIR / "predicted_gdt_cell_atlas_gold_fn_added_cells.csv.gz"
BY_SOURCE_CSV = TABLE_DIR / "predicted_gdt_cell_atlas_by_source.csv"
BY_TISSUE_CSV = TABLE_DIR / "predicted_gdt_cell_atlas_by_tissue.csv"
RUN_LOG = LOG_DIR / "promote_gdtai_v3_and_build_predicted_gdt_atlas.log"
MANIFEST_JSON = LOG_DIR / "predicted_gdt_cell_atlas_manifest.json"

PROMOTED_MODEL_NAME = "gdTAI_v3_model.pkl"
SELECTED_MODEL = "v3_round14_v2_score_trdc_gate_fixed_0p936"
SELECTED_THRESHOLD = 0.936

PREVIOUS_ATLAS_PATHS = [
    PROJECT_ROOT / "Integrated_dataset/gdT_atlas",
    PROJECT_ROOT / "Integrated_dataset/tables/gdT_atlas",
    PROJECT_ROOT / "Integrated_dataset/tables/gdT_atlas_curated_phenotypes",
    PROJECT_ROOT / "Integrated_dataset/tables/gdT_atlas_rigorous",
    PROJECT_ROOT / "Integrated_dataset/figures/gdT_atlas",
    PROJECT_ROOT / "Integrated_dataset/figures/gdT_atlas_curated_phenotypes",
    PROJECT_ROOT / "Integrated_dataset/figures/gdT_atlas_rigorous",
    PROJECT_ROOT / "Integrated_dataset/logs/gdT_atlas",
    PROJECT_ROOT / "Integrated_dataset/logs/gdT_atlas_curated_phenotypes",
    PROJECT_ROOT / "Integrated_dataset/logs/gdT_atlas_rigorous",
    PROJECT_ROOT / "Integrated_dataset/models/gdT_atlas_rigorous",
    PROJECT_ROOT / "gdT_atlas",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Promote gdTAI v3 and build predicted gdT cell-atlas H5AD.")
    parser.add_argument("--input-h5ad", type=Path, default=DEFAULT_INPUT_H5AD)
    parser.add_argument("--prediction-table", type=Path, default=PREDICTION_TABLE)
    parser.add_argument("--candidate-dir", type=Path, default=V3_CANDIDATE_DIR)
    parser.add_argument("--promoted-dir", type=Path, default=PROMOTED_DIR)
    parser.add_argument("--output-h5ad", type=Path, default=OUTPUT_H5AD)
    parser.add_argument("--skip-remove-previous", action="store_true")
    parser.add_argument("--skip-promote-model", action="store_true")
    parser.add_argument("--compression", default="gzip", choices=["gzip", "lzf", "none"])
    return parser.parse_args()


def setup_logging() -> None:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.FileHandler(RUN_LOG, mode="w", encoding="utf-8"), logging.StreamHandler()],
        force=True,
    )


def remove_previous_outputs() -> list[dict[str, Any]]:
    removed: list[dict[str, Any]] = []
    for path in PREVIOUS_ATLAS_PATHS:
        if not path.exists():
            continue
        kind = "directory" if path.is_dir() else "file"
        size_bytes = directory_size(path) if path.is_dir() else path.stat().st_size
        logging.info("Removing previous gdT atlas %s: %s", kind, path)
        if path.is_dir():
            shutil.rmtree(path)
        else:
            path.unlink()
        removed.append(
            {
                "path": str(path.relative_to(PROJECT_ROOT)),
                "kind": kind,
                "size_bytes": int(size_bytes),
            }
        )
    return removed


def directory_size(path: Path) -> int:
    total = 0
    for child in path.rglob("*"):
        if child.is_file() and not child.is_symlink():
            total += child.stat().st_size
    return total


def ensure_output_dirs() -> None:
    for path in [ATLAS_DIR, TABLE_DIR, LOG_DIR]:
        path.mkdir(parents=True, exist_ok=True)


def promote_model(candidate_dir: Path, promoted_dir: Path) -> dict[str, Any]:
    if not candidate_dir.exists():
        raise FileNotFoundError(candidate_dir)
    promoted_dir.mkdir(parents=True, exist_ok=True)
    copies = {
        "best_candidate_model.pkl": PROMOTED_MODEL_NAME,
        "README.md": "README.md",
        "METHODOLOGY.md": "METHODOLOGY.md",
        "USAGE_PROTOCOL.md": "USAGE_PROTOCOL.md",
        "feature_genes.csv": "feature_genes.csv",
        "mode_metrics.csv": "mode_metrics.csv",
        "candidate_target_selection.csv": "candidate_target_selection.csv",
        "external_false_positive_groups.csv": "external_false_positive_groups.csv",
        "external_recall_groups.csv": "external_recall_groups.csv",
    }
    copied: list[str] = []
    for source_name, target_name in copies.items():
        source = candidate_dir / source_name
        if not source.exists():
            logging.warning("Candidate package file missing, not copied: %s", source)
            continue
        target = promoted_dir / target_name
        shutil.copy2(source, target)
        copied.append(str(target.relative_to(PROJECT_ROOT)))
    examples_src = candidate_dir / "examples"
    examples_dst = promoted_dir / "examples"
    if examples_src.exists():
        if examples_dst.exists():
            shutil.rmtree(examples_dst)
        shutil.copytree(examples_src, examples_dst)
        copied.append(str(examples_dst.relative_to(PROJECT_ROOT)))

    manifest_source = candidate_dir / "model_manifest.json"
    manifest: dict[str, Any] = {}
    if manifest_source.exists():
        manifest = json.loads(manifest_source.read_text(encoding="utf-8"))
    manifest.update(
        {
            "version": "gdTAI_v3.0",
            "model_file": PROMOTED_MODEL_NAME,
            "accepted_for_promotion": True,
            "promotion_decision": "user_approved",
            "promotion_date_local": datetime.now().astimezone().strftime("%Y-%m-%d %Z"),
            "promotion_date_utc": datetime.now(timezone.utc).strftime("%Y-%m-%d"),
            "automated_gate_result_before_user_override": False,
            "promotion_note": (
                "Promoted by user decision after review. The prior automated gate failed only because "
                "external NK false positives were higher than v2 high-purity; this caveat remains documented."
            ),
        }
    )
    (promoted_dir / "model_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    (promoted_dir / "PROMOTION_DECISION.md").write_text(promotion_markdown(manifest), encoding="utf-8")
    copied.extend(
        [
            str((promoted_dir / "model_manifest.json").relative_to(PROJECT_ROOT)),
            str((promoted_dir / "PROMOTION_DECISION.md").relative_to(PROJECT_ROOT)),
        ]
    )
    return {"promoted_dir": str(promoted_dir.relative_to(PROJECT_ROOT)), "copied": copied, "manifest": manifest}


def promotion_markdown(manifest: dict[str, Any]) -> str:
    full = manifest.get("full_atlas", {})
    external = manifest.get("external_primary", {})
    return "\n".join(
        [
            "# gdTAI v3.0 Promotion Decision",
            "",
            "This directory contains the user-promoted gdTAI v3.0 model.",
            "",
            f"- Promoted model: `{manifest.get('model_name', SELECTED_MODEL)}`",
            f"- Model file: `{manifest.get('model_file', PROMOTED_MODEL_NAME)}`",
            "- Promotion decision: user-approved override after review",
            f"- Full-atlas predicted cells before curation: `{int(full.get('predicted_putative_gdT', 0)):,}`",
            f"- Full-atlas primary-gold recall before curation: `{float(full.get('primary_gold_recall', 0.0)):.6f}`",
            f"- Estimated FP fraction before curation: `{float(full.get('estimated_fp_fraction_of_predictions', 0.0)):.6f}`",
            f"- External precision/recall/F1: `{float(external.get('precision', 0.0)):.4f}` / `{float(external.get('recall', 0.0)):.4f}` / `{float(external.get('f1', 0.0)):.4f}`",
            "",
            "Caveat retained from review: external NK false positives were higher than v2 high-purity, but the user decided to promote this model.",
            "",
        ]
    )


def class_labels(class_code: np.ndarray) -> np.ndarray:
    labels = np.asarray(["unlabeled_or_ambiguous", "abT_gold", "gdT_gold", "gdT_silver"], dtype=object)
    clipped = np.clip(class_code.astype(int), 0, len(labels) - 1)
    return labels[clipped]


def build_curation_masks(input_h5ad: Path, prediction_table: Path) -> tuple[pd.DataFrame, dict[str, Any]]:
    predicted = pd.read_csv(prediction_table)
    if "obs_index" not in predicted:
        raise KeyError(f"`obs_index` missing from {prediction_table}")
    if predicted["obs_index"].duplicated().any():
        raise RuntimeError("Prediction table contains duplicated obs_index values.")

    with h5py.File(input_h5ad, "r") as handle:
        obs = build_obs_metadata(handle)
        n_obs = int(handle["obs"]["_index"].shape[0])
        obs_names = pd.Index(read_h5_string_dataset(handle["obs"]["_index"]), dtype="string").astype(str).to_numpy(dtype=object)

    predicted_rows = predicted["obs_index"].to_numpy(dtype=np.int64)
    if predicted_rows.min(initial=0) < 0 or predicted_rows.max(initial=-1) >= n_obs:
        raise RuntimeError("Prediction table obs_index is out of bounds for input H5AD.")

    pred_mask = np.zeros(n_obs, dtype=bool)
    pred_mask[predicted_rows] = True
    primary_gdt = obs.class_code == 2
    primary_abt = obs.class_code == 1
    paired_tcrab_no_gdt = (obs.has_TRA_TRB_paired | obs.has_any_ab_tcr) & (~primary_gdt)
    known_fp = pred_mask & (primary_abt | paired_tcrab_no_gdt)
    gold_fn_added = primary_gdt & (~pred_mask)
    final_mask = (pred_mask & (~known_fp)) | gold_fn_added

    score = np.full(n_obs, np.nan, dtype=np.float32)
    threshold = np.full(n_obs, np.nan, dtype=np.float32)
    quadrant = np.full(n_obs, "", dtype=object)
    pred_lookup = predicted.set_index("obs_index", drop=False)
    score[predicted_rows] = pred_lookup.loc[predicted_rows, "score"].to_numpy(dtype=np.float32)
    threshold[predicted_rows] = pred_lookup.loc[predicted_rows, "threshold"].to_numpy(dtype=np.float32)
    if "tcr_gene_quadrant" in predicted:
        quadrant[predicted_rows] = pred_lookup.loc[predicted_rows, "tcr_gene_quadrant"].astype(str).to_numpy(dtype=object)

    reason = np.full(n_obs, "not_in_final_atlas", dtype=object)
    reason[pred_mask & (~known_fp)] = "gdtai_v3_predicted_retained"
    reason[gold_fn_added] = "gold_gdT_FN_added"
    reason[known_fp] = "known_or_likely_FP_removed"

    curation = pd.DataFrame(
        {
            "obs_index": np.arange(n_obs, dtype=np.int64),
            "cell_id": obs_names,
            "source_gse_id": obs.source,
            "tissue": obs.tissue,
            "gdtai_v3_score": score,
            "gdtai_v3_threshold": threshold,
            "gdtai_v3_predicted": pred_mask,
            "gdtai_v3_tcr_gene_quadrant": quadrant,
            "gdtai_v3_primary_class": class_labels(obs.class_code),
            "gdtai_v3_class_code": obs.class_code.astype(np.int8),
            "gdtai_v3_gold_fn_added": gold_fn_added,
            "gdtai_v3_known_fp_removed": known_fp,
            "gdtai_v3_fp_primary_abt_gold": pred_mask & primary_abt,
            "gdtai_v3_fp_paired_TCRAB_no_gdT": pred_mask & paired_tcrab_no_gdt,
            "gdtai_v3_has_TRA_TRB_paired": obs.has_TRA_TRB_paired,
            "gdtai_v3_has_any_ab_tcr": obs.has_any_ab_tcr,
            "gdtai_v3_corrected_has_any_gd_tcr": obs.corrected_has_any_gd_tcr,
            "gdtai_v3_in_final_predicted_gdt_atlas": final_mask,
            "gdtai_v3_atlas_inclusion_reason": reason,
        }
    )
    metrics = {
        "n_input_cells": int(n_obs),
        "n_predicted_positive_raw": int(pred_mask.sum()),
        "n_known_or_likely_fp_removed": int(known_fp.sum()),
        "n_fp_primary_abt_gold_removed": int((pred_mask & primary_abt).sum()),
        "n_fp_paired_TCRAB_no_gdT_removed": int((pred_mask & paired_tcrab_no_gdt).sum()),
        "n_gold_gdT_fn_added": int(gold_fn_added.sum()),
        "n_final_atlas_cells": int(final_mask.sum()),
        "n_final_predicted_retained": int((pred_mask & (~known_fp)).sum()),
        "n_final_gold_fn_added": int(gold_fn_added.sum()),
        "n_final_primary_gdt_gold": int((final_mask & primary_gdt).sum()),
        "n_final_gdt_silver": int((final_mask & (obs.class_code == 3)).sum()),
        "n_final_unlabeled_or_ambiguous": int((final_mask & (obs.class_code == 0)).sum()),
        "n_final_primary_abt_gold": int((final_mask & primary_abt).sum()),
    }
    return curation, metrics


def read_h5_string_dataset(dataset: h5py.Dataset) -> np.ndarray:
    values = dataset[:]
    if values.dtype.kind in {"S", "O"}:
        return pd.Series(values, copy=False).map(lambda x: x.decode("utf-8") if isinstance(x, bytes) else str(x)).to_numpy(dtype=object)
    return values.astype(str)


def write_subset_h5ad(input_h5ad: Path, output_h5ad: Path, curation: pd.DataFrame, compression: str) -> None:
    final = curation["gdtai_v3_in_final_predicted_gdt_atlas"].to_numpy(dtype=bool)
    rows = curation.loc[final, "obs_index"].to_numpy(dtype=np.int64)
    logging.info("Reading final subset into memory: %s cells", f"{rows.size:,}")
    backed = ad.read_h5ad(input_h5ad, backed="r")
    try:
        subset = backed[rows, :].to_memory()
    finally:
        backed.file.close()

    curation_final = curation.loc[final].set_index("cell_id", drop=False)
    for column in [
        "obs_index",
        "gdtai_v3_score",
        "gdtai_v3_threshold",
        "gdtai_v3_predicted",
        "gdtai_v3_tcr_gene_quadrant",
        "gdtai_v3_primary_class",
        "gdtai_v3_class_code",
        "gdtai_v3_gold_fn_added",
        "gdtai_v3_known_fp_removed",
        "gdtai_v3_fp_primary_abt_gold",
        "gdtai_v3_fp_paired_TCRAB_no_gdT",
        "gdtai_v3_has_TRA_TRB_paired",
        "gdtai_v3_has_any_ab_tcr",
        "gdtai_v3_corrected_has_any_gd_tcr",
        "gdtai_v3_in_final_predicted_gdt_atlas",
        "gdtai_v3_atlas_inclusion_reason",
    ]:
        subset.obs[column] = curation_final.loc[subset.obs_names, column].to_numpy()

    subset.uns["predicted_gdt_cell_atlas"] = {
        "description": "Predicted gdT cell atlas from user-promoted gdTAI v3.0.",
        "source_h5ad": str(input_h5ad),
        "model_name": SELECTED_MODEL,
        "model_threshold": float(SELECTED_THRESHOLD),
        "selection_rule": "gdTAI v3 predicted positives, minus known/likely FP, plus primary-gold gdT FN add-back.",
        "created_utc": datetime.now(timezone.utc).isoformat(),
    }
    tmp = output_h5ad.with_suffix(output_h5ad.suffix + ".tmp")
    if tmp.exists():
        tmp.unlink()
    logging.info("Writing final H5AD: %s", output_h5ad)
    write_kwargs = {}
    if compression != "none":
        write_kwargs["compression"] = compression
    subset.write_h5ad(tmp, **write_kwargs)
    tmp.replace(output_h5ad)


def write_tables_and_handoff(
    curation: pd.DataFrame,
    metrics: dict[str, Any],
    removed_previous: list[dict[str, Any]],
    promoted: dict[str, Any] | None,
    input_h5ad: Path,
    output_h5ad: Path,
) -> None:
    final = curation["gdtai_v3_in_final_predicted_gdt_atlas"].to_numpy(dtype=bool)
    removed = curation["gdtai_v3_known_fp_removed"].to_numpy(dtype=bool)
    gold_fn = curation["gdtai_v3_gold_fn_added"].to_numpy(dtype=bool)
    summary = pd.DataFrame([metrics])
    summary.to_csv(SUMMARY_CSV, index=False)
    curation.loc[final].to_csv(METADATA_CSV, index=False, compression="gzip")
    curation.loc[removed].to_csv(REMOVED_FP_CSV, index=False, compression="gzip")
    curation.loc[gold_fn].to_csv(GOLD_FN_CSV, index=False, compression="gzip")

    by_source = (
        curation.loc[final]
        .groupby(["source_gse_id", "gdtai_v3_atlas_inclusion_reason"], dropna=False)
        .size()
        .reset_index(name="n_cells")
    )
    by_source.to_csv(BY_SOURCE_CSV, index=False)
    by_tissue = (
        curation.loc[final]
        .groupby(["tissue", "gdtai_v3_atlas_inclusion_reason"], dropna=False)
        .size()
        .reset_index(name="n_cells")
    )
    by_tissue.to_csv(BY_TISSUE_CSV, index=False)

    manifest = {
        "output_h5ad": str(output_h5ad.relative_to(PROJECT_ROOT)),
        "input_h5ad": str(input_h5ad.relative_to(PROJECT_ROOT)) if input_h5ad.is_relative_to(PROJECT_ROOT) else str(input_h5ad),
        "prediction_table": str(PREDICTION_TABLE.relative_to(PROJECT_ROOT)),
        "selected_model": SELECTED_MODEL,
        "selected_threshold": SELECTED_THRESHOLD,
        "metrics": metrics,
        "removed_previous_outputs": removed_previous,
        "promoted_model": promoted,
        "tables": {
            "summary": str(SUMMARY_CSV.relative_to(PROJECT_ROOT)),
            "metadata": str(METADATA_CSV.relative_to(PROJECT_ROOT)),
            "removed_fp": str(REMOVED_FP_CSV.relative_to(PROJECT_ROOT)),
            "gold_fn_added": str(GOLD_FN_CSV.relative_to(PROJECT_ROOT)),
            "by_source": str(BY_SOURCE_CSV.relative_to(PROJECT_ROOT)),
            "by_tissue": str(BY_TISSUE_CSV.relative_to(PROJECT_ROOT)),
        },
        "created_utc": datetime.now(timezone.utc).isoformat(),
    }
    MANIFEST_JSON.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    HANDOFF_MD.write_text(handoff_markdown(metrics, removed_previous, promoted, input_h5ad, output_h5ad), encoding="utf-8")


def handoff_markdown(
    metrics: dict[str, Any],
    removed_previous: list[dict[str, Any]],
    promoted: dict[str, Any] | None,
    input_h5ad: Path,
    output_h5ad: Path,
) -> str:
    promoted_dir = promoted["promoted_dir"] if promoted else str(PROMOTED_DIR.relative_to(PROJECT_ROOT))
    removed_lines = [f"- `{row['path']}` ({row['kind']}, {row['size_bytes']:,} bytes)" for row in removed_previous]
    if not removed_lines:
        removed_lines = ["- None found at run time."]
    return "\n".join(
        [
            "# Predicted gdT Cell Atlas H5AD",
            "",
            f"Output H5AD: `{output_h5ad.relative_to(PROJECT_ROOT)}`",
            "",
            "This H5AD is the handoff object for downstream gdT biology analysis. It is a subset of the full plus6 integrated atlas and contains cells selected by the user-promoted gdTAI v3.0 model, with primary-gold gdT false negatives added back and known/likely false positives removed.",
            "",
            "## Source",
            "",
            f"- Source H5AD: `{input_h5ad.relative_to(PROJECT_ROOT) if input_h5ad.is_relative_to(PROJECT_ROOT) else input_h5ad}`",
            f"- Prediction table: `{PREDICTION_TABLE.relative_to(PROJECT_ROOT)}`",
            f"- Promoted model directory: `{promoted_dir}`",
            f"- Promoted model: `{PROMOTED_MODEL_NAME}`",
            f"- Selected candidate: `{SELECTED_MODEL}`",
            f"- Threshold: `{SELECTED_THRESHOLD}`",
            "",
            "## Inclusion Rule",
            "",
            "1. Start with all full-atlas cells predicted positive by the selected gdTAI v3 candidate.",
            "2. Remove predicted positives that are known/likely false positives:",
            "   - primary abT-gold predicted positives (`class_code == 1`), and",
            "   - paired TCRAB / any abTCR cells that are not primary gdT gold.",
            "3. Add back primary-gold gdT cells that the model missed (`class_code == 2` and not predicted).",
            "",
            "Silver-only gdT-like cells are not added back unless they were already predicted positive.",
            "",
            "## Cell Counts",
            "",
            f"- Full input cells: `{metrics['n_input_cells']:,}`",
            f"- Raw gdTAI v3 predicted positives: `{metrics['n_predicted_positive_raw']:,}`",
            f"- Known/likely FP removed: `{metrics['n_known_or_likely_fp_removed']:,}`",
            f"- Primary abT-gold FP removed: `{metrics['n_fp_primary_abt_gold_removed']:,}`",
            f"- Paired TCRAB/no-gdT FP removed: `{metrics['n_fp_paired_TCRAB_no_gdT_removed']:,}`",
            f"- Primary-gold gdT FN added back: `{metrics['n_gold_gdT_fn_added']:,}`",
            f"- Final atlas cells: `{metrics['n_final_atlas_cells']:,}`",
            f"- Final predicted-retained cells: `{metrics['n_final_predicted_retained']:,}`",
            f"- Final gold-FN-added cells: `{metrics['n_final_gold_fn_added']:,}`",
            "",
            "## Added Obs Columns",
            "",
            "- `gdtai_v3_predicted`: whether the promoted model predicted the cell.",
            "- `gdtai_v3_score`: promoted model probability score; `NaN` for gold FN add-back cells.",
            "- `gdtai_v3_threshold`: operating threshold for predicted cells.",
            "- `gdtai_v3_tcr_gene_quadrant`: TRDC/TRDV expression quadrant when available from prediction output.",
            "- `gdtai_v3_primary_class`: `unlabeled_or_ambiguous`, `abT_gold`, `gdT_gold`, or `gdT_silver`.",
            "- `gdtai_v3_gold_fn_added`: true for primary-gold gdT cells added despite missed prediction.",
            "- `gdtai_v3_atlas_inclusion_reason`: `gdtai_v3_predicted_retained` or `gold_gdT_FN_added`.",
            "",
            "## Audit Tables",
            "",
            f"- Summary: `{SUMMARY_CSV.relative_to(PROJECT_ROOT)}`",
            f"- Final cell metadata: `{METADATA_CSV.relative_to(PROJECT_ROOT)}`",
            f"- Removed FP cells: `{REMOVED_FP_CSV.relative_to(PROJECT_ROOT)}`",
            f"- Gold FN add-back cells: `{GOLD_FN_CSV.relative_to(PROJECT_ROOT)}`",
            f"- By source: `{BY_SOURCE_CSV.relative_to(PROJECT_ROOT)}`",
            f"- By tissue: `{BY_TISSUE_CSV.relative_to(PROJECT_ROOT)}`",
            "",
            "## Previous Atlas Outputs Removed",
            "",
            *removed_lines,
            "",
            "## Downstream Use",
            "",
            "Use this H5AD as the starting object for downstream gdT-only analysis. The object has not been re-integrated or reclustered after subsetting; downstream agents should run their own gdT-only integration, clustering, annotation, and biological analyses from this file.",
            "",
        ]
    )


def main() -> None:
    args = parse_args()
    if not args.input_h5ad.exists():
        raise FileNotFoundError(args.input_h5ad)
    if not args.prediction_table.exists():
        raise FileNotFoundError(args.prediction_table)

    removed_previous = [] if args.skip_remove_previous else remove_previous_outputs()
    ensure_output_dirs()
    setup_logging()
    logging.info("Removed %s previous gdT atlas output paths", len(removed_previous))

    promoted = None if args.skip_promote_model else promote_model(args.candidate_dir, args.promoted_dir)
    curation, metrics = build_curation_masks(args.input_h5ad, args.prediction_table)
    write_subset_h5ad(args.input_h5ad, args.output_h5ad, curation, args.compression)
    write_tables_and_handoff(curation, metrics, removed_previous, promoted, args.input_h5ad, args.output_h5ad)
    logging.info("Final predicted gdT atlas cells: %s", f"{metrics['n_final_atlas_cells']:,}")
    logging.info("Wrote %s", args.output_h5ad)
    logging.info("Wrote %s", HANDOFF_MD)


if __name__ == "__main__":
    main()
