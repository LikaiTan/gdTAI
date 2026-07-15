#!/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python
"""Train and evaluate precision-first gdT classifiers from the plus6 object.

The script is read-only with respect to the H5AD. It reuses the project gdT
gold/silver label construction, trains candidate classifiers on a train split,
selects precision-first thresholds on a tune split, evaluates only on a held-out
gold validation split, and renders a static report.
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
import html
import json
import logging
import math
import pickle
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    confusion_matrix,
    fbeta_score,
    f1_score,
    matthews_corrcoef,
    precision_recall_curve,
    precision_score,
    recall_score,
    roc_auc_score,
    roc_curve,
)
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from xgboost import XGBClassifier

from run_gdt_prediction_package_evaluation import (
    DEFAULT_INPUT_H5AD,
    SORTED_GDT_SOURCES,
    TABLE_DIR,
    FIGURE_DIR,
    LOG_DIR,
    STATIC_ASSET_DIR,
    STATIC_DIR,
    build_corrected_tcr_evidence,
    build_sublibrary_labels,
    clean_group_values,
    dataframe_to_html,
    dataframe_to_markdown,
    make_truth_labels,
    normalize_strings,
    read_bool_obs,
    read_float_obs,
    read_nonempty_string_mask,
    read_obs_column,
    read_string_dataset,
    required_columns_present,
)


PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUTPUT_ROOT = PROJECT_ROOT / "Integrated_dataset"
CLASSIFIER_TABLE_DIR = TABLE_DIR / "classifier_training"
CLASSIFIER_FIGURE_DIR = FIGURE_DIR / "classifier_training"
CLASSIFIER_LOG_DIR = LOG_DIR / "classifier_training"
CLASSIFIER_MODEL_DIR = OUTPUT_ROOT / "models" / "gdT_prediction_classifier"

REPORT_MD = CLASSIFIER_LOG_DIR / "gdt_deg_tcr_classifier_report.md"
RUN_LOG = CLASSIFIER_LOG_DIR / "gdt_deg_tcr_classifier_training.log"
REPORT_HTML = STATIC_DIR / "gdt_deg_tcr_classifier_report.html"
REPORT_PDF = STATIC_DIR / "gdt_deg_tcr_classifier_report.pdf"
SUMMARY_JSON = CLASSIFIER_LOG_DIR / "gdt_deg_tcr_classifier_summary.json"
FEATURE_JSON = CLASSIFIER_LOG_DIR / "gdt_deg_tcr_classifier_features.json"

DEG_FEATURE_TABLE = (
    TABLE_DIR / "trd_gt0p1_high_no_ab_vs_expanded_negative_no_tcr_genes_deg_reproducibility.csv"
)

TARGET_SUM = 10_000.0
RANDOM_SEED = 20260501
MIN_SPECIFICITY = 0.995
F_BETA = 0.5
MAX_NEGATIVE_TRAIN = 250_000
MAX_DEG_FEATURES = 300
MAX_EXPRESSION_GENES = 300
MAX_TCR_EXPRESSION_GENES = 120
BOOTSTRAP_REPS = 500
FULL_APPLY_CHUNK = 50_000
DEATH_PENALTY_GENES = ["FOXP3", "CD4"]
PAN_CD3_GENES = ["CD3D", "CD3E", "CD3G"]
PREVALENCE_SCENARIOS = [0.001, 0.0025, 0.005, 0.01, 0.02, 0.05]

TCR_GENE_PREFIXES = (
    "TRAV",
    "TRAJ",
    "TRAC",
    "TRBV",
    "TRBJ",
    "TRBC",
    "TRGV",
    "TRGJ",
    "TRGC",
    "TRDV",
    "TRDD",
    "TRDJ",
    "TRDC",
)

TCR_FAMILIES = {
    "TRD": ("TRDC", "TRDV", "TRDD", "TRDJ"),
    "TRG": ("TRGC", "TRGV", "TRGJ"),
    "TRA": ("TRAC", "TRAV", "TRAJ"),
    "TRB": ("TRBC", "TRBV", "TRBJ"),
}

CORE_TCR_GENES = [
    "TRDC",
    "TRDV1",
    "TRDV2",
    "TRDV3",
    "TRDJ1",
    "TRDJ2",
    "TRDJ3",
    "TRDJ4",
    "TRDD1",
    "TRDD2",
    "TRDD3",
    "TRGC1",
    "TRGC2",
    "TRGV9",
    "TRGV4",
    "TRGV8",
    "TRGJ1",
    "TRGJ2",
    "TRAC",
    "TRBC1",
    "TRBC2",
]
SUBOPTIMAL_SORTED_GDT_SOURCES = {"GDTlung2023july_7p"}


@dataclass
class ObsMetadata:
    source: np.ndarray
    tissue: np.ndarray
    library_id: np.ndarray
    sample_id: np.ndarray
    has_TRA_TRB_paired: np.ndarray
    has_any_ab_tcr: np.ndarray
    has_TRG_TRD_paired_raw: np.ndarray
    has_any_gd_tcr_raw: np.ndarray
    tra_nonempty: np.ndarray
    trb_nonempty: np.ndarray
    trg_nonempty: np.ndarray
    trd_nonempty: np.ndarray
    corrected_has_any_gd_tcr: np.ndarray
    corrected_has_TRG_TRD_paired: np.ndarray
    is_gdtcr_sequenced_sublibrary: np.ndarray
    class_code: np.ndarray
    trd_score: np.ndarray
    trab_score: np.ndarray
    trd_minus_trab: np.ndarray
    sublibrary_summary: pd.DataFrame
    conflict_df: pd.DataFrame
    tcr_evidence_audit: pd.DataFrame


@dataclass
class FeatureSpec:
    feature_names: list[str]
    score_feature_cols: list[int]
    metadata_feature_cols: list[int]
    module_feature_cols: list[int]
    deg_feature_cols: list[int]
    specific_tcr_feature_cols: list[int]
    no_metadata_cols: list[int]
    compact_cols: list[int]
    expression_gene_count: int
    deg_genes: list[str]
    specific_tcr_genes: list[str]
    missing_deg_genes: list[str]
    gene_to_feature_col: np.ndarray
    gene_to_family_id: np.ndarray
    family_names: list[str]
    specific_gene_to_col: dict[int, int]
    penalty_gene_to_col: dict[int, int]
    pan_cd3_gene_indices: list[int]


@dataclass
class CandidateResult:
    model_name: str
    score_name: str
    score_tune: np.ndarray
    score_validation: np.ndarray
    score_silver: np.ndarray | None
    pred_tune: np.ndarray
    pred_validation: np.ndarray
    pred_silver: np.ndarray | None
    threshold: float
    tune_metrics: dict[str, Any]
    validation_metrics: dict[str, Any]
    silver_recall: float | None
    model_object: Any = None
    feature_cols: list[int] | None = None
    strategy_notes: str = ""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Train precision-first gdT classifiers on plus6 gold labels.")
    parser.add_argument("--input-h5ad", type=Path, default=DEFAULT_INPUT_H5AD)
    parser.add_argument("--seed", type=int, default=RANDOM_SEED)
    parser.add_argument("--max-negative-train", type=int, default=MAX_NEGATIVE_TRAIN)
    parser.add_argument("--max-deg-features", type=int, default=MAX_DEG_FEATURES)
    parser.add_argument("--max-expression-genes", type=int, default=MAX_EXPRESSION_GENES)
    parser.add_argument("--bootstrap-reps", type=int, default=BOOTSTRAP_REPS)
    parser.add_argument("--min-specificity", type=float, default=MIN_SPECIFICITY)
    parser.add_argument("--no-full-apply", action="store_true", help="Skip full plus6 application even if a model passes.")
    parser.add_argument("--no-pdf", action="store_true", help="Skip headless Chrome PDF rendering.")
    return parser.parse_args()


def ensure_dirs() -> None:
    for path in [
        CLASSIFIER_TABLE_DIR,
        CLASSIFIER_FIGURE_DIR,
        CLASSIFIER_LOG_DIR,
        CLASSIFIER_MODEL_DIR,
        STATIC_DIR,
        STATIC_ASSET_DIR,
    ]:
        path.mkdir(parents=True, exist_ok=True)


def setup_logging() -> None:
    ensure_dirs()
    handlers = [logging.FileHandler(RUN_LOG, mode="w", encoding="utf-8"), logging.StreamHandler()]
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=handlers,
        force=True,
    )


def is_tcr_gene(gene: str) -> bool:
    return str(gene).upper().startswith(TCR_GENE_PREFIXES)


def safe_div(numerator: float, denominator: float) -> float:
    return float(numerator / denominator) if denominator else float("nan")


def wilson_ci(successes: int, n: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if n <= 0:
        return float("nan"), float("nan")
    phat = successes / n
    denom = 1.0 + z * z / n
    center = (phat + z * z / (2.0 * n)) / denom
    half = z * math.sqrt((phat * (1.0 - phat) + z * z / (4.0 * n)) / n) / denom
    return max(0.0, center - half), min(1.0, center + half)


def ppv_at_prevalence(sensitivity: float, specificity: float, prevalence: float) -> float:
    numerator = sensitivity * prevalence
    denominator = numerator + (1.0 - specificity) * (1.0 - prevalence)
    return safe_div(numerator, denominator)


def json_ready(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): json_ready(v) for k, v in value.items()}
    if isinstance(value, list):
        return [json_ready(v) for v in value]
    if isinstance(value, tuple):
        return [json_ready(v) for v in value]
    if isinstance(value, np.ndarray):
        return json_ready(value.tolist())
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        val = float(value)
        return None if math.isnan(val) else val
    if isinstance(value, float) and math.isnan(value):
        return None
    return value


def build_obs_metadata(handle: h5py.File) -> ObsMetadata:
    required_columns_present(handle)
    source = clean_group_values(read_obs_column(handle, "source_gse_id"))
    library_id = normalize_strings(read_obs_column(handle, "library_id"))
    sample_id = normalize_strings(read_obs_column(handle, "sample_id"))
    tissue_key = "tissue_corrected" if "tissue_corrected" in handle["obs"] else "tissue"
    tissue = clean_group_values(read_obs_column(handle, tissue_key))
    has_TRA_TRB_paired = read_bool_obs(handle, "has_TRA_TRB_paired")
    has_any_ab_tcr = read_bool_obs(handle, "has_any_ab_tcr")
    has_TRG_TRD_paired_raw = read_bool_obs(handle, "has_TRG_TRD_paired")
    has_any_gd_tcr_raw = read_bool_obs(handle, "has_any_gd_tcr")
    n_obs = source.shape[0]
    tra_nonempty = read_nonempty_string_mask(handle, "TRA_cdr3") if "TRA_cdr3" in handle["obs"] else np.zeros(n_obs, dtype=bool)
    trb_nonempty = read_nonempty_string_mask(handle, "TRB_cdr3") if "TRB_cdr3" in handle["obs"] else np.zeros(n_obs, dtype=bool)
    trg_nonempty = read_nonempty_string_mask(handle, "TRG_cdr3")
    trd_nonempty = read_nonempty_string_mask(handle, "TRD_cdr3")
    trd_score = read_float_obs(handle, "phase4_trd_score")
    trab_score = read_float_obs(handle, "phase4_trab_score")
    trd_minus_trab = read_float_obs(handle, "phase4_trd_minus_trab")

    corrected_has_any_gd_tcr, corrected_has_TRG_TRD_paired, tcr_evidence_audit = build_corrected_tcr_evidence(
        source,
        has_TRA_TRB_paired,
        has_TRG_TRD_paired_raw,
        has_any_ab_tcr,
        has_any_gd_tcr_raw,
        trg_nonempty,
        trd_nonempty,
    )
    is_gdtcr_sequenced_sublibrary, sublibrary_summary = build_sublibrary_labels(
        source,
        library_id,
        sample_id,
        corrected_has_any_gd_tcr,
        trg_nonempty,
        trd_nonempty,
    )
    class_code, conflict_df = make_truth_labels(
        source,
        is_gdtcr_sequenced_sublibrary,
        has_TRA_TRB_paired,
        corrected_has_TRG_TRD_paired,
        has_any_ab_tcr,
        corrected_has_any_gd_tcr,
        trd_nonempty,
    )
    if int((class_code == 2).sum()) == 0 or int((class_code == 1).sum()) == 0:
        raise RuntimeError("Primary positive and negative gold labels must both be nonempty.")
    if int(conflict_df["n_cells"].max()) > 0:
        raise RuntimeError(f"Ground-truth label conflicts detected:\n{conflict_df}")
    return ObsMetadata(
        source=source,
        tissue=tissue,
        library_id=library_id,
        sample_id=sample_id,
        has_TRA_TRB_paired=has_TRA_TRB_paired,
        has_any_ab_tcr=has_any_ab_tcr,
        has_TRG_TRD_paired_raw=has_TRG_TRD_paired_raw,
        has_any_gd_tcr_raw=has_any_gd_tcr_raw,
        tra_nonempty=tra_nonempty,
        trb_nonempty=trb_nonempty,
        trg_nonempty=trg_nonempty,
        trd_nonempty=trd_nonempty,
        corrected_has_any_gd_tcr=corrected_has_any_gd_tcr,
        corrected_has_TRG_TRD_paired=corrected_has_TRG_TRD_paired,
        is_gdtcr_sequenced_sublibrary=is_gdtcr_sequenced_sublibrary,
        class_code=class_code,
        trd_score=trd_score,
        trab_score=trab_score,
        trd_minus_trab=trd_minus_trab,
        sublibrary_summary=sublibrary_summary,
        conflict_df=conflict_df,
        tcr_evidence_audit=tcr_evidence_audit,
    )


def select_deg_features(
    var_names: pd.Index,
    max_features: int,
    forced_genes: list[str] | None = None,
) -> tuple[list[str], list[str], pd.DataFrame]:
    if not DEG_FEATURE_TABLE.exists():
        raise FileNotFoundError(f"Missing DEG reproducibility feature table: {DEG_FEATURE_TABLE}")
    deg = pd.read_csv(DEG_FEATURE_TABLE)
    required = {"gene", "global_delta_log1p", "global_padj_bh", "reproducible_deg_strict"}
    missing = sorted(required - set(deg.columns))
    if missing:
        raise KeyError(f"Missing columns in DEG table: {missing}")
    forced_genes = forced_genes or []
    deg = deg.loc[~deg["gene"].astype(str).map(is_tcr_gene)].copy()
    deg["abs_delta"] = deg["global_delta_log1p"].abs()
    if deg["reproducible_deg_strict"].astype(bool).any():
        selected = deg.loc[deg["reproducible_deg_strict"].astype(bool)].copy()
    else:
        selected = deg.loc[deg.get("global_deg_pass", True).astype(bool)].copy()
    selected = selected.sort_values(
        ["abs_delta", "global_padj_bh", "gene"], ascending=[False, True, True]
    )
    forced = deg.loc[deg["gene"].astype(str).isin(forced_genes)].copy()
    forced["forced_penalty_gene"] = True
    selected["forced_penalty_gene"] = False
    selected = pd.concat([forced, selected], axis=0, ignore_index=True)
    selected = selected.drop_duplicates("gene", keep="first")
    selected = selected.sort_values(
        ["forced_penalty_gene", "abs_delta", "global_padj_bh", "gene"],
        ascending=[False, False, True, True],
    ).head(max_features)
    present = set(var_names.astype(str))
    deg_genes = [gene for gene in selected["gene"].astype(str).tolist() if gene in present]
    missing_genes = [gene for gene in selected["gene"].astype(str).tolist() if gene not in present]
    selected.assign(feature_present=selected["gene"].astype(str).isin(present)).to_csv(
        CLASSIFIER_TABLE_DIR / "selected_reproducible_deg_features.csv", index=False
    )
    return deg_genes, missing_genes, selected


def select_tcr_expression_genes(var_names: pd.Index, max_genes: int) -> list[str]:
    genes = [str(gene) for gene in var_names.astype(str) if is_tcr_gene(str(gene))]
    core_rank = {gene: rank for rank, gene in enumerate(CORE_TCR_GENES)}

    def priority(gene: str) -> tuple[int, int, str]:
        upper = gene.upper()
        if upper in core_rank:
            return (0, core_rank[upper], upper)
        if upper == "TRDC":
            return (1, 0, upper)
        if upper.startswith("TRDV"):
            return (2, 0, upper)
        if upper.startswith("TRDJ") or upper.startswith("TRDD"):
            return (3, 0, upper)
        if upper in {"TRAC", "TRBC1", "TRBC2"}:
            return (4, 0, upper)
        if upper in {"TRGC1", "TRGC2"}:
            return (5, 0, upper)
        if upper.startswith("TRGV") or upper.startswith("TRGJ"):
            return (6, 0, upper)
        if upper.startswith("TRBV") or upper.startswith("TRBJ") or upper.startswith("TRAV") or upper.startswith("TRAJ"):
            return (7, 0, upper)
        return (8, 0, upper)

    selected = sorted(genes, key=priority)[:max_genes]
    pd.DataFrame(
        {
            "gene": selected,
            "selection_rank": np.arange(1, len(selected) + 1, dtype=int),
            "note": "TRD genes prioritized; TRG V/J/C lower priority; TRA/TRB controls lowest priority.",
        }
    ).to_csv(CLASSIFIER_TABLE_DIR / "selected_tcr_expression_genes.csv", index=False)
    return selected


def build_feature_spec(handle: h5py.File, max_deg_features: int, max_expression_genes: int) -> FeatureSpec:
    var_names = pd.Index(read_string_dataset(handle["var"]["_index"]), dtype="string")
    var_upper = np.asarray([str(gene).upper() for gene in var_names], dtype=object)
    gene_lookup = {str(gene): idx for idx, gene in enumerate(var_names.astype(str))}
    tcr_budget = min(MAX_TCR_EXPRESSION_GENES, max_expression_genes)
    specific_tcr_genes = select_tcr_expression_genes(var_names, tcr_budget)
    penalty_genes = [gene for gene in DEATH_PENALTY_GENES if gene in gene_lookup]
    forced_control_genes = sorted(set(penalty_genes + [gene for gene in PAN_CD3_GENES if gene in gene_lookup]))
    pan_cd3_gene_indices = [gene_lookup[gene] for gene in PAN_CD3_GENES if gene in gene_lookup]
    remaining_expression_gene_budget = max(max_expression_genes - len(specific_tcr_genes), 0)
    deg_feature_budget = min(max_deg_features, remaining_expression_gene_budget)
    deg_genes, missing_deg_genes, _selected = select_deg_features(
        var_names, deg_feature_budget, forced_genes=forced_control_genes
    )

    feature_names: list[str] = []
    score_feature_cols: list[int] = []
    metadata_feature_cols: list[int] = []
    module_feature_cols: list[int] = []
    specific_tcr_feature_cols: list[int] = []
    deg_feature_cols: list[int] = []

    def add_feature(name: str, bucket: list[int]) -> int:
        idx = len(feature_names)
        feature_names.append(name)
        bucket.append(idx)
        return idx

    for name in [
        "phase4_trd_score",
        "phase4_trab_score",
        "phase4_trd_minus_trab",
        "neg_phase4_trab_score",
        "phase4_trd_score_gt_0",
        "phase4_trd_score_gt_0p1",
        "phase4_trd_minus_trab_gt_0p65",
        "phase4_trd_minus_trab_lt_0p35",
        "phase4_trab_score_gt_neg0p05",
        "pan_CD3_sum_low",
    ]:
        add_feature(name, score_feature_cols)

    for name in [
        "has_TRA_TRB_paired",
        "has_any_ab_tcr",
        "TRA_cdr3_nonempty",
        "TRB_cdr3_nonempty",
        "corrected_has_TRG_TRD_paired",
        "corrected_has_any_gd_tcr",
        "TRG_cdr3_nonempty",
        "TRD_cdr3_nonempty",
        "gdTCR_evidence_no_abTCR",
        "abTCR_evidence_no_gdTCR",
        "any_TCR_metadata_evidence",
        "sorted_gdT_source",
    ]:
        add_feature(name, metadata_feature_cols)

    family_names = list(TCR_FAMILIES.keys())
    for family in family_names:
        for stat in ["sum_log1p_cp10k", "max_log1p_cp10k", "n_detected"]:
            add_feature(f"{family}_{stat}", module_feature_cols)
    for family in ["gdT_family", "abT_family"]:
        for stat in ["sum_log1p_cp10k", "max_log1p_cp10k", "n_detected"]:
            add_feature(f"{family}_{stat}", module_feature_cols)

    specific_gene_to_col: dict[int, int] = {}
    for gene in specific_tcr_genes:
        col = add_feature(f"{gene}_log1p_cp10k", specific_tcr_feature_cols)
        specific_gene_to_col[int(gene_lookup[gene])] = col
    penalty_gene_to_col: dict[int, int] = {}

    gene_to_feature_col = np.full(len(var_names), -1, dtype=np.int32)
    for gene in deg_genes:
        col = add_feature(f"DEG_{gene}_log1p_cp10k", deg_feature_cols)
        gene_to_feature_col[gene_lookup[gene]] = col
        if gene in penalty_genes:
            penalty_gene_to_col[gene_lookup[gene]] = col

    gene_to_family_id = np.full(len(var_names), -1, dtype=np.int16)
    selected_tcr_gene_set = set(specific_tcr_genes)
    for family_id, family in enumerate(family_names):
        prefixes = TCR_FAMILIES[family]
        mask = np.zeros(len(var_names), dtype=bool)
        for prefix in prefixes:
            mask |= np.char.startswith(var_upper.astype(str), prefix)
        mask &= np.asarray([gene in selected_tcr_gene_set for gene in var_names.astype(str)], dtype=bool)
        gene_to_family_id[mask] = family_id

    no_metadata_cols = sorted(set(score_feature_cols + module_feature_cols + specific_tcr_feature_cols + deg_feature_cols))
    compact_cols = sorted(set(score_feature_cols + metadata_feature_cols + module_feature_cols + specific_tcr_feature_cols))
    spec = FeatureSpec(
        feature_names=feature_names,
        score_feature_cols=score_feature_cols,
        metadata_feature_cols=metadata_feature_cols,
        module_feature_cols=module_feature_cols,
        deg_feature_cols=deg_feature_cols,
        specific_tcr_feature_cols=specific_tcr_feature_cols,
        no_metadata_cols=no_metadata_cols,
        compact_cols=compact_cols,
        expression_gene_count=len(specific_tcr_feature_cols) + len(deg_feature_cols),
        deg_genes=deg_genes,
        specific_tcr_genes=specific_tcr_genes,
        missing_deg_genes=missing_deg_genes,
        gene_to_feature_col=gene_to_feature_col,
        gene_to_family_id=gene_to_family_id,
        family_names=family_names,
        specific_gene_to_col=specific_gene_to_col,
        penalty_gene_to_col=penalty_gene_to_col,
        pan_cd3_gene_indices=pan_cd3_gene_indices,
    )
    FEATURE_JSON.write_text(
        json.dumps(
            json_ready(
                {
                    "n_features": len(feature_names),
                    "n_score_features": len(score_feature_cols),
                    "n_metadata_features": len(metadata_feature_cols),
                    "n_module_features": len(module_feature_cols),
                    "n_specific_tcr_gene_features": len(specific_tcr_feature_cols),
                    "n_deg_features": len(deg_feature_cols),
                    "expression_gene_count": len(specific_tcr_feature_cols) + len(deg_feature_cols),
                    "max_expression_genes": max_expression_genes,
                    "max_tcr_expression_genes": tcr_budget,
                    "specific_tcr_gene_policy": "TRD constant/V/J/D genes prioritized; selected TRG V/J/C genes retained lower priority because TRG is less specific than TRD; TRAC/TRBC controls included.",
                    "specific_tcr_genes": specific_tcr_genes,
                    "death_penalty_genes": penalty_genes,
                    "pan_cd3_genes": [gene for gene in PAN_CD3_GENES if gene in gene_lookup],
                    "deg_genes": deg_genes,
                    "missing_deg_genes": missing_deg_genes,
                    "feature_names": feature_names,
                    "no_metadata_feature_names": [feature_names[idx] for idx in no_metadata_cols],
                    "compact_feature_names": [feature_names[idx] for idx in compact_cols],
                }
            ),
            indent=2,
        ),
        encoding="utf-8",
    )
    return spec


def fill_metadata_features(matrix: np.ndarray, rows: np.ndarray, obs: ObsMetadata, spec: FeatureSpec) -> None:
    name_to_col = {name: idx for idx, name in enumerate(spec.feature_names)}
    matrix[:, name_to_col["phase4_trd_score"]] = obs.trd_score[rows]
    matrix[:, name_to_col["phase4_trab_score"]] = obs.trab_score[rows]
    matrix[:, name_to_col["phase4_trd_minus_trab"]] = obs.trd_minus_trab[rows]
    matrix[:, name_to_col["neg_phase4_trab_score"]] = -obs.trab_score[rows]
    matrix[:, name_to_col["phase4_trd_score_gt_0"]] = (obs.trd_score[rows] > 0).astype(np.float32)
    matrix[:, name_to_col["phase4_trd_score_gt_0p1"]] = (obs.trd_score[rows] > 0.1).astype(np.float32)
    matrix[:, name_to_col["phase4_trd_minus_trab_gt_0p65"]] = (obs.trd_minus_trab[rows] > 0.65).astype(np.float32)
    matrix[:, name_to_col["phase4_trd_minus_trab_lt_0p35"]] = (obs.trd_minus_trab[rows] < 0.35).astype(np.float32)
    matrix[:, name_to_col["phase4_trab_score_gt_neg0p05"]] = (obs.trab_score[rows] > -0.05).astype(np.float32)
    matrix[:, name_to_col["pan_CD3_sum_low"]] = 0.0

    gd_no_ab = obs.corrected_has_any_gd_tcr[rows] & (~obs.has_any_ab_tcr[rows])
    ab_no_gd = obs.has_any_ab_tcr[rows] & (~obs.corrected_has_any_gd_tcr[rows])
    any_tcr = (
        obs.has_any_ab_tcr[rows]
        | obs.corrected_has_any_gd_tcr[rows]
        | obs.has_TRA_TRB_paired[rows]
        | obs.corrected_has_TRG_TRD_paired[rows]
        | obs.tra_nonempty[rows]
        | obs.trb_nonempty[rows]
        | obs.trg_nonempty[rows]
        | obs.trd_nonempty[rows]
    )
    matrix[:, name_to_col["has_TRA_TRB_paired"]] = obs.has_TRA_TRB_paired[rows].astype(np.float32)
    matrix[:, name_to_col["has_any_ab_tcr"]] = obs.has_any_ab_tcr[rows].astype(np.float32)
    matrix[:, name_to_col["TRA_cdr3_nonempty"]] = obs.tra_nonempty[rows].astype(np.float32)
    matrix[:, name_to_col["TRB_cdr3_nonempty"]] = obs.trb_nonempty[rows].astype(np.float32)
    matrix[:, name_to_col["corrected_has_TRG_TRD_paired"]] = obs.corrected_has_TRG_TRD_paired[rows].astype(np.float32)
    matrix[:, name_to_col["corrected_has_any_gd_tcr"]] = obs.corrected_has_any_gd_tcr[rows].astype(np.float32)
    matrix[:, name_to_col["TRG_cdr3_nonempty"]] = obs.trg_nonempty[rows].astype(np.float32)
    matrix[:, name_to_col["TRD_cdr3_nonempty"]] = obs.trd_nonempty[rows].astype(np.float32)
    matrix[:, name_to_col["gdTCR_evidence_no_abTCR"]] = gd_no_ab.astype(np.float32)
    matrix[:, name_to_col["abTCR_evidence_no_gdTCR"]] = ab_no_gd.astype(np.float32)
    matrix[:, name_to_col["any_TCR_metadata_evidence"]] = any_tcr.astype(np.float32)
    matrix[:, name_to_col["sorted_gdT_source"]] = np.isin(obs.source[rows], list(SORTED_GDT_SOURCES)).astype(np.float32)


def extract_feature_matrix(
    handle: h5py.File,
    rows: np.ndarray,
    obs: ObsMetadata,
    spec: FeatureSpec,
    *,
    label: str,
) -> np.ndarray:
    rows = np.asarray(rows, dtype=np.int64)
    matrix = np.zeros((rows.size, len(spec.feature_names)), dtype=np.float32)
    fill_metadata_features(matrix, rows, obs, spec)
    name_to_col = {name: idx for idx, name in enumerate(spec.feature_names)}
    family_sum_cols = [name_to_col[f"{family}_sum_log1p_cp10k"] for family in spec.family_names]
    family_max_cols = [name_to_col[f"{family}_max_log1p_cp10k"] for family in spec.family_names]
    family_n_cols = [name_to_col[f"{family}_n_detected"] for family in spec.family_names]
    gdt_sum_col = name_to_col["gdT_family_sum_log1p_cp10k"]
    gdt_max_col = name_to_col["gdT_family_max_log1p_cp10k"]
    gdt_n_col = name_to_col["gdT_family_n_detected"]
    abt_sum_col = name_to_col["abT_family_sum_log1p_cp10k"]
    abt_max_col = name_to_col["abT_family_max_log1p_cp10k"]
    abt_n_col = name_to_col["abT_family_n_detected"]

    x_group = handle["X"]
    indptr_ds = x_group["indptr"]
    indices_ds = x_group["indices"]
    data_ds = x_group["data"]
    for out_row, obs_idx in enumerate(rows):
        start = int(indptr_ds[obs_idx])
        end = int(indptr_ds[obs_idx + 1])
        if end <= start:
            continue
        row_indices = indices_ds[start:end].astype(np.int32, copy=False)
        row_data = data_ds[start:end].astype(np.float32, copy=False)
        row_sum = float(np.sum(row_data, dtype=np.float64))
        if row_sum <= 0:
            continue
        row_values = np.log1p(row_data * (TARGET_SUM / row_sum)).astype(np.float32, copy=False)

        feature_cols = spec.gene_to_feature_col[row_indices]
        has_feature = feature_cols >= 0
        if has_feature.any():
            matrix[out_row, feature_cols[has_feature]] = row_values[has_feature]

        for gene_idx, col in spec.specific_gene_to_col.items():
            pos = np.flatnonzero(row_indices == gene_idx)
            if pos.size:
                matrix[out_row, col] = row_values[pos[0]]
        if spec.pan_cd3_gene_indices:
            cd3_sum = 0.0
            for gene_idx in spec.pan_cd3_gene_indices:
                pos = np.flatnonzero(row_indices == gene_idx)
                if pos.size:
                    cd3_sum += float(row_values[pos[0]])
            matrix[out_row, name_to_col["pan_CD3_sum_low"]] = 1.0 if cd3_sum < 0.25 else 0.0

        family_ids = spec.gene_to_family_id[row_indices]
        valid_family = family_ids >= 0
        if valid_family.any():
            family_values = row_values[valid_family]
            row_family_ids = family_ids[valid_family].astype(np.int16, copy=False)
            sums = np.bincount(row_family_ids, weights=family_values, minlength=len(spec.family_names))
            counts = np.bincount(row_family_ids, minlength=len(spec.family_names))
            maxes = np.zeros(len(spec.family_names), dtype=np.float32)
            for family_id in range(len(spec.family_names)):
                vals = family_values[row_family_ids == family_id]
                if vals.size:
                    maxes[family_id] = float(np.max(vals))
            matrix[out_row, family_sum_cols] = sums[: len(spec.family_names)]
            matrix[out_row, family_n_cols] = counts[: len(spec.family_names)]
            matrix[out_row, family_max_cols] = maxes[: len(spec.family_names)]
            matrix[out_row, gdt_sum_col] = sums[0] + sums[1]
            matrix[out_row, gdt_max_col] = max(maxes[0], maxes[1])
            matrix[out_row, gdt_n_col] = counts[0] + counts[1]
            matrix[out_row, abt_sum_col] = sums[2] + sums[3]
            matrix[out_row, abt_max_col] = max(maxes[2], maxes[3])
            matrix[out_row, abt_n_col] = counts[2] + counts[3]

        if out_row and out_row % 50_000 == 0:
            logging.info("Extracted %s features for %s / %s rows", label, f"{out_row:,}", f"{rows.size:,}")
    return matrix


def make_gold_splits(obs: ObsMetadata, seed: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, pd.DataFrame]:
    primary_idx = np.flatnonzero((obs.class_code == 1) | (obs.class_code == 2)).astype(np.int64)
    rng = np.random.default_rng(seed)
    split_by_global: dict[int, str] = {}
    rare_idx: list[int] = []
    split_rows: list[dict[str, Any]] = []

    split_df = pd.DataFrame(
        {
            "obs_index": primary_idx,
            "source_gse_id": obs.source[primary_idx],
            "label": np.where(obs.class_code[primary_idx] == 2, "gdT_gold", "abT_gold"),
        }
    )
    for (source, label), group in split_df.groupby(["source_gse_id", "label"], sort=True):
        idx = group["obs_index"].to_numpy(dtype=np.int64)
        rng.shuffle(idx)
        n = idx.size
        if n < 5:
            rare_idx.extend(idx.tolist())
            continue
        n_train = max(1, int(round(n * 0.60)))
        n_tune = max(1, int(round(n * 0.20)))
        if n_train + n_tune >= n:
            n_train = max(1, n - 2)
            n_tune = 1
        train_part = idx[:n_train]
        tune_part = idx[n_train : n_train + n_tune]
        validation_part = idx[n_train + n_tune :]
        for value in train_part:
            split_by_global[int(value)] = "train"
        for value in tune_part:
            split_by_global[int(value)] = "tune"
        for value in validation_part:
            split_by_global[int(value)] = "validation"
        split_rows.append(
            {
                "source_gse_id": source,
                "label": label,
                "n_cells": n,
                "train": int(train_part.size),
                "tune": int(tune_part.size),
                "validation": int(validation_part.size),
                "rare_group": False,
            }
        )

    if rare_idx:
        rare = np.asarray(rare_idx, dtype=np.int64)
        rng.shuffle(rare)
        n_train = int(round(rare.size * 0.60))
        n_tune = int(round(rare.size * 0.20))
        for value in rare[:n_train]:
            split_by_global[int(value)] = "train"
        for value in rare[n_train : n_train + n_tune]:
            split_by_global[int(value)] = "tune"
        for value in rare[n_train + n_tune :]:
            split_by_global[int(value)] = "validation"
        split_rows.append(
            {
                "source_gse_id": "rare_combined",
                "label": "mixed",
                "n_cells": int(rare.size),
                "train": int(n_train),
                "tune": int(n_tune),
                "validation": int(rare.size - n_train - n_tune),
                "rare_group": True,
            }
        )

    split_labels = np.asarray([split_by_global[int(idx)] for idx in primary_idx], dtype=object)
    train_idx = primary_idx[split_labels == "train"]
    tune_idx = primary_idx[split_labels == "tune"]
    validation_idx = primary_idx[split_labels == "validation"]
    split_summary = pd.DataFrame(split_rows)
    split_summary.to_csv(CLASSIFIER_TABLE_DIR / "gold_split_by_source_label.csv", index=False)

    overall_rows = []
    for split_name, idx in [("train", train_idx), ("tune", tune_idx), ("validation", validation_idx)]:
        labels = obs.class_code[idx]
        overall_rows.append(
            {
                "split": split_name,
                "n_cells": int(idx.size),
                "gdT_gold": int((labels == 2).sum()),
                "abT_gold": int((labels == 1).sum()),
                "gdT_prevalence": safe_div(int((labels == 2).sum()), int(idx.size)),
            }
        )
    pd.DataFrame(overall_rows).to_csv(CLASSIFIER_TABLE_DIR / "gold_split_overall.csv", index=False)
    return train_idx, tune_idx, validation_idx, split_summary


def local_positions(eval_rows: np.ndarray, target_rows: np.ndarray) -> np.ndarray:
    lookup = pd.Series(np.arange(eval_rows.size, dtype=np.int64), index=eval_rows)
    return lookup.loc[target_rows].to_numpy(dtype=np.int64)


def sample_training_positions(pos_train: np.ndarray, y_eval: np.ndarray, max_negative: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    positives = pos_train[y_eval[pos_train] == 1]
    negatives = pos_train[y_eval[pos_train] == 0]
    max_neg = min(negatives.size, max(max_negative, positives.size * 4))
    if negatives.size > max_neg:
        negatives = rng.choice(negatives, size=max_neg, replace=False)
    sampled = np.concatenate([positives, negatives])
    rng.shuffle(sampled)
    return sampled.astype(np.int64)


def metric_dict(y_true: np.ndarray, y_pred: np.ndarray, score: np.ndarray | None = None) -> dict[str, Any]:
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()
    has_both = np.unique(y_true).size == 2
    out: dict[str, Any] = {
        "n_cells": int(y_true.size),
        "n_positive": int(y_true.sum()),
        "n_negative": int((y_true == 0).sum()),
        "predicted_positive": int(y_pred.sum()),
        "tp": int(tp),
        "fp": int(fp),
        "tn": int(tn),
        "fn": int(fn),
        "precision": float(precision_score(y_true, y_pred, zero_division=0)),
        "recall": float(recall_score(y_true, y_pred, zero_division=0)),
        "specificity": safe_div(tn, tn + fp),
        "f1": float(f1_score(y_true, y_pred, zero_division=0)),
        "f0.5": float(fbeta_score(y_true, y_pred, beta=0.5, zero_division=0)),
        "balanced_accuracy": float(balanced_accuracy_score(y_true, y_pred)) if has_both else float("nan"),
        "mcc": float(matthews_corrcoef(y_true, y_pred)) if has_both else float("nan"),
    }
    if score is not None and has_both:
        out["roc_auc"] = float(roc_auc_score(y_true, score))
        out["pr_auc"] = float(average_precision_score(y_true, score))
    else:
        out["roc_auc"] = float("nan")
        out["pr_auc"] = float("nan")
    return out


def choose_threshold(
    y_true: np.ndarray,
    score: np.ndarray,
    *,
    min_specificity: float,
    beta: float,
    min_threshold: float | None = None,
) -> tuple[float, dict[str, Any]]:
    fpr, tpr, thresholds = roc_curve(y_true, score)
    n_pos = int(y_true.sum())
    n_neg = int((y_true == 0).sum())
    tp = np.rint(tpr * n_pos).astype(int)
    fp = np.rint(fpr * n_neg).astype(int)
    fn = n_pos - tp
    tn = n_neg - fp
    specificity = np.divide(tn, tn + fp, out=np.zeros_like(tn, dtype=float), where=(tn + fp) > 0)
    recall = np.divide(tp, tp + fn, out=np.zeros_like(tp, dtype=float), where=(tp + fn) > 0)
    precision = np.divide(tp, tp + fp, out=np.zeros_like(tp, dtype=float), where=(tp + fp) > 0)
    beta2 = beta * beta
    fbeta = np.divide(
        (1.0 + beta2) * precision * recall,
        beta2 * precision + recall,
        out=np.zeros_like(precision, dtype=float),
        where=(beta2 * precision + recall) > 0,
    )
    threshold_ok = np.isfinite(thresholds)
    if min_threshold is not None:
        threshold_ok &= thresholds >= min_threshold
    valid = threshold_ok & (specificity >= min_specificity) & ((tp + fp) > 0)
    fallback_used = False
    if not valid.any():
        valid = threshold_ok & ((tp + fp) > 0)
        fallback_used = True
    if not valid.any():
        threshold = float(np.nanmax(score) + 1e-6)
        pred = score >= threshold
        metrics = metric_dict(y_true, pred.astype(np.int8), score)
        metrics["threshold_fallback_used"] = True
        return threshold, metrics
    valid_idx = np.flatnonzero(valid)
    order = np.lexsort(
        (
            thresholds[valid_idx],
            recall[valid_idx],
            precision[valid_idx],
            fbeta[valid_idx],
        )
    )
    best = int(valid_idx[order[-1]])
    threshold = float(thresholds[best])
    pred = (score >= threshold).astype(np.int8)
    metrics = metric_dict(y_true, pred, score)
    metrics["threshold_fallback_used"] = fallback_used
    return threshold, metrics


def choose_dual_threshold(
    y_true: np.ndarray,
    trd_score: np.ndarray,
    trab_score: np.ndarray,
    *,
    min_specificity: float,
    beta: float,
) -> tuple[tuple[float, float], dict[str, Any], np.ndarray]:
    trd_candidates = np.unique(np.quantile(trd_score, np.linspace(0.60, 0.995, 36)))
    trab_candidates = np.unique(np.quantile(trab_score, np.linspace(0.005, 0.45, 36)))
    best_metrics: dict[str, Any] | None = None
    best_thresholds = (float("nan"), float("nan"))
    best_pred = np.zeros(y_true.shape, dtype=np.int8)
    fallback_metrics: dict[str, Any] | None = None
    fallback_thresholds = (float("nan"), float("nan"))
    fallback_pred = np.zeros(y_true.shape, dtype=np.int8)
    pos_mask = y_true == 1
    neg_mask = ~pos_mask
    n_pos = int(pos_mask.sum())
    n_neg = int(neg_mask.sum())
    score = trd_score - trab_score
    roc_auc = float(roc_auc_score(y_true, score))
    pr_auc = float(average_precision_score(y_true, score))
    for trd_cut in trd_candidates:
        trd_mask = trd_score >= trd_cut
        for trab_cut in trab_candidates:
            pred = (trd_mask & (trab_score <= trab_cut)).astype(np.int8)
            predicted_positive = int(pred.sum())
            if predicted_positive == 0:
                continue
            tp = int(pred[pos_mask].sum())
            fp = int(pred[neg_mask].sum())
            fn = n_pos - tp
            tn = n_neg - fp
            precision = safe_div(tp, tp + fp)
            recall = safe_div(tp, tp + fn)
            specificity = safe_div(tn, tn + fp)
            beta2 = beta * beta
            fbeta = safe_div((1.0 + beta2) * precision * recall, (beta2 * precision) + recall)
            metrics = {
                "n_cells": int(y_true.size),
                "n_positive": n_pos,
                "n_negative": n_neg,
                "predicted_positive": predicted_positive,
                "tp": tp,
                "fp": fp,
                "tn": tn,
                "fn": fn,
                "precision": precision,
                "recall": recall,
                "specificity": specificity,
                "f1": safe_div(2.0 * precision * recall, precision + recall),
                "f0.5": fbeta,
                "balanced_accuracy": (recall + specificity) / 2.0,
                "mcc": float(matthews_corrcoef(y_true, pred)),
                "roc_auc": roc_auc,
                "pr_auc": pr_auc,
            }
            metrics["trd_threshold"] = float(trd_cut)
            metrics["trab_threshold"] = float(trab_cut)
            if fallback_metrics is None or (
                metrics["f0.5"],
                metrics["precision"],
                metrics["recall"],
            ) > (
                fallback_metrics["f0.5"],
                fallback_metrics["precision"],
                fallback_metrics["recall"],
            ):
                fallback_metrics = metrics
                fallback_thresholds = (float(trd_cut), float(trab_cut))
                fallback_pred = pred
            if metrics["specificity"] < min_specificity:
                continue
            if best_metrics is None:
                best_metrics = metrics
                best_thresholds = (float(trd_cut), float(trab_cut))
                best_pred = pred
                continue
            if (
                metrics["f0.5"],
                metrics["precision"],
                metrics["recall"],
            ) > (
                best_metrics["f0.5"],
                best_metrics["precision"],
                best_metrics["recall"],
            ):
                best_metrics = metrics
                best_thresholds = (float(trd_cut), float(trab_cut))
                best_pred = pred
    if best_metrics is None:
        if fallback_metrics is None:
            raise RuntimeError("No dual-rule threshold produced predicted positives.")
        best_thresholds = fallback_thresholds
        best_metrics = fallback_metrics
        best_metrics["threshold_fallback_used"] = True
        best_pred = fallback_pred
    else:
        best_metrics["threshold_fallback_used"] = False
    return best_thresholds, best_metrics, best_pred


def build_baseline_results(
    obs: ObsMetadata,
    tune_idx: np.ndarray,
    validation_idx: np.ndarray,
    silver_idx: np.ndarray,
    y_tune: np.ndarray,
    y_validation: np.ndarray,
    min_specificity: float,
) -> list[CandidateResult]:
    results: list[CandidateResult] = []
    baseline_specs = [
        ("baseline_TRD_score_high", "phase4_trd_score", obs.trd_score),
        ("baseline_TRAB_score_low", "neg_phase4_trab_score", -obs.trab_score),
        ("baseline_TRD_minus_TRAB_high", "phase4_trd_minus_trab", obs.trd_minus_trab),
    ]
    for name, score_name, score in baseline_specs:
        threshold, tune_metrics = choose_threshold(
            y_tune,
            score[tune_idx],
            min_specificity=min_specificity,
            beta=F_BETA,
        )
        tune_pred = (score[tune_idx] >= threshold).astype(np.int8)
        validation_score = score[validation_idx]
        validation_pred = (validation_score >= threshold).astype(np.int8)
        validation_metrics = metric_dict(y_validation, validation_pred, validation_score)
        silver_score = score[silver_idx] if silver_idx.size else None
        silver_pred = (silver_score >= threshold).astype(np.int8) if silver_score is not None and silver_score.size else None
        silver_recall = float(np.mean(silver_score >= threshold)) if silver_score is not None and silver_score.size else None
        results.append(
            CandidateResult(
                model_name=name,
                score_name=score_name,
                score_tune=score[tune_idx],
                score_validation=validation_score,
                score_silver=silver_score,
                pred_tune=tune_pred,
                pred_validation=validation_pred,
                pred_silver=silver_pred,
                threshold=threshold,
                tune_metrics=tune_metrics,
                validation_metrics=validation_metrics,
                silver_recall=silver_recall,
                strategy_notes="Pure score baseline; threshold tuned on tune split at specificity constraint.",
            )
        )

    dual_thresholds, tune_metrics, tune_pred = choose_dual_threshold(
        y_tune,
        obs.trd_score[tune_idx],
        obs.trab_score[tune_idx],
        min_specificity=min_specificity,
        beta=F_BETA,
    )
    trd_cut, trab_cut = dual_thresholds
    validation_pred = ((obs.trd_score[validation_idx] >= trd_cut) & (obs.trab_score[validation_idx] <= trab_cut)).astype(
        np.int8
    )
    validation_score = obs.trd_score[validation_idx] - obs.trab_score[validation_idx]
    validation_metrics = metric_dict(y_validation, validation_pred, validation_score)
    silver_score = obs.trd_score[silver_idx] - obs.trab_score[silver_idx] if silver_idx.size else None
    silver_pred = (
        ((obs.trd_score[silver_idx] >= trd_cut) & (obs.trab_score[silver_idx] <= trab_cut)).astype(np.int8)
        if silver_idx.size
        else None
    )
    silver_recall = (
        float(np.mean((obs.trd_score[silver_idx] >= trd_cut) & (obs.trab_score[silver_idx] <= trab_cut)))
        if silver_idx.size
        else None
    )
    results.append(
        CandidateResult(
            model_name="baseline_dual_TRD_high_TRAB_low",
            score_name="dual_rule",
            score_tune=obs.trd_score[tune_idx] - obs.trab_score[tune_idx],
            score_validation=validation_score,
            score_silver=silver_score,
            pred_tune=tune_pred,
            pred_validation=validation_pred,
            pred_silver=silver_pred,
            threshold=float("nan"),
            tune_metrics=tune_metrics,
            validation_metrics=validation_metrics,
            silver_recall=silver_recall,
            strategy_notes=f"Dual rule threshold: phase4_trd_score >= {trd_cut:.6g} and phase4_trab_score <= {trab_cut:.6g}.",
        )
    )
    return results


def fit_xgb_classifier(x_train: np.ndarray, y_train: np.ndarray, seed: int) -> XGBClassifier:
    n_pos = int(y_train.sum())
    n_neg = int((y_train == 0).sum())
    scale_pos_weight = max(n_neg / max(n_pos, 1), 1.0)
    model = XGBClassifier(
        n_estimators=420,
        max_depth=4,
        learning_rate=0.045,
        subsample=0.85,
        colsample_bytree=0.85,
        min_child_weight=2.0,
        reg_lambda=1.5,
        objective="binary:logistic",
        eval_metric="logloss",
        tree_method="hist",
        n_jobs=32,
        random_state=seed,
        scale_pos_weight=scale_pos_weight,
    )
    model.fit(x_train, y_train)
    return model


def train_candidate_models(
    x_eval: np.ndarray,
    y_eval: np.ndarray,
    pos_train: np.ndarray,
    pos_tune: np.ndarray,
    pos_validation: np.ndarray,
    pos_silver: np.ndarray,
    obs_eval_rows: np.ndarray,
    obs: ObsMetadata,
    spec: FeatureSpec,
    max_negative_train: int,
    min_specificity: float,
    seed: int,
) -> list[CandidateResult]:
    results: list[CandidateResult] = []
    y_tune = y_eval[pos_tune]
    y_validation = y_eval[pos_validation]
    train_sample = sample_training_positions(pos_train, y_eval, max_negative_train, seed)
    logging.info(
        "Training sampled set: %s cells, positives=%s, negatives=%s",
        f"{train_sample.size:,}",
        f"{int(y_eval[train_sample].sum()):,}",
        f"{int((y_eval[train_sample] == 0).sum()):,}",
    )

    xgb_specs = [
        (
            "strategy_A_score_RNA_TCRgene_XGB",
            spec.no_metadata_cols,
            "Scores + TCR gene-expression modules + reproducible non-TCR DEG genes; no TCR metadata features.",
        ),
        (
            "strategy_B_TCRaware_score_RNA_XGB",
            list(range(x_eval.shape[1])),
            "Scores + TCR gene-expression modules + reproducible DEG genes + TCR metadata evidence.",
        ),
    ]
    trained_models: dict[str, tuple[Any, list[int], str]] = {}
    for model_name, cols, note in xgb_specs:
        logging.info("Training %s with %s features", model_name, len(cols))
        model = fit_xgb_classifier(x_eval[train_sample][:, cols], y_eval[train_sample], seed)
        trained_models[model_name] = (model, cols, note)
        tune_score = model.predict_proba(x_eval[pos_tune][:, cols])[:, 1].astype(np.float32)
        validation_score = model.predict_proba(x_eval[pos_validation][:, cols])[:, 1].astype(np.float32)
        silver_score = (
            model.predict_proba(x_eval[pos_silver][:, cols])[:, 1].astype(np.float32) if pos_silver.size else None
        )
        threshold, tune_metrics = choose_threshold(
            y_tune,
            tune_score,
            min_specificity=min_specificity,
            beta=F_BETA,
        )
        tune_pred = (tune_score >= threshold).astype(np.int8)
        validation_pred = (validation_score >= threshold).astype(np.int8)
        validation_metrics = metric_dict(y_validation, validation_pred, validation_score)
        silver_pred = (silver_score >= threshold).astype(np.int8) if silver_score is not None and silver_score.size else None
        silver_recall = float(np.mean(silver_score >= threshold)) if silver_score is not None and silver_score.size else None
        results.append(
            CandidateResult(
                model_name=model_name,
                score_name="predicted_probability",
                score_tune=tune_score,
                score_validation=validation_score,
                score_silver=silver_score,
                pred_tune=tune_pred,
                pred_validation=validation_pred,
                pred_silver=silver_pred,
                threshold=threshold,
                tune_metrics=tune_metrics,
                validation_metrics=validation_metrics,
                silver_recall=silver_recall,
                model_object=model,
                feature_cols=cols,
                strategy_notes=note,
            )
        )

    base_a_model, base_a_cols, base_a_note = trained_models["strategy_A_score_RNA_TCRgene_XGB"]
    base_a_tune = base_a_model.predict_proba(x_eval[pos_tune][:, base_a_cols])[:, 1].astype(np.float32)
    base_a_validation = base_a_model.predict_proba(x_eval[pos_validation][:, base_a_cols])[:, 1].astype(np.float32)
    base_a_silver = (
        base_a_model.predict_proba(x_eval[pos_silver][:, base_a_cols])[:, 1].astype(np.float32)
        if pos_silver.size
        else None
    )
    penalty_a_tune = apply_expression_penalty_gates(base_a_tune, x_eval[pos_tune], spec)
    penalty_a_validation = apply_expression_penalty_gates(base_a_validation, x_eval[pos_validation], spec)
    penalty_a_silver = apply_expression_penalty_gates(base_a_silver, x_eval[pos_silver], spec) if base_a_silver is not None else None
    threshold, tune_metrics = choose_threshold(
        y_tune,
        penalty_a_tune,
        min_specificity=min_specificity,
        beta=F_BETA,
        min_threshold=0.5,
    )
    tune_pred = (penalty_a_tune >= threshold).astype(np.int8)
    validation_pred = (penalty_a_validation >= threshold).astype(np.int8)
    validation_metrics = metric_dict(y_validation, validation_pred, penalty_a_validation)
    silver_pred = (
        (penalty_a_silver >= threshold).astype(np.int8)
        if penalty_a_silver is not None and penalty_a_silver.size
        else None
    )
    silver_recall = (
        float(np.mean(penalty_a_silver >= threshold)) if penalty_a_silver is not None and penalty_a_silver.size else None
    )
    results.append(
        CandidateResult(
            model_name="strategy_A2_score_RNA_TCRgene_with_expression_death_penalties",
            score_name="penalty_gated_predicted_probability",
            score_tune=penalty_a_tune,
            score_validation=penalty_a_validation,
            score_silver=penalty_a_silver,
            pred_tune=tune_pred,
            pred_validation=validation_pred,
            pred_silver=silver_pred,
            threshold=threshold,
            tune_metrics=tune_metrics,
            validation_metrics=validation_metrics,
            silver_recall=silver_recall,
            model_object=base_a_model,
            feature_cols=base_a_cols,
            strategy_notes=(
                base_a_note
                + " Adds expression death penalties for FOXP3/CD4 positivity and low pan-CD3 evidence. "
                + "No TCR metadata fields are used."
            ),
        )
    )

    logging.info("Training strategy_C_compact_logistic")
    compact_cols = spec.compact_cols
    logistic = make_pipeline(
        StandardScaler(),
        LogisticRegression(
            solver="lbfgs",
            class_weight="balanced",
            max_iter=1000,
            random_state=seed,
        ),
    )
    logistic.fit(x_eval[train_sample][:, compact_cols], y_eval[train_sample])
    tune_score = logistic.predict_proba(x_eval[pos_tune][:, compact_cols])[:, 1].astype(np.float32)
    validation_score = logistic.predict_proba(x_eval[pos_validation][:, compact_cols])[:, 1].astype(np.float32)
    silver_score = logistic.predict_proba(x_eval[pos_silver][:, compact_cols])[:, 1].astype(np.float32) if pos_silver.size else None
    threshold, tune_metrics = choose_threshold(
        y_tune,
        tune_score,
        min_specificity=min_specificity,
        beta=F_BETA,
    )
    tune_pred = (tune_score >= threshold).astype(np.int8)
    validation_pred = (validation_score >= threshold).astype(np.int8)
    validation_metrics = metric_dict(y_validation, validation_pred, validation_score)
    silver_pred = (silver_score >= threshold).astype(np.int8) if silver_score is not None and silver_score.size else None
    silver_recall = float(np.mean(silver_score >= threshold)) if silver_score is not None and silver_score.size else None
    results.append(
        CandidateResult(
            model_name="strategy_C_compact_logistic",
            score_name="predicted_probability",
            score_tune=tune_score,
            score_validation=validation_score,
            score_silver=silver_score,
            pred_tune=tune_pred,
            pred_validation=validation_pred,
            pred_silver=silver_pred,
            threshold=threshold,
            tune_metrics=tune_metrics,
            validation_metrics=validation_metrics,
            silver_recall=silver_recall,
            model_object=logistic,
            feature_cols=compact_cols,
            strategy_notes="Compact score + TCR metadata + TCR-gene module logistic fallback strategy.",
        )
    )

    # Precision-gated probability adjustment using the full TCR-aware XGB as the
    # base estimator. This is intentionally reported separately from the model
    # score because it uses TCR metadata gates that overlap with the gold-label
    # definition.
    base_model, base_cols, base_note = trained_models["strategy_B_TCRaware_score_RNA_XGB"]
    for split_name, positions in [("tune", pos_tune), ("validation", pos_validation), ("silver", pos_silver)]:
        if positions.size == 0:
            continue
        logging.info("Scoring tiered TCR-aware gates on %s split", split_name)
    base_tune = base_model.predict_proba(x_eval[pos_tune][:, base_cols])[:, 1].astype(np.float32)
    base_validation = base_model.predict_proba(x_eval[pos_validation][:, base_cols])[:, 1].astype(np.float32)
    base_silver = (
        base_model.predict_proba(x_eval[pos_silver][:, base_cols])[:, 1].astype(np.float32) if pos_silver.size else None
    )
    tune_score = apply_tcr_precision_gates(base_tune, obs_eval_rows[pos_tune], obs)
    validation_score = apply_tcr_precision_gates(base_validation, obs_eval_rows[pos_validation], obs)
    silver_score = apply_tcr_precision_gates(base_silver, obs_eval_rows[pos_silver], obs) if base_silver is not None else None
    threshold, tune_metrics = choose_threshold(y_tune, tune_score, min_specificity=min_specificity, beta=F_BETA)
    tune_pred = (tune_score >= threshold).astype(np.int8)
    validation_pred = (validation_score >= threshold).astype(np.int8)
    validation_metrics = metric_dict(y_validation, validation_pred, validation_score)
    silver_pred = (silver_score >= threshold).astype(np.int8) if silver_score is not None and silver_score.size else None
    silver_recall = float(np.mean(silver_score >= threshold)) if silver_score is not None and silver_score.size else None
    results.append(
        CandidateResult(
            model_name="strategy_B2_tiered_TCRaware_precision_gated",
            score_name="gated_predicted_probability",
            score_tune=tune_score,
            score_validation=validation_score,
            score_silver=silver_score,
            pred_tune=tune_pred,
            pred_validation=validation_pred,
            pred_silver=silver_pred,
            threshold=threshold,
            tune_metrics=tune_metrics,
            validation_metrics=validation_metrics,
            silver_recall=silver_recall,
            model_object=base_model,
            feature_cols=base_cols,
            strategy_notes=(
                base_note
                + " Probability post-processing boosts paired gdTCR/no-abTCR evidence and suppresses abTCR-only evidence."
            ),
        )
    )
    penalty_tune = apply_expression_penalty_gates(tune_score, x_eval[pos_tune], spec)
    penalty_validation = apply_expression_penalty_gates(validation_score, x_eval[pos_validation], spec)
    penalty_silver = apply_expression_penalty_gates(silver_score, x_eval[pos_silver], spec) if silver_score is not None else None
    threshold, tune_metrics = choose_threshold(
        y_tune,
        penalty_tune,
        min_specificity=min_specificity,
        beta=F_BETA,
        min_threshold=0.5,
    )
    tune_pred = (penalty_tune >= threshold).astype(np.int8)
    validation_pred = (penalty_validation >= threshold).astype(np.int8)
    validation_metrics = metric_dict(y_validation, validation_pred, penalty_validation)
    silver_pred = (penalty_silver >= threshold).astype(np.int8) if penalty_silver is not None and penalty_silver.size else None
    silver_recall = float(np.mean(penalty_silver >= threshold)) if penalty_silver is not None and penalty_silver.size else None
    results.append(
        CandidateResult(
            model_name="strategy_B3_tiered_TCRaware_with_expression_death_penalties",
            score_name="penalty_gated_predicted_probability",
            score_tune=penalty_tune,
            score_validation=penalty_validation,
            score_silver=penalty_silver,
            pred_tune=tune_pred,
            pred_validation=validation_pred,
            pred_silver=silver_pred,
            threshold=threshold,
            tune_metrics=tune_metrics,
            validation_metrics=validation_metrics,
            silver_recall=silver_recall,
            model_object=base_model,
            feature_cols=base_cols,
            strategy_notes=(
                base_note
                + " Adds expression death penalties for FOXP3/CD4 positivity and low pan-CD3 evidence."
            ),
        )
    )
    return results


def apply_tcr_precision_gates(probability: np.ndarray | None, rows: np.ndarray, obs: ObsMetadata) -> np.ndarray | None:
    if probability is None:
        return None
    out = probability.astype(np.float32, copy=True)
    paired_gd_no_ab = obs.corrected_has_TRG_TRD_paired[rows] & (~obs.has_any_ab_tcr[rows])
    any_gd_no_ab = obs.corrected_has_any_gd_tcr[rows] & (~obs.has_any_ab_tcr[rows])
    ab_only = obs.has_any_ab_tcr[rows] & (~obs.corrected_has_any_gd_tcr[rows])
    low_trdminus_ab = ab_only & (obs.trd_minus_trab[rows] < 0.35)
    out[any_gd_no_ab] = np.maximum(out[any_gd_no_ab], 0.97)
    out[paired_gd_no_ab] = np.maximum(out[paired_gd_no_ab], 0.995)
    out[ab_only] = np.minimum(out[ab_only], out[ab_only] * 0.20)
    out[low_trdminus_ab] = np.minimum(out[low_trdminus_ab], 0.02)
    return out


def apply_expression_penalty_gates(
    probability: np.ndarray | None,
    x_values: np.ndarray,
    spec: FeatureSpec,
) -> np.ndarray | None:
    if probability is None:
        return None
    out = probability.astype(np.float32, copy=True)
    name_to_col = {name: idx for idx, name in enumerate(spec.feature_names)}
    penalty_mask = np.zeros(out.shape[0], dtype=bool)
    foxp3_col = next((col for gene_idx, col in spec.penalty_gene_to_col.items() if spec.feature_names[col] == "DEG_FOXP3_log1p_cp10k"), None)
    cd4_col = next((col for gene_idx, col in spec.penalty_gene_to_col.items() if spec.feature_names[col] == "DEG_CD4_log1p_cp10k"), None)
    if foxp3_col is not None:
        penalty_mask |= x_values[:, foxp3_col] > 0.25
    if cd4_col is not None:
        penalty_mask |= x_values[:, cd4_col] > 0.75
    penalty_mask |= x_values[:, name_to_col["pan_CD3_sum_low"]] > 0.5
    if penalty_mask.any():
        out[penalty_mask] = np.minimum(out[penalty_mask], 0.03)
    return out


def results_to_frame(results: list[CandidateResult], split: str) -> pd.DataFrame:
    rows = []
    for result in results:
        metrics = result.tune_metrics if split == "tune" else result.validation_metrics
        row = {
            "model": result.model_name,
            "score_name": result.score_name,
            "threshold": result.threshold,
            "silver_recall_at_threshold": result.silver_recall,
            "strategy_notes": result.strategy_notes,
        }
        row.update(metrics)
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["f0.5", "precision", "recall"], ascending=False).reset_index(drop=True)


def bootstrap_delta_fbeta(
    y_true: np.ndarray,
    model_pred: np.ndarray,
    baseline_pred: np.ndarray,
    reps: int,
    seed: int,
) -> dict[str, Any]:
    rng = np.random.default_rng(seed)
    deltas = np.zeros(reps, dtype=np.float64)
    n = y_true.size
    for rep in range(reps):
        idx = rng.integers(0, n, size=n)
        deltas[rep] = fbeta_score(y_true[idx], model_pred[idx], beta=F_BETA, zero_division=0) - fbeta_score(
            y_true[idx], baseline_pred[idx], beta=F_BETA, zero_division=0
        )
    return {
        "bootstrap_reps": reps,
        "delta_f0.5_mean": float(np.mean(deltas)),
        "delta_f0.5_ci_low": float(np.quantile(deltas, 0.025)),
        "delta_f0.5_ci_high": float(np.quantile(deltas, 0.975)),
        "p_delta_le_0": float(np.mean(deltas <= 0)),
    }


def model_selection_key(result: CandidateResult) -> tuple[float, float, float, int]:
    penalty_preference = 1 if "death_penalties" in result.model_name else 0
    return (
        float(result.validation_metrics["f0.5"]),
        float(result.validation_metrics["precision"]),
        float(result.validation_metrics["recall"]),
        penalty_preference,
    )


def choose_deployable_model(candidate_results: list[CandidateResult], accepted_names: set[str]) -> CandidateResult | None:
    accepted = [result for result in candidate_results if result.model_name in accepted_names]
    if not accepted:
        return None
    no_metadata_penalty = [result for result in accepted if result.model_name.startswith("strategy_A2")]
    if no_metadata_penalty:
        return max(no_metadata_penalty, key=model_selection_key)
    no_metadata = [result for result in accepted if result.model_name.startswith("strategy_A")]
    pool = no_metadata if no_metadata else accepted
    return max(pool, key=model_selection_key)


def add_acceptance_columns(
    validation_df: pd.DataFrame,
    results: list[CandidateResult],
    y_validation: np.ndarray,
    min_specificity: float,
    bootstrap_reps: int,
    seed: int,
) -> tuple[pd.DataFrame, CandidateResult | None, CandidateResult]:
    baseline_results = [result for result in results if result.model_name.startswith("baseline_")]
    candidate_results = [result for result in results if not result.model_name.startswith("baseline_")]
    best_baseline = max(
        baseline_results,
        key=lambda item: (
            item.validation_metrics["f0.5"],
            item.validation_metrics["precision"],
            item.validation_metrics["recall"],
        ),
    )
    baseline_pred = best_baseline.pred_validation

    rows = []
    accepted_model: CandidateResult | None = None
    accepted_names: set[str] = set()
    best_candidate: CandidateResult | None = None
    for result in candidate_results:
        pred = result.pred_validation
        boot = bootstrap_delta_fbeta(y_validation, pred, baseline_pred, bootstrap_reps, seed)
        metrics = result.validation_metrics
        recall_delta = float(metrics["recall"] - best_baseline.validation_metrics["recall"])
        fbeta_delta = float(metrics["f0.5"] - best_baseline.validation_metrics["f0.5"])
        accepted = (
            metrics["specificity"] >= min_specificity
            and metrics["precision"] >= 0.90
            and fbeta_delta > 0
            and recall_delta >= 0.05
            and boot["delta_f0.5_ci_low"] > 0
        )
        row = {
            "model": result.model_name,
            "best_baseline_model": best_baseline.model_name,
            "validation_f0.5": metrics["f0.5"],
            "best_baseline_f0.5": best_baseline.validation_metrics["f0.5"],
            "delta_f0.5": fbeta_delta,
            "validation_precision": metrics["precision"],
            "validation_recall": metrics["recall"],
            "best_baseline_recall": best_baseline.validation_metrics["recall"],
            "delta_recall": recall_delta,
            "validation_specificity": metrics["specificity"],
            "accepted": accepted,
        }
        row.update(boot)
        rows.append(row)
        if best_candidate is None or model_selection_key(result) > model_selection_key(best_candidate):
            best_candidate = result
        if accepted and (accepted_model is None or model_selection_key(result) > model_selection_key(accepted_model)):
            accepted_model = result
        if accepted:
            accepted_names.add(result.model_name)

    acceptance_df = pd.DataFrame(rows)
    acceptance_df.to_csv(CLASSIFIER_TABLE_DIR / "classifier_acceptance_vs_best_baseline.csv", index=False)
    validation_df = validation_df.merge(
        acceptance_df[["model", "delta_f0.5", "delta_f0.5_ci_low", "delta_f0.5_ci_high", "p_delta_le_0", "accepted"]],
        on="model",
        how="left",
    )
    deployable_model = choose_deployable_model(candidate_results, accepted_names)
    return validation_df, deployable_model, best_baseline


def write_group_metrics(
    result: CandidateResult,
    y_validation: np.ndarray,
    validation_idx: np.ndarray,
    obs: ObsMetadata,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    pred = result.pred_validation
    rows_source: list[dict[str, Any]] = []
    rows_tissue: list[dict[str, Any]] = []
    for key_name, groups, out_rows in [
        ("source_gse_id", obs.source[validation_idx], rows_source),
        ("tissue", obs.tissue[validation_idx], rows_tissue),
    ]:
        for group in sorted(pd.Series(groups).astype(str).unique()):
            mask = groups.astype(str) == group
            if int(mask.sum()) == 0:
                continue
            metrics = metric_dict(y_validation[mask], pred[mask], result.score_validation[mask])
            out_rows.append({key_name: group, **metrics})
    source_df = pd.DataFrame(rows_source).sort_values(["n_cells", "source_gse_id"], ascending=[False, True])
    tissue_df = pd.DataFrame(rows_tissue).sort_values(["n_cells", "tissue"], ascending=[False, True])
    source_df.to_csv(CLASSIFIER_TABLE_DIR / "selected_model_validation_metrics_by_source.csv", index=False)
    tissue_df.to_csv(CLASSIFIER_TABLE_DIR / "selected_model_validation_metrics_by_tissue.csv", index=False)
    return source_df, tissue_df


def build_source_tcr_mode_table(obs: ObsMetadata) -> pd.DataFrame:
    metadata_ab = obs.has_any_ab_tcr | obs.has_TRA_TRB_paired | obs.tra_nonempty | obs.trb_nonempty
    metadata_gd = (
        obs.corrected_has_any_gd_tcr
        | obs.corrected_has_TRG_TRD_paired
        | obs.trg_nonempty
        | obs.trd_nonempty
    )
    rows: list[dict[str, Any]] = []
    for source in sorted(pd.Series(obs.source).astype(str).unique()):
        mask = obs.source.astype(str) == source
        has_ab = bool(metadata_ab[mask].any())
        has_gd = bool(metadata_gd[mask].any())
        if has_ab and has_gd:
            mode = "TRAB_and_gdTCR_metadata"
        elif has_ab:
            mode = "TRAB_only_metadata"
        elif has_gd:
            mode = "gdTCR_only_metadata"
        else:
            mode = "no_TCRseq_metadata"
        rows.append(
            {
                "source_gse_id": source,
                "n_cells": int(mask.sum()),
                "has_any_TRAB_metadata": has_ab,
                "has_any_gdTCR_metadata": has_gd,
                "tcr_metadata_mode": mode,
            }
        )
    out = pd.DataFrame(rows)
    out.to_csv(CLASSIFIER_TABLE_DIR / "source_tcr_metadata_modes.csv", index=False)
    return out


def write_trab_only_false_positive_tables(
    results: list[CandidateResult],
    validation_idx: np.ndarray,
    y_validation: np.ndarray,
    obs: ObsMetadata,
    source_mode: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    mode_map = source_mode.set_index("source_gse_id")["tcr_metadata_mode"].to_dict()
    source_values = pd.Series(obs.source[validation_idx]).astype(str).to_numpy(dtype=object)
    trab_only_mask = np.asarray([mode_map.get(str(source), "") == "TRAB_only_metadata" for source in source_values], dtype=bool)
    negative_mask = y_validation == 0
    rows_overall: list[dict[str, Any]] = []
    rows_source: list[dict[str, Any]] = []
    for result in results:
        pred = result.pred_validation.astype(bool)
        fp_mask = trab_only_mask & negative_mask & pred
        denom_mask = trab_only_mask & negative_mask
        rows_overall.append(
            {
                "model": result.model_name,
                "n_validation_cells_in_TRAB_only_sources": int(trab_only_mask.sum()),
                "n_validation_abT_gold_in_TRAB_only_sources": int(denom_mask.sum()),
                "false_positive_in_TRAB_only_sources": int(fp_mask.sum()),
                "false_positive_rate_in_TRAB_only_sources": safe_div(int(fp_mask.sum()), int(denom_mask.sum())),
                "predicted_positive_in_TRAB_only_sources": int((trab_only_mask & pred).sum()),
            }
        )
        for source in sorted(pd.Series(source_values[trab_only_mask]).astype(str).unique()):
            source_mask = trab_only_mask & (source_values == source)
            source_denom = source_mask & negative_mask
            source_fp = source_denom & pred
            rows_source.append(
                {
                    "model": result.model_name,
                    "source_gse_id": source,
                    "n_validation_cells": int(source_mask.sum()),
                    "n_validation_abT_gold": int(source_denom.sum()),
                    "false_positive": int(source_fp.sum()),
                    "false_positive_rate": safe_div(int(source_fp.sum()), int(source_denom.sum())),
                    "predicted_positive": int((source_mask & pred).sum()),
                }
            )
    overall_df = pd.DataFrame(rows_overall).sort_values("false_positive_in_TRAB_only_sources", ascending=False)
    by_source_df = pd.DataFrame(rows_source)
    if not by_source_df.empty:
        by_source_df = by_source_df.sort_values(["false_positive", "source_gse_id"], ascending=[False, True])
    overall_df.to_csv(CLASSIFIER_TABLE_DIR / "trab_only_source_false_positives_validation_overall.csv", index=False)
    by_source_df.to_csv(CLASSIFIER_TABLE_DIR / "trab_only_source_false_positives_validation_by_source.csv", index=False)
    return overall_df, by_source_df


def write_quality_sensitivity_metrics(
    results: list[CandidateResult],
    validation_idx: np.ndarray,
    y_validation: np.ndarray,
    obs: ObsMetadata,
) -> pd.DataFrame:
    source_values = pd.Series(obs.source[validation_idx]).astype(str).to_numpy(dtype=object)
    low_quality_mask = np.isin(source_values, list(SUBOPTIMAL_SORTED_GDT_SOURCES))
    rows: list[dict[str, Any]] = []
    for result in results:
        for subset_name, mask in [
            ("all_validation", np.ones(y_validation.shape[0], dtype=bool)),
            ("exclude_suboptimal_GDTlung2023july_7p", ~low_quality_mask),
            ("suboptimal_GDTlung2023july_7p_only", low_quality_mask),
        ]:
            if int(mask.sum()) == 0:
                continue
            metrics = metric_dict(y_validation[mask], result.pred_validation[mask], result.score_validation[mask])
            rows.append({"model": result.model_name, "validation_subset": subset_name, **metrics})
    out = pd.DataFrame(rows)
    out.to_csv(CLASSIFIER_TABLE_DIR / "suboptimal_gdtlung_quality_sensitivity_metrics.csv", index=False)
    return out


def write_prevalence_aware_metrics(
    results: list[CandidateResult],
    validation_idx: np.ndarray,
    y_validation: np.ndarray,
    obs: ObsMetadata,
) -> pd.DataFrame:
    gold_prevalence = float(y_validation.mean())
    all_gold_prevalence = safe_div(int((obs.class_code == 2).sum()), int(((obs.class_code == 1) | (obs.class_code == 2)).sum()))
    plus6_gold_fraction = safe_div(int((obs.class_code == 2).sum()), int(obs.class_code.size))
    scenario_values = sorted(set(PREVALENCE_SCENARIOS + [plus6_gold_fraction, all_gold_prevalence, gold_prevalence]))
    rows: list[dict[str, Any]] = []
    for result in results:
        metrics = result.validation_metrics
        tp = int(metrics["tp"])
        fn = int(metrics["fn"])
        tn = int(metrics["tn"])
        fp = int(metrics["fp"])
        sensitivity_low, sensitivity_high = wilson_ci(tp, tp + fn)
        specificity_low, specificity_high = wilson_ci(tn, tn + fp)
        for prevalence in scenario_values:
            rows.append(
                {
                    "model": result.model_name,
                    "prevalence": prevalence,
                    "prevalence_percent": prevalence * 100.0,
                    "validation_recall_observed": metrics["recall"],
                    "validation_specificity_observed": metrics["specificity"],
                    "validation_specificity_wilson95_low": specificity_low,
                    "validation_recall_wilson95_low": sensitivity_low,
                    "ppv_observed_at_prevalence": ppv_at_prevalence(metrics["recall"], metrics["specificity"], prevalence),
                    "ppv_conservative_wilson_at_prevalence": ppv_at_prevalence(
                        sensitivity_low,
                        specificity_low,
                        prevalence,
                    ),
                    "expected_predicted_positive_per_million_observed": (
                        metrics["recall"] * prevalence + (1.0 - metrics["specificity"]) * (1.0 - prevalence)
                    )
                    * 1_000_000.0,
                    "expected_false_positive_per_million_observed": (1.0 - metrics["specificity"])
                    * (1.0 - prevalence)
                    * 1_000_000.0,
                    "expected_false_positive_per_million_conservative_wilson": (1.0 - specificity_low)
                    * (1.0 - prevalence)
                    * 1_000_000.0,
                    "validation_tp": tp,
                    "validation_fp": fp,
                    "validation_tn": tn,
                    "validation_fn": fn,
                    "validation_negative_cells": tn + fp,
                    "prevalence_scenario_note": (
                        "plus6_gdT_gold_fraction"
                        if abs(prevalence - plus6_gold_fraction) < 1e-12
                        else "primary_gold_validation_prevalence"
                        if abs(prevalence - gold_prevalence) < 1e-12
                        else "primary_gold_overall_prevalence"
                        if abs(prevalence - all_gold_prevalence) < 1e-12
                        else "fixed_low_prevalence_scenario"
                    ),
                }
            )
    out = pd.DataFrame(rows).sort_values(["model", "prevalence"])
    out.to_csv(CLASSIFIER_TABLE_DIR / "prevalence_aware_ppv_scenarios.csv", index=False)
    return out


def save_model_artifact(result: CandidateResult, spec: FeatureSpec) -> Path:
    path = CLASSIFIER_MODEL_DIR / "selected_gdt_classifier.pkl"
    payload = {
        "model_name": result.model_name,
        "threshold": result.threshold,
        "feature_cols": result.feature_cols,
        "feature_names": spec.feature_names,
        "selected_feature_names": [spec.feature_names[idx] for idx in (result.feature_cols or [])],
        "strategy_notes": result.strategy_notes,
        "model_object": result.model_object,
    }
    with path.open("wb") as handle:
        pickle.dump(payload, handle)
    return path


def apply_model_to_full_dataset(
    handle: h5py.File,
    result: CandidateResult,
    obs: ObsMetadata,
    spec: FeatureSpec,
    source_mode: pd.DataFrame,
    *,
    chunk_size: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if result.model_object is None or result.feature_cols is None:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    n_obs = int(handle["obs"]["_index"].shape[0])
    total_pred = 0
    source_counts: dict[str, int] = {}
    tissue_counts: dict[str, int] = {}
    mode_map = source_mode.set_index("source_gse_id")["tcr_metadata_mode"].to_dict()
    overlap = {
        "all_cells": n_obs,
        "predicted_putative_gdT": 0,
        "predicted_in_TRAB_only_metadata_sources": 0,
        "predicted_and_sorted_gdT_source": 0,
        "predicted_and_corrected_paired_gdTCR_no_abTCR": 0,
        "predicted_and_corrected_any_gdTCR_no_abTCR": 0,
        "predicted_and_legacy_TRD_minus_TRAB_gt_0p4": 0,
    }
    for start in range(0, n_obs, chunk_size):
        end = min(start + chunk_size, n_obs)
        rows = np.arange(start, end, dtype=np.int64)
        x_chunk = extract_feature_matrix(handle, rows, obs, spec, label=f"full_apply_{start}_{end}")
        score = result.model_object.predict_proba(x_chunk[:, result.feature_cols])[:, 1].astype(np.float32)
        if result.model_name in {
            "strategy_B2_tiered_TCRaware_precision_gated",
            "strategy_B3_tiered_TCRaware_with_expression_death_penalties",
        }:
            score = apply_tcr_precision_gates(score, rows, obs)
        if result.model_name in {
            "strategy_A2_score_RNA_TCRgene_with_expression_death_penalties",
            "strategy_B3_tiered_TCRaware_with_expression_death_penalties",
        }:
            score = apply_expression_penalty_gates(score, x_chunk, spec)
        pred = score >= result.threshold
        pred_count = int(pred.sum())
        total_pred += pred_count
        if pred_count:
            source_series = pd.Series(obs.source[rows][pred]).astype(str).value_counts()
            tissue_series = pd.Series(obs.tissue[rows][pred]).astype(str).value_counts()
            for key, value in source_series.items():
                source_counts[key] = source_counts.get(key, 0) + int(value)
            for key, value in tissue_series.items():
                tissue_counts[key] = tissue_counts.get(key, 0) + int(value)
            pred_rows = rows[pred]
            pred_source_modes = np.asarray(
                [mode_map.get(str(source), "") for source in obs.source[pred_rows]], dtype=object
            )
            overlap["predicted_in_TRAB_only_metadata_sources"] += int(
                (pred_source_modes == "TRAB_only_metadata").sum()
            )
            overlap["predicted_and_sorted_gdT_source"] += int(np.isin(obs.source[pred_rows], list(SORTED_GDT_SOURCES)).sum())
            overlap["predicted_and_corrected_paired_gdTCR_no_abTCR"] += int(
                (obs.corrected_has_TRG_TRD_paired[pred_rows] & (~obs.has_any_ab_tcr[pred_rows])).sum()
            )
            overlap["predicted_and_corrected_any_gdTCR_no_abTCR"] += int(
                (obs.corrected_has_any_gd_tcr[pred_rows] & (~obs.has_any_ab_tcr[pred_rows])).sum()
            )
            overlap["predicted_and_legacy_TRD_minus_TRAB_gt_0p4"] += int((obs.trd_minus_trab[pred_rows] > 0.4).sum())
        if start and start % 500_000 == 0:
            logging.info("Applied selected classifier to %s / %s cells", f"{start:,}", f"{n_obs:,}")
    overlap["predicted_putative_gdT"] = total_pred
    overall_df = pd.DataFrame(
        [
            {
                "model": result.model_name,
                "threshold": result.threshold,
                "n_cells": n_obs,
                "predicted_putative_gdT": total_pred,
                "predicted_fraction": safe_div(total_pred, n_obs),
            }
        ]
    )
    source_df = (
        pd.DataFrame({"source_gse_id": list(source_counts.keys()), "predicted_putative_gdT": list(source_counts.values())})
        .sort_values("predicted_putative_gdT", ascending=False)
        .reset_index(drop=True)
    )
    tissue_df = (
        pd.DataFrame({"tissue": list(tissue_counts.keys()), "predicted_putative_gdT": list(tissue_counts.values())})
        .sort_values("predicted_putative_gdT", ascending=False)
        .reset_index(drop=True)
    )
    overlap_df = pd.DataFrame([overlap])
    overall_df.to_csv(CLASSIFIER_TABLE_DIR / "selected_model_full_dataset_prediction_overall.csv", index=False)
    source_df.to_csv(CLASSIFIER_TABLE_DIR / "selected_model_full_dataset_prediction_by_source.csv", index=False)
    tissue_df.to_csv(CLASSIFIER_TABLE_DIR / "selected_model_full_dataset_prediction_by_tissue.csv", index=False)
    overlap_df.to_csv(CLASSIFIER_TABLE_DIR / "selected_model_full_dataset_prediction_overlap.csv", index=False)
    return overall_df, source_df, tissue_df


def plot_roc_pr(results: list[CandidateResult], y_validation: np.ndarray) -> tuple[Path, Path]:
    roc_path = CLASSIFIER_FIGURE_DIR / "classifier_validation_roc.png"
    pr_path = CLASSIFIER_FIGURE_DIR / "classifier_validation_pr.png"
    fig, ax = plt.subplots(figsize=(6.5, 5.5), constrained_layout=True)
    for result in results:
        if np.unique(y_validation).size < 2:
            continue
        fpr, tpr, _ = roc_curve(y_validation, result.score_validation)
        ax.plot(fpr, tpr, lw=1.4, label=f"{result.model_name} ({result.validation_metrics['roc_auc']:.3f})")
    ax.plot([0, 1], [0, 1], color="#777777", lw=0.8, linestyle="--")
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.set_title("Held-out gold validation ROC")
    ax.legend(fontsize=7, frameon=False, loc="lower right")
    fig.savefig(roc_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.5, 5.5), constrained_layout=True)
    for result in results:
        precision, recall, _ = precision_recall_curve(y_validation, result.score_validation)
        ax.plot(recall, precision, lw=1.4, label=f"{result.model_name} ({result.validation_metrics['pr_auc']:.3f})")
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.set_title("Held-out gold validation precision-recall")
    ax.legend(fontsize=7, frameon=False, loc="lower left")
    fig.savefig(pr_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return roc_path, pr_path


def plot_confusion(result: CandidateResult, y_validation: np.ndarray) -> Path:
    path = CLASSIFIER_FIGURE_DIR / "selected_model_validation_confusion_matrix.png"
    pred = result.pred_validation
    cm = confusion_matrix(y_validation, pred, labels=[0, 1])
    fig, ax = plt.subplots(figsize=(4.8, 4.2), constrained_layout=True)
    im = ax.imshow(cm, cmap="Blues")
    ax.set_xticks([0, 1], labels=["Pred abT", "Pred gdT"])
    ax.set_yticks([0, 1], labels=["True abT", "True gdT"])
    for i in range(2):
        for j in range(2):
            ax.text(j, i, f"{cm[i, j]:,}", ha="center", va="center", color="#111111", fontsize=13)
    ax.set_title(f"{result.model_name}\nHeld-out validation confusion matrix")
    fig.colorbar(im, ax=ax, shrink=0.8)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_score_distribution(results: list[CandidateResult], y_validation: np.ndarray) -> Path:
    path = CLASSIFIER_FIGURE_DIR / "classifier_validation_score_distributions.png"
    top = sorted(results, key=lambda r: r.validation_metrics["f0.5"], reverse=True)[:5]
    fig, axes = plt.subplots(len(top), 1, figsize=(7, max(3, 1.75 * len(top))), constrained_layout=True)
    if len(top) == 1:
        axes = [axes]
    for ax, result in zip(axes, top):
        neg_scores = result.score_validation[y_validation == 0]
        pos_scores = result.score_validation[y_validation == 1]
        finite_scores = np.concatenate([neg_scores[np.isfinite(neg_scores)], pos_scores[np.isfinite(pos_scores)]])
        hist_kwargs: dict[str, Any] = {"bins": 80, "range": (0.0, 1.0), "density": True}
        if finite_scores.size:
            low = float(np.nanmin(finite_scores))
            high = float(np.nanmax(finite_scores))
            if high - low < 1e-8:
                low -= 0.5
                high += 0.5
                hist_kwargs["bins"] = 10
            else:
                pad = 0.02 * (high - low)
                low -= pad
                high += pad
            hist_kwargs["range"] = (low, high)
        ax.hist(neg_scores, alpha=0.55, color="#2d6f9f", label="abT_gold", **hist_kwargs)
        ax.hist(pos_scores, alpha=0.55, color="#b51f2a", label="gdT_gold", **hist_kwargs)
        if np.isfinite(result.threshold):
            ax.axvline(result.threshold, color="black", linestyle="--", lw=1.0)
        ax.set_title(result.model_name, fontsize=9)
        ax.set_ylabel("Density")
    axes[-1].set_xlabel("Score/probability")
    axes[0].legend(frameon=False, fontsize=8)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_feature_importance(result: CandidateResult, spec: FeatureSpec) -> Path | None:
    if result.model_object is None or result.feature_cols is None or not hasattr(result.model_object, "feature_importances_"):
        return None
    path = CLASSIFIER_FIGURE_DIR / "selected_model_feature_importance.png"
    importances = np.asarray(result.model_object.feature_importances_, dtype=np.float64)
    names = np.asarray([spec.feature_names[idx] for idx in result.feature_cols], dtype=object)
    order = np.argsort(importances)[::-1][:35]
    fig, ax = plt.subplots(figsize=(8, max(5, 0.22 * len(order) + 1.5)), constrained_layout=True)
    y = np.arange(order.size)
    ax.barh(y, importances[order], color="#4a6fa5")
    ax.set_yticks(y, labels=names[order])
    ax.invert_yaxis()
    ax.set_xlabel("XGBoost feature importance")
    ax.set_title(f"Selected model feature importance: {result.model_name}")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def copy_assets(paths: list[Path | None]) -> list[str]:
    copied: list[str] = []
    for path in paths:
        if path is None or not path.exists():
            continue
        target = STATIC_ASSET_DIR / path.name
        shutil.copyfile(path, target)
        copied.append(path.name)
    return copied


def render_pdf(no_pdf: bool) -> None:
    if no_pdf:
        return
    chrome_candidates = [
        Path("/usr/bin/google-chrome"),
        Path("/usr/bin/google-chrome-stable"),
        shutil.which("google-chrome"),
        shutil.which("google-chrome-stable"),
    ]
    chrome = next((Path(path) for path in chrome_candidates if path and Path(path).exists()), None)
    if chrome is None:
        logging.warning("Skipping PDF export because google-chrome was not found")
        return
    subprocess.run(
        [
            str(chrome),
            "--headless",
            "--disable-gpu",
            "--no-sandbox",
            "--print-to-pdf-no-header",
            f"--print-to-pdf={REPORT_PDF}",
            str(REPORT_HTML),
        ],
        check=True,
    )


def write_reports(
    input_h5ad: Path,
    n_obs: int,
    obs: ObsMetadata,
    split_summary: pd.DataFrame,
    tune_df: pd.DataFrame,
    validation_df: pd.DataFrame,
    acceptance_df: pd.DataFrame,
    selected_result: CandidateResult | None,
    best_baseline: CandidateResult,
    source_metrics: pd.DataFrame,
    tissue_metrics: pd.DataFrame,
    trab_fp_overall: pd.DataFrame,
    trab_fp_by_source: pd.DataFrame,
    quality_sensitivity_df: pd.DataFrame,
    prevalence_df: pd.DataFrame,
    full_overall: pd.DataFrame,
    full_source: pd.DataFrame,
    spec: FeatureSpec,
    copied_assets: list[str],
    model_artifact: Path | None,
) -> None:
    selected = selected_result
    selected_status = "accepted" if selected is not None else "not accepted"
    best_candidate_name = validation_df.loc[~validation_df["model"].str.startswith("baseline_"), "model"].iloc[0]
    primary_counts = pd.DataFrame(
        [
            {
                "class": "gdT_gold",
                "n_cells": int((obs.class_code == 2).sum()),
            },
            {
                "class": "abT_gold",
                "n_cells": int((obs.class_code == 1).sum()),
            },
            {
                "class": "gdT_silver",
                "n_cells": int((obs.class_code == 3).sum()),
            },
        ]
    )
    top_features = pd.DataFrame({"feature": spec.deg_genes[:40]})
    selected_prevalence = (
        prevalence_df.loc[prevalence_df["model"] == selected.model_name].copy() if selected is not None else pd.DataFrame()
    )
    if selected_prevalence.empty:
        selected_prevalence = prevalence_df.loc[prevalence_df["model"] == best_candidate_name].copy()
    selected_prevalence_display = selected_prevalence[
        [
            "prevalence_percent",
            "prevalence_scenario_note",
            "ppv_observed_at_prevalence",
            "ppv_conservative_wilson_at_prevalence",
            "expected_false_positive_per_million_observed",
            "expected_false_positive_per_million_conservative_wilson",
        ]
    ].head(20)
    selected_block = []
    if selected is not None:
        selected_block = [
            f"- Selected accepted classifier: `{selected.model_name}`",
            f"- Validation precision: `{selected.validation_metrics['precision']:.4f}`",
            f"- Validation recall: `{selected.validation_metrics['recall']:.4f}`",
            f"- Validation specificity: `{selected.validation_metrics['specificity']:.4f}`",
            f"- Validation F0.5: `{selected.validation_metrics['f0.5']:.4f}`",
            f"- Silver gdT recall: `{selected.silver_recall}`",
            f"- Model artifact: `{model_artifact}`",
        ]
    else:
        selected_block = [
            "- No candidate satisfied all acceptance rules.",
            f"- Best non-baseline candidate by validation F0.5: `{best_candidate_name}`",
            "- No production full-dataset classifier artifact was published.",
        ]
    full_block = []
    if not full_overall.empty:
        row = full_overall.iloc[0]
        full_block = [
            f"- Full plus6 predicted putative gdT cells: `{int(row['predicted_putative_gdT']):,}` / `{int(row['n_cells']):,}` (`{row['predicted_fraction']:.4%}`).",
            f"- Full-dataset count table: `{CLASSIFIER_TABLE_DIR / 'selected_model_full_dataset_prediction_overall.csv'}`",
        ]
    else:
        full_block = ["- Full-dataset application was skipped because no classifier passed the acceptance gate or because `--no-full-apply` was used."]

    lines = [
        "# gdT DEG/TCR Classifier Training Report",
        "",
        "## Scope",
        "",
        f"- Input H5AD: `{input_h5ad}`",
        f"- H5AD cells: `{n_obs:,}`",
        "- Analysis is read-only. The H5AD was opened in read mode and checked for unchanged size and mtime after execution.",
        "- Gold labels are the same project gdT/abT labels used by the package-style evaluation, including the local GSE144469 TCR evidence correction.",
        "- `GDTlung2023july_7p` remains part of the sorted gdT gold label definition, but it is flagged as suboptimal library quality and reported in a separate sensitivity audit.",
        "- Silver gdT cells are never used for training or threshold selection; they are only used for sensitivity recall reporting.",
        "- Thresholds are selected on the tune split with a precision-first constraint: specificity >= 0.995 and F0.5 maximization.",
        "- Final claims are made only on the held-out gold validation split.",
        "- The real dataset is extremely imbalanced: gdT cells are a small minority, so report interpretation emphasizes prevalence-aware PPV and false-positive burden, not only validation precision.",
        "",
        "## Ground Truth Counts",
        "",
        dataframe_to_markdown(primary_counts),
        "",
        "## Feature Design",
        "",
        f"- Score features: `{len(spec.score_feature_cols)}`",
        f"- TCR metadata features: `{len(spec.metadata_feature_cols)}`",
        f"- TCR gene-expression module features: `{len(spec.module_feature_cols)}`",
        f"- Specific TCR gene-expression features: `{len(spec.specific_tcr_feature_cols)}`",
        f"- Reproducible non-TCR DEG features: `{len(spec.deg_feature_cols)}`",
        f"- Total explicit expression genes used: `{spec.expression_gene_count}` (cap: `{MAX_EXPRESSION_GENES}`).",
        "- TCR genes are excluded from the DEG-derived feature list, but TCR gene-expression modules are explicit classifier features.",
        "- Explicit TCR gene features prioritize TRD constant/V/J/D genes; selected TRG V/J/C genes are retained but lower priority because TRG expression is less specific than TRD; TRAC/TRBC controls are included.",
        "- FOXP3, CD4, and low pan-CD3 evidence are treated as expression-level penalty signals in the tiered penalty strategy; additional negative-enriched reproducible genes can enter through the DEG feature set.",
        "- TCR metadata-aware strategies are retained as diagnostics. Because gold labels are partly defined from TCR metadata and the dataset is highly imbalanced, the deployable selected classifier is chosen from the score/RNA/TCR-gene strategies when one passes the acceptance gate.",
        "",
        "Top non-TCR DEG features:",
        "",
        dataframe_to_markdown(top_features),
        "",
        "## Gold Split",
        "",
        dataframe_to_markdown(pd.read_csv(CLASSIFIER_TABLE_DIR / "gold_split_overall.csv")),
        "",
        "## Validation Outcome",
        "",
        f"- Acceptance status: `{selected_status}`",
        f"- Best pure score baseline on held-out validation: `{best_baseline.model_name}` with F0.5 `{best_baseline.validation_metrics['f0.5']:.4f}` and recall `{best_baseline.validation_metrics['recall']:.4f}`.",
        *selected_block,
        "",
        "Acceptance rules:",
        "",
        "- validation specificity >= 0.995",
        "- validation precision >= 0.90",
        "- validation F0.5 greater than the best pure score baseline",
        "- validation recall at least 0.05 higher than the best pure score baseline",
        "- paired bootstrap 95% CI for delta F0.5 entirely above zero",
        "",
        "## Tune Metrics",
        "",
        dataframe_to_markdown(tune_df, max_rows=20),
        "",
        "## Held-Out Validation Metrics",
        "",
        dataframe_to_markdown(validation_df, max_rows=20),
        "",
        "## Acceptance Against Best Baseline",
        "",
        dataframe_to_markdown(acceptance_df, max_rows=20),
        "",
        "## Prevalence-Aware Imbalance Check",
        "",
        "Held-out gold validation is enriched for gdT compared with deployment on the whole plus6 object. The table below recalculates PPV under low-prevalence scenarios and also gives a conservative PPV using Wilson 95% lower bounds for recall/specificity. This is important because even a very small false-positive rate can dominate when gdT prevalence is below 1%.",
        "",
        dataframe_to_markdown(selected_prevalence_display, max_rows=20),
        "",
        "## Selected Model Strata",
        "",
        "By source:",
        "",
        dataframe_to_markdown(source_metrics, max_rows=40),
        "",
        "By tissue:",
        "",
        dataframe_to_markdown(tissue_metrics, max_rows=40),
        "",
        "## TRAB-Only Dataset False Positive Audit",
        "",
        "These rows count false positives in validation cells from sources with TRAB metadata but no gdTCR metadata.",
        "",
        dataframe_to_markdown(trab_fp_overall, max_rows=20),
        "",
        dataframe_to_markdown(trab_fp_by_source, max_rows=40),
        "",
        "## GDTlung Library Quality Sensitivity",
        "",
        "`GDTlung2023july_7p` is sorted gdT but marked as suboptimal library quality, so performance is also shown after excluding it from validation.",
        "",
        dataframe_to_markdown(quality_sensitivity_df, max_rows=30),
        "",
        "## Full Dataset Application",
        "",
        *full_block,
        "",
        "Top full-dataset predicted sources:",
        "",
        dataframe_to_markdown(full_source.head(30) if not full_source.empty else full_source),
        "",
        "## Output Files",
        "",
        f"- Tune metrics: `{CLASSIFIER_TABLE_DIR / 'classifier_metrics_tune.csv'}`",
        f"- Validation metrics: `{CLASSIFIER_TABLE_DIR / 'classifier_metrics_validation.csv'}`",
        f"- Acceptance table: `{CLASSIFIER_TABLE_DIR / 'classifier_acceptance_vs_best_baseline.csv'}`",
        f"- Prevalence-aware PPV table: `{CLASSIFIER_TABLE_DIR / 'prevalence_aware_ppv_scenarios.csv'}`",
        f"- TRAB-only FP audit: `{CLASSIFIER_TABLE_DIR / 'trab_only_source_false_positives_validation_overall.csv'}`",
        f"- GDTlung sensitivity audit: `{CLASSIFIER_TABLE_DIR / 'suboptimal_gdtlung_quality_sensitivity_metrics.csv'}`",
        f"- Feature manifest: `{FEATURE_JSON}`",
        f"- HTML report: `{REPORT_HTML}`",
        f"- PDF report: `{REPORT_PDF}`",
        "",
    ]
    REPORT_MD.write_text("\n".join(lines), encoding="utf-8")

    metric_cards = []
    if selected is not None:
        metric_cards.extend(
            [
                ("Selected classifier", selected.model_name),
                ("Precision", f"{selected.validation_metrics['precision']:.3f}"),
                ("Recall", f"{selected.validation_metrics['recall']:.3f}"),
                ("Specificity", f"{selected.validation_metrics['specificity']:.3f}"),
                ("F0.5", f"{selected.validation_metrics['f0.5']:.3f}"),
                ("Silver recall", "NA" if selected.silver_recall is None else f"{selected.silver_recall:.3f}"),
            ]
        )
    else:
        metric_cards.extend(
            [
                ("Acceptance", "No model passed"),
                ("Best baseline", best_baseline.model_name),
                ("Baseline F0.5", f"{best_baseline.validation_metrics['f0.5']:.3f}"),
                ("Best candidate", best_candidate_name),
            ]
        )
    figure_titles = {
        "classifier_validation_roc.png": "Validation ROC curves",
        "classifier_validation_pr.png": "Validation precision-recall curves",
        "selected_model_validation_confusion_matrix.png": "Selected model confusion matrix",
        "classifier_validation_score_distributions.png": "Validation score distributions",
        "selected_model_feature_importance.png": "Selected model feature importance",
    }
    figure_html = "\n".join(
        f"<section class='figure'><h3>{html.escape(figure_titles.get(name, name))}</h3>"
        f"<img src='assets/{html.escape(name)}' alt='{html.escape(figure_titles.get(name, name))}'></section>"
        for name in copied_assets
    )
    css = """
    :root{--bg:#f5f6f7;--paper:#fff;--ink:#1d252c;--muted:#64707d;--line:#d9dee4;--accent:#9f2731;--blue:#2f6f8f}
    *{box-sizing:border-box} body{margin:0;background:var(--bg);color:var(--ink);font-family:Arial,Helvetica,sans-serif;line-height:1.5}
    main{width:min(1240px,calc(100vw - 32px));margin:22px auto 46px}
    section{background:var(--paper);border:1px solid var(--line);padding:22px;margin-top:16px}.hero{border-top:6px solid var(--accent)}
    h1{font-size:31px;line-height:1.15;margin:0 0 10px;letter-spacing:0} h2{font-size:22px;margin:0 0 12px} h3{font-size:15px;margin:0 0 8px}
    p{margin:0 0 12px;color:var(--muted)} code{background:#eef1f4;padding:1px 5px;border-radius:4px}
    .metrics{display:grid;grid-template-columns:repeat(auto-fit,minmax(170px,1fr));gap:10px;margin-top:14px}.metric{background:#fbfcfd;border:1px solid var(--line);padding:10px 12px}.metric b{display:block;font-size:20px}.metric span{display:block;font-size:11px;text-transform:uppercase;color:var(--muted)}
    table{border-collapse:collapse;width:100%;font-size:12px;margin:8px 0 14px} th,td{border:1px solid var(--line);padding:5px 7px;text-align:left;vertical-align:top} th{background:#eef1f4}
    .figures{display:grid;grid-template-columns:repeat(auto-fit,minmax(430px,1fr));gap:14px}.figure{padding:12px;margin:0;background:#fbfcfd}.figure img{width:100%;height:auto;border:1px solid var(--line);background:white}
    """
    html_parts = [
        "<!DOCTYPE html><html lang='en'><head><meta charset='utf-8'>",
        "<meta name='viewport' content='width=device-width, initial-scale=1'>",
        "<title>gdT DEG/TCR Classifier Training</title>",
        f"<style>{css}</style></head><body><main>",
        "<section class='hero'>",
        "<h1>gdT DEG/TCR Classifier Training</h1>",
        f"<p>Read-only classifier training and held-out validation from <code>{html.escape(str(input_h5ad))}</code>.</p>",
        "<div class='metrics'>",
        *[
            f"<div class='metric'><b>{html.escape(value)}</b><span>{html.escape(label)}</span></div>"
            for label, value in metric_cards
        ],
        "</div></section>",
        "<section><h2>Interpretation</h2>",
        f"<p>Acceptance status: <code>{html.escape(selected_status)}</code>. The best pure score baseline was <code>{html.escape(best_baseline.model_name)}</code>. Thresholds were selected on the tune split at specificity >= {MIN_SPECIFICITY}; metrics below are held-out validation metrics.</p>",
        "<p>TCR-aware metadata strategies use evidence that overlaps with gold-label construction. They are retained as diagnostics, but deployment selection prefers score/RNA/TCR-gene strategies when they pass because gdT cells are rare and validation labels are partly TCR-defined.</p>",
        "</section>",
        f"<section><h2>Figures</h2><div class='figures'>{figure_html}</div></section>",
        "<section><h2>Held-Out Validation Metrics</h2>",
        dataframe_to_html(validation_df, max_rows=20),
        "</section>",
        "<section><h2>Acceptance Against Baseline</h2>",
        dataframe_to_html(acceptance_df, max_rows=20),
        "</section>",
        "<section><h2>Prevalence-Aware Imbalance Check</h2>",
        "<p>Because gdT cells are rare in deployment, PPV is shown under low-prevalence scenarios. Conservative PPV uses Wilson 95% lower bounds for recall and specificity.</p>",
        dataframe_to_html(selected_prevalence_display, max_rows=20),
        "</section>",
        "<section><h2>Feature Design</h2>",
        f"<p>Score features: <code>{len(spec.score_feature_cols)}</code>; TCR metadata: <code>{len(spec.metadata_feature_cols)}</code>; TCR gene modules: <code>{len(spec.module_feature_cols)}</code>; explicit expression genes: <code>{spec.expression_gene_count}</code> / {MAX_EXPRESSION_GENES} cap; non-TCR DEG genes: <code>{len(spec.deg_feature_cols)}</code>.</p>",
        "<p>Explicit TCR genes prioritize TRD constant/V/J/D genes. Selected TRG genes are retained but lower priority because TRG expression is less specific than TRD; TRAC/TRBC controls are included.</p>",
        "<p>FOXP3, CD4, and low pan-CD3 evidence are expression-level penalty signals in the tiered penalty strategy.</p>",
        dataframe_to_html(top_features),
        "</section>",
        "<section><h2>Selected Model Strata</h2><h3>By Source</h3>",
        dataframe_to_html(source_metrics, max_rows=40),
        "<h3>By Tissue</h3>",
        dataframe_to_html(tissue_metrics, max_rows=40),
        "</section>",
        "<section><h2>TRAB-Only Dataset False Positives</h2>",
        "<p>False positives are counted on held-out validation cells from sources with TRAB metadata but no gdTCR metadata.</p>",
        dataframe_to_html(trab_fp_overall, max_rows=20),
        dataframe_to_html(trab_fp_by_source, max_rows=40),
        "</section>",
        "<section><h2>GDTlung Quality Sensitivity</h2>",
        "<p>GDTlung2023july_7p is sorted gdT but flagged as suboptimal library quality; this table shows validation metrics with and without that source.</p>",
        dataframe_to_html(quality_sensitivity_df, max_rows=30),
        "</section>",
        "<section><h2>Full Dataset Application</h2>",
        dataframe_to_html(full_overall if not full_overall.empty else pd.DataFrame({"status": [full_block[0]]})),
        dataframe_to_html(full_source.head(30) if not full_source.empty else full_source),
        "</section>",
        "</main></body></html>",
    ]
    REPORT_HTML.write_text("\n".join(html_parts), encoding="utf-8")


def main() -> None:
    args = parse_args()
    setup_logging()
    input_h5ad = args.input_h5ad.resolve()
    stat_before = input_h5ad.stat()
    with h5py.File(input_h5ad, "r") as handle:
        n_obs = int(handle["obs"]["_index"].shape[0])
        logging.info("Reading labels and metadata for %s cells", f"{n_obs:,}")
        obs = build_obs_metadata(handle)
        source_mode = build_source_tcr_mode_table(obs)
        logging.info(
            "Gold/silver labels: gdT_gold=%s abT_gold=%s gdT_silver=%s",
            f"{int((obs.class_code == 2).sum()):,}",
            f"{int((obs.class_code == 1).sum()):,}",
            f"{int((obs.class_code == 3).sum()):,}",
        )
        spec = build_feature_spec(handle, args.max_deg_features, args.max_expression_genes)
        logging.info("Classifier feature count: %s", len(spec.feature_names))
        train_idx, tune_idx, validation_idx, split_summary = make_gold_splits(obs, args.seed)
        primary_eval_rows = np.concatenate([train_idx, tune_idx, validation_idx])
        silver_idx = np.flatnonzero(obs.class_code == 3).astype(np.int64)
        eval_rows = np.concatenate([primary_eval_rows, silver_idx])
        y_eval = (obs.class_code[eval_rows] == 2).astype(np.int8)
        pos_train = local_positions(eval_rows, train_idx)
        pos_tune = local_positions(eval_rows, tune_idx)
        pos_validation = local_positions(eval_rows, validation_idx)
        pos_silver = local_positions(eval_rows, silver_idx) if silver_idx.size else np.array([], dtype=np.int64)

        logging.info("Extracting feature matrix for gold+silver cells: %s rows", f"{eval_rows.size:,}")
        x_eval = extract_feature_matrix(handle, eval_rows, obs, spec, label="gold_silver")

        y_tune = (obs.class_code[tune_idx] == 2).astype(np.int8)
        y_validation = (obs.class_code[validation_idx] == 2).astype(np.int8)
        baseline_results = build_baseline_results(
            obs,
            tune_idx,
            validation_idx,
            silver_idx,
            y_tune,
            y_validation,
            args.min_specificity,
        )
        candidate_results = train_candidate_models(
            x_eval,
            y_eval,
            pos_train,
            pos_tune,
            pos_validation,
            pos_silver,
            eval_rows,
            obs,
            spec,
            args.max_negative_train,
            args.min_specificity,
            args.seed,
        )
        all_results = baseline_results + candidate_results

        tune_df = results_to_frame(all_results, "tune")
        validation_df = results_to_frame(all_results, "validation")
        validation_df, accepted_model, best_baseline = add_acceptance_columns(
            validation_df,
            all_results,
            y_validation,
            args.min_specificity,
            args.bootstrap_reps,
            args.seed,
        )
        acceptance_df = pd.read_csv(CLASSIFIER_TABLE_DIR / "classifier_acceptance_vs_best_baseline.csv")
        trab_fp_overall, trab_fp_by_source = write_trab_only_false_positive_tables(
            all_results, validation_idx, y_validation, obs, source_mode
        )
        quality_sensitivity_df = write_quality_sensitivity_metrics(all_results, validation_idx, y_validation, obs)
        prevalence_df = write_prevalence_aware_metrics(all_results, validation_idx, y_validation, obs)
        tune_df.to_csv(CLASSIFIER_TABLE_DIR / "classifier_metrics_tune.csv", index=False)
        validation_df.to_csv(CLASSIFIER_TABLE_DIR / "classifier_metrics_validation.csv", index=False)

        selected_result = accepted_model
        if selected_result is None:
            selected_result = max(candidate_results, key=model_selection_key)
            logging.warning("No model passed acceptance; reporting best candidate but not saving production artifact.")
            model_artifact = None
        else:
            logging.info("Accepted selected model: %s", selected_result.model_name)
            model_artifact = save_model_artifact(selected_result, spec)

        source_metrics, tissue_metrics = write_group_metrics(selected_result, y_validation, validation_idx, obs)
        full_overall = pd.DataFrame()
        full_source = pd.DataFrame()
        full_tissue = pd.DataFrame()
        if accepted_model is not None and not args.no_full_apply:
            logging.info("Applying accepted model to the full plus6 dataset")
            full_overall, full_source, full_tissue = apply_model_to_full_dataset(
                handle,
                accepted_model,
                obs,
                spec,
                source_mode,
                chunk_size=FULL_APPLY_CHUNK,
            )

    stat_after = input_h5ad.stat()
    if (stat_before.st_size != stat_after.st_size) or (stat_before.st_mtime_ns != stat_after.st_mtime_ns):
        raise RuntimeError("Input H5AD changed during read-only classifier analysis.")

    roc_path, pr_path = plot_roc_pr(all_results, y_validation)
    confusion_path = plot_confusion(selected_result, y_validation)
    distribution_path = plot_score_distribution(all_results, y_validation)
    importance_path = plot_feature_importance(selected_result, spec)
    copied_assets = copy_assets([roc_path, pr_path, confusion_path, distribution_path, importance_path])
    write_reports(
        input_h5ad,
        n_obs,
        obs,
        split_summary,
        tune_df,
        validation_df,
        acceptance_df,
        accepted_model,
        best_baseline,
        source_metrics,
        tissue_metrics,
        trab_fp_overall,
        trab_fp_by_source,
        quality_sensitivity_df,
        prevalence_df,
        full_overall,
        full_source,
        spec,
        copied_assets,
        model_artifact,
    )
    render_pdf(args.no_pdf)

    SUMMARY_JSON.write_text(
        json.dumps(
            json_ready(
                {
                    "input_h5ad": str(input_h5ad),
                    "n_obs": n_obs,
                    "accepted_model": None if accepted_model is None else accepted_model.model_name,
                    "reported_selected_model": selected_result.model_name,
                    "best_baseline": best_baseline.model_name,
                    "validation_metrics_path": str(CLASSIFIER_TABLE_DIR / "classifier_metrics_validation.csv"),
                    "acceptance_path": str(CLASSIFIER_TABLE_DIR / "classifier_acceptance_vs_best_baseline.csv"),
                    "report_md": str(REPORT_MD),
                    "report_html": str(REPORT_HTML),
                    "report_pdf": str(REPORT_PDF),
                    "model_artifact": None if model_artifact is None else str(model_artifact),
                }
            ),
            indent=2,
        ),
        encoding="utf-8",
    )
    logging.info("Saved classifier report: %s", REPORT_MD)
    logging.info("Saved classifier HTML: %s", REPORT_HTML)
    if REPORT_PDF.exists():
        logging.info("Saved classifier PDF: %s", REPORT_PDF)


if __name__ == "__main__":
    main()
