#!/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python
"""Compare the four frozen gdTAI operating profiles on extension cohorts.

The workflow has two intentionally separate stages. ``select`` excludes sealed
holdouts and writes a checksum-pinned decision. ``holdout`` verifies that
decision before opening the sealed cohorts and never reselects a runner-up.
Input H5AD files are opened read-only and must contain raw count-like CSR X.
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
    _TNK_PROJECT_ROOT / "workflows" / "gdtai",
):
    _tnk_value = str(_tnk_path)
    if _tnk_value not in _tnk_sys.path:
        _tnk_sys.path.insert(0, _tnk_value)

import argparse
import gzip
import hashlib
import html
import json
import logging
import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

os.environ.setdefault("MPLCONFIGDIR", "/tmp/gdtai_matplotlib_cache")

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    f1_score,
    matthews_corrcoef,
    precision_score,
    recall_score,
    roc_auc_score,
)
from scipy.stats import binomtest

from tnk_atlas.model_io import load_trusted_pickle
from tnk_atlas.provenance import sha256_file
from run_gdt_prediction_package_evaluation import read_obs_column, read_string_dataset
from run_gdtai_v3_trdc_nk_guard_classifier import (
    FeatureSpec,
    append_engineered_features,
    extract_gene_features,
    h5ad_shape,
)


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = PROJECT_ROOT / "configs/models/gdtai/extension_evaluation.json"
DEFAULT_INPUT_MANIFEST = PROJECT_ROOT / "data/interim/extension_intake/evaluation_inputs.csv"
MODEL_REGISTRY = PROJECT_ROOT / "configs/models/gdtai/model_registry.csv"
OUTPUT_ROOT = PROJECT_ROOT / "Integrated_dataset"
PREFIX = "gdtai_extension_evaluation"
TABLE_DIR = OUTPUT_ROOT / "tables/gdT_prediction" / PREFIX
FIGURE_DIR = OUTPUT_ROOT / "figures/gdT_prediction" / PREFIX
LOG_DIR = OUTPUT_ROOT / "logs/gdT_prediction" / PREFIX
STATIC_DIR = PROJECT_ROOT / "gdT_prediction" / PREFIX
DECISION_PATH = LOG_DIR / "selection_decision.json"
TARGET_SUM = 10_000.0

CANONICAL_OBS_COLUMNS = (
    "source_gse_id",
    "dataset_id",
    "sample_id",
    "library_id",
    "donor_id",
    "patient_id",
    "donor",
    "tissue_corrected",
    "tissue",
    "evaluation_annotation",
    "simple_annotation_plus6",
    "author_annotation",
    "cell_type",
    "annot_final",
    "major_cluster",
    "subset",
    "group",
    "evaluation_lineage",
    "ground_truth_class",
    "has_TRA_TRB_paired",
    "has_TRG_TRD_paired",
    "has_any_ab_tcr",
    "has_any_gd_tcr",
    "has_TRA",
    "has_TRB",
    "has_TRG",
    "has_TRD",
    "TRA_cdr3",
    "TRB_cdr3",
    "TRG_cdr3",
    "TRD_cdr3",
    "is_doublet",
    "gdT_silver",
    "Sorted_gdT",
    "targeted_cohort",
    "sample",
    "batch",
)


@dataclass(frozen=True)
class Profile:
    profile_id: str
    model_id: str
    mode: str
    artifact: Path
    sha256: str
    payload: dict[str, Any]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("select", "holdout"), required=True)
    parser.add_argument("--input-manifest", type=Path, default=DEFAULT_INPUT_MANIFEST)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--decision", type=Path, default=DECISION_PATH)
    parser.add_argument("--chunk-size", type=int, default=50_000)
    parser.add_argument("--bootstrap-iterations", type=int, default=None)
    parser.add_argument("--no-pdf", action="store_true")
    return parser.parse_args()


def stable_json_bytes(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("utf-8")


def sha256_json(value: Any) -> str:
    return hashlib.sha256(stable_json_bytes(value)).hexdigest()


def as_bool(values: Iterable[Any]) -> np.ndarray:
    series = pd.Series(values, copy=False).astype("string").fillna("").str.strip().str.lower()
    return series.isin({"true", "1", "yes", "y", "t"}).to_numpy(dtype=bool)


def clean_strings(values: Iterable[Any], default: str = "") -> np.ndarray:
    return pd.Series(values, copy=False).astype("string").fillna(default).astype(str).to_numpy(dtype=object)


def chain_present(obs: dict[str, np.ndarray], chain: str, n: int) -> np.ndarray:
    """Return productive-chain evidence from canonical flags or CDR3 metadata."""
    flag = as_bool(obs.get(f"has_{chain}", np.zeros(n, dtype=bool)))
    cdr3 = clean_strings(obs.get(f"{chain}_cdr3", np.full(n, "", dtype=object)))
    text = pd.Series(cdr3, copy=False).astype(str).str.strip()
    missing = text.str.lower().isin({"", "na", "nan", "none", "null", "<na>"}).to_numpy()
    return flag | ~missing


def read_obs(handle: h5py.File, column: str, default: Any = "") -> np.ndarray:
    n_obs = n_obs_from_handle(handle)
    if column not in handle["obs"]:
        return np.full(n_obs, default, dtype=object)
    return np.asarray(read_obs_column(handle, column))


def n_obs_from_handle(handle: h5py.File) -> int:
    obs = handle["obs"]
    index_key = obs.attrs.get("_index", "_index")
    if isinstance(index_key, bytes):
        index_key = index_key.decode("utf-8")
    if str(index_key) in obs:
        return int(obs[str(index_key)].shape[0])
    if "_index" in obs:
        return int(obs["_index"].shape[0])
    return h5ad_shape(handle, "X")[0]


def obs_index(handle: h5py.File) -> np.ndarray:
    obs = handle["obs"]
    index_key = obs.attrs.get("_index", "_index")
    if isinstance(index_key, bytes):
        index_key = index_key.decode("utf-8")
    for key in (str(index_key), "_index", "cell_id", "barcode"):
        if key in obs:
            obj = obs[key]
            if isinstance(obj, h5py.Dataset):
                return read_string_dataset(obj)
            return clean_strings(read_obs_column(handle, key))
    return np.asarray([str(i) for i in range(n_obs_from_handle(handle))], dtype=object)


def matrix_encoding(handle: h5py.File, matrix_path: str = "X") -> str:
    obj: h5py.Group | h5py.Dataset = handle[matrix_path]
    if isinstance(obj, h5py.Dataset):
        return "dense"
    value = obj.attrs.get("encoding-type", "")
    if isinstance(value, bytes):
        value = value.decode("utf-8")
    if value:
        return str(value)
    if {"data", "indices", "indptr"}.issubset(obj.keys()):
        return "csr_matrix"
    return "unknown"


def validate_count_matrix(handle: h5py.File, matrix_path: str = "X") -> dict[str, Any]:
    encoding = matrix_encoding(handle, matrix_path)
    if encoding not in {"csr_matrix", "dense"}:
        raise TypeError(f"{matrix_path} must be CSR or dense raw counts; found {encoding!r}")
    obj = handle[matrix_path]
    values = obj["data"] if isinstance(obj, h5py.Group) else obj
    size = int(np.prod(values.shape))
    if size:
        take = min(size, 100_000)
        sample = np.asarray(values[:take], dtype=np.float64).ravel()
        if np.any(~np.isfinite(sample)) or np.any(sample < 0):
            raise ValueError("Expression matrix contains non-finite or negative sampled values.")
        if np.any(np.abs(sample - np.rint(sample)) > 1e-6):
            raise ValueError("Expression matrix is not raw integer-like counts.")
    return {"matrix_path": matrix_path, "encoding": encoding, "shape": list(h5ad_shape(handle, matrix_path))}


def load_profiles(config: dict[str, Any]) -> list[Profile]:
    registry_df = pd.read_csv(MODEL_REGISTRY, dtype=str).fillna("")
    registry = {str(row["model_id"]): row for _, row in registry_df.iterrows()}
    out: list[Profile] = []
    for item in config["profiles"]:
        model_id = str(item["model_id"])
        if model_id not in registry:
            raise KeyError(f"Unknown model_id in evaluation config: {model_id}")
        row = registry[model_id]
        artifact = PROJECT_ROOT / str(row["artifact_path"])
        observed = sha256_file(artifact)
        accepted = {str(row["sha256"])}
        if str(row.get("allowed_workspace_sha256", "")):
            accepted.add(str(row["allowed_workspace_sha256"]))
        if observed not in accepted:
            raise RuntimeError(f"Checksum mismatch for {model_id}: {observed}")
        payload = load_trusted_pickle(artifact, PROJECT_ROOT)
        out.append(
            Profile(
                profile_id=str(item["profile_id"]),
                model_id=model_id,
                mode=str(item["mode"]),
                artifact=artifact,
                sha256=observed,
                payload=payload,
            )
        )
    return out


def feature_schema(profiles: list[Profile]) -> tuple[list[str], dict[str, int]]:
    genes: list[str] = []
    for profile in profiles:
        payload = profile.payload["base_model"] if profile.model_id == "gdtai_v2" else profile.payload
        for gene in payload["gene_names"]:
            gene = str(gene)
            if gene not in genes:
                genes.append(gene)
    return genes, {gene: i for i, gene in enumerate(genes)}


def build_union_spec(handle: h5py.File, genes: list[str]) -> tuple[FeatureSpec, pd.DataFrame]:
    var = handle["var"]
    index_key = var.attrs.get("_index", "_index")
    if isinstance(index_key, bytes):
        index_key = index_key.decode("utf-8")
    if str(index_key) not in var:
        raise KeyError("H5AD var index is missing.")
    var_names = [str(x) for x in read_string_dataset(var[str(index_key)])]
    lookup = {gene: i for i, gene in enumerate(var_names)}
    available_genes = [gene for gene in genes if gene in lookup]
    availability = pd.DataFrame(
        {
            "gene": genes,
            "available_in_h5ad": [gene in lookup for gene in genes],
            "h5ad_gene_index": [lookup.get(gene, "") for gene in genes],
        }
    )
    engineered = [
        "any_TRDV", "any_TRDJ", "any_TRG", "any_ab_TCR_gene", "TRDC_only",
        "TRDC_plus_TRDV", "TRDC_plus_TRDJ", "CD3_score", "NK_score",
        "gdT_TCR_score", "abT_TCR_score", "NK_minus_CD3_score", "TRDC_log1p",
        "TRDV_score", "TRDJ_score", "TRG_score",
    ]
    spec = FeatureSpec(
        gene_names=available_genes,
        gene_indices=np.asarray([lookup[g] for g in available_genes], dtype=np.int32),
        gene_feature_names=[f"{g}_log1p_cp10k" for g in available_genes],
        engineered_feature_names=engineered,
        model_feature_names=[*[f"{g}_log1p_cp10k" for g in available_genes], *engineered],
        gene_to_col={gene: i for i, gene in enumerate(available_genes)},
        engineered_to_col={name: len(available_genes) + i for i, name in enumerate(engineered)},
    )
    return spec, availability


def remap_gene_matrix(x_union: np.ndarray, union_spec: FeatureSpec, target_genes: list[str]) -> np.ndarray:
    out = np.zeros((x_union.shape[0], len(target_genes)), dtype=np.float32)
    for target_col, gene in enumerate(target_genes):
        source_col = union_spec.gene_to_col.get(str(gene))
        if source_col is not None:
            out[:, target_col] = x_union[:, source_col]
    return out


def normalize_annotation(obs: dict[str, np.ndarray], x_union: np.ndarray, spec: FeatureSpec) -> tuple[np.ndarray, np.ndarray]:
    n = x_union.shape[0]
    source = np.full(n, "derived_expression", dtype=object)
    annotation = np.full(n, "OTHER", dtype=object)
    for column in (
        "evaluation_annotation",
        "simple_annotation_plus6",
        "author_annotation",
        "annot_final",
        "cell_type",
        "subset",
        "major_cluster",
        "group",
    ):
        values = clean_strings(obs.get(column, np.full(n, "", dtype=object)))
        present = np.char.str_len(values.astype(str)) > 0
        if not present.any():
            continue
        upper = pd.Series(values).str.upper().to_numpy(dtype=object)
        mapped = np.full(n, "OTHER", dtype=object)
        mapped[pd.Series(upper).str.contains("TREG|REGULATORY", regex=True).to_numpy()] = "TREG"
        mapped[pd.Series(upper).str.contains("GDT|GD T|GAMMA", regex=True).to_numpy()] = "GDT_CELL"
        mapped[pd.Series(upper).str.contains("CD8|CYTOTOXIC", regex=True).to_numpy()] = "CD8_T"
        mapped[pd.Series(upper).str.contains("CD4|HELPER|TH17|TFH", regex=True).to_numpy()] = "CD4_T"
        nk = pd.Series(upper).str.contains(r"(?:^|[^A-Z])NK(?:[^A-Z]|$)", regex=True).to_numpy()
        nkt = pd.Series(upper).str.contains("NKT", regex=False).to_numpy()
        mapped[nk & ~nkt] = "NK_CELL"
        recognized = mapped != "OTHER"
        assign = (source == "derived_expression") & present & recognized
        annotation[assign] = mapped[assign]
        source[assign] = column

    def gene(name: str) -> np.ndarray:
        col = spec.gene_to_col.get(name)
        return np.zeros(n, dtype=np.float32) if col is None else x_union[:, col]

    missing = source == "derived_expression"
    cd3 = gene("CD3D") + gene("CD3E") + gene("CD3G")
    nk_score = gene("TYROBP") + gene("FCER1G") + gene("KLRD1") + gene("NCAM1")
    strict_nk = missing & (cd3 < 0.25) & (nk_score >= 1.5)
    annotation[strict_nk] = "NK_CELL"
    annotation[missing & ~strict_nk & ((gene("CD8A") + gene("CD8B")) > 0.5)] = "CD8_T"
    annotation[missing & ~strict_nk & (gene("FOXP3") > 0.25)] = "TREG"
    annotation[missing & ~strict_nk & (gene("CD4") > 0.75)] = "CD4_T"
    return annotation, source


def v2_thresholds(profile: Profile, annotation: np.ndarray) -> np.ndarray:
    mode = profile.payload["operating_modes"][profile.mode]
    if profile.mode == "high_f1":
        return np.full(annotation.size, float(mode["threshold"]), dtype=np.float32)
    rules = mode["annotation_thresholds"]
    out = np.full(annotation.size, float(rules["other_threshold"]), dtype=np.float32)
    mapping = {
        "GDT_CELL": "gdt_threshold", "CD8_T": "cd8_threshold", "CD4_T": "cd4_threshold",
        "TREG": "treg_threshold", "NK_CELL": "nk_threshold",
    }
    for label, key in mapping.items():
        value = rules[key]
        out[annotation == label] = np.inf if str(value) == "disabled" else float(value)
    return out


def predict_profiles(
    profiles: list[Profile],
    x_union: np.ndarray,
    union_spec: FeatureSpec,
    annotation: np.ndarray,
) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    out: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for profile in profiles:
        if profile.model_id == "gdtai_v2":
            base = profile.payload["base_model"]
            x = remap_gene_matrix(x_union, union_spec, [str(g) for g in base["gene_names"]])
            score = base["model_object"].predict_proba(x)[:, 1].astype(np.float32)
            threshold = v2_thresholds(profile, annotation)
        else:
            payload = profile.payload
            genes = [str(g) for g in payload["gene_names"]]
            x_gene = remap_gene_matrix(x_union, union_spec, genes)
            engineered = [str(x) for x in payload["engineered_feature_names"]]
            target_spec = FeatureSpec(
                gene_names=genes,
                gene_indices=np.arange(len(genes), dtype=np.int32),
                gene_feature_names=[f"{g}_log1p_cp10k" for g in genes],
                engineered_feature_names=engineered,
                model_feature_names=[str(x) for x in payload["feature_names"]],
                gene_to_col={gene: i for i, gene in enumerate(genes)},
                engineered_to_col={name: len(genes) + i for i, name in enumerate(engineered)},
            )
            x = append_engineered_features(x_gene, target_spec)
            score = payload["model_object"].predict_proba(x)[:, 1].astype(np.float32)
            threshold = np.full(score.size, float(payload["threshold"]), dtype=np.float32)
        out[profile.profile_id] = (score, score >= threshold)
    return out


def truth_frame(
    obs: dict[str, np.ndarray],
    annotation: np.ndarray,
    annotation_source: np.ndarray | None = None,
) -> pd.DataFrame:
    n = annotation.size
    has_tra = chain_present(obs, "TRA", n)
    has_trb = chain_present(obs, "TRB", n)
    has_trg = chain_present(obs, "TRG", n)
    has_trd = chain_present(obs, "TRD", n)
    has_ab = as_bool(obs.get("has_any_ab_tcr", np.zeros(n, dtype=bool))) | has_tra | has_trb
    paired_ab = as_bool(obs.get("has_TRA_TRB_paired", np.zeros(n, dtype=bool))) | (has_tra & has_trb)
    has_gd = as_bool(obs.get("has_any_gd_tcr", np.zeros(n, dtype=bool))) | has_trg | has_trd
    paired_gd = as_bool(obs.get("has_TRG_TRD_paired", np.zeros(n, dtype=bool))) | (has_trg & has_trd)
    doublet = as_bool(obs.get("is_doublet", np.zeros(n, dtype=bool)))
    silver = derive_silver_mask(obs, n_rows=n)
    sorted_gdt = as_bool(obs.get("Sorted_gdT", np.zeros(n, dtype=bool)))
    explicit = clean_strings(obs.get("ground_truth_class", np.full(n, "", dtype=object)))
    gd_gold = ((paired_gd & ~has_ab) | (explicit == "gdT_gold") | sorted_gdt) & ~doublet
    ab_gold = ((has_ab & ~has_gd) | (explicit == "abT_gold")) & ~doublet & ~gd_gold
    lineage = clean_strings(obs.get("evaluation_lineage", np.full(n, "", dtype=object)))
    strict_nk = ((lineage == "strict_nk") | (annotation == "NK_CELL")) & ~has_ab & ~has_gd & ~doublet
    if annotation_source is None:
        author_nk = annotation == "NK_CELL"
    else:
        author_nk = (annotation == "NK_CELL") & (annotation_source != "derived_expression")
    label = np.full(n, "unlabeled_or_ambiguous", dtype=object)
    label[ab_gold] = "abT_gold"
    label[gd_gold] = "gdT_gold"
    label[silver & ~gd_gold & ~ab_gold] = "gdT_silver"
    label[strict_nk & ~gd_gold & ~ab_gold] = "strict_NK_negative"
    return pd.DataFrame(
        {
            "truth_class": label,
            "gdT_gold": gd_gold,
            "abT_gold": ab_gold,
            "paired_abT": paired_ab & ~has_gd,
            "strict_NK": strict_nk,
            "author_NK": author_nk,
            "gdT_silver": silver,
            "doublet": doublet,
        }
    )


def derive_silver_mask(obs: dict[str, np.ndarray], *, n_rows: int | None = None) -> np.ndarray:
    """Derive productive TRD-only silver positives within gdT-sequenced libraries."""
    n = n_rows if n_rows is not None else (len(next(iter(obs.values()))) if obs else 0)
    explicit = as_bool(obs.get("gdT_silver", np.zeros(n, dtype=bool)))
    has_tra = chain_present(obs, "TRA", n)
    has_trb = chain_present(obs, "TRB", n)
    has_trg = chain_present(obs, "TRG", n)
    has_trd = chain_present(obs, "TRD", n)
    has_ab = as_bool(obs.get("has_any_ab_tcr", np.zeros(n, dtype=bool))) | has_tra | has_trb
    has_gd = as_bool(obs.get("has_any_gd_tcr", np.zeros(n, dtype=bool))) | has_trg | has_trd
    paired_gd = as_bool(obs.get("has_TRG_TRD_paired", np.zeros(n, dtype=bool))) | (has_trg & has_trd)
    library = first_nonempty(obs, ("library_id", "sample_id", "sample", "batch"))
    valid_library = pd.Series(library, copy=False).astype(str).str.len().to_numpy() > 0
    if valid_library.any():
        frame = pd.DataFrame({"library": library, "has_gd": has_gd})
        library_has_gd = frame.groupby("library", sort=False)["has_gd"].transform("any").to_numpy(
            dtype=bool,
            copy=True,
        )
        library_has_gd &= valid_library
    else:
        library_has_gd = np.zeros(n, dtype=bool)
    derived = library_has_gd & has_trd & ~paired_gd & ~has_ab
    return explicit | derived


def author_nk_mask(obs: dict[str, np.ndarray]) -> np.ndarray:
    """Identify explicit author-level NK labels without treating NKT as NK."""
    n = len(next(iter(obs.values()))) if obs else 0
    mask = np.zeros(n, dtype=bool)
    for column in (
        "evaluation_annotation",
        "simple_annotation_plus6",
        "author_annotation",
        "annot_final",
        "cell_type",
        "subset",
        "major_cluster",
        "group",
    ):
        values = clean_strings(obs.get(column, np.full(n, "", dtype=object)))
        upper = pd.Series(values, copy=False).astype(str).str.upper()
        is_nk = upper.str.contains(r"(?:^|[^A-Z])NK(?:[^A-Z]|$)", regex=True).to_numpy()
        is_nkt = upper.str.contains("NKT", regex=False).to_numpy()
        mask |= is_nk & ~is_nkt
    return mask


def infer_one_dataset(
    row: pd.Series,
    profiles: list[Profile],
    config: dict[str, Any],
    chunk_size: int,
    output_path: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    h5ad_path = Path(str(row["h5ad_path"])).expanduser()
    if not h5ad_path.is_absolute():
        h5ad_path = PROJECT_ROOT / h5ad_path
    dataset_id = str(row["dataset_id"])
    union_genes, _ = feature_schema(profiles)
    summary_rows: list[dict[str, Any]] = []
    availability_rows: list[pd.DataFrame] = []
    wrote = False
    with h5py.File(h5ad_path, "r") as handle, gzip.open(output_path, "wt", encoding="utf-8", newline="") as sink:
        matrix_path = str(row.get("matrix_path", "X") or "X").strip()
        matrix_info = validate_count_matrix(handle, matrix_path)
        spec, availability = build_union_spec(handle, union_genes)
        availability.insert(0, "dataset_id", dataset_id)
        availability_rows.append(availability)
        available_fraction = float(availability["available_in_h5ad"].mean())
        missing_critical = sorted(set(config["critical_genes"]) - set(spec.gene_names))
        schema_primary_eligible = available_fraction >= float(config["minimum_model_gene_fraction"])
        schema_warning = available_fraction < float(config["model_gene_fraction_warning_below"])
        if missing_critical:
            raise RuntimeError(
                f"{dataset_id} fails model-gene compatibility: fraction={available_fraction:.3f}, "
                f"missing critical={missing_critical}"
            )
        n_obs = n_obs_from_handle(handle)
        ids = obs_index(handle)
        obs_all = {column: read_obs(handle, column) for column in CANONICAL_OBS_COLUMNS}
        obs_all["gdT_silver"] = derive_silver_mask(obs_all)
        row_filter = str(row.get("row_filter", "all") or "all").strip().lower()
        selected_rows = select_input_rows(
            obs_all,
            row_filter,
            include_sources=parse_source_list(row.get("include_source_gse_id", "")),
            exclude_sources=parse_source_list(row.get("exclude_source_gse_id", "")),
        )
        if selected_rows.size == 0:
            raise RuntimeError(f"{dataset_id} has no cells after row_filter={row_filter!r}.")
        for chunk_start in range(0, selected_rows.size, chunk_size):
            chunk_end = min(chunk_start + chunk_size, selected_rows.size)
            rows = selected_rows[chunk_start:chunk_end]
            x_union, row_sum, n_detected = extract_gene_features(
                handle,
                matrix_path,
                rows,
                spec,
                label=f"{dataset_id}_{chunk_start}_{chunk_end}",
            )
            obs = {key: value[rows] for key, value in obs_all.items()}
            annotation, annotation_source = normalize_annotation(obs, x_union, spec)
            truth = truth_frame(obs, annotation, annotation_source)
            predictions = predict_profiles(profiles, x_union, spec, annotation)
            sample_id = first_nonempty(obs, ("sample_id", "sample"))
            library_id = first_nonempty(obs, ("library_id", "batch", "sample"))
            donor_id = first_nonempty(obs, ("donor_id", "patient_id", "donor"))
            data: dict[str, Any] = {
                "dataset_id": dataset_id,
                "cell_id": ids[rows],
                "sample_id": sample_id,
                "library_id": library_id,
                "donor_id": donor_id,
                "evaluation_annotation": annotation,
                "evaluation_annotation_source": annotation_source,
                "row_sum_counts": row_sum,
                "n_detected_genes": n_detected,
                "model_schema_primary_eligible": schema_primary_eligible,
                "model_schema_warning": schema_warning,
            }
            for column in truth.columns:
                data[column] = truth[column].to_numpy()
            for profile_id, (score, pred) in predictions.items():
                data[f"{profile_id}_score"] = score
                data[f"{profile_id}_pred"] = pred
            frame = pd.DataFrame(data)
            frame.to_csv(sink, index=False, header=not wrote)
            wrote = True
        summary_rows.append(
            {
                "dataset_id": dataset_id,
                "h5ad_path": str(h5ad_path),
                "n_cells_total": n_obs,
                "n_cells_evaluated": int(selected_rows.size),
                "row_filter": row_filter,
                "matrix_path": matrix_path,
                "matrix_encoding": matrix_info["encoding"],
                "model_gene_fraction": available_fraction,
                "model_schema_primary_eligible": schema_primary_eligible,
                "model_schema_warning": schema_warning,
                "missing_critical_genes": ";".join(missing_critical),
            }
        )
    return pd.DataFrame(summary_rows), pd.concat(availability_rows, ignore_index=True)


def parse_source_list(value: Any) -> set[str]:
    text = str(value or "").strip()
    if not text:
        return set()
    return {item.strip() for item in text.replace("|", ";").split(";") if item.strip()}


def first_nonempty(obs: dict[str, np.ndarray], columns: tuple[str, ...]) -> np.ndarray:
    n = len(next(iter(obs.values()))) if obs else 0
    out = np.full(n, "", dtype=object)
    for column in columns:
        values = clean_strings(obs.get(column, np.full(n, "", dtype=object)))
        mask = (pd.Series(out).astype(str).str.len() == 0).to_numpy() & (pd.Series(values).astype(str).str.len() > 0).to_numpy()
        out[mask] = values[mask]
    return out


def select_input_rows(
    obs: dict[str, np.ndarray],
    row_filter: str,
    *,
    include_sources: set[str] | None = None,
    exclude_sources: set[str] | None = None,
) -> np.ndarray:
    n = len(next(iter(obs.values()))) if obs else 0
    include_sources = include_sources or set()
    exclude_sources = exclude_sources or set()
    has_tra = chain_present(obs, "TRA", n)
    has_trb = chain_present(obs, "TRB", n)
    has_trg = chain_present(obs, "TRG", n)
    has_trd = chain_present(obs, "TRD", n)
    has_ab = as_bool(obs.get("has_any_ab_tcr", np.zeros(n, dtype=bool))) | has_tra | has_trb
    has_gd = as_bool(obs.get("has_any_gd_tcr", np.zeros(n, dtype=bool))) | has_trg | has_trd
    paired_gd = as_bool(obs.get("has_TRG_TRD_paired", np.zeros(n, dtype=bool))) | (has_trg & has_trd)
    sorted_gdt = as_bool(obs.get("Sorted_gdT", np.zeros(n, dtype=bool)))
    explicit = clean_strings(obs.get("ground_truth_class", np.full(n, "", dtype=object)))
    gd_gold = (paired_gd & ~has_ab) | sorted_gdt | (explicit == "gdT_gold")
    ab_gold = (has_ab & ~has_gd) | (explicit == "abT_gold")
    silver = derive_silver_mask(obs) | (explicit == "gdT_silver")
    if row_filter in {"", "all"}:
        mask = np.ones(n, dtype=bool)
    elif row_filter == "gold_only":
        mask = gd_gold
    elif row_filter == "primary_labeled":
        mask = gd_gold | ab_gold
    elif row_filter == "labeled_with_silver":
        mask = gd_gold | ab_gold | silver
    elif row_filter == "evaluation_labeled":
        mask = gd_gold | ab_gold | silver | author_nk_mask(obs)
    else:
        raise ValueError(f"Unsupported row_filter: {row_filter!r}")
    source = clean_strings(obs.get("source_gse_id", np.full(n, "", dtype=object)))
    if include_sources:
        mask &= np.isin(source, list(include_sources))
    if exclude_sources:
        mask &= ~np.isin(source, list(exclude_sources))
    return np.flatnonzero(mask).astype(np.int64)


def confusion_counts(y: np.ndarray, pred: np.ndarray) -> tuple[int, int, int, int]:
    y = np.asarray(y, dtype=bool)
    pred = np.asarray(pred, dtype=bool)
    tp = int(np.sum(y & pred))
    fp = int(np.sum(~y & pred))
    tn = int(np.sum(~y & ~pred))
    fn = int(np.sum(y & ~pred))
    return tp, fp, tn, fn


def metric_row(name: str, y: np.ndarray, score: np.ndarray, pred: np.ndarray) -> dict[str, Any]:
    tp, fp, tn, fn = confusion_counts(y, pred)
    has_both = np.unique(y).size == 2
    specificity = tn / (tn + fp) if (tn + fp) else np.nan
    return {
        "profile_id": name,
        "n_cells": int(y.size), "n_positive": int(y.sum()), "n_negative": int((~y).sum()),
        "predicted_positive": int(pred.sum()), "tp": tp, "fp": fp, "tn": tn, "fn": fn,
        "precision": precision_score(y, pred, zero_division=0),
        "recall": recall_score(y, pred, zero_division=0),
        "specificity": specificity,
        "fpr": 1.0 - specificity if np.isfinite(specificity) else np.nan,
        "f1": f1_score(y, pred, zero_division=0),
        "balanced_accuracy": balanced_accuracy_score(y, pred) if has_both else np.nan,
        "mcc": matthews_corrcoef(y, pred) if has_both else np.nan,
        "roc_auc": roc_auc_score(y, score) if has_both else np.nan,
        "pr_auc": average_precision_score(y, score) if has_both else np.nan,
    }


def evaluate_predictions(predictions: pd.DataFrame, profiles: list[Profile]) -> tuple[pd.DataFrame, pd.DataFrame]:
    metric_rows: list[dict[str, Any]] = []
    strata_rows: list[dict[str, Any]] = []
    schema_ok = predictions["model_schema_primary_eligible"].astype(bool).to_numpy()
    primary = predictions["truth_class"].isin(["gdT_gold", "abT_gold"]).to_numpy() & schema_ok
    y = predictions.loc[primary, "gdT_gold"].to_numpy(dtype=bool)
    for profile in profiles:
        score = predictions.loc[primary, f"{profile.profile_id}_score"].to_numpy(dtype=float)
        pred = predictions.loc[primary, f"{profile.profile_id}_pred"].to_numpy(dtype=bool)
        row = metric_row(profile.profile_id, y, score, pred)
        row["evaluation"] = "primary_gold"
        row.update(macro_f1_values(predictions.loc[primary].copy(), profile.profile_id))
        metric_rows.append(row)
        for dataset_id, group in predictions.groupby("dataset_id", sort=True):
            for stratum, mask in {
                "abT_gold": group["abT_gold"].to_numpy(dtype=bool),
                "paired_abT": group["paired_abT"].to_numpy(dtype=bool),
                "strict_NK": group["strict_NK"].to_numpy(dtype=bool),
                "author_NK": group["author_NK"].to_numpy(dtype=bool),
                "gdT_gold": group["gdT_gold"].to_numpy(dtype=bool),
                "gdT_silver": group["gdT_silver"].to_numpy(dtype=bool),
            }.items():
                n = int(mask.sum())
                if n == 0:
                    continue
                called = group.loc[mask, f"{profile.profile_id}_pred"].to_numpy(dtype=bool)
                strata_rows.append(
                    {
                        "profile_id": profile.profile_id, "dataset_id": dataset_id,
                        "stratum": stratum, "n_cells": n, "predicted_gdT": int(called.sum()),
                        "call_rate": float(called.mean()),
                        "model_schema_primary_eligible": bool(group["model_schema_primary_eligible"].astype(bool).all()),
                    }
                )
    return pd.DataFrame(metric_rows), pd.DataFrame(strata_rows)


def macro_f1_values(primary: pd.DataFrame, profile_id: str) -> dict[str, float]:
    values: dict[str, float] = {}
    for level, columns in {
        "dataset_macro_f1": ["dataset_id"],
        "donor_macro_f1": ["dataset_id", "donor_id"],
    }.items():
        group_scores: list[float] = []
        frame = primary.copy()
        if "donor_id" in columns:
            donor = frame["donor_id"].astype("string").fillna("").astype(str)
            sample = frame["sample_id"].astype("string").fillna("").astype(str)
            library = frame["library_id"].astype("string").fillna("").astype(str)
            fallback = sample.where(sample.ne(""), library)
            frame["donor_id"] = donor.where(donor.ne(""), fallback)
        for _, group in frame.groupby(columns, sort=True, dropna=False):
            y = group["gdT_gold"].to_numpy(dtype=bool)
            if np.unique(y).size < 2:
                continue
            pred = group[f"{profile_id}_pred"].to_numpy(dtype=bool)
            group_scores.append(float(f1_score(y, pred, zero_division=0)))
        values[level] = float(np.mean(group_scores)) if group_scores else np.nan
    candidates = [values["dataset_macro_f1"], values["donor_macro_f1"]]
    finite = [value for value in candidates if np.isfinite(value)]
    values["selection_macro_f1"] = float(np.mean(finite)) if finite else np.nan
    return values


def donor_cluster_bootstrap(
    predictions: pd.DataFrame,
    profiles: list[Profile],
    iterations: int,
    seed: int,
) -> pd.DataFrame:
    primary = predictions[
        predictions["truth_class"].isin(["gdT_gold", "abT_gold"])
        & predictions["model_schema_primary_eligible"].astype(bool)
    ].copy()
    if primary.empty or iterations <= 0:
        return pd.DataFrame()
    donor = primary["donor_id"].astype("string").fillna("").astype(str)
    sample = primary["sample_id"].astype("string").fillna("").astype(str)
    library = primary["library_id"].astype("string").fillna("").astype(str)
    fallback = sample.where(sample.ne(""), library)
    fallback = fallback.where(fallback.ne(""), primary.index.astype(str))
    primary["cluster_key"] = primary["dataset_id"].astype(str) + "|" + donor.where(donor.ne(""), fallback)
    keys = primary["cluster_key"].drop_duplicates().to_numpy(dtype=object)
    rng = np.random.default_rng(seed)
    rows: list[dict[str, Any]] = []
    for profile in profiles:
        records = []
        for _, group in primary.groupby("cluster_key", sort=True):
            y = group["gdT_gold"].to_numpy(dtype=bool)
            pred = group[f"{profile.profile_id}_pred"].to_numpy(dtype=bool)
            records.append(confusion_counts(y, pred))
        counts = np.asarray(records, dtype=np.int64)
        draws = rng.integers(0, len(keys), size=(iterations, len(keys)))
        sampled = counts[draws].sum(axis=1)
        tp, fp, tn, fn = (sampled[:, i].astype(float) for i in range(4))
        with np.errstate(divide="ignore", invalid="ignore"):
            values = {
                "precision": tp / (tp + fp),
                "recall": tp / (tp + fn),
                "specificity": tn / (tn + fp),
                "f1": 2 * tp / (2 * tp + fp + fn),
            }
        for metric, array in values.items():
            finite = array[np.isfinite(array)]
            rows.append(
                {
                    "profile_id": profile.profile_id,
                    "metric": metric,
                    "bootstrap_iterations": iterations,
                    "n_clusters": len(keys),
                    "ci_low": float(np.quantile(finite, 0.025)) if finite.size else np.nan,
                    "ci_high": float(np.quantile(finite, 0.975)) if finite.size else np.nan,
                    "bootstrap_median": float(np.median(finite)) if finite.size else np.nan,
                }
            )
    return pd.DataFrame(rows)


def paired_mcnemar(predictions: pd.DataFrame, profiles: list[Profile]) -> pd.DataFrame:
    primary = predictions[
        predictions["truth_class"].isin(["gdT_gold", "abT_gold"])
        & predictions["model_schema_primary_eligible"].astype(bool)
    ]
    y = primary["gdT_gold"].to_numpy(dtype=bool)
    rows: list[dict[str, Any]] = []
    for i, left in enumerate(profiles):
        left_error = primary[f"{left.profile_id}_pred"].to_numpy(dtype=bool) != y
        for right in profiles[i + 1 :]:
            right_error = primary[f"{right.profile_id}_pred"].to_numpy(dtype=bool) != y
            left_only = int(np.sum(left_error & ~right_error))
            right_only = int(np.sum(~left_error & right_error))
            discordant = left_only + right_only
            pvalue = float(binomtest(min(left_only, right_only), discordant, 0.5).pvalue) if discordant else 1.0
            rows.append(
                {
                    "profile_left": left.profile_id,
                    "profile_right": right.profile_id,
                    "left_error_right_correct": left_only,
                    "left_correct_right_error": right_only,
                    "discordant": discordant,
                    "exact_two_sided_pvalue": pvalue,
                }
            )
    return pd.DataFrame(rows)


def prevalence_aware_ppv(metrics: pd.DataFrame, prevalence: list[float]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for _, metric in metrics.iterrows():
        sensitivity = float(metric["recall"])
        specificity = float(metric["specificity"])
        for value in prevalence:
            denominator = sensitivity * value + (1.0 - specificity) * (1.0 - value)
            rows.append(
                {
                    "profile_id": metric["profile_id"],
                    "prevalence": float(value),
                    "expected_ppv": float(sensitivity * value / denominator) if denominator > 0 else np.nan,
                    "predicted_positive_per_million": float(1_000_000 * denominator),
                    "false_positive_per_million": float(1_000_000 * (1.0 - specificity) * (1.0 - value)),
                }
            )
    return pd.DataFrame(rows)


def select_profile(metrics: pd.DataFrame, strata: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    guard = config["selection_guardrails"]
    rows = []
    schema_ok = (
        strata["model_schema_primary_eligible"].astype(bool)
        if "model_schema_primary_eligible" in strata
        else pd.Series(True, index=strata.index)
    )
    for _, metric in metrics.loc[metrics["evaluation"] == "primary_gold"].iterrows():
        profile = str(metric["profile_id"])
        subset = strata[
            (strata["profile_id"] == profile)
            & schema_ok
        ]
        ab = subset[subset["stratum"] == "abT_gold"]
        nk = subset[subset["stratum"] == "strict_NK"]
        pooled_ab_fpr = float(ab["predicted_gdT"].sum() / ab["n_cells"].sum()) if ab["n_cells"].sum() else np.nan
        pooled_nk_fpr = float(nk["predicted_gdT"].sum() / nk["n_cells"].sum()) if nk["n_cells"].sum() else np.nan
        large = subset[(subset["stratum"].isin(["abT_gold", "strict_NK"])) & (subset["n_cells"] >= int(guard["large_cohort_minimum_negatives"]))]
        max_large = float(large["call_rate"].max()) if not large.empty else np.nan
        predicted = int(metric["predicted_positive"])
        fp_fraction = float(metric["fp"] / predicted) if predicted else 0.0
        checks = {
            "gold_recall_pass": float(metric["recall"]) >= float(guard["minimum_gold_recall"]),
            "pooled_abt_fpr_pass": np.isfinite(pooled_ab_fpr) and pooled_ab_fpr <= float(guard["maximum_pooled_abt_fpr"]),
            "pooled_strict_nk_fpr_pass": np.isfinite(pooled_nk_fpr) and pooled_nk_fpr <= float(guard["maximum_pooled_strict_nk_fpr"]),
            "large_cohort_fpr_pass": (not np.isfinite(max_large)) or max_large <= float(guard["maximum_large_cohort_fpr"]),
            "labeled_fp_fraction_pass": fp_fraction <= float(guard["maximum_labeled_fp_fraction"]),
        }
        rows.append(
            {
                "profile_id": profile, "primary_f1": float(metric["f1"]),
                "dataset_macro_f1": float(metric.get("dataset_macro_f1", np.nan)),
                "donor_macro_f1": float(metric.get("donor_macro_f1", np.nan)),
                "selection_macro_f1": float(metric.get("selection_macro_f1", np.nan)),
                "primary_precision": float(metric["precision"]), "gold_recall": float(metric["recall"]),
                "pooled_abt_fpr": pooled_ab_fpr, "pooled_strict_nk_fpr": pooled_nk_fpr,
                "maximum_large_cohort_fpr": max_large, "labeled_fp_fraction": fp_fraction,
                **checks, "eligible": all(checks.values()),
            }
        )
    ranking = pd.DataFrame(rows)
    ranking["ranking_f1"] = ranking["selection_macro_f1"].where(
        ranking["selection_macro_f1"].notna(), ranking["primary_f1"]
    )
    ranking = ranking.sort_values(["eligible", "ranking_f1", "primary_f1"], ascending=[False, False, False]).reset_index(drop=True)
    ranking["selected"] = False
    if bool(ranking["eligible"].any()):
        ranking.loc[ranking.index[0], "selected"] = True
    return ranking


def load_prediction_files(paths: list[Path]) -> pd.DataFrame:
    frames = [pd.read_csv(path, low_memory=False) for path in paths]
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def profile_identity(profiles: list[Profile]) -> dict[str, dict[str, str]]:
    return {p.profile_id: {"model_id": p.model_id, "mode": p.mode, "artifact": str(p.artifact), "sha256": p.sha256} for p in profiles}


def write_decision(
    path: Path,
    ranking: pd.DataFrame,
    profiles: list[Profile],
    config: dict[str, Any],
    config_path: Path,
    manifest_path: Path,
    input_rows: pd.DataFrame,
) -> dict[str, Any]:
    selected = ranking.loc[ranking["selected"], "profile_id"].astype(str).tolist()
    payload = {
        "schema_version": 1,
        "status": "selected" if selected else "no_profile_passed_guardrails",
        "selected_profile": selected[0] if selected else None,
        "profile_identity": profile_identity(profiles),
        "config_path": str(config_path.resolve()),
        "config_sha256": sha256_file(config_path),
        "input_manifest_path": str(manifest_path.resolve()),
        "input_manifest_sha256": sha256_file(manifest_path),
        "input_manifest_rows": input_rows[["dataset_id", "h5ad_path", "sealed_holdout"]].to_dict(orient="records"),
        "selection_ranking": ranking.replace({np.nan: None}).to_dict(orient="records"),
    }
    payload["decision_sha256"] = sha256_json(payload)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return payload


def verify_decision(
    path: Path,
    profiles: list[Profile],
    config_path: Path,
    manifest_path: Path,
) -> dict[str, Any]:
    decision = json.loads(path.read_text(encoding="utf-8"))
    digest = decision.pop("decision_sha256")
    if sha256_json(decision) != digest:
        raise RuntimeError("Selection decision digest is invalid; holdout remains sealed.")
    decision["decision_sha256"] = digest
    if decision.get("status") != "selected" or not decision.get("selected_profile"):
        raise RuntimeError("No eligible profile was selected; holdout remains sealed.")
    if decision.get("profile_identity") != profile_identity(profiles):
        raise RuntimeError("Model identity changed after selection; holdout remains sealed.")
    if decision.get("config_sha256") != sha256_file(config_path):
        raise RuntimeError("Evaluation config changed after selection; holdout remains sealed.")
    if decision.get("input_manifest_sha256") != sha256_file(manifest_path):
        raise RuntimeError("Input manifest changed after selection; holdout remains sealed.")
    return decision


def plot_false_positive_rates(strata: pd.DataFrame, path: Path, title: str) -> None:
    view = strata[strata["stratum"].isin(["abT_gold", "strict_NK"])].copy()
    if view.empty:
        return
    agg = view.groupby(["profile_id", "stratum"], as_index=False).agg(predicted_gdT=("predicted_gdT", "sum"), n_cells=("n_cells", "sum"))
    agg["rate"] = agg["predicted_gdT"] / agg["n_cells"]
    profiles = agg["profile_id"].drop_duplicates().tolist()
    x = np.arange(len(profiles))
    width = 0.36
    fig, ax = plt.subplots(figsize=(9, 5), constrained_layout=True)
    for offset, stratum, color in [(-width / 2, "abT_gold", "#287271"), (width / 2, "strict_NK", "#C75146")]:
        values = agg[agg["stratum"] == stratum].set_index("profile_id")["rate"].reindex(profiles).fillna(0)
        ax.bar(x + offset, values, width, label=stratum, color=color)
    ax.set_xticks(x, profiles, rotation=20, ha="right")
    ax.set_ylabel("False-positive rate")
    ax.set_title(title)
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.25)
    fig.savefig(path, dpi=220)
    plt.close(fig)


def dataframe_html(df: pd.DataFrame) -> str:
    return df.to_html(index=False, escape=True, border=0, classes="dataframe")


def render_report(
    stage: str,
    metrics: pd.DataFrame,
    strata: pd.DataFrame,
    ranking: pd.DataFrame | None,
    decision: dict[str, Any] | None,
    bootstrap: pd.DataFrame,
    mcnemar: pd.DataFrame,
    prevalence: pd.DataFrame,
    figure_path: Path,
) -> Path:
    STATIC_DIR.mkdir(parents=True, exist_ok=True)
    title = "gdTAI frozen-profile selection" if stage == "select" else "gdTAI sealed-holdout evaluation"
    decision_text = ""
    if decision:
        decision_text = f"<p><strong>Selected before holdout:</strong> <code>{html.escape(str(decision.get('selected_profile')))}</code></p>"
    rank_html = "" if ranking is None else f"<h2>Purity-constrained selection</h2>{dataframe_html(ranking)}"
    figure_html = ""
    if figure_path.exists():
        relative = Path(os.path.relpath(figure_path, STATIC_DIR)).as_posix()
        figure_html = f"<figure><img src='{html.escape(relative)}'><figcaption>False-positive rates by negative population.</figcaption></figure>"
    body = f"""<!doctype html><html><head><meta charset='utf-8'><title>{html.escape(title)}</title>
<style>
@page{{size:A4 landscape;margin:12mm}}body{{font-family:Arial,sans-serif;max-width:1400px;margin:24px auto;color:#172126;line-height:1.45}}
h1,h2{{letter-spacing:0}}table{{border-collapse:collapse;width:100%;font-size:8.5px;table-layout:fixed;page-break-inside:auto}}tr{{page-break-inside:avoid}}th,td{{border:1px solid #d4dadd;padding:3px 4px;text-align:right;white-space:normal;overflow-wrap:anywhere}}th:first-child,td:first-child{{text-align:left}}th{{background:#eef2f3}}code{{background:#eef2f3;padding:2px 4px}}img{{max-width:100%;max-height:155mm}}figure{{page-break-inside:avoid}}
</style></head><body><h1>{html.escape(title)}</h1>
<p>The four artifacts were checksum-pinned. Alpha-beta T cells and strict NK cells are evaluated as separate false-positive populations. Sealed holdouts cannot influence profile selection.</p>
{decision_text}<h2>Primary metrics</h2>{dataframe_html(metrics)}
{figure_html}<h2>Per-dataset and lineage audit</h2>{dataframe_html(strata)}{rank_html}
<h2>Donor-cluster bootstrap intervals</h2>{dataframe_html(bootstrap)}
<h2>Paired profile comparisons</h2>{dataframe_html(mcnemar)}
<h2>Prevalence-aware projections</h2>{dataframe_html(prevalence)}
</body></html>"""
    path = STATIC_DIR / ("selection.html" if stage == "select" else "index.html")
    path.write_text(body, encoding="utf-8")
    return path


def render_pdf(report_html: Path, output_pdf: Path) -> bool:
    chrome = shutil.which("google-chrome") or shutil.which("chromium") or shutil.which("chromium-browser")
    if chrome is None:
        logging.warning("Chrome is unavailable; PDF was not rendered.")
        return False
    output_pdf.parent.mkdir(parents=True, exist_ok=True)
    command = [
        chrome,
        "--headless",
        "--no-sandbox",
        "--disable-gpu",
        f"--print-to-pdf={output_pdf.resolve()}",
        report_html.resolve().as_uri(),
    ]
    result = subprocess.run(command, check=False, capture_output=True, text=True)
    if result.returncode != 0 or not output_pdf.exists() or output_pdf.stat().st_size == 0:
        logging.warning("PDF rendering failed: %s", result.stderr.strip())
        return False
    return True


def main() -> None:
    args = parse_args()
    for path in (TABLE_DIR, FIGURE_DIR, LOG_DIR, STATIC_DIR):
        path.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
    config = json.loads(args.config.read_text(encoding="utf-8"))
    profiles = load_profiles(config)
    manifest = pd.read_csv(args.input_manifest, dtype=str).fillna("")
    required = {"dataset_id", "h5ad_path", "sealed_holdout"}
    missing = required - set(manifest.columns)
    if missing:
        raise KeyError(f"Input manifest is missing columns: {sorted(missing)}")
    manifest["sealed_holdout"] = manifest["sealed_holdout"].str.lower().isin({"true", "1", "yes"})
    if args.stage == "select":
        rows = manifest.loc[~manifest["sealed_holdout"]].copy()
        decision = None
    else:
        decision = verify_decision(args.decision, profiles, args.config, args.input_manifest)
        rows = manifest.loc[manifest["sealed_holdout"]].copy()
    if rows.empty:
        raise RuntimeError(f"No datasets are available for stage {args.stage!r}.")

    prediction_paths: list[Path] = []
    dataset_summaries: list[pd.DataFrame] = []
    availability: list[pd.DataFrame] = []
    for _, row in rows.iterrows():
        dataset_id = str(row["dataset_id"])
        output_path = TABLE_DIR / f"{args.stage}_{dataset_id}_predictions.csv.gz"
        logging.info("Evaluating %s", dataset_id)
        summary, genes = infer_one_dataset(row, profiles, config, args.chunk_size, output_path)
        prediction_paths.append(output_path)
        dataset_summaries.append(summary)
        availability.append(genes)

    predictions = load_prediction_files(prediction_paths)
    metrics, strata = evaluate_predictions(predictions, profiles)
    bootstrap_iterations = args.bootstrap_iterations
    if bootstrap_iterations is None:
        bootstrap_iterations = int(config["bootstrap_iterations"])
    bootstrap = donor_cluster_bootstrap(predictions, profiles, bootstrap_iterations, int(config["random_seed"]))
    mcnemar = paired_mcnemar(predictions, profiles)
    prevalence = prevalence_aware_ppv(metrics, [float(x) for x in config["prevalence_scenarios"]])
    metrics.to_csv(TABLE_DIR / f"{args.stage}_metrics.csv", index=False)
    strata.to_csv(TABLE_DIR / f"{args.stage}_strata.csv", index=False)
    bootstrap.to_csv(TABLE_DIR / f"{args.stage}_donor_cluster_bootstrap.csv", index=False)
    mcnemar.to_csv(TABLE_DIR / f"{args.stage}_paired_mcnemar.csv", index=False)
    prevalence.to_csv(TABLE_DIR / f"{args.stage}_prevalence_aware_ppv.csv", index=False)
    pd.concat(dataset_summaries, ignore_index=True).to_csv(TABLE_DIR / f"{args.stage}_dataset_qc.csv", index=False)
    pd.concat(availability, ignore_index=True).to_csv(TABLE_DIR / f"{args.stage}_gene_availability.csv", index=False)

    ranking: pd.DataFrame | None = None
    if args.stage == "select":
        ranking = select_profile(metrics, strata, config)
        ranking.to_csv(TABLE_DIR / "selection_ranking.csv", index=False)
        decision = write_decision(args.decision, ranking, profiles, config, args.config, args.input_manifest, rows)
    else:
        selected = str(decision["selected_profile"])
        selected_strata = strata[
            (strata["profile_id"] == selected)
            & strata["model_schema_primary_eligible"].astype(bool)
        ]
        guard = config["selection_guardrails"]
        holdout_dataset_qc = pd.concat(dataset_summaries, ignore_index=True)
        schema_failure = not bool(holdout_dataset_qc["model_schema_primary_eligible"].astype(bool).all())
        veto = schema_failure
        for stratum, limit in (("abT_gold", "maximum_pooled_abt_fpr"), ("strict_NK", "maximum_pooled_strict_nk_fpr")):
            view = selected_strata[selected_strata["stratum"] == stratum]
            rate = float(view["predicted_gdT"].sum() / view["n_cells"].sum()) if view["n_cells"].sum() else np.nan
            if not np.isfinite(rate) or rate > float(guard[limit]):
                veto = True
        holdout_status = {
            "selected_profile": selected,
            "selection_decision_sha256": decision["decision_sha256"],
            "holdout_schema_failure": schema_failure,
            "holdout_veto": veto,
            "promotion_status": "vetoed_no_runner_up_selected" if veto else "holdout_passed",
        }
        (LOG_DIR / "holdout_status.json").write_text(json.dumps(holdout_status, indent=2) + "\n", encoding="utf-8")

    figure_path = FIGURE_DIR / f"{args.stage}_false_positive_rates.png"
    plot_false_positive_rates(strata, figure_path, f"{args.stage.title()} false-positive rates")
    report = render_report(args.stage, metrics, strata, ranking, decision, bootstrap, mcnemar, prevalence, figure_path)
    pdf_path = STATIC_DIR / ("selection_report.pdf" if args.stage == "select" else "gdtai_extension_evaluation_report.pdf")
    pdf_rendered = False if args.no_pdf else render_pdf(report, pdf_path)
    summary = {
        "stage": args.stage,
        "report": str(report),
        "profile_identity": profile_identity(profiles),
        "datasets": rows["dataset_id"].tolist(),
        "prediction_files": [str(path) for path in prediction_paths],
        "pdf": str(pdf_path) if pdf_rendered else None,
    }
    (LOG_DIR / f"{args.stage}_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    logging.info("Wrote %s report: %s", args.stage, report)


if __name__ == "__main__":
    main()
