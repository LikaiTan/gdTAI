#!/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python
"""Compare the four frozen gdTAI operating profiles on extension cohorts.

``select`` and ``holdout`` implement a guarded positive/negative evaluation.
``screen`` applies every frozen profile to alpha-beta-TCR extension cohorts and
reports negative-control false-positive rates without selecting a model. Input
H5AD files are opened read-only and must contain raw count-like CSR X.
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
from itertools import combinations
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
DEFAULT_SCREEN_MANIFEST = PROJECT_ROOT / "data/interim/extension_intake/tnk_filtered_h5ads_manifest.csv"
MODEL_REGISTRY = PROJECT_ROOT / "configs/models/gdtai/model_registry.csv"
OUTPUT_ROOT = PROJECT_ROOT / "Integrated_dataset"
PREFIX = "gdtai_extension_evaluation"
TABLE_DIR = OUTPUT_ROOT / "tables/gdT_prediction" / PREFIX
FIGURE_DIR = OUTPUT_ROOT / "figures/gdT_prediction" / PREFIX
LOG_DIR = OUTPUT_ROOT / "logs/gdT_prediction" / PREFIX
STATIC_DIR = PROJECT_ROOT / "gdT_prediction" / PREFIX
DECISION_PATH = LOG_DIR / "selection_decision.json"
TARGET_SUM = 10_000.0

AUXILIARY_AUDIT_GENES = (
    "CD4",
    "CD8A",
    "CD8B",
    "FOXP3",
    "TRDC",
    "TRDV1",
    "TRDV2",
    "TRDV3",
    "CD3D",
    "CD3E",
    "CD3G",
    "NKG7",
    "GNLY",
    "KLRD1",
    "TYROBP",
    "FCER1G",
    "NCAM1",
)

CANONICAL_OBS_COLUMNS = (
    "source_gse_id",
    "dataset_id",
    "extension_cohort_id",
    "sample_id",
    "library_id",
    "donor_id",
    "patient_id",
    "donor",
    "tissue_corrected",
    "tissue_harmonized",
    "tissue",
    "specimen_context",
    "specimen_context_source",
    "disease",
    "disease_status",
    "diagnosis",
    "cancer_type",
    "cancer_type_refined",
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
    parser.add_argument("--stage", choices=("select", "holdout", "screen"), required=True)
    parser.add_argument(
        "--input-manifest",
        type=Path,
        default=None,
        help="Defaults to the TNK-filter manifest for screen and the guarded manifest otherwise.",
    )
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--decision", type=Path, default=DECISION_PATH)
    parser.add_argument("--chunk-size", type=int, default=50_000)
    parser.add_argument("--bootstrap-iterations", type=int, default=None)
    parser.add_argument("--no-pdf", action="store_true")
    return parser.parse_args()


def resolve_input_manifest(stage: str, requested: Path | None) -> Path:
    if requested is not None:
        return requested
    return DEFAULT_SCREEN_MANIFEST if stage == "screen" else DEFAULT_INPUT_MANIFEST


def load_input_manifest(path: Path, stage: str) -> pd.DataFrame:
    """Load either the guarded evaluation manifest or the TNK-filter manifest."""
    manifest = pd.read_csv(path, dtype=str).fillna("")
    if stage == "screen" and {"cohort_id", "output_h5ad"}.issubset(manifest.columns):
        manifest = manifest.rename(
            columns={
                "cohort_id": "dataset_id",
                "output_h5ad": "h5ad_path",
                "output_sha256": "expected_h5ad_sha256",
            }
        )
        manifest["sealed_holdout"] = False
        manifest["matrix_path"] = "X"
        manifest["row_filter"] = "all"
    required = {"dataset_id", "h5ad_path", "sealed_holdout"}
    missing = required - set(manifest.columns)
    if missing:
        raise KeyError(f"Input manifest is missing columns: {sorted(missing)}")
    manifest["sealed_holdout"] = (
        manifest["sealed_holdout"].astype("string").fillna("").str.lower().isin({"true", "1", "yes"})
    )
    if manifest["dataset_id"].duplicated().any():
        duplicates = manifest.loc[manifest["dataset_id"].duplicated(False), "dataset_id"].tolist()
        raise ValueError(f"Input manifest has duplicate dataset_id values: {duplicates}")
    return manifest


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
    paired_ab_negative = paired_ab & ~has_gd & ~doublet
    single_ab_negative = has_ab & ~paired_ab & ~has_gd & ~doublet
    no_productive_tcr = ~has_ab & ~has_gd & ~doublet
    known_negative = ab_gold | strict_nk
    return pd.DataFrame(
        {
            "truth_class": label,
            "gdT_gold": gd_gold,
            "abT_gold": ab_gold,
            "paired_abT": paired_ab_negative,
            "single_abT": single_ab_negative,
            "strict_NK": strict_nk,
            "author_NK": author_nk,
            "gdT_silver": silver,
            "doublet": doublet,
            "has_TRA": has_tra,
            "has_TRB": has_trb,
            "has_TRG": has_trg,
            "has_TRD": has_trd,
            "has_any_abT": has_ab,
            "has_any_gdT": has_gd,
            "no_productive_tcr": no_productive_tcr,
            "known_negative_union": known_negative,
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
    model_genes, _ = feature_schema(profiles)
    union_genes = [*model_genes, *[gene for gene in AUXILIARY_AUDIT_GENES if gene not in model_genes]]
    summary_rows: list[dict[str, Any]] = []
    availability_rows: list[pd.DataFrame] = []
    wrote = False
    temporary_path = output_path.with_name(f".{output_path.name}.tmp")
    if temporary_path.exists():
        temporary_path.unlink()
    try:
        with h5py.File(h5ad_path, "r") as handle, gzip.open(temporary_path, "wt", encoding="utf-8", newline="") as sink:
            matrix_path = str(row.get("matrix_path", "X") or "X").strip()
            matrix_info = validate_count_matrix(handle, matrix_path)
            spec, availability = build_union_spec(handle, union_genes)
            availability.insert(0, "dataset_id", dataset_id)
            availability["used_by_model"] = availability["gene"].isin(model_genes)
            availability_rows.append(availability)
            model_availability = availability.loc[availability["used_by_model"]]
            auxiliary_availability = availability.loc[~availability["used_by_model"]]
            available_fraction = float(model_availability["available_in_h5ad"].mean())
            auxiliary_fraction = float(auxiliary_availability["available_in_h5ad"].mean())
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
                source_gse_id = first_nonempty(obs, ("source_gse_id", "dataset_id"))
                extension_cohort_id = first_nonempty(obs, ("extension_cohort_id",))
                extension_cohort_id[extension_cohort_id == ""] = dataset_id
                tissue = first_nonempty(obs, ("tissue_corrected", "tissue_harmonized", "tissue"))
                specimen_context = first_nonempty(obs, ("specimen_context",))

                def gene_values(gene: str) -> np.ndarray:
                    column = spec.gene_to_col.get(gene)
                    if column is None:
                        return np.zeros(rows.size, dtype=np.float32)
                    return x_union[:, column]

                trdc_positive = gene_values("TRDC") > 0
                trdv_positive = np.logical_or.reduce(
                    [gene_values("TRDV1") > 0, gene_values("TRDV2") > 0, gene_values("TRDV3") > 0]
                )
                tcr_gene_quadrant = np.select(
                    [
                        trdc_positive & trdv_positive,
                        trdc_positive & ~trdv_positive,
                        ~trdc_positive & trdv_positive,
                    ],
                    ["TRDC+TRDV+", "TRDC+TRDV-", "TRDC-TRDV+"],
                    default="TRDC-TRDV-",
                )
                data: dict[str, Any] = {
                    "dataset_id": dataset_id,
                    "extension_cohort_id": extension_cohort_id,
                    "source_gse_id": source_gse_id,
                    "cell_id": ids[rows],
                    "sample_id": sample_id,
                    "library_id": library_id,
                    "donor_id": donor_id,
                    "tissue": tissue,
                    "specimen_context": specimen_context,
                    "evaluation_annotation": annotation,
                    "evaluation_annotation_source": annotation_source,
                    "tcr_gene_quadrant": tcr_gene_quadrant,
                    "TRDC_expressed": trdc_positive,
                    "TRDV_expressed": trdv_positive,
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
                    "auxiliary_audit_gene_fraction": auxiliary_fraction,
                    "model_schema_primary_eligible": schema_primary_eligible,
                    "model_schema_warning": schema_warning,
                    "missing_critical_genes": ";".join(missing_critical),
                }
            )
        os.replace(temporary_path, output_path)
    except Exception:
        temporary_path.unlink(missing_ok=True)
        raise
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


def column_bool(frame: pd.DataFrame, column: str) -> np.ndarray:
    if column not in frame:
        return np.zeros(len(frame), dtype=bool)
    values = frame[column]
    if pd.api.types.is_bool_dtype(values.dtype):
        return values.to_numpy(dtype=bool)
    return as_bool(values)


def wilson_interval(successes: int, total: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if total <= 0:
        return np.nan, np.nan
    proportion = successes / total
    denominator = 1.0 + z * z / total
    center = (proportion + z * z / (2.0 * total)) / denominator
    half_width = z * np.sqrt(
        proportion * (1.0 - proportion) / total + z * z / (4.0 * total * total)
    ) / denominator
    return max(0.0, center - half_width), min(1.0, center + half_width)


def screen_masks(frame: pd.DataFrame) -> dict[str, tuple[np.ndarray, str, str]]:
    annotation = frame["evaluation_annotation"].astype(str).to_numpy(dtype=object)
    return {
        "all_cells": (np.ones(len(frame), dtype=bool), "candidate burden", "All TNK-filtered cells"),
        "abT_any_chain": (
            column_bool(frame, "abT_gold"),
            "negative control",
            "Productive TRA or TRB, no productive TRG/TRD, non-doublet",
        ),
        "paired_abT": (
            column_bool(frame, "paired_abT"),
            "negative control",
            "Productive TRA and TRB, no productive TRG/TRD, non-doublet",
        ),
        "single_abT": (
            column_bool(frame, "single_abT"),
            "negative control",
            "Productive TRA or TRB but not both, no productive TRG/TRD, non-doublet",
        ),
        "strict_NK": (
            column_bool(frame, "strict_NK"),
            "negative control",
            "NK annotation/expression, no productive alpha-beta or gamma-delta TCR, non-doublet",
        ),
        "author_NK": (
            column_bool(frame, "author_NK"),
            "annotation control",
            "Explicit author-level NK label; NKT labels excluded",
        ),
        "known_negative_union": (
            column_bool(frame, "known_negative_union"),
            "negative control",
            "Union of alpha-beta TCR and strict-NK negative controls",
        ),
        "no_productive_tcr": (
            column_bool(frame, "no_productive_tcr"),
            "unevaluable candidate",
            "No productive TRA/TRB/TRG/TRD evidence; not a positive or negative truth label",
        ),
        "annotation_CD4": (annotation == "CD4_T", "annotation control", "Normalized CD4 T annotation"),
        "annotation_Treg": (annotation == "TREG", "annotation control", "Normalized Treg annotation"),
        "annotation_CD8": (annotation == "CD8_T", "annotation control", "Normalized CD8 T annotation"),
        "annotation_NK": (annotation == "NK_CELL", "annotation control", "Normalized NK annotation"),
    }


def screen_strata(predictions: pd.DataFrame, profiles: list[Profile]) -> pd.DataFrame:
    groups: list[tuple[str, str, str, pd.DataFrame]] = [
        ("pooled", "__POOLED__", "__POOLED__", predictions)
    ]
    for dataset_id, group in predictions.groupby("dataset_id", sort=True):
        groups.append(("cohort", str(dataset_id), "__ALL_SOURCE_GSE__", group))
        for source_gse_id, source_group in group.groupby("source_gse_id", sort=True, dropna=False):
            source = str(source_gse_id).strip() or "unresolved"
            groups.append(("source_gse", str(dataset_id), source, source_group))
    rows: list[dict[str, Any]] = []
    for group_type, dataset_id, source_gse_id, group in groups:
        for stratum, (mask, interpretation, definition) in screen_masks(group).items():
            n_cells = int(mask.sum())
            if n_cells == 0:
                continue
            for profile in profiles:
                called = column_bool(group.loc[mask], f"{profile.profile_id}_pred")
                n_called = int(called.sum())
                lower, upper = wilson_interval(n_called, n_cells)
                rows.append(
                    {
                        "group_type": group_type,
                        "dataset_id": dataset_id,
                        "source_gse_id": source_gse_id,
                        "profile_id": profile.profile_id,
                        "stratum": stratum,
                        "interpretation": interpretation,
                        "definition": definition,
                        "n_cells": n_cells,
                        "predicted_gdT": n_called,
                        "call_rate": n_called / n_cells,
                        "wilson_95_low": lower,
                        "wilson_95_high": upper,
                    }
                )
    return pd.DataFrame(rows)


def screen_grouped_negative_controls(predictions: pd.DataFrame, profiles: list[Profile]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for group_level in ("library_id", "donor_id"):
        for keys, group in predictions.groupby(["dataset_id", "source_gse_id", group_level], sort=True, dropna=False):
            dataset_id, source_gse_id, group_value = (str(value) for value in keys)
            if not group_value.strip():
                group_value = "unresolved"
            masks = screen_masks(group)
            for stratum in ("abT_any_chain", "paired_abT", "single_abT", "strict_NK", "known_negative_union"):
                mask = masks[stratum][0]
                n_cells = int(mask.sum())
                if n_cells == 0:
                    continue
                for profile in profiles:
                    n_called = int(column_bool(group.loc[mask], f"{profile.profile_id}_pred").sum())
                    lower, upper = wilson_interval(n_called, n_cells)
                    rows.append(
                        {
                            "dataset_id": dataset_id,
                            "source_gse_id": source_gse_id,
                            "group_level": group_level,
                            "group_value": group_value,
                            "profile_id": profile.profile_id,
                            "stratum": stratum,
                            "n_cells": n_cells,
                            "predicted_gdT": n_called,
                            "call_rate": n_called / n_cells,
                            "wilson_95_low": lower,
                            "wilson_95_high": upper,
                        }
                    )
    return pd.DataFrame(rows)


def screen_pairwise(predictions: pd.DataFrame, profiles: list[Profile]) -> pd.DataFrame:
    groups: list[tuple[str, pd.DataFrame]] = [("__POOLED__", predictions)]
    groups.extend((str(key), group) for key, group in predictions.groupby("dataset_id", sort=True))
    rows: list[dict[str, Any]] = []
    for dataset_id, group in groups:
        populations = {
            "all_cells": np.ones(len(group), dtype=bool),
            "known_negative_union": column_bool(group, "known_negative_union"),
            "abT_any_chain": column_bool(group, "abT_gold"),
            "strict_NK": column_bool(group, "strict_NK"),
        }
        for population, mask in populations.items():
            n_cells = int(mask.sum())
            if n_cells == 0:
                continue
            for left, right in combinations(profiles, 2):
                left_pred = column_bool(group.loc[mask], f"{left.profile_id}_pred")
                right_pred = column_bool(group.loc[mask], f"{right.profile_id}_pred")
                both = int((left_pred & right_pred).sum())
                left_only = int((left_pred & ~right_pred).sum())
                right_only = int((~left_pred & right_pred).sum())
                union = both + left_only + right_only
                discordant = left_only + right_only
                pvalue = (
                    float(binomtest(min(left_only, right_only), discordant, 0.5).pvalue)
                    if population != "all_cells" and discordant
                    else (1.0 if population != "all_cells" else np.nan)
                )
                rows.append(
                    {
                        "dataset_id": dataset_id,
                        "population": population,
                        "profile_left": left.profile_id,
                        "profile_right": right.profile_id,
                        "n_cells": n_cells,
                        "left_calls": int(left_pred.sum()),
                        "right_calls": int(right_pred.sum()),
                        "both": both,
                        "left_only": left_only,
                        "right_only": right_only,
                        "union": union,
                        "jaccard": both / union if union else 1.0,
                        "discordant": discordant,
                        "exact_mcnemar_p": pvalue,
                        "note": "Cell-level exact paired test is descriptive because cells within donors are correlated.",
                    }
                )
    return pd.DataFrame(rows)


def screen_schema_sensitivity(
    predictions: pd.DataFrame,
    profiles: list[Profile],
    schema_warning_cohorts: set[str],
) -> pd.DataFrame:
    dataset = predictions["dataset_id"].astype(str)
    scopes: list[tuple[str, pd.DataFrame]] = [("all_cohorts", predictions)]
    if schema_warning_cohorts:
        scopes.extend(
            [
                ("complete_schema_sensitivity", predictions.loc[~dataset.isin(schema_warning_cohorts)]),
                ("schema_warning_cohorts", predictions.loc[dataset.isin(schema_warning_cohorts)]),
            ]
        )
    rows: list[dict[str, Any]] = []
    keep_strata = {
        "all_cells",
        "abT_any_chain",
        "paired_abT",
        "single_abT",
        "strict_NK",
        "known_negative_union",
        "no_productive_tcr",
    }
    for scope, frame in scopes:
        if frame.empty:
            continue
        pooled = screen_strata(frame, profiles)
        pooled = pooled[(pooled["group_type"] == "pooled") & pooled["stratum"].isin(keep_strata)].copy()
        pooled.insert(0, "scope", scope)
        pooled.insert(1, "n_cohorts", int(frame["dataset_id"].nunique()))
        rows.extend(pooled.to_dict(orient="records"))
    return pd.DataFrame(rows)


def screen_tcr_quadrants(predictions: pd.DataFrame, profiles: list[Profile]) -> pd.DataFrame:
    groups: list[tuple[str, pd.DataFrame]] = [("__POOLED__", predictions)]
    groups.extend((str(key), group) for key, group in predictions.groupby("dataset_id", sort=True))
    quadrants = ("TRDC+TRDV+", "TRDC+TRDV-", "TRDC-TRDV+", "TRDC-TRDV-")
    rows: list[dict[str, Any]] = []
    for dataset_id, group in groups:
        populations = {
            "all_predictions": np.ones(len(group), dtype=bool),
            "abT_false_positives": column_bool(group, "abT_gold"),
            "strict_NK_false_positives": column_bool(group, "strict_NK"),
            "no_TCR_candidates": column_bool(group, "no_productive_tcr"),
        }
        quadrant_values = group["tcr_gene_quadrant"].astype(str).to_numpy(dtype=object)
        for profile in profiles:
            predicted = column_bool(group, f"{profile.profile_id}_pred")
            for population, population_mask in populations.items():
                selected = predicted & population_mask
                total = int(selected.sum())
                for quadrant in quadrants:
                    count = int((selected & (quadrant_values == quadrant)).sum())
                    rows.append(
                        {
                            "dataset_id": dataset_id,
                            "profile_id": profile.profile_id,
                            "population": population,
                            "tcr_gene_quadrant": quadrant,
                            "n_selected": total,
                            "n_cells": count,
                            "fraction_of_selected": count / total if total else np.nan,
                        }
                    )
    return pd.DataFrame(rows)


def screen_negative_control_outliers(grouped: pd.DataFrame, minimum_n: int = 100) -> pd.DataFrame:
    eligible = grouped[
        (grouped["n_cells"] >= minimum_n)
        & grouped["stratum"].isin(["abT_any_chain", "strict_NK", "known_negative_union"])
    ].copy()
    if eligible.empty:
        return eligible
    eligible = eligible.sort_values(
        ["profile_id", "stratum", "group_level", "call_rate", "n_cells"],
        ascending=[True, True, True, False, False],
    )
    return eligible.groupby(["profile_id", "stratum", "group_level"], sort=True).head(5).reset_index(drop=True)


def historical_benchmarks() -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    v2_path = OUTPUT_ROOT / "models/gdT_prediction_classifier/gdTAI_v2.0/mode_metrics.csv"
    for _, item in pd.read_csv(v2_path).iterrows():
        profile_id = "v2_high_f1" if item["mode"] == "high_f1" else "v2_high_purity"
        rows.append(
            {
                "profile_id": profile_id,
                "evaluation_frame": "packaged V2 validation",
                "precision": item["validation_precision"],
                "recall": item["validation_recall"],
                "specificity": item["validation_specificity"],
                "f1": item["validation_f1"],
                "source_file": str(v2_path.relative_to(PROJECT_ROOT)),
            }
        )
    v3_path = OUTPUT_ROOT / "tables/gdT_prediction/gdtai_v3_round12_vs_round14/full_atlas_performance.csv"
    mapping = {"round12": "v3_round12_high_purity", "round14": "v3_round14_balanced"}
    for _, item in pd.read_csv(v3_path).iterrows():
        rows.append(
            {
                "profile_id": mapping[str(item["model"])],
                "evaluation_frame": "historical full-atlas primary gold",
                "precision": item["precision"],
                "recall": item["recall"],
                "specificity": item["specificity"],
                "f1": item["f1"],
                "source_file": str(v3_path.relative_to(PROJECT_ROOT)),
            }
        )
    return pd.DataFrame(rows)


def profile_manifest(profiles: list[Profile]) -> pd.DataFrame:
    rows = []
    for profile in profiles:
        payload = profile.payload["base_model"] if profile.model_id == "gdtai_v2" else profile.payload
        threshold: Any
        if profile.model_id == "gdtai_v2" and profile.mode == "high_purity":
            threshold = "annotation-specific"
        elif profile.model_id == "gdtai_v2":
            threshold = profile.payload["operating_modes"][profile.mode]["threshold"]
        else:
            threshold = payload["threshold"]
        rows.append(
            {
                "profile_id": profile.profile_id,
                "model_id": profile.model_id,
                "mode": profile.mode,
                "threshold": threshold,
                "gene_features": len(payload["gene_names"]),
                "artifact_sha256": profile.sha256,
                "artifact_path": str(profile.artifact.relative_to(PROJECT_ROOT)),
            }
        )
    return pd.DataFrame(rows)


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


def screen_cohort_table(strata: pd.DataFrame) -> pd.DataFrame:
    view = strata[strata["group_type"] == "cohort"]
    rows: list[dict[str, Any]] = []
    for (dataset_id, profile_id), group in view.groupby(["dataset_id", "profile_id"], sort=True):
        lookup = group.set_index("stratum")

        def value(stratum: str, column: str, default: Any = np.nan) -> Any:
            return lookup.at[stratum, column] if stratum in lookup.index else default

        rows.append(
            {
                "dataset_id": dataset_id,
                "profile_id": profile_id,
                "total_cells": value("all_cells", "n_cells", 0),
                "predicted_gdT": value("all_cells", "predicted_gdT", 0),
                "predicted_rate": value("all_cells", "call_rate"),
                "abT_n": value("abT_any_chain", "n_cells", 0),
                "abT_fp": value("abT_any_chain", "predicted_gdT", 0),
                "abT_fpr": value("abT_any_chain", "call_rate"),
                "paired_abT_n": value("paired_abT", "n_cells", 0),
                "paired_abT_fp": value("paired_abT", "predicted_gdT", 0),
                "paired_abT_fpr": value("paired_abT", "call_rate"),
                "single_abT_n": value("single_abT", "n_cells", 0),
                "single_abT_fp": value("single_abT", "predicted_gdT", 0),
                "single_abT_fpr": value("single_abT", "call_rate"),
                "strict_NK_n": value("strict_NK", "n_cells", 0),
                "strict_NK_fp": value("strict_NK", "predicted_gdT", 0),
                "strict_NK_fpr": value("strict_NK", "call_rate"),
                "no_TCR_n": value("no_productive_tcr", "n_cells", 0),
                "no_TCR_calls": value("no_productive_tcr", "predicted_gdT", 0),
            }
        )
    return pd.DataFrame(rows)


def profile_display_name(profile_id: str) -> str:
    return {
        "v2_high_f1": "V2 high F1",
        "v2_high_purity": "V2 high purity",
        "v3_round14_balanced": "V3 R14 balanced",
        "v3_round12_high_purity": "V3 R12 high purity",
    }.get(profile_id, profile_id)


def plot_screen_overview(
    strata: pd.DataFrame,
    predictions: pd.DataFrame,
    profiles: list[Profile],
    path: Path,
) -> None:
    profile_ids = [profile.profile_id for profile in profiles]
    profile_labels = [profile_display_name(profile_id) for profile_id in profile_ids]
    cohort_ids = sorted(predictions["dataset_id"].astype(str).unique())
    cohort = strata[strata["group_type"] == "cohort"]

    def matrix_for(stratum: str) -> np.ndarray:
        view = cohort[cohort["stratum"] == stratum]
        pivot = view.pivot(index="dataset_id", columns="profile_id", values="call_rate")
        return pivot.reindex(index=cohort_ids, columns=profile_ids).fillna(0).to_numpy(dtype=float)

    all_rates = matrix_for("all_cells")
    ab_rates = matrix_for("abT_any_chain")
    nk_rates = matrix_for("strict_NK")
    colors = ["#266B8E", "#4A8F5C", "#C47A2C", "#A44942"]
    fig, axes = plt.subplots(2, 2, figsize=(14, 10.5), constrained_layout=True)
    x = np.arange(len(cohort_ids))
    width = 0.2
    for index, profile_id in enumerate(profile_ids):
        axes[0, 0].bar(
            x + (index - 1.5) * width,
            all_rates[:, index] * 100,
            width,
            color=colors[index],
            label=profile_display_name(profile_id),
        )
    axes[0, 0].set_xticks(x, cohort_ids, rotation=28, ha="right")
    axes[0, 0].set_ylabel("Predicted gdT cells (%)")
    axes[0, 0].set_title("A. Total call burden")
    axes[0, 0].grid(axis="y", alpha=0.2)
    axes[0, 0].legend(frameon=False, fontsize=8, ncol=2)

    for ax, values, title in (
        (axes[0, 1], ab_rates * 100, "B. Alpha-beta TCR negative-control FPR (%)"),
        (axes[1, 0], nk_rates * 100, "C. Strict-NK negative-control FPR (%)"),
    ):
        image_obj = ax.imshow(values, aspect="auto", cmap="YlOrRd", vmin=0)
        ax.set_xticks(np.arange(len(profile_ids)), profile_labels, rotation=20, ha="right")
        ax.set_yticks(np.arange(len(cohort_ids)), cohort_ids)
        ax.set_title(title)
        for row_index in range(values.shape[0]):
            for column_index in range(values.shape[1]):
                value = values[row_index, column_index]
                ax.text(column_index, row_index, f"{value:.2f}", ha="center", va="center", fontsize=7)
        fig.colorbar(image_obj, ax=ax, fraction=0.03, pad=0.02)

    quadrant_order = ["TRDC+TRDV+", "TRDC+TRDV-", "TRDC-TRDV+", "TRDC-TRDV-"]
    quadrant_colors = ["#2B7A78", "#D49A32", "#6A7FDB", "#9AA2A6"]
    bottom = np.zeros(len(profile_ids), dtype=float)
    for quadrant, color in zip(quadrant_order, quadrant_colors, strict=True):
        values = []
        for profile_id in profile_ids:
            called = column_bool(predictions, f"{profile_id}_pred")
            total = int(called.sum())
            count = int((called & predictions["tcr_gene_quadrant"].eq(quadrant).to_numpy()).sum())
            values.append(100 * count / total if total else 0.0)
        axes[1, 1].bar(profile_labels, values, bottom=bottom, color=color, label=quadrant)
        bottom += np.asarray(values)
    axes[1, 1].set_ylim(0, 100)
    axes[1, 1].set_ylabel("Share of predicted cells (%)")
    axes[1, 1].set_title("D. TCR-gene expression composition of calls")
    axes[1, 1].tick_params(axis="x", rotation=25)
    axes[1, 1].legend(frameon=False, fontsize=8, ncol=2)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=240)
    plt.close(fig)


def plot_screen_overlap(pairwise: pd.DataFrame, profiles: list[Profile], path: Path) -> None:
    profile_ids = [profile.profile_id for profile in profiles]
    profile_labels = [profile_display_name(profile_id) for profile_id in profile_ids]
    matrix = np.eye(len(profile_ids), dtype=float)
    pooled = pairwise[
        (pairwise["dataset_id"] == "__POOLED__") & (pairwise["population"] == "all_cells")
    ]
    for _, row in pooled.iterrows():
        left = profile_ids.index(str(row["profile_left"]))
        right = profile_ids.index(str(row["profile_right"]))
        matrix[left, right] = matrix[right, left] = float(row["jaccard"])
    fig, ax = plt.subplots(figsize=(8.2, 6.6), constrained_layout=True)
    image_obj = ax.imshow(matrix, vmin=0, vmax=1, cmap="viridis")
    ax.set_xticks(np.arange(len(profile_ids)), profile_labels, rotation=20, ha="right")
    ax.set_yticks(np.arange(len(profile_ids)), profile_labels)
    ax.set_title("Pooled prediction overlap (Jaccard index)")
    for row_index in range(matrix.shape[0]):
        for column_index in range(matrix.shape[1]):
            color = "white" if matrix[row_index, column_index] < 0.55 else "black"
            ax.text(column_index, row_index, f"{matrix[row_index, column_index]:.3f}", ha="center", va="center", color=color)
    fig.colorbar(image_obj, ax=ax, fraction=0.045, pad=0.03)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=240)
    plt.close(fig)


def report_table(frame: pd.DataFrame, percent_columns: tuple[str, ...] = ()) -> str:
    view = frame.copy()
    for column in percent_columns:
        if column in view:
            view[column] = view[column].map(lambda value: "" if pd.isna(value) else f"{100 * float(value):.3f}%")
    for column in view.select_dtypes(include=["float"]).columns:
        if column not in percent_columns:
            view[column] = view[column].map(lambda value: "" if pd.isna(value) else f"{float(value):.4g}")
    return dataframe_html(view)


def render_screen_report(
    strata: pd.DataFrame,
    cohort_table: pd.DataFrame,
    pairwise: pd.DataFrame,
    sensitivity: pd.DataFrame,
    quadrants: pd.DataFrame,
    outliers: pd.DataFrame,
    historical: pd.DataFrame,
    dataset_qc: pd.DataFrame,
    artifacts: pd.DataFrame,
    predictions: pd.DataFrame,
    overview_figure: Path,
    overlap_figure: Path,
) -> tuple[Path, dict[str, Any]]:
    STATIC_DIR.mkdir(parents=True, exist_ok=True)
    pooled = strata[(strata["group_type"] == "pooled")]

    def pooled_row(profile_id: str, stratum: str) -> pd.Series:
        view = pooled[(pooled["profile_id"] == profile_id) & (pooled["stratum"] == stratum)]
        if view.empty:
            return pd.Series(dtype=object)
        return view.iloc[0]

    profile_ids = artifacts["profile_id"].tolist()
    pooled_rows = []
    for profile_id in profile_ids:
        total = pooled_row(profile_id, "all_cells")
        ab = pooled_row(profile_id, "abT_any_chain")
        paired = pooled_row(profile_id, "paired_abT")
        single = pooled_row(profile_id, "single_abT")
        nk = pooled_row(profile_id, "strict_NK")
        known = pooled_row(profile_id, "known_negative_union")
        no_tcr = pooled_row(profile_id, "no_productive_tcr")
        pooled_rows.append(
            {
                "profile_id": profile_id,
                "all_calls": int(total.get("predicted_gdT", 0)),
                "all_call_rate": total.get("call_rate", np.nan),
                "abT_FP": int(ab.get("predicted_gdT", 0)),
                "abT_FPR": ab.get("call_rate", np.nan),
                "paired_abT_FP": int(paired.get("predicted_gdT", 0)),
                "paired_abT_FPR": paired.get("call_rate", np.nan),
                "single_abT_FP": int(single.get("predicted_gdT", 0)),
                "single_abT_FPR": single.get("call_rate", np.nan),
                "strict_NK_FP": int(nk.get("predicted_gdT", 0)),
                "strict_NK_FPR": nk.get("call_rate", np.nan),
                "known_negative_FP": int(known.get("predicted_gdT", 0)),
                "known_negative_FPR": known.get("call_rate", np.nan),
                "no_TCR_calls": int(no_tcr.get("predicted_gdT", 0)),
            }
        )
    pooled_summary = pd.DataFrame(pooled_rows)
    best_known = pooled_summary.sort_values(["known_negative_FPR", "all_calls"]).iloc[0]
    best_ab = pooled_summary.sort_values(["abT_FPR", "all_calls"]).iloc[0]
    best_nk = pooled_summary.sort_values(["strict_NK_FPR", "all_calls"]).iloc[0]
    compatible_known = sensitivity[
        (sensitivity["scope"] == "complete_schema_sensitivity")
        & (sensitivity["stratum"] == "known_negative_union")
    ].sort_values(["call_rate", "predicted_gdT"])
    compatible_best = compatible_known.iloc[0] if not compatible_known.empty else pd.Series(dtype=object)
    observed_gold = int(column_bool(predictions, "gdT_gold").sum())
    observed_silver = int(column_bool(predictions, "gdT_silver").sum())
    warned = dataset_qc[dataset_qc["model_schema_warning"].astype(bool)]
    summary = {
        "n_cells": int(len(predictions)),
        "n_cohorts": int(predictions["dataset_id"].nunique()),
        "n_source_gse": int(predictions["source_gse_id"].replace("", np.nan).nunique()),
        "observed_gdT_gold": observed_gold,
        "observed_gdT_silver": observed_silver,
        "lowest_known_negative_fpr_profile": str(best_known["profile_id"]),
        "lowest_known_negative_fpr": float(best_known["known_negative_FPR"]),
        "lowest_abt_fpr_profile": str(best_ab["profile_id"]),
        "lowest_abt_fpr": float(best_ab["abT_FPR"]),
        "lowest_strict_nk_fpr_profile": str(best_nk["profile_id"]),
        "lowest_strict_nk_fpr": float(best_nk["strict_NK_FPR"]),
        "complete_schema_lowest_known_negative_fpr_profile": (
            str(compatible_best.get("profile_id", "not_applicable"))
        ),
        "complete_schema_lowest_known_negative_fpr": (
            float(compatible_best.get("call_rate", np.nan))
        ),
        "schema_warning_cohorts": warned["dataset_id"].astype(str).tolist(),
        "decision": "No model selected or promoted: the new cohorts do not contain an unbiased gdT positive truth set.",
    }

    def image_tag(path: Path, caption: str) -> str:
        relative = Path(os.path.relpath(path, STATIC_DIR)).as_posix()
        return f"<figure><img src='{html.escape(relative)}'><figcaption>{html.escape(caption)}</figcaption></figure>"

    compact_cohort = cohort_table[
        [
            "dataset_id", "profile_id", "total_cells", "predicted_gdT", "predicted_rate",
            "abT_n", "abT_fp", "abT_fpr", "strict_NK_n", "strict_NK_fp", "strict_NK_fpr",
        ]
    ]
    feature_qc = dataset_qc[
        [
            "dataset_id", "n_cells_evaluated", "matrix_encoding", "model_gene_fraction",
            "auxiliary_audit_gene_fraction", "model_schema_warning", "missing_critical_genes",
        ]
    ]
    sensitivity_report = sensitivity[
        sensitivity["stratum"].isin(["abT_any_chain", "strict_NK", "known_negative_union"])
    ][["scope", "n_cohorts", "profile_id", "stratum", "n_cells", "predicted_gdT", "call_rate", "wilson_95_low", "wilson_95_high"]]
    quadrant_report = quadrants[
        (quadrants["dataset_id"] == "__POOLED__")
        & quadrants["population"].isin(["abT_false_positives", "strict_NK_false_positives"])
    ][["profile_id", "population", "tcr_gene_quadrant", "n_selected", "n_cells", "fraction_of_selected"]]
    outlier_report = outliers[
        (outliers["group_level"] == "library_id")
        & outliers["stratum"].isin(["abT_any_chain", "strict_NK"])
    ].groupby(["profile_id", "stratum"], sort=True).head(2)
    outlier_report = outlier_report[
        ["dataset_id", "source_gse_id", "group_value", "profile_id", "stratum", "n_cells", "predicted_gdT", "call_rate", "wilson_95_low", "wilson_95_high"]
    ]
    title = "Frozen gdTAI V2/V3 extension-cohort negative-control evaluation"
    body = f"""<!doctype html><html><head><meta charset='utf-8'><title>{html.escape(title)}</title>
<style>
@page{{size:A4 landscape;margin:10mm}}*{{box-sizing:border-box}}body{{font-family:Arial,sans-serif;max-width:1380px;margin:0 auto;padding:24px;color:#172126;line-height:1.42;background:#fff}}
h1{{font-size:27px;margin:0 0 6px}}h2{{font-size:19px;margin:28px 0 8px;border-bottom:2px solid #314a56;padding-bottom:4px}}h3{{font-size:15px;margin:18px 0 6px}}p,li{{font-size:12.5px}}.lede{{font-size:14px;max-width:1050px}}.status{{border-left:5px solid #a44942;background:#f5f1ed;padding:10px 14px;margin:16px 0;page-break-inside:avoid}}
.metrics{{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:10px;margin:16px 0}}.metric{{border:1px solid #cdd5d8;padding:10px;background:#f7f9f9}}.metric b{{display:block;font-size:20px;color:#174f68}}.metric span{{font-size:11px;color:#4d5d64}}
table{{border-collapse:collapse;width:100%;font-size:7.4px;table-layout:auto;page-break-inside:auto}}tr{{page-break-inside:avoid}}th,td{{border:1px solid #d3dadd;padding:3px 4px;text-align:right;white-space:normal;overflow-wrap:anywhere}}th:first-child,td:first-child{{text-align:left}}th{{background:#eaf0f2;color:#172126}}code{{background:#eef2f3;padding:2px 4px}}figure{{margin:12px 0;page-break-inside:avoid}}img{{display:block;max-width:100%;max-height:172mm;margin:auto}}figcaption{{font-size:10px;color:#526169;margin-top:4px}}.page-break{{break-before:page;page-break-before:always}}.note{{font-size:11px;color:#48585f}}
@media(max-width:800px){{.metrics{{grid-template-columns:repeat(2,minmax(0,1fr))}}body{{padding:14px;overflow-x:auto}}}}
</style></head><body>
<h1>{html.escape(title)}</h1>
<p class='lede'>All four registered operating profiles were applied directly to the eight independently TNK-filtered cohorts using <code>log1p(raw counts per 10,000)</code>. Model files were checksum-pinned and every H5AD was opened read-only.</p>
<div class='status'><strong>Interpretation boundary:</strong> these studies used alpha-beta TCR assays and contributed {observed_gold} gdT-gold and {observed_silver} gdT-silver cells under the project TCR rules. Therefore this report estimates false-positive behavior and prediction burden only. It does not estimate recall, precision, F1, ROC-AUC, or PR-AUC on these cohorts, and it does not promote a model.</div>
<div class='metrics'><div class='metric'><b>{len(predictions):,}</b><span>TNK-filtered cells screened</span></div><div class='metric'><b>{predictions['dataset_id'].nunique()}</b><span>extension cohorts</span></div><div class='metric'><b>{predictions['source_gse_id'].replace('', np.nan).nunique()}</b><span>represented RNA GSE accessions</span></div><div class='metric'><b>{html.escape(str(best_known['profile_id']))}</b><span>lowest pooled known-negative FPR ({100 * float(best_known['known_negative_FPR']):.3f}%)</span></div></div>
<h2>Main Results</h2>
{image_tag(overview_figure, 'Prediction burden, alpha-beta T-cell and strict-NK negative controls, and TCR-gene expression composition.')}
<h3>Pooled negative-control results</h3>
{report_table(pooled_summary, ('all_call_rate','abT_FPR','paired_abT_FPR','single_abT_FPR','strict_NK_FPR','known_negative_FPR'))}
<div class='status'><strong>Trade-off:</strong> <code>{html.escape(str(best_ab['profile_id']))}</code> has the lowest pooled alpha-beta TCR FPR ({100 * float(best_ab['abT_FPR']):.3f}%), whereas <code>{html.escape(str(best_nk['profile_id']))}</code> has the lowest strict-NK FPR ({100 * float(best_nk['strict_NK_FPR']):.3f}%). After excluding the reduced-feature GSE169246 sensitivity case, the lowest known-negative union FPR is <code>{html.escape(str(compatible_best.get('profile_id', 'not applicable')))}</code>{'' if compatible_best.empty else f" ({100 * float(compatible_best['call_rate']):.3f}%)"}. Historical positive sensitivity still favors V3 Round 14 over Round 12. These data therefore do not support replacing the current operating profiles with a single new winner.</div>
<p class='note'>A productive single TRA or single TRB chain is accepted as alpha-beta negative evidence, as requested. Strict NK controls require no productive alpha-beta or gamma-delta chain evidence. Calls among cells with no productive TCR are candidates, not validated true positives.</p>
<h3>Feature-coverage sensitivity analysis</h3>
{report_table(sensitivity_report, ('call_rate','wilson_95_low','wilson_95_high'))}
<h2>Agreement Between Profiles</h2>
{image_tag(overlap_figure, 'Pooled cell-level Jaccard overlap. Detailed cohort and negative-control discordances are in screen_pairwise_overlap.csv.')}
<h2 class='page-break'>Cohort-by-Cohort Audit</h2>
{report_table(compact_cohort, ('predicted_rate','abT_fpr','strict_NK_fpr'))}
<p class='note'>Wilson 95% confidence intervals, source-GSE rows, single-chain versus paired-chain strata, and donor/library results are provided in the canonical CSV tables.</p>
<h2>Failure-Mode Inspection</h2>
<h3>TRDC/TRDV expression among known false positives</h3>
{report_table(quadrant_report, ('fraction_of_selected',))}
<h3>Highest-rate libraries with at least 100 controls</h3>
{report_table(outlier_report, ('call_rate','wilson_95_low','wilson_95_high'))}
<h2>Historical Positive/Negative Benchmarks</h2>
<p>These values come from the frozen models' prior validation frames and are shown only to retain the sensitivity context missing from the new studies. They are not pooled with the extension results because the cohorts and truth definitions differ.</p>
{report_table(historical, ('precision','recall','specificity','f1'))}
<h2>Model and Input Integrity</h2>
<h3>Frozen artifacts</h3>{report_table(artifacts)}
<h3>Expression schema</h3>{report_table(feature_qc, ('model_gene_fraction','auxiliary_audit_gene_fraction'))}
<p class='note'>A schema warning does not mean critical genes are absent. It means fewer than 80% of the full model feature genes were available. Such cohorts must be interpreted separately because zero-filled missing genes can shift predictions.</p>
<h2>Methods and Definitions</h2>
<p><strong>Input:</strong> cohort-separated Phase 1 TNK H5ADs. <strong>Normalization:</strong> raw integer-like <code>X</code> counts were normalized cell-wise to 10,000 and transformed with <code>log1p</code> in chunks. <strong>Models:</strong> V2 high-F1, V2 annotation-specific high-purity, V3 Round 14 balanced, and V3 Round 12 high-purity were loaded only after SHA256 verification.</p>
<p><strong>Alpha-beta negative:</strong> productive TRA or TRB, no productive TRG/TRD, and not marked doublet. Paired and single-chain cells are reported separately. <strong>Strict NK negative:</strong> explicit/derived NK classification with no productive TRA/TRB/TRG/TRD and not marked doublet. <strong>Uncertainty:</strong> Wilson 95% binomial intervals. <strong>Paired comparisons:</strong> exact cell-level McNemar tests are descriptive because cells within donors are correlated.</p>
<p><strong>Decision:</strong> {html.escape(summary['decision'])}</p>
</body></html>"""
    path = STATIC_DIR / "index.html"
    path.write_text(body, encoding="utf-8")
    return path, summary


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
    manifest_path = resolve_input_manifest(args.stage, args.input_manifest)
    manifest = load_input_manifest(manifest_path, args.stage)
    if args.stage == "select":
        rows = manifest.loc[~manifest["sealed_holdout"]].copy()
        decision = None
    elif args.stage == "holdout":
        decision = verify_decision(args.decision, profiles, args.config, manifest_path)
        rows = manifest.loc[manifest["sealed_holdout"]].copy()
    else:
        decision = None
        rows = manifest.copy()
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
        expected_sha256 = str(row.get("expected_h5ad_sha256", "")).strip()
        if expected_sha256:
            h5ad_path = Path(str(row["h5ad_path"])).expanduser()
            if not h5ad_path.is_absolute():
                h5ad_path = PROJECT_ROOT / h5ad_path
            observed_sha256 = sha256_file(h5ad_path)
            if observed_sha256 != expected_sha256:
                raise RuntimeError(
                    f"Read-only source checksum mismatch for {dataset_id}: "
                    f"expected {expected_sha256}, observed {observed_sha256}"
                )
            summary["expected_h5ad_sha256"] = expected_sha256
            summary["observed_h5ad_sha256_after_inference"] = observed_sha256
            summary["source_checksum_unchanged"] = True
        prediction_paths.append(output_path)
        dataset_summaries.append(summary)
        availability.append(genes)

    predictions = load_prediction_files(prediction_paths)
    dataset_qc = pd.concat(dataset_summaries, ignore_index=True)
    gene_availability = pd.concat(availability, ignore_index=True)
    dataset_qc.to_csv(TABLE_DIR / f"{args.stage}_dataset_qc.csv", index=False)
    gene_availability.to_csv(TABLE_DIR / f"{args.stage}_gene_availability.csv", index=False)
    if args.stage == "screen":
        strata = screen_strata(predictions, profiles)
        cohort_table = screen_cohort_table(strata)
        grouped = screen_grouped_negative_controls(predictions, profiles)
        pairwise = screen_pairwise(predictions, profiles)
        warning_cohorts = set(
            dataset_qc.loc[dataset_qc["model_schema_warning"].astype(bool), "dataset_id"].astype(str)
        )
        sensitivity = screen_schema_sensitivity(predictions, profiles, warning_cohorts)
        quadrants = screen_tcr_quadrants(predictions, profiles)
        outliers = screen_negative_control_outliers(grouped)
        historical = historical_benchmarks()
        artifacts = profile_manifest(profiles)
        strata.to_csv(TABLE_DIR / "screen_strata.csv", index=False)
        cohort_table.to_csv(TABLE_DIR / "screen_cohort_summary.csv", index=False)
        grouped.to_csv(TABLE_DIR / "screen_library_donor_negative_controls.csv", index=False)
        pairwise.to_csv(TABLE_DIR / "screen_pairwise_overlap.csv", index=False)
        sensitivity.to_csv(TABLE_DIR / "screen_schema_sensitivity.csv", index=False)
        quadrants.to_csv(TABLE_DIR / "screen_tcr_gene_quadrants.csv", index=False)
        outliers.to_csv(TABLE_DIR / "screen_negative_control_outliers.csv", index=False)
        historical.to_csv(TABLE_DIR / "screen_historical_benchmarks.csv", index=False)
        artifacts.to_csv(TABLE_DIR / "screen_model_artifacts.csv", index=False)
        rows.to_csv(TABLE_DIR / "screen_input_manifest.csv", index=False)

        overview_figure = FIGURE_DIR / "screen_overview.png"
        overlap_figure = FIGURE_DIR / "screen_prediction_overlap.png"
        plot_screen_overview(strata, predictions, profiles, overview_figure)
        plot_screen_overlap(pairwise, profiles, overlap_figure)
        served_figure_dir = STATIC_DIR / "figures"
        served_figure_dir.mkdir(parents=True, exist_ok=True)
        served_overview_figure = served_figure_dir / overview_figure.name
        served_overlap_figure = served_figure_dir / overlap_figure.name
        shutil.copy2(overview_figure, served_overview_figure)
        shutil.copy2(overlap_figure, served_overlap_figure)
        report, scientific_summary = render_screen_report(
            strata,
            cohort_table,
            pairwise,
            sensitivity,
            quadrants,
            outliers,
            historical,
            dataset_qc,
            artifacts,
            predictions,
            served_overview_figure,
            served_overlap_figure,
        )
        pdf_path = STATIC_DIR / "gdtai_extension_evaluation_report.pdf"
        pdf_rendered = False if args.no_pdf else render_pdf(report, pdf_path)
        summary = {
            "stage": args.stage,
            "report": str(report),
            "pdf": str(pdf_path) if pdf_rendered else None,
            "input_manifest": str(manifest_path),
            "input_manifest_sha256": sha256_file(manifest_path),
            "profile_identity": profile_identity(profiles),
            "datasets": rows["dataset_id"].tolist(),
            "prediction_files": [str(path) for path in prediction_paths],
            "scientific_summary": scientific_summary,
            "source_checksums_unchanged": bool(dataset_qc.get("source_checksum_unchanged", pd.Series(False)).all()),
        }
        (LOG_DIR / "screen_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
        markdown = f"""# Frozen gdTAI extension-cohort screen

- Cells screened: `{scientific_summary['n_cells']:,}` across `{scientific_summary['n_cohorts']}` cohorts.
- Positive truth observed: `{scientific_summary['observed_gdT_gold']}` gold and `{scientific_summary['observed_gdT_silver']}` silver cells.
- Lowest pooled known-negative FPR: `{scientific_summary['lowest_known_negative_fpr_profile']}` at `{scientific_summary['lowest_known_negative_fpr']:.6%}`.
- Lowest pooled alpha-beta TCR FPR: `{scientific_summary['lowest_abt_fpr_profile']}` at `{scientific_summary['lowest_abt_fpr']:.6%}`.
- Lowest pooled strict-NK FPR: `{scientific_summary['lowest_strict_nk_fpr_profile']}` at `{scientific_summary['lowest_strict_nk_fpr']:.6%}`.
- Complete-schema lowest known-negative FPR: `{scientific_summary['complete_schema_lowest_known_negative_fpr_profile']}` at `{scientific_summary['complete_schema_lowest_known_negative_fpr']:.6%}`.
- Decision: {scientific_summary['decision']}
- Source H5AD checksums unchanged: `{summary['source_checksums_unchanged']}`.
- HTML: `{report}`
- PDF: `{pdf_path if pdf_rendered else 'not rendered'}`
- Canonical tables: `{TABLE_DIR}`
- Canonical figures: `{FIGURE_DIR}`
"""
        (LOG_DIR / "gdtai_extension_evaluation_summary.md").write_text(markdown, encoding="utf-8")
        logging.info("Wrote screen report: %s", report)
        return

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

    ranking: pd.DataFrame | None = None
    if args.stage == "select":
        ranking = select_profile(metrics, strata, config)
        ranking.to_csv(TABLE_DIR / "selection_ranking.csv", index=False)
        decision = write_decision(args.decision, ranking, profiles, config, args.config, manifest_path, rows)
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
