#!/usr/bin/env python3
"""Repair the V4.2 NK definition with within-lane T anchors and two evidences.

This workflow is read-only with respect to every H5AD. It restores productive
TRA/TRB anchors omitted from the existing-atlas sidecar, calibrates a latent
nearest-anchor rule on held-out known cells, and requires an independent,
conservative transcriptome rule before assigning a low-weight pseudo-NK label.
"""

from __future__ import annotations

import html
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

_mps = os.environ.get("CUDA_MPS_PIPE_DIRECTORY", "")
if _mps == "/tmp/nvidia-mps" or (
    _mps and not _mps.startswith("/tmp/gdtai-v42-nk-repair-mps-")
):
    os.environ.pop("CUDA_MPS_PIPE_DIRECTORY", None)
os.environ.setdefault(
    "CUDA_MPS_PIPE_DIRECTORY",
    f"/tmp/gdtai-v42-nk-repair-mps-{os.getuid()}-{os.getpid()}",
)
os.environ.setdefault("CUDA_VISIBLE_DEVICES", "0")
os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")
os.environ.setdefault("PYTHONHASHSEED", "20260817")
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/gdtai-v42-nk-repair-numba")

import torch  # noqa: E402  # Load the CUDA runtime before RAPIDS.
import anndata as ad  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from cuml.neighbors import NearestNeighbors  # noqa: E402

from workflows.gdtai.gdtai_v4_2_integration_core import (  # noqa: E402
    h5ad_obs_frame,
    read_json,
    read_sparse_rows_genes,
    resolve,
    sha256_file,
)


CONFIG_PATH = PROJECT_ROOT / "configs/models/gdtai/v4_2_integration_execution.json"
SCRIPT_PATH = Path(__file__).resolve()
BASE_NAME = "gdtai_v4_2_nk_definition_repair"
TABLE_DIR = PROJECT_ROOT / f"Integrated_dataset/tables/gdT_prediction/{BASE_NAME}"
FIGURE_DIR = PROJECT_ROOT / f"Integrated_dataset/figures/gdT_prediction/{BASE_NAME}"
LOG_DIR = PROJECT_ROOT / f"Integrated_dataset/logs/gdT_prediction/{BASE_NAME}"
STATIC_DIR = PROJECT_ROOT / f"gdT_prediction/{BASE_NAME}"
RANDOM_SEED = 20260817
K_NEIGHBORS = 50
T_REFERENCE_CAP_PER_SOURCE = 2_000
T_VALIDATION_CAP_PER_SOURCE = 2_000
NK_VALIDATION_FRACTION = 0.20
MINIMUM_LATENT_NK_FRACTION = 0.90
LATENT_T_FPR_QUANTILE = 0.999
MAXIMUM_HELDOUT_T_LATENT_FPR = 1.0 - LATENT_T_FPR_QUANTILE
LATENT_NK_DISTANCE_QUANTILE = 0.99
MAXIMUM_TRAINING_SOURCE_FRACTION = 0.70

CORE_T_GENES = ["CD3D", "CD3E", "CD3G", "TRAT1", "BCL11B"]
NK_ADAPTOR_GENES = ["FCER1G", "TYROBP"]
NK_RECEPTOR_GENES = ["KLRD1", "NCR1", "FCGR3A", "KLRC1"]
CYTOTOXIC_AUDIT_GENES = ["NKG7", "GNLY", "PRF1", "GZMB", "XCL1", "XCL2"]
MARKER_GENES = (
    CORE_T_GENES + NK_ADAPTOR_GENES + NK_RECEPTOR_GENES + CYTOTOXIC_AUDIT_GENES
)
LOCKED_GSE169246 = (
    PROJECT_ROOT / "data/interim/extension_intake/tnk_filtered_h5ads/GSE169246.h5ad"
)

COHORT_LABELS = {
    "current_atlas": "Existing atlas reference",
    "GSE114724": "GSE114724",
    "GSE159251": "GSE159251",
    "GSE292700": "GSE292700",
    "GSE294273_GSE294274": "GSE294273/4",
    "GSE296954": "GSE296954",
}
SOURCE_COLORS = {
    "GSE114724": "#16697a",
    "GSE159251": "#d98f39",
    "GSE292700": "#8f5da2",
    "GSE294273": "#2a9d6f",
    "GSE294274": "#2a9d6f",
    "GSE296954": "#c64b44",
}


def as_bool(values: pd.Series) -> np.ndarray:
    if pd.api.types.is_bool_dtype(values) or pd.api.types.is_numeric_dtype(values):
        return values.fillna(False).astype(bool).to_numpy()
    return (
        values.astype("string")
        .fillna("")
        .str.strip()
        .str.lower()
        .isin(["true", "1", "yes", "y", "t"])
        .to_numpy()
    )


def nonempty(values: pd.Series) -> np.ndarray:
    return (
        ~values.astype("string")
        .fillna("")
        .str.strip()
        .str.lower()
        .isin(["", "nan", "none", "na"])
        .to_numpy()
    )


def configure_plotting() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.titlesize": 13,
            "axes.labelsize": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )


def finish_figure(fig: plt.Figure, path: Path) -> None:
    fig.tight_layout()
    fig.savefig(path, dpi=260, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def source_balanced_indices(
    mask: np.ndarray, sources: np.ndarray, cap: int, seed: int
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    selected: list[np.ndarray] = []
    for offset, source in enumerate(np.sort(np.unique(sources[mask]))):
        rows = np.flatnonzero(mask & (sources == source))
        if rows.size > cap:
            rows = np.sort(
                np.random.default_rng(seed + 104729 * (offset + 1)).choice(
                    rows, cap, replace=False
                )
            )
        else:
            rows = np.sort(rows)
        selected.append(rows)
    if not selected:
        return np.empty(0, dtype=np.int64)
    result = np.sort(np.concatenate(selected))
    rng.shuffle(result)
    return np.sort(result)


def split_nk_reference_validation(
    nk_mask: np.ndarray, sources: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    reference: list[np.ndarray] = []
    validation: list[np.ndarray] = []
    for offset, source in enumerate(np.sort(np.unique(sources[nk_mask]))):
        rows = np.flatnonzero(nk_mask & (sources == source))
        shuffled = np.random.default_rng(RANDOM_SEED + 7919 * (offset + 1)).permutation(
            rows
        )
        n_validation = max(1, int(round(rows.size * NK_VALIDATION_FRACTION)))
        validation.append(np.sort(shuffled[:n_validation]))
        reference.append(np.sort(shuffled[n_validation:]))
    return np.sort(np.concatenate(reference)), np.sort(np.concatenate(validation))


def split_t_reference_validation(
    t_mask: np.ndarray, sources: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    reference: list[np.ndarray] = []
    validation: list[np.ndarray] = []
    for offset, source in enumerate(np.sort(np.unique(sources[t_mask]))):
        rows = np.flatnonzero(t_mask & (sources == source))
        shuffled = np.random.default_rng(
            RANDOM_SEED + 15485863 * (offset + 1)
        ).permutation(rows)
        reference.append(
            np.sort(shuffled[: min(T_REFERENCE_CAP_PER_SOURCE, shuffled.size)])
        )
        remaining = shuffled[min(T_REFERENCE_CAP_PER_SOURCE, shuffled.size) :]
        validation.append(
            np.sort(remaining[: min(T_VALIDATION_CAP_PER_SOURCE, remaining.size)])
        )
    return np.sort(np.concatenate(reference)), np.sort(np.concatenate(validation))


def gpu_knn_score(
    reference: np.ndarray,
    reference_is_nk: np.ndarray,
    query: np.ndarray,
    n_neighbors: int = K_NEIGHBORS,
) -> tuple[np.ndarray, np.ndarray]:
    if reference.shape[0] < n_neighbors:
        raise ValueError("Reference has fewer cells than requested neighbors")
    model = NearestNeighbors(
        n_neighbors=n_neighbors,
        algorithm="brute",
        metric="euclidean",
        output_type="numpy",
    )
    model.fit(np.asarray(reference, dtype=np.float32, order="C"))
    distances, indices = model.kneighbors(
        np.asarray(query, dtype=np.float32, order="C")
    )
    distances = np.asarray(distances, dtype=np.float32)
    indices = np.asarray(indices, dtype=np.int64)
    labels = np.asarray(reference_is_nk, dtype=bool)[indices]
    nk_fraction = labels.mean(axis=1, dtype=np.float64).astype(np.float32)
    nk_distances = np.where(labels, distances, np.inf)
    nearest_nk_distance = nk_distances.min(axis=1).astype(np.float32)
    return nk_fraction, nearest_nk_distance


def select_latent_threshold(t_scores: np.ndarray) -> float:
    scores = np.asarray(t_scores, dtype=np.float32)
    candidates = np.unique(scores[scores >= MINIMUM_LATENT_NK_FRACTION])
    for threshold in np.sort(candidates):
        if float((scores >= threshold).mean()) <= MAXIMUM_HELDOUT_T_LATENT_FPR:
            return float(threshold)
    return float(np.nextafter(scores.max(), np.float32(np.inf)))


def restore_anchor_masks(
    obs: pd.DataFrame, current_path: Path
) -> tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    primary_nk = obs["primary_nk_anchor"].astype(bool).to_numpy()
    t_anchor = obs["productive_tcr_anchor"].astype(bool).to_numpy()
    doublet = obs["doublet_flag_effective"].astype(bool).to_numpy()
    current = h5ad_obs_frame(current_path, ["source_gse_id", "TRA_cdr3", "TRB_cdr3"])
    current_productive = nonempty(current["TRA_cdr3"]) | nonempty(current["TRB_cdr3"])
    current_rows = np.flatnonzero(
        obs["input_cohort_id"].astype(str).eq("current_atlas").to_numpy()
    )
    input_rows = obs.iloc[current_rows]["input_row"].to_numpy(np.int64)
    restored = current_productive[input_rows]
    t_anchor[current_rows] |= restored
    t_anchor &= ~primary_nk & ~doublet
    sources = obs["source_gse_id"].astype(str).to_numpy()
    rows = []
    for source in np.sort(np.unique(sources[primary_nk | t_anchor])):
        rows.append(
            {
                "source_gse_id": source,
                "n_primary_nk": int((primary_nk & (sources == source)).sum()),
                "n_productive_t": int((t_anchor & (sources == source)).sum()),
                "input_lane": "existing_atlas"
                if (
                    obs.loc[sources == source, "input_cohort_id"].astype(str)
                    == "current_atlas"
                ).all()
                else "new_cohort",
            }
        )
    return primary_nk, t_anchor, pd.DataFrame(rows)


def strict_expression_evidence(detected: np.ndarray) -> dict[str, np.ndarray]:
    lookup = {gene: index for index, gene in enumerate(MARKER_GENES)}
    core_t = detected[:, [lookup[g] for g in CORE_T_GENES]].sum(axis=1)
    adaptor = detected[:, [lookup[g] for g in NK_ADAPTOR_GENES]].sum(axis=1)
    receptor = detected[:, [lookup[g] for g in NK_RECEPTOR_GENES]].sum(axis=1)
    cytotoxic = detected[:, [lookup[g] for g in CYTOTOXIC_AUDIT_GENES]].sum(axis=1)
    return {
        "core_t_detected": np.asarray(core_t, dtype=np.int8),
        "nk_adaptor_detected": np.asarray(adaptor, dtype=np.int8),
        "nk_receptor_detected": np.asarray(receptor, dtype=np.int8),
        "cytotoxic_detected": np.asarray(cytotoxic, dtype=np.int8),
        "strict_expression_nk": np.asarray(
            (core_t == 0) & (adaptor == 2) & (receptor >= 1), dtype=bool
        ),
        "cytotoxic_only_rule": np.asarray(cytotoxic >= 2, dtype=bool),
    }


def extension_expression_evidence(
    obs: pd.DataFrame,
    source_contract: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    candidate_rows = np.flatnonzero(obs["candidate_eligible"].astype(bool).to_numpy())
    candidate_position = {int(row): pos for pos, row in enumerate(candidate_rows)}
    candidate = pd.DataFrame(
        {
            "integration_cell_id": obs.index[candidate_rows].astype(str),
            "global_row": candidate_rows,
            "input_cohort_id": obs.iloc[candidate_rows]["input_cohort_id"]
            .astype(str)
            .to_numpy(),
            "source_gse_id": obs.iloc[candidate_rows]["source_gse_id"]
            .astype(str)
            .to_numpy(),
            "donor_id": obs.iloc[candidate_rows]["donor_id"].astype(str).to_numpy(),
            "sample_id": obs.iloc[candidate_rows]["sample_id"].astype(str).to_numpy(),
            "library_id": obs.iloc[candidate_rows]["library_id"].astype(str).to_numpy(),
        }
    )
    for column in [
        "core_t_detected",
        "nk_adaptor_detected",
        "nk_receptor_detected",
        "cytotoxic_detected",
    ]:
        candidate[column] = np.zeros(candidate.shape[0], dtype=np.int8)
    candidate["strict_expression_nk"] = False
    candidate["cytotoxic_only_rule"] = False
    productive_rows: list[dict[str, Any]] = []

    paths = source_contract.set_index("cohort_id")["path"].to_dict()
    for cohort in sorted(set(candidate["input_cohort_id"]) - {"current_atlas"}):
        path = Path(paths[cohort])
        matrix = read_sparse_rows_genes(
            path, "X", MARKER_GENES, rows=None, row_chunk_size=25_000
        )
        detected = matrix.toarray() > 0
        evidence = strict_expression_evidence(detected)
        global_rows = np.flatnonzero(
            obs["input_cohort_id"].astype(str).eq(cohort).to_numpy()
        )
        input_rows = obs.iloc[global_rows]["input_row"].to_numpy(np.int64)
        cohort_candidate = (
            obs.iloc[global_rows]["candidate_eligible"].astype(bool).to_numpy()
        )
        for local_index, global_row in enumerate(global_rows[cohort_candidate]):
            position = candidate_position[int(global_row)]
            source_row = int(input_rows[cohort_candidate][local_index])
            for column, values in evidence.items():
                candidate.at[position, column] = values[source_row]
        productive = (
            obs.iloc[global_rows]["productive_tcr_anchor"].astype(bool).to_numpy()
        )
        if productive.any():
            selected = input_rows[productive]
            strict_rate = float(evidence["strict_expression_nk"][selected].mean())
            cytotoxic_rate = float(evidence["cytotoxic_only_rule"][selected].mean())
        else:
            strict_rate = float("nan")
            cytotoxic_rate = float("nan")
        productive_rows.append(
            {
                "cohort_id": cohort,
                "n_productive_t": int(productive.sum()),
                "strict_expression_nk_rate": strict_rate,
                "cytotoxic_only_rate": cytotoxic_rate,
            }
        )
        del matrix, detected
    return candidate, pd.DataFrame(productive_rows)


def locked_gse169246_validation() -> tuple[pd.DataFrame, pd.DataFrame]:
    columns = [
        "major_cluster",
        "has_TRA_TRB_paired",
        "has_any_ab_tcr",
        "has_any_gd_tcr",
        "doublet_flag",
        "predicted_doublet",
    ]
    obs = h5ad_obs_frame(LOCKED_GSE169246, columns)
    detected = (
        read_sparse_rows_genes(
            LOCKED_GSE169246, "X", MARKER_GENES, rows=None, row_chunk_size=25_000
        ).toarray()
        > 0
    )
    evidence = strict_expression_evidence(detected)
    author_nk = (
        obs["major_cluster"]
        .astype("string")
        .fillna("")
        .str.strip()
        .str.lower()
        .eq("nkcell")
        .to_numpy()
    )
    any_tcr = as_bool(obs["has_any_ab_tcr"]) | as_bool(obs["has_any_gd_tcr"])
    doublet = as_bool(obs["doublet_flag"]) | as_bool(obs["predicted_doublet"])
    author_nk_no_tcr = author_nk & ~any_tcr & ~doublet
    paired_ab = (
        as_bool(obs["has_TRA_TRB_paired"]) & ~as_bool(obs["has_any_gd_tcr"]) & ~doublet
    )
    groups = {
        "Author NK, no productive TCR": author_nk_no_tcr,
        "Paired alpha-beta T": paired_ab,
    }
    rows: list[dict[str, Any]] = []
    for group, mask in groups.items():
        rows.append(
            {
                "group": group,
                "n_cells": int(mask.sum()),
                "n_strict_expression_nk": int(
                    (mask & evidence["strict_expression_nk"]).sum()
                ),
                "strict_expression_nk_rate": float(
                    evidence["strict_expression_nk"][mask].mean()
                )
                if mask.any()
                else float("nan"),
                "n_cytotoxic_only": int((mask & evidence["cytotoxic_only_rule"]).sum()),
                "cytotoxic_only_rate": float(
                    evidence["cytotoxic_only_rule"][mask].mean()
                )
                if mask.any()
                else float("nan"),
            }
        )
    cell = pd.DataFrame(
        {
            "cell_id": obs.index.astype(str),
            "author_nk_no_productive_tcr": author_nk_no_tcr,
            "paired_ab_t": paired_ab,
            **evidence,
        }
    )
    return pd.DataFrame(rows), cell


def source_balanced_training_weights(candidate: pd.DataFrame) -> np.ndarray:
    confident = candidate["nk_confident"].to_numpy(bool)
    sources = candidate["source_gse_id"].astype(str).to_numpy()
    levels, counts = np.unique(sources[confident], return_counts=True)
    weights = np.zeros(candidate.shape[0], dtype=np.float32)
    if levels.size < 2:
        return weights
    weights[confident] = 0.25
    effective = counts.astype(np.float64) * 0.25
    dominant_index = int(np.argmax(effective))
    other_total = float(effective.sum() - effective[dominant_index])
    maximum_dominant = (
        MAXIMUM_TRAINING_SOURCE_FRACTION
        * other_total
        / (1.0 - MAXIMUM_TRAINING_SOURCE_FRACTION)
    )
    if effective[dominant_index] > maximum_dominant:
        dominant = levels[dominant_index]
        weights[confident & (sources == dominant)] = (
            maximum_dominant / counts[dominant_index]
        )
    return weights


def make_figures(
    anchors: pd.DataFrame,
    latent_validation: pd.DataFrame,
    candidate: pd.DataFrame,
    locked: pd.DataFrame,
    marker_summary: pd.DataFrame,
) -> None:
    configure_plotting()

    ordered = anchors.sort_values(["input_lane", "source_gse_id"])
    fig, ax = plt.subplots(figsize=(9.2, 5.6))
    y = np.arange(ordered.shape[0])
    ax.barh(y, ordered["n_productive_t"], color="#247ba0", label="Productive T anchors")
    ax.barh(
        y,
        ordered["n_primary_nk"],
        left=ordered["n_productive_t"],
        color="#2a7f62",
        label="Primary NK anchors",
    )
    ax.set_yticks(y, ordered["source_gse_id"])
    ax.set_xscale("symlog", linthresh=10)
    ax.set_xlabel("Anchor cells, symlog scale")
    ax.set_title("Restored anchor classes now coexist in the existing-atlas lane")
    ax.legend(frameon=False)
    finish_figure(fig, FIGURE_DIR / "restored_anchor_counts_by_source.png")

    fig, ax = plt.subplots(figsize=(8.6, 4.8))
    colors = {
        "Held-out primary NK": "#2a7f62",
        "Held-out productive T": "#247ba0",
        "Eligible candidate": "#d98f39",
    }
    for group in ["Held-out productive T", "Held-out primary NK", "Eligible candidate"]:
        values = latent_validation.loc[
            latent_validation["group"].eq(group), "latent_nk_fraction"
        ].to_numpy()
        if values.size > 50_000:
            values = np.random.default_rng(RANDOM_SEED).choice(
                values, 50_000, replace=False
            )
        ax.hist(
            values,
            bins=np.linspace(0, 1, 51),
            density=True,
            histtype="step",
            linewidth=2,
            color=colors[group],
            label=f"{group} ({values.size:,})",
        )
    threshold = float(candidate["latent_threshold"].iloc[0])
    ax.axvline(
        threshold,
        color="#b3261e",
        linestyle="--",
        linewidth=1.5,
        label=f"Frozen threshold {threshold:.3f}",
    )
    ax.set_yscale("log")
    ax.set_xlabel("NK fraction among 50 nearest balanced anchors")
    ax.set_ylabel("Density, log scale")
    ax.set_title("Latent NK evidence calibrated on held-out known cells")
    ax.legend(frameon=False, fontsize=8)
    finish_figure(fig, FIGURE_DIR / "latent_anchor_score_distributions.png")

    funnel = pd.DataFrame(
        {
            "Step": [
                "Eligible no-TCR candidates",
                "Strict expression NK",
                "Latent NK evidence",
                "Both evidences",
                "Low-weight future-training rows",
            ],
            "Cells": [
                candidate.shape[0],
                int(candidate.strict_expression_nk.sum()),
                int(candidate.latent_pass.sum()),
                int(candidate.nk_confident.sum()),
                int(candidate.selected_for_future_training.sum()),
            ],
        }
    )
    fig, ax = plt.subplots(figsize=(8.4, 4.4))
    bars = ax.barh(
        np.arange(funnel.shape[0]),
        funnel["Cells"],
        color=["#8c979e", "#6b9e78", "#d98f39", "#2a7f62", "#17324d"],
    )
    ax.set_yticks(np.arange(funnel.shape[0]), funnel["Step"])
    ax.invert_yaxis()
    ax.bar_label(
        bars, labels=[f"{value:,}" for value in funnel["Cells"]], padding=4, fontsize=9
    )
    ax.set_xlabel("Cells")
    ax.set_title("Conservative NK evidence funnel")
    finish_figure(fig, FIGURE_DIR / "candidate_evidence_funnel.png")

    source = (
        candidate.groupby("source_gse_id", observed=True)
        .agg(
            n_candidates=("integration_cell_id", "size"),
            n_expression=("strict_expression_nk", "sum"),
            n_latent=("latent_pass", "sum"),
            n_confident=("nk_confident", "sum"),
            n_training=("selected_for_future_training", "sum"),
        )
        .reset_index()
    )
    fig, ax = plt.subplots(figsize=(9.0, 4.7))
    x = np.arange(source.shape[0])
    width = 0.20
    for offset, column, label, color in [
        (-1.5, "n_candidates", "Eligible", "#b8c1c7"),
        (-0.5, "n_expression", "Expression", "#6b9e78"),
        (0.5, "n_latent", "Latent", "#d98f39"),
        (1.5, "n_confident", "Both", "#2a7f62"),
    ]:
        ax.bar(
            x + offset * width, source[column], width=width, label=label, color=color
        )
    ax.set_xticks(x, source["source_gse_id"], rotation=25, ha="right")
    ax.set_ylabel("Cells")
    ax.set_title("NK evidence by development source")
    ax.set_yscale("symlog", linthresh=10)
    ax.legend(frameon=False, ncol=4)
    finish_figure(fig, FIGURE_DIR / "candidate_counts_by_source.png")

    fig, ax = plt.subplots(figsize=(7.5, 4.4))
    x = np.arange(locked.shape[0])
    width = 0.34
    ax.bar(
        x - width / 2,
        100 * locked["strict_expression_nk_rate"],
        width=width,
        color="#2a7f62",
        label="Strict lineage rule",
    )
    ax.bar(
        x + width / 2,
        100 * locked["cytotoxic_only_rate"],
        width=width,
        color="#c64b44",
        label="Cytotoxic-only rule",
    )
    ax.set_xticks(x, locked["group"], rotation=12, ha="right")
    ax.set_ylabel("Cells passing rule (%)")
    ax.set_title("Locked GSE169246: lineage evidence versus cytotoxicity")
    ax.legend(frameon=False)
    ax.grid(axis="y", color="#e1e7ea", linewidth=0.6)
    finish_figure(fig, FIGURE_DIR / "locked_gse169246_expression_validation.png")

    display = marker_summary.set_index("group")
    fig, ax = plt.subplots(figsize=(8.8, 4.8))
    matrix = (
        100
        * display[
            [
                "core_t_zero_rate",
                "both_adaptors_rate",
                "nk_receptor_rate",
                "cytotoxic_two_gene_rate",
            ]
        ].to_numpy()
    )
    image = ax.imshow(matrix, aspect="auto", cmap="viridis", vmin=0, vmax=100)
    ax.set_yticks(np.arange(display.shape[0]), display.index)
    ax.set_xticks(
        np.arange(4),
        ["No core T genes", "Both adaptors", "NK receptor", "2+ cytotoxic genes"],
        rotation=20,
        ha="right",
    )
    for row in range(matrix.shape[0]):
        for column in range(matrix.shape[1]):
            ax.text(
                column,
                row,
                f"{matrix[row, column]:.1f}%",
                ha="center",
                va="center",
                color="white" if matrix[row, column] < 55 else "#17212b",
                fontsize=8,
            )
    fig.colorbar(image, ax=ax, label="Cells (%)")
    ax.set_title("Expression components supporting the conservative label")
    finish_figure(fig, FIGURE_DIR / "marker_component_heatmap.png")

    diagnostic_path = (
        PROJECT_ROOT
        / "Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_nk_reference_diagnostics/diagnostic_visualization_sample.csv.gz"
    )
    if diagnostic_path.exists():
        sample = pd.read_csv(
            diagnostic_path,
            usecols=["integration_cell_id", "UMAP1", "UMAP2", "diagnostic_role"],
        )
        tiers = candidate.set_index("integration_cell_id")["nk_definition"].to_dict()
        sample["nk_definition"] = (
            sample["integration_cell_id"].map(tiers).fillna("Not an eligible candidate")
        )
        fig, ax = plt.subplots(figsize=(8.0, 6.4))
        background = ~sample["nk_definition"].isin(
            ["NK_CONFIDENT", "NK_EXPRESSION_ONLY", "NK_LATENT_ONLY"]
        )
        ax.scatter(
            sample.loc[background, "UMAP1"],
            sample.loc[background, "UMAP2"],
            s=0.45,
            c="#d5dadd",
            alpha=0.10,
            edgecolors="none",
            rasterized=True,
        )
        for label, color, size in [
            ("NK_LATENT_ONLY", "#d98f39", 1.2),
            ("NK_EXPRESSION_ONLY", "#6b9e78", 1.4),
            ("NK_CONFIDENT", "#17324d", 2.0),
        ]:
            mask = sample["nk_definition"].eq(label)
            ax.scatter(
                sample.loc[mask, "UMAP1"],
                sample.loc[mask, "UMAP2"],
                s=size,
                c=color,
                alpha=0.72,
                edgecolors="none",
                rasterized=True,
                label=f"{label} ({mask.sum():,} sampled)",
            )
        ax.set_title("Repaired NK evidence on the existing diagnostic UMAP")
        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.spines[["top", "right", "bottom", "left"]].set_visible(False)
        ax.legend(frameon=False, fontsize=8, markerscale=5)
        finish_figure(fig, FIGURE_DIR / "nk_definition_umap.png")


def table_html(frame: pd.DataFrame, columns: list[str] | None = None) -> str:
    shown = frame if columns is None else frame[columns]
    return shown.to_html(index=False, border=0, classes="data-table", escape=True)


def render_report(
    summary: dict[str, Any],
    anchor_summary: pd.DataFrame,
    source_summary: pd.DataFrame,
    locked_summary: pd.DataFrame,
    productive_summary: pd.DataFrame,
) -> Path:
    prefix = f"../../Integrated_dataset/figures/gdT_prediction/{BASE_NAME}"
    source_display = source_summary.copy()
    for column in [
        "expression_rate",
        "latent_rate",
        "confident_rate",
        "training_fraction",
        "effective_training_source_fraction",
    ]:
        source_display[column] = source_display[column].map(lambda x: f"{x:.2%}")
    locked_display = locked_summary.copy()
    for column in ["strict_expression_nk_rate", "cytotoxic_only_rate"]:
        locked_display[column] = locked_display[column].map(lambda x: f"{x:.2%}")
    productive_display = productive_summary.copy()
    for column in ["strict_expression_nk_rate", "cytotoxic_only_rate"]:
        productive_display[column] = productive_display[column].map(
            lambda x: f"{x:.2%}"
        )
    training_status = (
        "ready for low-weight development use with capped source influence"
        if summary["training_subset_ready"]
        else "not yet training-ready"
    )
    document = f"""<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>gdTAI V4.2 NK definition repair</title><style>
:root{{--ink:#17212b;--navy:#17324d;--teal:#197278;--red:#b3261e;--amber:#b86e16;--line:#ccd6dc;--soft:#f3f6f7}}
*{{box-sizing:border-box}}body{{font-family:Arial,Helvetica,sans-serif;color:var(--ink);line-height:1.48;max-width:1120px;margin:0 auto;padding:30px 42px;background:white}}h1{{font-size:29px;color:var(--navy);border-bottom:4px solid var(--teal);padding-bottom:12px}}h2{{font-size:20px;color:var(--navy);margin-top:31px}}p,li{{font-size:14px}}.status{{border-left:6px solid var(--teal);background:#eef8f7;padding:14px 18px;margin:18px 0}}.warning{{border-left:6px solid var(--amber);background:#fff9ef;padding:13px 17px;margin:18px 0}}.cards{{display:grid;grid-template-columns:repeat(4,1fr);gap:10px;margin:18px 0}}.card{{border:1px solid var(--line);padding:12px;background:var(--soft)}}.card b{{display:block;font-size:21px;color:var(--navy)}}.card span{{font-size:11px;color:#52636e}}.grid{{display:grid;grid-template-columns:1fr 1fr;gap:18px}}figure{{margin:18px 0;break-inside:avoid}}img{{display:block;width:100%;height:auto;margin:0 auto}}figcaption{{font-size:11px;color:#52636e;text-align:center;margin-top:7px}}.table-wrap{{overflow-x:auto;margin:14px 0 24px}}table{{border-collapse:collapse;width:100%;font-size:9px}}th{{background:var(--navy);color:white;text-align:left}}th,td{{padding:6px;border:1px solid var(--line);vertical-align:top;overflow-wrap:anywhere}}tr:nth-child(even){{background:var(--soft)}}code{{background:#edf1f4;padding:1px 4px}}
@media(max-width:800px){{body{{padding:18px}}.cards,.grid{{grid-template-columns:1fr}}}}@media print{{@page{{size:A4;margin:9mm}}body{{max-width:none;padding:0}}h1{{font-size:23px}}h2{{font-size:17px;break-after:avoid}}p,li{{font-size:10px}}.cards{{grid-template-columns:repeat(4,1fr)}}.card b{{font-size:16px}}.card span{{font-size:8px}}.grid{{display:block}}figure{{break-inside:avoid}}figure img{{max-height:145mm;object-fit:contain}}table{{font-size:7px}}tr{{break-inside:avoid}}.table-wrap{{overflow:visible}}}}
</style></head><body>
<h1>gdTAI V4.2 NK definition repair</h1>
<p>This read-only repair addresses the cohort-confounded pseudo-NK failure. It does not alter the canonical atlas, refit scVI, train gdTAI, or use locked cohorts for threshold selection.</p>
<div class="status"><b>Result: {html.escape(summary["result"])}.</b> The repaired definition identifies {summary["n_confident_nk"]:,} cells with both latent and strict lineage evidence. The source-balanced subset is {training_status}.</div>
<div class="cards"><div class="card"><b>{summary["n_restored_existing_t_anchors"]:,}</b><span>restored existing-atlas T anchors</span></div><div class="card"><b>{summary["n_primary_nk"]:,}</b><span>dual-annotation NK anchors</span></div><div class="card"><b>{summary["n_confident_nk"]:,}</b><span>two-evidence pseudo-NK cells</span></div><div class="card"><b>{summary["n_training_subset"]:,}</b><span>capped-weight training rows</span></div></div>

<h2>1. Defect repaired</h2><p>The prior sidecar mapped no productive-T anchors into the existing-atlas lane even though its source object contains {summary["n_restored_existing_t_anchors"]:,} cells with productive TRA or TRB evidence. Restoring these cells means NK and T anchors now coexist within the same old-atlas lane; productive-T anchors also remain present in every new cohort.</p><figure><img src="{prefix}/restored_anchor_counts_by_source.png" alt="Restored anchors"><figcaption>Anchor counts use expression-independent TCR metadata and the frozen dual-annotation NK labels.</figcaption></figure><div class="table-wrap">{table_html(anchor_summary)}</div>

<h2>2. NK definition</h2><p>A development cell is called <b>NK_CONFIDENT</b> only when it has no productive TCR, no doublet flag, passes the latent nearest-anchor threshold, remains within the known-NK distance envelope, expresses both <code>FCER1G</code> and <code>TYROBP</code> plus at least one of <code>KLRD1/NCR1/FCGR3A/KLRC1</code>, and has no detected <code>CD3D/CD3E/CD3G/TRAT1/BCL11B</code>. Cytotoxic genes are audit-only and cannot establish NK identity.</p><div class="grid"><figure><img src="{prefix}/latent_anchor_score_distributions.png" alt="Latent calibration"><figcaption>The threshold is the lowest discrete 50-neighbor score whose empirical held-out productive-T FPR is at most 0.1%; score ties fail closed.</figcaption></figure><figure><img src="{prefix}/candidate_evidence_funnel.png" alt="Evidence funnel"><figcaption>Only the intersection of independent latent and lineage-expression evidence is accepted.</figcaption></figure></div>

<h2>3. Source-level result</h2><div class="grid"><figure><img src="{prefix}/candidate_counts_by_source.png" alt="Counts by source"><figcaption>All development sources are reported separately; source imbalance is controlled by capped effective weights for any future low-weight training use.</figcaption></figure><figure><img src="{prefix}/nk_definition_umap.png" alt="NK definition UMAP"><figcaption>Diagnostic UMAP overlay; the full label table, not the visualization sample, defines counts.</figcaption></figure></div><div class="table-wrap">{table_html(source_display)}</div>

<h2>4. Locked validation and specificity</h2><p>GSE169246 was not integrated and did not control any threshold. Its author-defined NK cells and paired alpha-beta T cells test only the fixed strict expression rule. The cytotoxic-only comparator demonstrates why NKG7/GNLY-like activity is not accepted as NK identity.</p><figure><img src="{prefix}/locked_gse169246_expression_validation.png" alt="Locked validation"><figcaption>Fixed-rule passage in locked author-NK and paired-abT populations.</figcaption></figure><div class="table-wrap">{table_html(locked_display)}</div><h3>Development productive-T controls</h3><div class="table-wrap">{table_html(productive_display)}</div>

<h2>5. Marker components</h2><figure><img src="{prefix}/marker_component_heatmap.png" alt="Marker components"><figcaption>Cytotoxicity is frequent in T cells and is shown only as a warning diagnostic.</figcaption></figure>

<div class="warning"><b>Scope:</b> these are conservative NK reference labels, not a finished gdTAI V4 model. No classifier was fitted, no threshold was tuned on GSE169246 or the candidate counts, no H5AD was changed, and nothing was pushed to GitHub.</div>
<h2>Reproducibility</h2><ul><li>Latent SHA-256: <code>{html.escape(summary["latent_sha256"])}</code></li><li>Script SHA-256: <code>{html.escape(summary["script_sha256"])}</code></li><li>Random seed: <code>{RANDOM_SEED}</code></li><li>Neighbor search: <code>{html.escape(summary["knn_algorithm"])}</code></li><li>Bit-identical exact-GPU repeat: <b>{str(summary["deterministic_exact_repeat_pass"]).lower()}</b></li><li>GPU: <code>{html.escape(summary["gpu"])}</code>; CPU fallback: <b>no</b></li></ul>
</body></html>"""
    output = STATIC_DIR / "index.html"
    output.write_text(document, encoding="utf-8")
    return output


def main() -> None:
    if not torch.cuda.is_available():
        raise RuntimeError(
            "CUDA is required for latent nearest-anchor scoring; CPU fallback is forbidden"
        )
    for directory in [TABLE_DIR, FIGURE_DIR, LOG_DIR, STATIC_DIR]:
        directory.mkdir(parents=True, exist_ok=True)

    config = read_json(CONFIG_PATH)
    ssd_root = Path(config["resources"]["ssd_root"])
    staged_path = ssd_root / "development_hvg_counts.h5ad"
    latent_path = ssd_root / "X_scVI.npy"
    base_log_dir = resolve(config["outputs"]["log_dir"])
    base_table_dir = resolve(config["outputs"]["table_dir"])
    fit_summary = read_json(base_log_dir / "fit_summary.json")
    if sha256_file(latent_path) != fit_summary["latent_sha256"]:
        raise RuntimeError("Saved latent checksum differs from the completed V4.2 fit")

    source_contract = pd.read_csv(base_table_dir / "development_input_contract.csv")
    before_state = {
        str(path): (Path(path).stat().st_size, Path(path).stat().st_mtime_ns)
        for path in source_contract["path"]
    }
    before_state[str(LOCKED_GSE169246)] = (
        LOCKED_GSE169246.stat().st_size,
        LOCKED_GSE169246.stat().st_mtime_ns,
    )

    backed = ad.read_h5ad(staged_path, backed="r")
    obs = backed.obs.copy()
    backed.file.close()
    locked_ids = {"GSE169246", "GSE315928", "GSE121636_GSE121637"}
    if obs["input_cohort_id"].astype(str).isin(locked_ids).any():
        raise RuntimeError("Locked cohort leaked into the development sidecar")
    latent = np.load(latent_path, mmap_mode="r")
    if (
        latent.shape != (obs.shape[0], 30)
        or not np.isfinite(np.asarray(latent[:1000])).all()
    ):
        raise RuntimeError("Saved latent is malformed")

    current_path = Path(
        source_contract.loc[
            source_contract["cohort_id"].eq("current_atlas"), "path"
        ].iloc[0]
    )
    primary_nk, productive_t, anchor_summary = restore_anchor_masks(obs, current_path)
    sources = obs["source_gse_id"].astype(str).to_numpy()
    nk_reference, nk_validation = split_nk_reference_validation(primary_nk, sources)
    t_reference, t_validation = split_t_reference_validation(productive_t, sources)
    if (
        nk_reference.size == 0
        or nk_validation.size == 0
        or t_reference.size == 0
        or t_validation.size == 0
    ):
        raise RuntimeError(
            "Anchor split produced an empty reference or validation class"
        )

    reference_rows = np.concatenate([nk_reference, t_reference])
    reference_labels = np.concatenate(
        [np.ones(nk_reference.size, dtype=bool), np.zeros(t_reference.size, dtype=bool)]
    )
    validation_rows = np.concatenate([nk_validation, t_validation])
    validation_groups = np.concatenate(
        [
            np.repeat("Held-out primary NK", nk_validation.size),
            np.repeat("Held-out productive T", t_validation.size),
        ]
    )
    validation_score, validation_distance = gpu_knn_score(
        np.asarray(latent[reference_rows]),
        reference_labels,
        np.asarray(latent[validation_rows]),
    )
    validation_score_repeat, validation_distance_repeat = gpu_knn_score(
        np.asarray(latent[reference_rows]),
        reference_labels,
        np.asarray(latent[validation_rows]),
    )
    if not (
        np.array_equal(validation_score, validation_score_repeat)
        and np.array_equal(validation_distance, validation_distance_repeat)
    ):
        raise RuntimeError("Exact GPU validation-neighbor repeat was not bit-identical")
    t_scores = validation_score[validation_groups == "Held-out productive T"]
    threshold = select_latent_threshold(t_scores)
    nk_distances = validation_distance[validation_groups == "Held-out primary NK"]
    finite_nk_distances = nk_distances[np.isfinite(nk_distances)]
    if finite_nk_distances.size == 0:
        raise RuntimeError("No held-out primary NK cell had an NK reference neighbor")
    distance_cap = float(
        np.quantile(finite_nk_distances, LATENT_NK_DISTANCE_QUANTILE, method="higher")
    )

    candidate, productive_expression = extension_expression_evidence(
        obs, source_contract
    )
    candidate_rows = candidate["global_row"].to_numpy(np.int64)
    final_t_reference = source_balanced_indices(
        productive_t, sources, T_REFERENCE_CAP_PER_SOURCE, RANDOM_SEED
    )
    final_reference_rows = np.concatenate(
        [np.flatnonzero(primary_nk), final_t_reference]
    )
    final_reference_labels = np.concatenate(
        [
            np.ones(int(primary_nk.sum()), dtype=bool),
            np.zeros(final_t_reference.size, dtype=bool),
        ]
    )
    candidate_score, candidate_distance = gpu_knn_score(
        np.asarray(latent[final_reference_rows]),
        final_reference_labels,
        np.asarray(latent[candidate_rows]),
    )
    candidate_score_repeat, candidate_distance_repeat = gpu_knn_score(
        np.asarray(latent[final_reference_rows]),
        final_reference_labels,
        np.asarray(latent[candidate_rows]),
    )
    if not (
        np.array_equal(candidate_score, candidate_score_repeat)
        and np.array_equal(candidate_distance, candidate_distance_repeat)
    ):
        raise RuntimeError("Exact GPU candidate-neighbor repeat was not bit-identical")
    candidate["latent_nk_fraction"] = candidate_score
    candidate["nearest_nk_distance"] = candidate_distance
    candidate["latent_threshold"] = threshold
    candidate["nk_distance_cap"] = distance_cap
    candidate["latent_pass"] = (candidate_score >= threshold) & (
        candidate_distance <= distance_cap
    )
    candidate["nk_confident"] = candidate["strict_expression_nk"].astype(
        bool
    ) & candidate["latent_pass"].astype(bool)
    candidate["nk_definition"] = "UNRESOLVED_NO_TCR"
    candidate.loc[
        candidate["latent_pass"] & ~candidate["strict_expression_nk"], "nk_definition"
    ] = "NK_LATENT_ONLY"
    candidate.loc[
        candidate["strict_expression_nk"] & ~candidate["latent_pass"], "nk_definition"
    ] = "NK_EXPRESSION_ONLY"
    candidate.loc[candidate["nk_confident"], "nk_definition"] = "NK_CONFIDENT"
    candidate["future_training_weight"] = source_balanced_training_weights(candidate)
    candidate["selected_for_future_training"] = candidate["future_training_weight"] > 0

    candidate_validation = pd.DataFrame(
        {
            "group": "Eligible candidate",
            "source_gse_id": candidate["source_gse_id"],
            "latent_nk_fraction": candidate_score,
            "nearest_nk_distance": candidate_distance,
        }
    )
    latent_validation = pd.DataFrame(
        {
            "group": validation_groups,
            "source_gse_id": sources[validation_rows],
            "latent_nk_fraction": validation_score,
            "nearest_nk_distance": validation_distance,
        }
    )
    latent_validation = pd.concat(
        [latent_validation, candidate_validation], ignore_index=True
    )

    locked_summary, locked_cells = locked_gse169246_validation()
    source_summary = (
        candidate.groupby("source_gse_id", observed=True)
        .agg(
            n_candidates=("integration_cell_id", "size"),
            n_expression=("strict_expression_nk", "sum"),
            n_latent=("latent_pass", "sum"),
            n_confident=("nk_confident", "sum"),
            n_training=("selected_for_future_training", "sum"),
            effective_training_weight=("future_training_weight", "sum"),
            n_donors=("donor_id", "nunique"),
            n_libraries=("library_id", "nunique"),
        )
        .reset_index()
    )
    source_summary["expression_rate"] = (
        source_summary["n_expression"] / source_summary["n_candidates"]
    )
    source_summary["latent_rate"] = (
        source_summary["n_latent"] / source_summary["n_candidates"]
    )
    source_summary["confident_rate"] = (
        source_summary["n_confident"] / source_summary["n_candidates"]
    )
    source_summary["training_fraction"] = (
        source_summary["n_training"] / source_summary["n_candidates"]
    )
    source_summary["effective_training_source_fraction"] = (
        source_summary["effective_training_weight"]
        / source_summary["effective_training_weight"].sum()
    )

    marker_groups = {
        "All eligible candidates": np.ones(candidate.shape[0], dtype=bool),
        "NK_CONFIDENT": candidate["nk_confident"].to_numpy(bool),
        "NK_EXPRESSION_ONLY": candidate["nk_definition"]
        .eq("NK_EXPRESSION_ONLY")
        .to_numpy(),
        "NK_LATENT_ONLY": candidate["nk_definition"].eq("NK_LATENT_ONLY").to_numpy(),
    }
    marker_rows = []
    for group, mask in marker_groups.items():
        marker_rows.append(
            {
                "group": group,
                "n_cells": int(mask.sum()),
                "core_t_zero_rate": float(
                    (candidate.loc[mask, "core_t_detected"] == 0).mean()
                )
                if mask.any()
                else float("nan"),
                "both_adaptors_rate": float(
                    (candidate.loc[mask, "nk_adaptor_detected"] == 2).mean()
                )
                if mask.any()
                else float("nan"),
                "nk_receptor_rate": float(
                    (candidate.loc[mask, "nk_receptor_detected"] >= 1).mean()
                )
                if mask.any()
                else float("nan"),
                "cytotoxic_two_gene_rate": float(
                    (candidate.loc[mask, "cytotoxic_detected"] >= 2).mean()
                )
                if mask.any()
                else float("nan"),
            }
        )
    marker_summary = pd.DataFrame(marker_rows)

    t_validation_mask = latent_validation["group"].eq("Held-out productive T")
    nk_validation_mask = latent_validation["group"].eq("Held-out primary NK")
    latent_t_fpr = float(
        (
            latent_validation.loc[t_validation_mask, "latent_nk_fraction"] >= threshold
        ).mean()
    )
    latent_nk_recall = float(
        (
            (
                latent_validation.loc[nk_validation_mask, "latent_nk_fraction"]
                >= threshold
            )
            & (
                latent_validation.loc[nk_validation_mask, "nearest_nk_distance"]
                <= distance_cap
            )
        ).mean()
    )
    locked_nk_rate = float(
        locked_summary.loc[
            locked_summary["group"].eq("Author NK, no productive TCR"),
            "strict_expression_nk_rate",
        ].iloc[0]
    )
    locked_ab_rate = float(
        locked_summary.loc[
            locked_summary["group"].eq("Paired alpha-beta T"),
            "strict_expression_nk_rate",
        ].iloc[0]
    )
    confident_sources = int((source_summary["n_confident"] > 0).sum())
    training_sources = int((source_summary["n_training"] > 0).sum())
    training_n = int(candidate["selected_for_future_training"].sum())
    training_dominance = float(
        source_summary["effective_training_source_fraction"].max()
    )
    label_gate_pass = (
        latent_t_fpr <= MAXIMUM_HELDOUT_T_LATENT_FPR + 1e-12
        and locked_ab_rate <= 0.01
        and locked_nk_rate >= 0.20
        and int(candidate["nk_confident"].sum()) >= 100
    )
    training_ready = (
        label_gate_pass
        and training_n >= 100
        and training_sources >= 2
        and training_dominance <= MAXIMUM_TRAINING_SOURCE_FRACTION + 1e-6
    )
    result = (
        "PASS_NK_REFERENCE_READY"
        if training_ready
        else (
            "PASS_CONSERVATIVE_NK_LABELS_ONLY"
            if label_gate_pass
            else "FAIL_NK_DEFINITION_GUARDRAILS"
        )
    )

    candidate.to_csv(
        TABLE_DIR / "candidate_nk_evidence.csv.gz", index=False, compression="gzip"
    )
    candidate.loc[candidate["selected_for_future_training"]].to_csv(
        TABLE_DIR / "source_balanced_pseudo_nk_for_future_training.csv.gz",
        index=False,
        compression="gzip",
    )
    anchor_summary.to_csv(
        TABLE_DIR / "restored_anchor_counts_by_source.csv", index=False
    )
    source_summary.to_csv(TABLE_DIR / "candidate_nk_counts_by_source.csv", index=False)
    latent_validation.to_csv(
        TABLE_DIR / "latent_anchor_validation.csv.gz", index=False, compression="gzip"
    )
    locked_summary.to_csv(
        TABLE_DIR / "locked_gse169246_validation_summary.csv", index=False
    )
    locked_cells.to_csv(
        TABLE_DIR / "locked_gse169246_validation_cells.csv.gz",
        index=False,
        compression="gzip",
    )
    productive_expression.to_csv(
        TABLE_DIR / "development_productive_t_expression_specificity.csv", index=False
    )
    marker_summary.to_csv(TABLE_DIR / "marker_component_summary.csv", index=False)

    make_figures(
        anchor_summary, latent_validation, candidate, locked_summary, marker_summary
    )
    after_state = {
        path: (Path(path).stat().st_size, Path(path).stat().st_mtime_ns)
        for path in before_state
    }
    if before_state != after_state:
        raise RuntimeError(
            "A source H5AD size or modification time changed during the read-only repair"
        )

    summary = {
        "result": result,
        "n_cells_sidecar": int(obs.shape[0]),
        "n_primary_nk": int(primary_nk.sum()),
        "n_productive_t_anchors_total": int(productive_t.sum()),
        "n_restored_existing_t_anchors": int(
            productive_t[
                obs["input_cohort_id"].astype(str).eq("current_atlas").to_numpy()
            ].sum()
        ),
        "n_eligible_candidates": int(candidate.shape[0]),
        "n_expression_nk": int(candidate["strict_expression_nk"].sum()),
        "n_latent_nk": int(candidate["latent_pass"].sum()),
        "n_confident_nk": int(candidate["nk_confident"].sum()),
        "n_confident_sources": confident_sources,
        "n_training_subset": training_n,
        "n_training_sources": training_sources,
        "training_maximum_source_fraction": training_dominance,
        "training_effective_weight": float(candidate["future_training_weight"].sum()),
        "training_subset_ready": training_ready,
        "latent_threshold": threshold,
        "latent_nk_distance_cap": distance_cap,
        "heldout_productive_t_latent_fpr": latent_t_fpr,
        "heldout_primary_nk_latent_recall": latent_nk_recall,
        "locked_author_nk_expression_recall": locked_nk_rate,
        "locked_paired_abt_expression_fpr": locked_ab_rate,
        "label_gate_pass": label_gate_pass,
        "locked_cohorts_used_for_threshold_selection": False,
        "classifier_fitting_performed": False,
        "source_h5ad_mutation_performed": False,
        "github_push_performed": False,
        "gpu": torch.cuda.get_device_name(0),
        "knn_algorithm": "exact_gpu_brute",
        "deterministic_exact_repeat_pass": True,
        "cpu_fallback": False,
        "latent_sha256": fit_summary["latent_sha256"],
        "script_sha256": sha256_file(SCRIPT_PATH),
        "random_seed": RANDOM_SEED,
        "strict_expression_rule": {
            "no_detected_core_t_genes": CORE_T_GENES,
            "both_detected_nk_adaptors": NK_ADAPTOR_GENES,
            "at_least_one_detected_nk_receptor": NK_RECEPTOR_GENES,
            "cytotoxic_audit_only": CYTOTOXIC_AUDIT_GENES,
        },
    }
    (LOG_DIR / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    thresholds = {
        "k_neighbors": K_NEIGHBORS,
        "minimum_latent_nk_fraction": MINIMUM_LATENT_NK_FRACTION,
        "productive_t_fpr_quantile": LATENT_T_FPR_QUANTILE,
        "maximum_heldout_t_latent_fpr": MAXIMUM_HELDOUT_T_LATENT_FPR,
        "selected_latent_threshold": threshold,
        "known_nk_distance_quantile": LATENT_NK_DISTANCE_QUANTILE,
        "selected_nk_distance_cap": distance_cap,
        "maximum_training_source_fraction": MAXIMUM_TRAINING_SOURCE_FRACTION,
        "t_reference_cap_per_source": T_REFERENCE_CAP_PER_SOURCE,
        "t_validation_cap_per_source": T_VALIDATION_CAP_PER_SOURCE,
        "knn_algorithm": "exact_gpu_brute",
    }
    (LOG_DIR / "thresholds.json").write_text(
        json.dumps(thresholds, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    report_html = render_report(
        summary, anchor_summary, source_summary, locked_summary, productive_expression
    )
    chrome = shutil.which("google-chrome") or shutil.which("google-chrome-stable")
    if chrome is None:
        raise RuntimeError("Google Chrome is required for PDF export")
    profile = Path("/tmp/gdtai-v42-nk-definition-repair-chrome")
    profile.mkdir(parents=True, exist_ok=True)
    pdf_path = STATIC_DIR / f"{BASE_NAME}_report.pdf"
    subprocess.run(
        [
            chrome,
            "--headless",
            "--no-sandbox",
            "--disable-gpu",
            "--disable-dev-shm-usage",
            "--disable-breakpad",
            "--disable-crash-reporter",
            "--allow-file-access-from-files",
            "--no-pdf-header-footer",
            f"--user-data-dir={profile}",
            f"--print-to-pdf={pdf_path}",
            report_html.resolve().as_uri(),
        ],
        cwd=PROJECT_ROOT,
        check=True,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
