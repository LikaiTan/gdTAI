#!/usr/bin/env python3
"""Refresh the two report-critical raw Phase 4 scatter panels.

These report panels should reflect the full integrated dataset rather than a
tiny speed-oriented subset:
- raw TRAB-versus-TRD colored by raw TRD-minus-TRAB
- raw TRAB-versus-TRD with paired TRA/TRB in red and all other cells in blue
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


from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
FIGURE_DIR = PROJECT_ROOT / "Integrated_dataset" / "figures"
INTEGRATED_H5AD = PROJECT_ROOT / "high_speed_temp" / "Integrated_dataset" / "integrated.h5ad"

FIGURE_DPI = 300
POINT_SIZE = 0.35
ALPHA = 0.40


def read_full_obs_frame() -> pd.DataFrame:
    """Read raw scores and paired-TRA/TRB status from the integrated milestone."""

    def decode_string_array(values: np.ndarray) -> np.ndarray:
        return np.asarray(
            [value.decode("utf-8") if isinstance(value, bytes) else str(value) for value in values],
            dtype=object,
        )

    with h5py.File(INTEGRATED_H5AD, "r") as handle:
        obs = handle["obs"]
        trab = np.asarray(obs["phase4_trab_score"][:], dtype=np.float32)
        trd = np.asarray(obs["phase4_trd_score"][:], dtype=np.float32)
        trd_minus_trab = np.asarray(obs["phase4_trd_minus_trab"][:], dtype=np.float32)
        tra_cdr3 = decode_string_array(obs["TRA_cdr3"][:]) if "TRA_cdr3" in obs else np.repeat("", trab.shape[0])
        trb_cdr3 = decode_string_array(obs["TRB_cdr3"][:]) if "TRB_cdr3" in obs else np.repeat("", trab.shape[0])
    paired_mask = (pd.Series(tra_cdr3).str.strip() != "") & (pd.Series(trb_cdr3).str.strip() != "")
    return pd.DataFrame(
        {
            "phase4_trab_score": trab,
            "phase4_trd_score": trd,
            "phase4_trd_minus_trab": trd_minus_trab,
            "tcr_pairing_group": np.where(paired_mask.to_numpy(), "Paired TRA/TRB", "Not paired"),
        }
    )


def write_raw_trab_trd_scatter(frame: pd.DataFrame) -> None:
    """Write the raw TRAB-versus-TRD scatter colored by raw TRD-minus-TRAB."""
    fig, ax = plt.subplots(figsize=(6.8, 5.8))
    points = ax.scatter(
        frame["phase4_trab_score"],
        frame["phase4_trd_score"],
        c=frame["phase4_trd_minus_trab"],
        cmap="coolwarm",
        s=POINT_SIZE,
        linewidths=0,
        alpha=ALPHA,
        rasterized=True,
    )
    ax.set_title("Raw TRAB-versus-TRD Score Space", fontsize=15)
    ax.set_xlabel("Raw TRAB score", fontsize=11)
    ax.set_ylabel("Raw TRD score", fontsize=11)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.colorbar(points, ax=ax, fraction=0.046, pad=0.03, label="Raw TRD - TRAB")
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "phase4_trab_vs_trd_raw_only.png", dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)


def write_raw_paired_tcr_scatter(frame: pd.DataFrame) -> None:
    """Write the raw TRAB-versus-TRD panel for paired TRA/TRB versus all others."""
    palette = {"Not paired": "#0077B6", "Paired TRA/TRB": "#D1495B"}
    fig, ax = plt.subplots(figsize=(6.8, 5.8))
    for group_name in ["Not paired", "Paired TRA/TRB"]:
        group_df = frame.loc[frame["tcr_pairing_group"] == group_name]
        ax.scatter(
            group_df["phase4_trab_score"],
            group_df["phase4_trd_score"],
            s=POINT_SIZE,
            linewidths=0,
            alpha=ALPHA,
            rasterized=True,
            color=palette[group_name],
            label=group_name,
        )
    ax.set_title("Raw TRAB-versus-TRD Space with Paired TRA/TRB Context", fontsize=15)
    ax.set_xlabel("Raw TRAB score", fontsize=11)
    ax.set_ylabel("Raw TRD score", fontsize=11)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(title="Pairing context", frameon=False, fontsize=9, title_fontsize=10)
    fig.tight_layout()
    fig.savefig(
        FIGURE_DIR / "phase4_trab_vs_trd_paired_tratrb_vs_no_tcr_raw_only.png",
        dpi=FIGURE_DPI,
        bbox_inches="tight",
    )
    plt.close(fig)


def main() -> None:
    """Refresh the two raw scatter panels used in the report."""
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    frame = read_full_obs_frame()
    write_raw_trab_trd_scatter(frame)
    write_raw_paired_tcr_scatter(frame)


if __name__ == "__main__":
    main()
