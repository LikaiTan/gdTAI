#!/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python
"""Audit productive chain evidence in staged extension-cohort TCR sources."""

from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
import re
from typing import Any

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_COHORTS = PROJECT_ROOT / "configs/datasets/extension_cohorts.csv"
DEFAULT_STAGE_ROOT = PROJECT_ROOT / "data/interim/extension_intake/staged"
DEFAULT_OUTPUT = (
    PROJECT_ROOT
    / "Integrated_dataset/tables/extension_intake/extension_tcr_source_chain_audit.csv"
)
DEFAULT_LOG = (
    PROJECT_ROOT
    / "Integrated_dataset/logs/extension_intake/extension_tcr_source_chain_audit.md"
)
CHAINS = ("TRA", "TRB", "TRG", "TRD")
TRUE_VALUES = {"true", "1", "yes", "y", "t"}
MISSING_VALUES = {"", "na", "nan", "none", "null", "<na>"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohorts", type=Path, default=DEFAULT_COHORTS)
    parser.add_argument("--stage-root", type=Path, default=DEFAULT_STAGE_ROOT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--log", type=Path, default=DEFAULT_LOG)
    return parser.parse_args()


def normalized_column(columns: list[str], aliases: tuple[str, ...]) -> str | None:
    lookup = {
        re.sub(r"[^a-z0-9]+", "", str(column).casefold()): column
        for column in columns
    }
    for alias in aliases:
        key = re.sub(r"[^a-z0-9]+", "", alias.casefold())
        if key in lookup:
            return lookup[key]
    return None


def bool_mask(
    series: pd.Series | None, index: pd.Index, *, default: bool
) -> pd.Series:
    if series is None:
        return pd.Series(default, index=index, dtype=bool)
    return series.astype("string").fillna("").str.strip().str.casefold().isin(TRUE_VALUES)


def summarize_contig_frame(frame: pd.DataFrame) -> dict[str, Counter[str]]:
    chain_column = normalized_column(list(frame.columns), ("chain", "locus"))
    if chain_column is None:
        raise ValueError("Contig table lacks a chain/locus column")
    productive_column = normalized_column(
        list(frame.columns), ("productive", "is_productive")
    )
    cdr3_column = normalized_column(list(frame.columns), ("cdr3", "junction_aa"))
    cell_column = normalized_column(list(frame.columns), ("is_cell",))
    confidence_column = normalized_column(list(frame.columns), ("high_confidence",))

    chain = frame[chain_column].astype("string").fillna("").str.strip().str.upper()
    productive = bool_mask(
        frame[productive_column] if productive_column else None,
        frame.index,
        default=False,
    )
    is_cell = bool_mask(
        frame[cell_column] if cell_column else None,
        frame.index,
        default=True,
    )
    high_confidence = bool_mask(
        frame[confidence_column] if confidence_column else None,
        frame.index,
        default=True,
    )
    if cdr3_column is None:
        cdr3 = pd.Series(False, index=frame.index, dtype=bool)
    else:
        cdr3 = ~(
            frame[cdr3_column]
            .astype("string")
            .fillna("")
            .str.strip()
            .str.casefold()
            .isin(MISSING_VALUES)
        )
    result: dict[str, Counter[str]] = {
        "total": Counter(),
        "cell_high_confidence": Counter(),
        "productive_flag": Counter(),
        "eligible_productive_cdr3": Counter(),
    }
    for chain_name in CHAINS:
        selected = chain.eq(chain_name)
        result["total"][chain_name] += int(selected.sum())
        result["cell_high_confidence"][chain_name] += int(
            (selected & is_cell & high_confidence).sum()
        )
        result["productive_flag"][chain_name] += int((selected & productive).sum())
        result["eligible_productive_cdr3"][chain_name] += int(
            (selected & productive & is_cell & high_confidence & cdr3).sum()
        )
    return result


def staged_contig_paths(stage_dir: Path) -> list[Path]:
    manifest_path = stage_dir / "stage_manifest.json"
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    paths = []
    for source in payload.get("sources", []):
        path = stage_dir / str(source.get("staged_path", ""))
        if path.is_file() and "contig" in path.name.casefold() and ".csv" in path.name.casefold():
            paths.append(path)
    return sorted(set(paths), key=str)


def audit_cohort(cohort: dict[str, str], stage_root: Path) -> list[dict[str, Any]]:
    cohort_id = cohort["cohort_id"]
    paths = staged_contig_paths(stage_root / cohort_id)
    totals = {
        metric: Counter({chain: 0 for chain in CHAINS})
        for metric in (
            "total",
            "cell_high_confidence",
            "productive_flag",
            "eligible_productive_cdr3",
        )
    }
    for path in paths:
        frame = pd.read_csv(path, low_memory=False)
        summary = summarize_contig_frame(frame)
        for metric, counts in summary.items():
            totals[metric].update(counts)
    return [
        {
            "cohort_id": cohort_id,
            "tcr_schema": cohort["tcr_schema"],
            "source_contig_files": len(paths),
            "chain": chain,
            "source_contig_rows": totals["total"][chain],
            "cell_high_confidence_rows": totals["cell_high_confidence"][chain],
            "productive_flag_true_rows": totals["productive_flag"][chain],
            "eligible_productive_cdr3_rows": totals["eligible_productive_cdr3"][chain],
        }
        for chain in CHAINS
    ]


def markdown_table(frame: pd.DataFrame) -> list[str]:
    columns = list(frame.columns)
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for row in frame.itertuples(index=False, name=None):
        lines.append("| " + " | ".join(str(value) for value in row) + " |")
    return lines


def main() -> None:
    args = parse_args()
    cohorts = pd.read_csv(args.cohorts, dtype=str, keep_default_na=False)
    cohorts = cohorts[
        cohorts["build_enabled"].str.casefold().isin({"true", "1", "yes"})
    ]
    rows: list[dict[str, Any]] = []
    for cohort in cohorts.to_dict(orient="records"):
        rows.extend(audit_cohort(cohort, args.stage_root))
    result = pd.DataFrame(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, index=False)

    gd = result[result["chain"].isin(["TRG", "TRD"])].copy()
    view = gd[
        [
            "cohort_id",
            "chain",
            "source_contig_rows",
            "productive_flag_true_rows",
            "eligible_productive_cdr3_rows",
        ]
    ]
    source_rows = int(view["source_contig_rows"].sum())
    productive_rows = int(view["eligible_productive_cdr3_rows"].sum())
    lines = [
        "# Extension TCR Source-Chain Audit",
        "",
        "- Mode: read-only staged-source audit",
        f"- TRG/TRD source contig rows: `{source_rows:,}`",
        f"- Eligible productive TRG/TRD CDR3 rows: `{productive_rows:,}`",
        "",
        *markdown_table(view),
        "",
        "The builder requires productive, cell-associated, high-confidence contigs with a nonempty CDR3. "
        "TRG/TRD fragments that lack productive CDR3 evidence are not ground truth.",
        "",
    ]
    args.log.parent.mkdir(parents=True, exist_ok=True)
    args.log.write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
