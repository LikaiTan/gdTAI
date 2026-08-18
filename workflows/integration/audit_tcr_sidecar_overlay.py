#!/usr/bin/env python3
"""Audit whether the validated TCR sidecar can overlay the rebuilt full atlas.

This workflow is read-only. It validates immutable identities and join coverage
before any metadata replacement is considered.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_ATLAS = PROJECT_ROOT / "Integrated_dataset/integrated_full_atlas.h5ad"
DEFAULT_SIDECAR = (
    PROJECT_ROOT
    / "Integrated_dataset/tables/tcr_join_rebuild/validated_tcr_replacement_sidecar.parquet"
)
DEFAULT_EXCLUSIONS = PROJECT_ROOT / "Integrated_dataset/tables/phase3_input_sanitization.csv"
DEFAULT_TABLE_DIR = PROJECT_ROOT / "Integrated_dataset/tables/tcr_sidecar_overlay_preflight"
DEFAULT_LOG_DIR = PROJECT_ROOT / "Integrated_dataset/logs/tcr_sidecar_overlay_preflight"
EXPECTED_SIDECAR_SHA256 = "3114e70719301d693ae1a2bc2c63bac6c8bd57e3e8ac73a88c24320eaabfc2f0"
EXPECTED_SOURCES = 14
EXPECTED_ROWS = 3_041_871


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atlas", type=Path, default=DEFAULT_ATLAS)
    parser.add_argument("--sidecar", type=Path, default=DEFAULT_SIDECAR)
    parser.add_argument("--exclusions", type=Path, default=DEFAULT_EXCLUSIONS)
    parser.add_argument("--table-dir", type=Path, default=DEFAULT_TABLE_DIR)
    parser.add_argument("--log-dir", type=Path, default=DEFAULT_LOG_DIR)
    return parser.parse_args()


def sha256_file(path: Path, chunk_size: int = 16 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(chunk_size):
            digest.update(block)
    return digest.hexdigest()


def decode_strings(values: np.ndarray) -> np.ndarray:
    if values.dtype.kind in {"O", "S"}:
        return np.asarray(
            [value.decode("utf-8") if isinstance(value, bytes) else str(value) for value in values],
            dtype=object,
        )
    return values.astype(str).astype(object)


class CategoricalColumn:
    def __init__(self, obs: h5py.Group, name: str) -> None:
        group = obs[name]
        if not isinstance(group, h5py.Group) or not {"categories", "codes"}.issubset(group.keys()):
            raise ValueError(f"obs[{name!r}] is not an AnnData categorical column")
        self.categories = decode_strings(group["categories"][:])
        self.codes = group["codes"][:]

    def values(self, positions: np.ndarray) -> np.ndarray:
        codes = self.codes[positions]
        if np.any(codes < 0):
            raise ValueError("Unexpected missing category code in overlay key column")
        return self.categories[codes]


def write_json(path: Path, payload: dict[str, object]) -> None:
    temporary = path.with_suffix(path.suffix + ".partial")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def main() -> None:
    args = parse_args()
    atlas = args.atlas.resolve()
    sidecar_path = args.sidecar.resolve()
    exclusions_path = args.exclusions.resolve()
    args.table_dir.mkdir(parents=True, exist_ok=True)
    args.log_dir.mkdir(parents=True, exist_ok=True)

    atlas_stat_before = (atlas.stat().st_size, atlas.stat().st_mtime_ns)
    sidecar_sha256 = sha256_file(sidecar_path)
    columns = [
        "source_gse_id",
        "source_obs_name",
        "join_sample_id",
        "barcode_core",
        "join_status",
        "join_reason",
        "tcr_assignment_eligible",
        "replacement_eligible",
        "has_any_ab_tcr_rebuilt",
        "has_any_gd_tcr_rebuilt",
        "has_TRA_TRB_paired_rebuilt",
        "has_TRG_TRD_paired_rebuilt",
    ]
    sidecar = pd.read_parquet(sidecar_path, columns=columns)
    for column in ("source_gse_id", "source_obs_name", "join_sample_id", "barcode_core"):
        sidecar[column] = sidecar[column].fillna("").astype(str)

    exclusions = pd.read_csv(exclusions_path, usecols=["obs_name", "source_gse_id"], dtype=str)
    excluded_by_source = {}
    for source, group in exclusions.groupby("source_gse_id", sort=False):
        prefix = f"{source}__"
        excluded_by_source[source] = {
            value[len(prefix) :] if value.startswith(prefix) else value
            for value in group["obs_name"].astype(str)
        }

    per_source: list[dict[str, object]] = []
    missing_rows: list[pd.DataFrame] = []
    with h5py.File(atlas, "r") as handle:
        obs = handle["obs"]
        source_column = CategoricalColumn(obs, "source_gse_id")
        original_cell_column = CategoricalColumn(obs, "original_cell_id")
        sample_column = CategoricalColumn(obs, "sample_id")
        barcode_column = CategoricalColumn(obs, "barcode_core")
        atlas_shape = tuple(int(value) for value in handle["X"].attrs["shape"])
        matrix_nnz = int(handle["X"]["data"].shape[0])
        atlas_old_any_ab_total = int(
            np.asarray(obs["has_any_ab_tcr"][:], dtype=bool).sum()
        )
        atlas_old_paired_ab_total = int(
            np.asarray(obs["has_TRA_TRB_paired"][:], dtype=bool).sum()
        )

        for source in sorted(sidecar["source_gse_id"].unique()):
            source_codes = np.flatnonzero(source_column.categories == source)
            if len(source_codes) != 1:
                raise RuntimeError(f"Atlas source category is not unique for {source}: {source_codes}")
            atlas_positions = np.flatnonzero(source_column.codes == source_codes[0])
            atlas_names = original_cell_column.values(atlas_positions).astype(str)
            atlas_samples = sample_column.values(atlas_positions).astype(str)
            atlas_barcodes = barcode_column.values(atlas_positions).astype(str)
            atlas_frame = pd.DataFrame(
                {
                    "source_obs_name": atlas_names,
                    "atlas_sample_id": atlas_samples,
                    "atlas_barcode_core": atlas_barcodes,
                    "atlas_old_any_ab": np.asarray(obs["has_any_ab_tcr"][atlas_positions], dtype=bool),
                    "atlas_old_any_gd": np.asarray(obs["has_any_gd_tcr"][atlas_positions], dtype=bool),
                    "atlas_old_paired_ab": np.asarray(
                        obs["has_TRA_TRB_paired"][atlas_positions], dtype=bool
                    ),
                    "atlas_old_paired_gd": np.asarray(
                        obs["has_TRG_TRD_paired"][atlas_positions], dtype=bool
                    ),
                }
            )
            atlas_unique = bool(atlas_frame["source_obs_name"].is_unique)

            source_sidecar = sidecar.loc[sidecar["source_gse_id"].eq(source)].copy()
            sidecar_unique = bool(source_sidecar["source_obs_name"].is_unique)
            merged = source_sidecar.merge(
                atlas_frame,
                on="source_obs_name",
                how="left",
                validate="one_to_one" if atlas_unique and sidecar_unique else None,
                indicator=True,
            )
            matched = merged["_merge"].eq("both")
            missing = merged.loc[~matched, ["source_gse_id", "source_obs_name", "join_status", "join_reason"]].copy()
            if not missing.empty:
                allowed_exclusions = excluded_by_source.get(source, set())
                missing["listed_phase3_exclusion"] = missing["source_obs_name"].isin(allowed_exclusions)
                missing_rows.append(missing)

            matched_frame = merged.loc[matched]
            sidecar_only_frame = merged.loc[~matched]
            sample_comparable = matched_frame["atlas_sample_id"].ne("")
            barcode_comparable = matched_frame["atlas_barcode_core"].ne("")
            sample_match = matched_frame["join_sample_id"].eq(matched_frame["atlas_sample_id"])
            barcode_match = matched_frame["barcode_core"].eq(matched_frame["atlas_barcode_core"])
            old_any = matched_frame["atlas_old_any_ab"] | matched_frame["atlas_old_any_gd"]
            rebuilt_any = (
                matched_frame["has_any_ab_tcr_rebuilt"]
                | matched_frame["has_any_gd_tcr_rebuilt"]
            )
            atlas_without_sidecar = int((~atlas_frame["source_obs_name"].isin(source_sidecar["source_obs_name"])).sum())
            per_source.append(
                {
                    "source_gse_id": source,
                    "n_atlas_rows": int(len(atlas_frame)),
                    "n_sidecar_rows": int(len(source_sidecar)),
                    "n_exact_source_obs_matches": int(matched.sum()),
                    "n_sidecar_rows_missing_from_atlas": int((~matched).sum()),
                    "n_missing_listed_phase3_exclusion": int(
                        missing.get("listed_phase3_exclusion", pd.Series(dtype=bool)).sum()
                    ),
                    "n_atlas_rows_without_sidecar": atlas_without_sidecar,
                    "n_sample_id_comparable": int(sample_comparable.sum()),
                    "n_sample_id_unavailable_in_atlas": int((~sample_comparable).sum()),
                    "n_sample_id_mismatches_where_available": int((sample_comparable & ~sample_match).sum()),
                    "n_barcode_core_comparable": int(barcode_comparable.sum()),
                    "n_barcode_core_unavailable_in_atlas": int((~barcode_comparable).sum()),
                    "n_barcode_core_mismatches_where_available": int((barcode_comparable & ~barcode_match).sum()),
                    "n_tcr_assignment_ineligible": int((~source_sidecar["tcr_assignment_eligible"]).sum()),
                    "n_matched_tcr_assignment_ineligible": int(
                        (~matched_frame["tcr_assignment_eligible"]).sum()
                    ),
                    "n_replacement_eligible": int(source_sidecar["replacement_eligible"].sum()),
                    "n_rebuilt_any_ab": int(source_sidecar["has_any_ab_tcr_rebuilt"].sum()),
                    "n_rebuilt_any_gd": int(source_sidecar["has_any_gd_tcr_rebuilt"].sum()),
                    "n_matched_old_any_tcr": int(old_any.sum()),
                    "n_matched_old_any_ab": int(
                        matched_frame["atlas_old_any_ab"].sum()
                    ),
                    "n_matched_old_paired_ab": int(
                        matched_frame["atlas_old_paired_ab"].sum()
                    ),
                    "n_matched_rebuilt_any_tcr": int(rebuilt_any.sum()),
                    "n_matched_rebuilt_any_ab": int(
                        matched_frame["has_any_ab_tcr_rebuilt"].sum()
                    ),
                    "n_matched_rebuilt_paired_ab": int(
                        matched_frame["has_TRA_TRB_paired_rebuilt"].sum()
                    ),
                    "n_matched_rebuilt_any_gd": int(
                        matched_frame["has_any_gd_tcr_rebuilt"].sum()
                    ),
                    "n_matched_rebuilt_paired_gd": int(
                        matched_frame["has_TRG_TRD_paired_rebuilt"].sum()
                    ),
                    "n_matched_old_only_any_tcr": int((old_any & ~rebuilt_any).sum()),
                    "n_matched_rebuilt_only_any_tcr": int((~old_any & rebuilt_any).sum()),
                    "n_matched_any_tcr_agreement": int((old_any == rebuilt_any).sum()),
                    "n_sidecar_only_rebuilt_any_tcr": int(
                        (
                            sidecar_only_frame["has_any_ab_tcr_rebuilt"]
                            | sidecar_only_frame["has_any_gd_tcr_rebuilt"]
                        ).sum()
                    ),
                    "n_sidecar_only_rebuilt_any_ab": int(
                        sidecar_only_frame["has_any_ab_tcr_rebuilt"].sum()
                    ),
                    "n_sidecar_only_rebuilt_any_gd": int(
                        sidecar_only_frame["has_any_gd_tcr_rebuilt"].sum()
                    ),
                    "n_sidecar_only_rebuilt_paired_ab": int(
                        sidecar_only_frame["has_TRA_TRB_paired_rebuilt"].sum()
                    ),
                    "n_sidecar_only_rebuilt_paired_gd": int(
                        sidecar_only_frame["has_TRG_TRD_paired_rebuilt"].sum()
                    ),
                    "atlas_source_obs_unique": atlas_unique,
                    "sidecar_source_obs_unique": sidecar_unique,
                }
            )

    source_table = pd.DataFrame(per_source)
    missing_table = (
        pd.concat(missing_rows, ignore_index=True)
        if missing_rows
        else pd.DataFrame(
            columns=[
                "source_gse_id",
                "source_obs_name",
                "join_status",
                "join_reason",
                "listed_phase3_exclusion",
            ]
        )
    )
    source_table.to_csv(args.table_dir / "overlay_join_by_source.csv", index=False)
    missing_table.to_csv(args.table_dir / "sidecar_rows_missing_from_atlas.csv", index=False)

    checks = {
        "sidecar_sha256": sidecar_sha256 == EXPECTED_SIDECAR_SHA256,
        "sidecar_row_count": len(sidecar) == EXPECTED_ROWS,
        "sidecar_source_count": sidecar["source_gse_id"].nunique() == EXPECTED_SOURCES,
        "sidecar_keys_unique": not sidecar.duplicated(["source_gse_id", "source_obs_name"]).any(),
        "atlas_keys_unique_within_source": bool(source_table["atlas_source_obs_unique"].all()),
        "no_atlas_rows_without_sidecar_in_affected_sources": bool(
            source_table["n_atlas_rows_without_sidecar"].eq(0).all()
        ),
        "ambiguous_rows_fail_closed": int((~sidecar["tcr_assignment_eligible"]).sum()) == 110,
        "atlas_file_unchanged": atlas_stat_before == (atlas.stat().st_size, atlas.stat().st_mtime_ns),
    }
    status = "PASS_OVERLAY_READY" if all(checks.values()) else "FAIL_OVERLAY_NOT_READY"
    summary = {
        "status": status,
        "atlas": str(atlas),
        "atlas_shape": list(atlas_shape),
        "atlas_matrix_nnz": matrix_nnz,
        "sidecar": str(sidecar_path),
        "sidecar_sha256": sidecar_sha256,
        "n_sidecar_rows": int(len(sidecar)),
        "n_sidecar_sources": int(sidecar["source_gse_id"].nunique()),
        "n_exact_matches": int(source_table["n_exact_source_obs_matches"].sum()),
        "n_missing_from_atlas": int(source_table["n_sidecar_rows_missing_from_atlas"].sum()),
        "n_missing_listed_phase3_exclusion": int(source_table["n_missing_listed_phase3_exclusion"].sum()),
        "n_atlas_rows_without_sidecar": int(source_table["n_atlas_rows_without_sidecar"].sum()),
        "n_tcr_assignment_ineligible": int((~sidecar["tcr_assignment_eligible"]).sum()),
        "n_matched_tcr_assignment_ineligible": int(
            source_table["n_matched_tcr_assignment_ineligible"].sum()
        ),
        "n_sample_id_mismatches_where_available": int(
            source_table["n_sample_id_mismatches_where_available"].sum()
        ),
        "n_barcode_core_mismatches_where_available": int(
            source_table["n_barcode_core_mismatches_where_available"].sum()
        ),
        "n_sidecar_only_rebuilt_any_tcr": int(
            source_table["n_sidecar_only_rebuilt_any_tcr"].sum()
        ),
        "n_sidecar_only_rebuilt_any_ab": int(
            source_table["n_sidecar_only_rebuilt_any_ab"].sum()
        ),
        "n_sidecar_only_rebuilt_any_gd": int(
            source_table["n_sidecar_only_rebuilt_any_gd"].sum()
        ),
        "n_sidecar_only_rebuilt_paired_ab": int(
            source_table["n_sidecar_only_rebuilt_paired_ab"].sum()
        ),
        "n_sidecar_only_rebuilt_paired_gd": int(
            source_table["n_sidecar_only_rebuilt_paired_gd"].sum()
        ),
        "n_matched_old_any_tcr": int(source_table["n_matched_old_any_tcr"].sum()),
        "n_matched_old_any_ab": int(source_table["n_matched_old_any_ab"].sum()),
        "n_matched_old_paired_ab": int(
            source_table["n_matched_old_paired_ab"].sum()
        ),
        "n_matched_rebuilt_any_tcr": int(
            source_table["n_matched_rebuilt_any_tcr"].sum()
        ),
        "n_matched_rebuilt_any_ab": int(
            source_table["n_matched_rebuilt_any_ab"].sum()
        ),
        "n_matched_rebuilt_paired_ab": int(
            source_table["n_matched_rebuilt_paired_ab"].sum()
        ),
        "n_matched_rebuilt_any_gd": int(
            source_table["n_matched_rebuilt_any_gd"].sum()
        ),
        "n_matched_rebuilt_paired_gd": int(
            source_table["n_matched_rebuilt_paired_gd"].sum()
        ),
        "n_matched_old_only_any_tcr": int(
            source_table["n_matched_old_only_any_tcr"].sum()
        ),
        "n_matched_rebuilt_only_any_tcr": int(
            source_table["n_matched_rebuilt_only_any_tcr"].sum()
        ),
        "checks": checks,
        "h5ad_write_performed": False,
    }
    summary["n_atlas_overlay_any_ab"] = (
        atlas_old_any_ab_total
        - summary["n_matched_old_any_ab"]
        + summary["n_matched_rebuilt_any_ab"]
    )
    summary["n_atlas_overlay_paired_ab"] = (
        atlas_old_paired_ab_total
        - summary["n_matched_old_paired_ab"]
        + summary["n_matched_rebuilt_paired_ab"]
    )
    summary["sidecar_only_productive_ab_pct_of_available"] = 100.0 * (
        summary["n_sidecar_only_rebuilt_any_ab"]
        / (
            summary["n_matched_rebuilt_any_ab"]
            + summary["n_sidecar_only_rebuilt_any_ab"]
        )
    )
    summary["sidecar_only_paired_ab_pct_of_available"] = 100.0 * (
        summary["n_sidecar_only_rebuilt_paired_ab"]
        / (
            summary["n_matched_rebuilt_paired_ab"]
            + summary["n_sidecar_only_rebuilt_paired_ab"]
        )
    )
    summary["sidecar_only_productive_ab_pct_of_atlas"] = 100.0 * (
        summary["n_sidecar_only_rebuilt_any_ab"] / atlas_shape[0]
    )
    summary["sidecar_only_productive_ab_pct_of_whole_available"] = 100.0 * (
        summary["n_sidecar_only_rebuilt_any_ab"]
        / (
            summary["n_atlas_overlay_any_ab"]
            + summary["n_sidecar_only_rebuilt_any_ab"]
        )
    )
    summary["sidecar_only_paired_ab_pct_of_whole_available"] = 100.0 * (
        summary["n_sidecar_only_rebuilt_paired_ab"]
        / (
            summary["n_atlas_overlay_paired_ab"]
            + summary["n_sidecar_only_rebuilt_paired_ab"]
        )
    )
    write_json(args.log_dir / "overlay_preflight_summary.json", summary)

    check_lines = [
        f"| {name} | {'PASS' if passed else 'FAIL'} |"
        for name, passed in checks.items()
    ]
    markdown = "\n".join(
        [
            "# TCR sidecar overlay preflight",
            "",
            f"Status: **{status}**",
            "",
            "This audit is read-only. No H5AD was modified.",
            "",
            "## Coverage",
            "",
            f"- Sidecar rows: {len(sidecar):,}",
            f"- Exact atlas matches: {summary['n_exact_matches']:,}",
            f"- Sidecar rows outside the T/NK-filtered atlas: {summary['n_missing_from_atlas']:,}",
            f"- Of those, rows explained by frozen Phase 3 exclusions: {summary['n_missing_listed_phase3_exclusion']:,}",
            f"- Atlas rows without a sidecar in affected sources: {summary['n_atlas_rows_without_sidecar']:,}",
            f"- Sidecar-only cells with rebuilt productive TCR: {summary['n_sidecar_only_rebuilt_any_tcr']:,}",
            f"- Sidecar-only cells with rebuilt paired TRA/TRB: {summary['n_sidecar_only_rebuilt_paired_ab']:,}",
            f"- Sidecar-only cells with rebuilt paired TRG/TRD: {summary['n_sidecar_only_rebuilt_paired_gd']:,}",
            f"- Retained atlas cells with rebuilt productive alpha-beta TCR: {summary['n_matched_rebuilt_any_ab']:,}",
            f"- Retained atlas cells with rebuilt paired TRA/TRB: {summary['n_matched_rebuilt_paired_ab']:,}",
            f"- Estimated whole-atlas productive-alpha-beta count after overlay: {summary['n_atlas_overlay_any_ab']:,}",
            f"- Estimated whole-atlas paired-TRA/TRB count after overlay: {summary['n_atlas_overlay_paired_ab']:,}",
            f"- Omitted productive-alpha-beta fraction within the 14 repaired sources: {summary['sidecar_only_productive_ab_pct_of_available']:.3f}%",
            f"- Omitted paired-alpha-beta fraction within the 14 repaired sources: {summary['sidecar_only_paired_ab_pct_of_available']:.3f}%",
            f"- Omitted productive-alpha-beta fraction of the whole available atlas-plus-omitted pool: {summary['sidecar_only_productive_ab_pct_of_whole_available']:.3f}%",
            f"- Omitted paired-alpha-beta fraction of the whole available atlas-plus-omitted pool: {summary['sidecar_only_paired_ab_pct_of_whole_available']:.3f}%",
            f"- Omitted productive-alpha-beta cells as a fraction of the atlas: {summary['sidecar_only_productive_ab_pct_of_atlas']:.3f}%",
            f"- Current atlas cells losing all old TCR evidence after repair: {summary['n_matched_old_only_any_tcr']:,}",
            f"- Current atlas cells gaining rebuilt TCR evidence: {summary['n_matched_rebuilt_only_any_tcr']:,}",
            f"- Fail-closed TCR-ineligible rows: {summary['n_tcr_assignment_ineligible']:,}",
            f"- Fail-closed TCR-ineligible rows present in the atlas: {summary['n_matched_tcr_assignment_ineligible']:,}",
            f"- Atlas display-sample differences from normalized join keys: {summary['n_sample_id_mismatches_where_available']:,} (diagnostic only; not an overlay key)",
            "",
            "## Checks",
            "",
            "| Check | Result |",
            "|---|---:|",
            *check_lines,
            "",
            "## Proposed write boundary",
            "",
            "Keep the exact 5,933,312-cell atlas membership. Join the immutable sidecar by source_gse_id plus atlas original_cell_id, apply it to a temporary copy of the rebuilt atlas, preserve legacy TCR fields, recompute canonical flags, validate expression and embedding identities, and atomically switch the canonical symlink only after all post-write checks pass. The 20,875 productive-alpha-beta sidecar-only rows, including 14,495 paired cells and no productive gamma-delta cells, are intentionally not rescued in this repair.",
            "",
        ]
    )
    (args.log_dir / "overlay_preflight_summary.md").write_text(markdown, encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))
    if status != "PASS_OVERLAY_READY":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
