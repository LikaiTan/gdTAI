#!/usr/bin/env python3
"""Build a row-keyed additive metadata overlay without modifying an H5AD."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path
from typing import Any

import h5py
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = ROOT / "configs/metadata/full_atlas_metadata_harmonization_v2.json"
DEFAULT_ATLAS = ROOT / "Integrated_dataset/integrated_full_atlas.h5ad"
DEFAULT_TCR_SIDECAR = (
    ROOT
    / "Integrated_dataset/tables/tcr_join_rebuild/validated_tcr_replacement_sidecar.parquet"
)
DEFAULT_OUTPUT = ROOT / "Integrated_dataset/tables/metadata_harmonization/full_atlas_v2"
EXPECTED_ATLAS_CELLS = 5_933_312
EXPECTED_TCR_IDENTITY_ROWS = 2_155_409
EXPECTED_TCR_SIDECAR_SHA256 = "3114e70719301d693ae1a2bc2c63bac6c8bd57e3e8ac73a88c24320eaabfc2f0"
FORBIDDEN_CONTEXT_VALUES = {"cd4", "cd8", "nk", "treg", "naive", "proliferating"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--atlas", type=Path, default=DEFAULT_ATLAS)
    parser.add_argument("--tcr-sidecar", type=Path, default=DEFAULT_TCR_SIDECAR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def sha256_file(path: Path, chunk_size: int = 16 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(chunk_size):
            digest.update(block)
    return digest.hexdigest()


def file_state(path: Path) -> dict[str, int | str]:
    stat = path.stat()
    return {
        "path": str(path.resolve()),
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "inode": stat.st_ino,
    }


def decode(values: np.ndarray) -> np.ndarray:
    return np.asarray(
        [value.decode() if isinstance(value, bytes) else str(value) for value in values],
        dtype=object,
    )


def values_from_obs(obs: h5py.Group, name: str) -> np.ndarray:
    node = obs[name]
    if isinstance(node, h5py.Group) and {"categories", "codes"}.issubset(node):
        categories = decode(node["categories"][:])
        codes = node["codes"][:]
        result = np.full(len(codes), "", dtype=object)
        present = codes >= 0
        result[present] = categories[codes[present]]
        return result
    if isinstance(node, h5py.Dataset):
        return decode(node[:])
    raise ValueError(f"Unsupported obs encoding for {name!r}")


class CategoricalColumn:
    def __init__(self, obs: h5py.Group, name: str) -> None:
        group = obs[name]
        if not isinstance(group, h5py.Group) or not {"categories", "codes"}.issubset(group):
            raise ValueError(f"obs[{name!r}] is not categorical")
        self.categories = decode(group["categories"][:])
        self.codes = group["codes"][:]

    def values(self, positions: np.ndarray) -> np.ndarray:
        codes = self.codes[positions]
        result = np.full(len(codes), "", dtype=object)
        present = codes >= 0
        result[present] = self.categories[codes[present]]
        return result


def normalize(value: Any) -> str:
    return " ".join(str(value).strip().casefold().replace("_", " ").split())


def normalize_token(value: str) -> str:
    return "_".join(value.strip().casefold().replace("/", "_").split())


def nonempty(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip().ne("")


def source_obs_index(obs: h5py.Group) -> np.ndarray:
    index_name = obs.attrs.get("_index", "_index")
    if isinstance(index_name, bytes):
        index_name = index_name.decode()
    return values_from_obs(obs, str(index_name))


def load_source_overlays(config: dict[str, Any]) -> tuple[dict[str, pd.DataFrame], dict[str, dict[str, int | str]]]:
    overlays: dict[str, pd.DataFrame] = {}
    states: dict[str, dict[str, int | str]] = {}
    for rule in config["source_metadata_overlays"]:
        source = rule["source"]
        path = (ROOT / rule["path"]).resolve()
        states[source] = file_state(path)
        with h5py.File(path, "r") as handle:
            obs = handle["obs"]
            index = source_obs_index(obs)
            tissue_raw = values_from_obs(obs, rule["tissue_field"])
            if rule.get("sample_rule") == "donor__tissue":
                donor = values_from_obs(obs, rule["donor_field"])
                sample = np.asarray(
                    [f"{normalize_token(d)}__{normalize_token(t)}" for d, t in zip(donor, tissue_raw)],
                    dtype=object,
                )
                library = sample.copy()
            else:
                sample = values_from_obs(obs, rule["sample_field"])
                library = values_from_obs(obs, rule["library_field"])
                donor = values_from_obs(obs, rule["donor_field"])
            tissue = np.asarray(
                [rule["tissue_map"].get(value, "unresolved") for value in tissue_raw],
                dtype=object,
            )
            context = np.asarray(
                [rule["context_map"].get(value, "unresolved") for value in tissue_raw],
                dtype=object,
            )
            frame = pd.DataFrame(
                {
                    "original_cell_id": index,
                    "overlay_sample_id": sample,
                    "overlay_library_id": library,
                    "overlay_donor_id": donor,
                    "overlay_tissue": tissue,
                    "overlay_context": context,
                    "overlay_rule_id": f"source_overlay:{source}",
                }
            )
        if states[source] != file_state(path):
            raise RuntimeError(f"Source H5AD changed during overlay build: {source}")
        if not frame["original_cell_id"].is_unique:
            raise RuntimeError(f"Source H5AD index is not unique: {source}")
        overlays[source] = frame
    return overlays, states


def load_exact_tumor_maps(config: dict[str, Any]) -> dict[str, dict[str, str]]:
    result: dict[str, dict[str, str]] = {}
    for imported in config.get("imports", []):
        source_rules = json.loads((ROOT / imported["path"]).read_text(encoding="utf-8"))
        target_id = imported["selector"]["sample_value_map_id"]
        source_map = next(row for row in source_rules["sample_value_maps"] if row["id"] == target_id)
        mapping: dict[str, str] = {}
        for group in source_map["value_groups"]:
            output = imported["value_map"][group["value"]]
            for key in group["keys"]:
                mapping[str(key)] = output
        result[imported["source"]] = mapping
    return result


def matches_any_field(frame: pd.DataFrame, fields: list[str], pattern: str) -> pd.Series:
    result = pd.Series(False, index=frame.index, dtype=bool)
    for field in fields:
        result |= frame[field].fillna("").astype("string").str.contains(
            pattern, case=False, regex=True, na=False
        )
    return result


def main() -> None:
    args = parse_args()
    config = json.loads(args.config.resolve().read_text(encoding="utf-8"))
    atlas = args.atlas.resolve()
    sidecar_path = args.tcr_sidecar.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "metadata_overlay_v2.parquet"
    temporary_path = output_path.with_suffix(".parquet.partial")
    if temporary_path.exists():
        temporary_path.unlink()

    atlas_before = file_state(atlas)
    source_overlays, source_states = load_source_overlays(config)
    exact_tumor_maps = load_exact_tumor_maps(config)
    legacy_rule = config["legacy_refined_metadata_overlay"]
    legacy_path = (ROOT / legacy_rule["path"]).resolve()
    legacy_before = file_state(legacy_path)
    legacy_sha = sha256_file(legacy_path)
    legacy = pd.read_csv(
        legacy_path,
        usecols=["source_gse_id", "obs_name", "tissue", "sample_source_refined"],
        dtype=str,
    ).fillna("")
    legacy["legacy_tissue_v2"] = legacy["tissue"].map(legacy_rule["tissue_map"]).fillna(
        "unresolved"
    )
    legacy["legacy_context_v2"] = legacy["tissue"].map(legacy_rule["context_map"]).fillna(
        "unresolved"
    )
    legacy_by_source = {
        source: group.drop(columns="source_gse_id").copy()
        for source, group in legacy.groupby("source_gse_id", sort=False)
    }
    defaults = {row["source"]: row for row in config["source_defaults"]}
    scoped_rules: dict[str, list[dict[str, Any]]] = {}
    for rule in config["source_scoped_rules"]:
        scoped_rules.setdefault(rule["source"], []).append(rule)
    synonyms = {
        normalize(raw): values for raw, values in config["global_rebuilt_tissue_synonyms"].items()
    }

    sidecar_sha = sha256_file(sidecar_path)
    tcr_identity = pd.read_parquet(
        sidecar_path,
        columns=["source_gse_id", "source_obs_name", "join_sample_id"],
    ).fillna("")
    for column in tcr_identity.columns:
        tcr_identity[column] = tcr_identity[column].astype(str)
    tcr_sources = set(tcr_identity["source_gse_id"].unique())

    summary_rows: list[dict[str, object]] = []
    tissue_counts: dict[str, int] = {}
    context_counts: dict[str, int] = {}
    tumor_counts: dict[str, int] = {}
    writer: pq.ParquetWriter | None = None
    total_rows = 0
    total_tcr_identity = 0
    total_sample_unresolved = 0
    overlay_source_unmatched = 0
    tcr_identity_unmatched = 0
    legacy_expected_rows = 0
    legacy_unmatched_rows = 0

    with h5py.File(atlas, "r") as handle:
        obs = handle["obs"]
        source_col = CategoricalColumn(obs, "source_gse_id")
        columns = {
            name: CategoricalColumn(obs, name)
            for name in (
                "source_obs_name",
                "original_cell_id",
                "sample_id",
                "library_id",
                "donor_id",
                "tissue_original",
                "tissue_harmonized",
            )
        }

        for source_code, source in enumerate(source_col.categories):
            positions = np.flatnonzero(source_col.codes == source_code)
            frame = pd.DataFrame(
                {
                    "source_gse_id": source,
                    **{name: column.values(positions) for name, column in columns.items()},
                }
            )
            if frame["original_cell_id"].duplicated().any():
                raise RuntimeError(f"Duplicate original_cell_id within source {source}")

            raw_tissue = frame["tissue_harmonized"].where(
                nonempty(frame["tissue_harmonized"]), frame["tissue_original"]
            )
            normalized_tissue = raw_tissue.map(normalize)
            frame["tissue_harmonized_v2"] = normalized_tissue.map(
                lambda value: synonyms.get(value, {}).get("tissue", "unresolved")
            )
            frame["specimen_context_harmonized_v2"] = normalized_tissue.map(
                lambda value: synonyms.get(value, {}).get("context", "unresolved")
            )
            frame["metadata_rule_id_v2"] = np.where(
                frame["tissue_harmonized_v2"].ne("unresolved"),
                "global_tissue_synonym",
                "fail_closed_unresolved",
            )

            default = defaults[source]
            for field, config_key in (
                ("tissue_harmonized_v2", "tissue"),
                ("specimen_context_harmonized_v2", "context"),
            ):
                if config_key in default:
                    missing = frame[field].eq("unresolved")
                    frame.loc[missing, field] = default[config_key]
                    frame.loc[missing, "metadata_rule_id_v2"] = f"source_default:{source}"
            frame["tumor_type_harmonized_v2"] = default.get("tumor_type", "unresolved")
            frame["metadata_evidence_url_v2"] = default.get("evidence_url", "")
            frame["metadata_evidence_level_v2"] = default.get("evidence_level", "none")
            frame["metadata_status_v2"] = default.get("status", "unresolved")

            for rule in scoped_rules.get(source, []):
                matcher = rule["any_field_regex"]
                matched = matches_any_field(frame, matcher["fields"], matcher["pattern"])
                for output_name, value in rule["outputs"].items():
                    target = {
                        "tissue": "tissue_harmonized_v2",
                        "context": "specimen_context_harmonized_v2",
                        "tumor_type": "tumor_type_harmonized_v2",
                    }[output_name]
                    frame.loc[matched, target] = value
                frame.loc[matched, "metadata_rule_id_v2"] = f"scoped_rule:{rule['id']}"

            if source in exact_tumor_maps:
                tumor_map = exact_tumor_maps[source]
                exact_value = frame["sample_id"].map(tumor_map).fillna("")
                exact_value = exact_value.where(nonempty(exact_value), frame["library_id"].map(tumor_map).fillna(""))
                matched = nonempty(exact_value)
                frame.loc[matched, "tumor_type_harmonized_v2"] = exact_value[matched]
                frame.loc[matched, "metadata_rule_id_v2"] = "exact_sample_tumor_map"

            if source in legacy_by_source:
                frame = frame.merge(
                    legacy_by_source[source],
                    left_on="source_obs_name",
                    right_on="obs_name",
                    how="left",
                    validate="one_to_one",
                    indicator="_legacy_merge",
                )
                matched = frame["_legacy_merge"].eq("both")
                legacy_expected_rows += len(frame)
                legacy_unmatched_rows += int((~matched).sum())
                tissue_known = matched & frame["legacy_tissue_v2"].ne("unresolved")
                context_known = matched & frame["legacy_context_v2"].ne("unresolved")
                frame.loc[tissue_known, "tissue_harmonized_v2"] = frame.loc[
                    tissue_known, "legacy_tissue_v2"
                ]
                frame.loc[context_known, "specimen_context_harmonized_v2"] = frame.loc[
                    context_known, "legacy_context_v2"
                ]
                frame.loc[tissue_known | context_known, "metadata_rule_id_v2"] = (
                    "legacy_refined_metadata_overlay"
                )

            frame["source_accession_harmonized_v2"] = config["source_accession_aliases"].get(
                source, source
            )
            frame["sample_id_harmonized_v2"] = frame["sample_id"]
            frame["library_id_harmonized_v2"] = frame["library_id"]
            frame["donor_id_harmonized_v2"] = frame["donor_id"]
            frame["sample_identity_rule_v2"] = "preserve_current_identity"
            frame["tcr_library_join_id_v2"] = ""

            if source in tcr_sources:
                source_tcr = tcr_identity.loc[tcr_identity["source_gse_id"].eq(source)].drop(
                    columns="source_gse_id"
                )
                frame = frame.merge(
                    source_tcr,
                    left_on="original_cell_id",
                    right_on="source_obs_name",
                    how="left",
                    validate="one_to_one",
                    indicator="_tcr_merge",
                    suffixes=("", "_tcr"),
                )
                matched = frame["_tcr_merge"].eq("both")
                tcr_identity_unmatched += int((~matched).sum())
                frame["tcr_library_join_id_v2"] = frame["join_sample_id"].fillna("")
                total_tcr_identity += int(nonempty(frame["tcr_library_join_id_v2"]).sum())

            if source in source_overlays:
                frame = frame.merge(
                    source_overlays[source],
                    on="original_cell_id",
                    how="left",
                    validate="one_to_one",
                    indicator="_source_overlay_merge",
                )
                matched = frame["_source_overlay_merge"].eq("both")
                overlay_source_unmatched += int((~matched).sum())
                frame.loc[matched, "sample_id_harmonized_v2"] = frame.loc[
                    matched, "overlay_sample_id"
                ]
                frame.loc[matched, "library_id_harmonized_v2"] = frame.loc[
                    matched, "overlay_library_id"
                ]
                if source == "GSE125527":
                    frame.loc[matched, "donor_id_harmonized_v2"] = frame.loc[
                        matched, "overlay_donor_id"
                    ]
                frame.loc[matched, "tissue_harmonized_v2"] = frame.loc[
                    matched, "overlay_tissue"
                ]
                frame.loc[matched, "specimen_context_harmonized_v2"] = frame.loc[
                    matched, "overlay_context"
                ]
                frame.loc[matched, "sample_identity_rule_v2"] = frame.loc[
                    matched, "overlay_rule_id"
                ]
                frame.loc[matched, "metadata_rule_id_v2"] = frame.loc[
                    matched, "overlay_rule_id"
                ]

            if source == "GSE228597":
                sample_missing = ~nonempty(frame["sample_id_harmonized_v2"])
                library_missing = ~nonempty(frame["library_id_harmonized_v2"])
                frame.loc[sample_missing, "sample_id_harmonized_v2"] = "unresolved"
                frame.loc[library_missing, "library_id_harmonized_v2"] = frame.loc[
                    library_missing, "tcr_library_join_id_v2"
                ]
                frame.loc[sample_missing, "sample_identity_rule_v2"] = (
                    "unresolved_biological_sample_known_tcr_library"
                )

            output_columns = [
                "source_gse_id",
                "original_cell_id",
                "source_accession_harmonized_v2",
                "tissue_harmonized_v2",
                "specimen_context_harmonized_v2",
                "tumor_type_harmonized_v2",
                "sample_id_harmonized_v2",
                "library_id_harmonized_v2",
                "donor_id_harmonized_v2",
                "sample_identity_rule_v2",
                "tcr_library_join_id_v2",
                "metadata_rule_id_v2",
                "metadata_evidence_url_v2",
                "metadata_evidence_level_v2",
                "metadata_status_v2",
            ]
            output_frame = frame[output_columns].fillna("").astype(str)
            table = pa.Table.from_pandas(output_frame, preserve_index=False)
            if writer is None:
                writer = pq.ParquetWriter(temporary_path, table.schema, compression="zstd")
            writer.write_table(table, row_group_size=250_000)

            for value, count in output_frame["tissue_harmonized_v2"].value_counts().items():
                tissue_counts[value] = tissue_counts.get(value, 0) + int(count)
            for value, count in output_frame["specimen_context_harmonized_v2"].value_counts().items():
                context_counts[value] = context_counts.get(value, 0) + int(count)
            for value, count in output_frame["tumor_type_harmonized_v2"].value_counts().items():
                tumor_counts[value] = tumor_counts.get(value, 0) + int(count)

            source_summary = {
                "source_gse_id": source,
                "n_cells": len(output_frame),
                "n_tissue_unresolved": int(output_frame["tissue_harmonized_v2"].eq("unresolved").sum()),
                "n_context_unresolved": int(output_frame["specimen_context_harmonized_v2"].eq("unresolved").sum()),
                "n_tumor_type_unresolved": int(
                    output_frame["tumor_type_harmonized_v2"].str.startswith("unresolved").sum()
                ),
                "n_sample_unresolved": int(output_frame["sample_id_harmonized_v2"].eq("unresolved").sum()),
                "n_sample_blank": int((~nonempty(output_frame["sample_id_harmonized_v2"])).sum()),
                "n_library_blank": int((~nonempty(output_frame["library_id_harmonized_v2"])).sum()),
                "n_donor_blank": int((~nonempty(output_frame["donor_id_harmonized_v2"])).sum()),
                "n_tcr_library_join_id": int(nonempty(output_frame["tcr_library_join_id_v2"]).sum()),
            }
            summary_rows.append(source_summary)
            total_rows += len(output_frame)
            total_sample_unresolved += source_summary["n_sample_unresolved"]

    if writer is None:
        raise RuntimeError("No metadata rows were written")
    writer.close()
    temporary_path.replace(output_path)

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(output_dir / "metadata_overlay_by_source.csv", index=False)
    for name, counts in (
        ("tissue_harmonized_v2", tissue_counts),
        ("specimen_context_harmonized_v2", context_counts),
        ("tumor_type_harmonized_v2", tumor_counts),
    ):
        pd.DataFrame(
            [{name: value, "n_cells": count} for value, count in counts.items()]
        ).sort_values("n_cells", ascending=False).to_csv(
            output_dir / f"{name}_counts.csv", index=False
        )

    atlas_after = file_state(atlas)
    source_files_unchanged = all(
        state == file_state((ROOT / next(
            row["path"] for row in config["source_metadata_overlays"] if row["source"] == source
        )).resolve())
        for source, state in source_states.items()
    )
    legacy_unchanged = legacy_before == file_state(legacy_path)
    checks = {
        "atlas_unchanged": atlas_before == atlas_after,
        "source_files_unchanged": source_files_unchanged,
        "legacy_review_unchanged": legacy_unchanged,
        "legacy_review_sha256": legacy_sha == legacy_rule["sha256"],
        "legacy_review_row_count": len(legacy) == int(legacy_rule["rows"]),
        "legacy_review_exact_atlas_coverage": (
            legacy_expected_rows == int(legacy_rule["rows"]) and legacy_unmatched_rows == 0
        ),
        "tcr_sidecar_sha256": sidecar_sha == EXPECTED_TCR_SIDECAR_SHA256,
        "row_count": total_rows == EXPECTED_ATLAS_CELLS,
        "all_source_overlay_rows_match": overlay_source_unmatched == 0,
        "all_affected_rows_have_tcr_library_identity": tcr_identity_unmatched == 0,
        "expected_tcr_library_identity_rows": total_tcr_identity == EXPECTED_TCR_IDENTITY_ROWS,
        "sample_nonblank": int(summary["n_sample_blank"].sum()) == 0,
        "library_nonblank": int(summary["n_library_blank"].sum()) == 0,
        "only_expected_sample_unresolved": total_sample_unresolved == 4_611,
        "no_cell_type_specimen_values": not (set(context_counts) & FORBIDDEN_CONTEXT_VALUES),
        "tissue_controlled_vocabulary": set(tissue_counts) <= set(
            config["controlled_vocabularies"]["tissue_harmonized_v2"]
        ),
        "context_controlled_vocabulary": set(context_counts) <= set(
            config["controlled_vocabularies"]["specimen_context_harmonized_v2"]
        ),
        "tumor_controlled_vocabulary": set(tumor_counts) <= set(
            config["controlled_vocabularies"]["tumor_type_harmonized_v2"]
        ),
    }
    manifest = {
        "status": "PASS_METADATA_OVERLAY_READY" if all(checks.values()) else "FAIL_METADATA_OVERLAY",
        "h5ad_write_performed": False,
        "tcr_chain_calls_applied": False,
        "atlas": str(atlas),
        "overlay": str(output_path),
        "overlay_sha256": sha256_file(output_path),
        "overlay_size_bytes": output_path.stat().st_size,
        "n_rows": total_rows,
        "n_sources": len(summary),
        "n_tissue_unresolved": int(summary["n_tissue_unresolved"].sum()),
        "n_context_unresolved": int(summary["n_context_unresolved"].sum()),
        "n_tumor_type_unresolved": int(summary["n_tumor_type_unresolved"].sum()),
        "n_sample_unresolved": total_sample_unresolved,
        "n_donor_blank": int(summary["n_donor_blank"].sum()),
        "n_tcr_library_join_id": total_tcr_identity,
        "checks": checks,
    }
    (output_dir / "metadata_overlay_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    if not all(checks.values()):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
