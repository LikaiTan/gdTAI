#!/usr/bin/env python3
"""Audit biological sample identities before repaired TCR calls are applied.

The TCR sidecar's join-sample field is retained as a technical VDJ-library
identity. It is never used wholesale as a biological sample identifier.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_ATLAS = ROOT / "Integrated_dataset/integrated_full_atlas.h5ad"
DEFAULT_TCR_SIDECAR = (
    ROOT
    / "Integrated_dataset/tables/tcr_join_rebuild/validated_tcr_replacement_sidecar.parquet"
)
DEFAULT_OUTPUT = ROOT / "Integrated_dataset/tables/metadata_harmonization/full_atlas_v2"
SOURCE_IDENTITY_FILES = {
    "GSE125527": ROOT / "data/datasets/GSE125527/processed/current.h5ad",
    "GSE254249": ROOT / "data/datasets/GSE254249/processed/current.h5ad",
}
EXPECTED_ATLAS_CELLS = 5_933_312
EXPECTED_SIDECAR_ROWS = 3_041_871
EXPECTED_SIDECAR_SOURCES = 14
EXPECTED_SIDECAR_SHA256 = "3114e70719301d693ae1a2bc2c63bac6c8bd57e3e8ac73a88c24320eaabfc2f0"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
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


def file_state(path: Path) -> dict[str, int | str]:
    stat = path.stat()
    return {
        "path": str(path.resolve()),
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "inode": stat.st_ino,
    }


def nonempty(values: pd.Series) -> pd.Series:
    return values.fillna("").astype(str).str.strip().ne("")


def normalize_token(value: str) -> str:
    return "_".join(value.strip().casefold().replace("/", "_").split())


def source_obs_index(obs: h5py.Group) -> np.ndarray:
    index_name = obs.attrs.get("_index", "_index")
    if isinstance(index_name, bytes):
        index_name = index_name.decode()
    return values_from_obs(obs, str(index_name))


def load_source_identity(source: str, path: Path) -> tuple[pd.DataFrame, dict[str, int | str]]:
    resolved = path.resolve()
    before = file_state(resolved)
    with h5py.File(resolved, "r") as handle:
        obs = handle["obs"]
        index = source_obs_index(obs)
        if source == "GSE125527":
            donor = values_from_obs(obs, "donor")
            tissue = values_from_obs(obs, "tissue")
            sample = np.asarray(
                [f"{normalize_token(d)}__{normalize_token(t)}" for d, t in zip(donor, tissue)],
                dtype=object,
            )
            frame = pd.DataFrame(
                {
                    "source_obs_name": index,
                    "source_sample_id": sample,
                    "source_library_id": sample,
                    "source_donor_id": donor,
                    "source_identity_rule": "source_donor_by_tissue",
                }
            )
        elif source == "GSE254249":
            sample = values_from_obs(obs, "Ident")
            frame = pd.DataFrame(
                {
                    "source_obs_name": index,
                    "source_sample_id": sample,
                    "source_library_id": sample,
                    "source_donor_id": values_from_obs(obs, "PatientID"),
                    "source_identity_rule": "source_ident",
                }
            )
        else:
            raise ValueError(f"No source identity rule for {source}")
    if before != file_state(resolved):
        raise RuntimeError(f"Source H5AD changed during read-only identity audit: {source}")
    if not frame["source_obs_name"].is_unique:
        raise RuntimeError(f"Source H5AD index is not unique: {source}")
    return frame, before


def mapping_structure(
    frame: pd.DataFrame,
    source: str,
    current_column: str,
    proposed_column: str,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    usable = frame.loc[nonempty(frame[current_column]) & nonempty(frame[proposed_column])]
    for direction, left, right in (
        ("current_value_splits", current_column, proposed_column),
        ("proposed_value_merges", proposed_column, current_column),
    ):
        counts = usable.groupby(left, observed=True)[right].nunique()
        for value, n_targets in counts[counts > 1].items():
            subset = usable.loc[usable[left].eq(value)]
            rows.append(
                {
                    "source_gse_id": source,
                    "field": current_column.removeprefix("current_"),
                    "direction": direction,
                    "value": value,
                    "n_mapped_values": int(n_targets),
                    "n_cells": int(len(subset)),
                    "mapped_values": ";".join(sorted(subset[right].astype(str).unique())),
                }
            )
    return rows


def main() -> None:
    args = parse_args()
    atlas = args.atlas.resolve()
    sidecar_path = args.tcr_sidecar.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)

    atlas_before = file_state(atlas)
    source_identity = {}
    source_states = {}
    for source, path in SOURCE_IDENTITY_FILES.items():
        source_identity[source], source_states[source] = load_source_identity(source, path)

    sidecar_sha = sha256_file(sidecar_path)
    sidecar = pd.read_parquet(
        sidecar_path,
        columns=["source_gse_id", "source_obs_name", "join_sample_id"],
    ).fillna("")
    for column in sidecar.columns:
        sidecar[column] = sidecar[column].astype(str)
    sidecar_sources = sorted(sidecar["source_gse_id"].unique())

    source_rows: list[dict[str, object]] = []
    structure_rows: list[dict[str, object]] = []
    exact_tcr_identity_matches = 0
    affected_atlas_rows = 0
    affected_without_tcr_identity = 0
    source_recovery_rows = 0
    source_recovery_unmatched = 0

    with h5py.File(atlas, "r") as handle:
        obs = handle["obs"]
        source_col = CategoricalColumn(obs, "source_gse_id")
        original_col = CategoricalColumn(obs, "original_cell_id")
        sample_col = CategoricalColumn(obs, "sample_id")
        library_col = CategoricalColumn(obs, "library_id")
        donor_col = CategoricalColumn(obs, "donor_id")
        n_cells = len(source_col.codes)

        for source_code, source in enumerate(source_col.categories):
            positions = np.flatnonzero(source_col.codes == source_code)
            current = pd.DataFrame(
                {
                    "source_obs_name": original_col.values(positions),
                    "current_sample_id": sample_col.values(positions),
                    "current_library_id": library_col.values(positions),
                    "current_donor_id": donor_col.values(positions),
                }
            )
            proposed = current.copy()
            proposed["proposed_sample_id"] = current["current_sample_id"]
            proposed["proposed_library_id"] = current["current_library_id"]
            proposed["proposed_donor_id"] = current["current_donor_id"]
            proposed["sample_identity_rule"] = "preserve_current_identity"
            proposed["tcr_library_join_id"] = ""

            if source in sidecar_sources:
                source_sidecar = sidecar.loc[sidecar["source_gse_id"].eq(source)].drop(
                    columns="source_gse_id"
                )
                proposed = proposed.merge(
                    source_sidecar,
                    on="source_obs_name",
                    how="left",
                    validate="one_to_one",
                    indicator="_tcr_identity_merge",
                )
                matched = proposed["_tcr_identity_merge"].eq("both")
                exact_tcr_identity_matches += int(matched.sum())
                affected_atlas_rows += len(proposed)
                affected_without_tcr_identity += int((~matched).sum())
                proposed["tcr_library_join_id"] = proposed["join_sample_id"].fillna("")

            if source in source_identity:
                proposed = proposed.merge(
                    source_identity[source],
                    on="source_obs_name",
                    how="left",
                    validate="one_to_one",
                    indicator="_source_identity_merge",
                )
                recovered = proposed["_source_identity_merge"].eq("both")
                source_recovery_rows += len(proposed)
                source_recovery_unmatched += int((~recovered).sum())
                proposed.loc[recovered, "proposed_sample_id"] = proposed.loc[
                    recovered, "source_sample_id"
                ]
                proposed.loc[recovered, "proposed_library_id"] = proposed.loc[
                    recovered, "source_library_id"
                ]
                if source == "GSE125527":
                    proposed.loc[recovered, "proposed_donor_id"] = proposed.loc[
                        recovered, "source_donor_id"
                    ]
                proposed.loc[recovered, "sample_identity_rule"] = proposed.loc[
                    recovered, "source_identity_rule"
                ]

            if source == "GSE228597":
                sample_missing = ~nonempty(proposed["proposed_sample_id"])
                library_missing = ~nonempty(proposed["proposed_library_id"])
                proposed.loc[sample_missing, "proposed_sample_id"] = "unresolved"
                proposed.loc[library_missing, "proposed_library_id"] = proposed.loc[
                    library_missing, "tcr_library_join_id"
                ]
                proposed.loc[sample_missing, "sample_identity_rule"] = (
                    "unresolved_biological_sample_known_tcr_library"
                )

            structure_rows.extend(
                mapping_structure(
                    proposed,
                    source,
                    "current_sample_id",
                    "proposed_sample_id",
                )
            )
            structure_rows.extend(
                mapping_structure(
                    proposed,
                    source,
                    "current_library_id",
                    "proposed_library_id",
                )
            )

            source_rows.append(
                {
                    "source_gse_id": source,
                    "n_cells": int(len(proposed)),
                    "n_current_sample_blank": int((~nonempty(proposed["current_sample_id"])).sum()),
                    "n_proposed_sample_blank": int((~nonempty(proposed["proposed_sample_id"])).sum()),
                    "n_proposed_sample_unresolved": int(
                        proposed["proposed_sample_id"].eq("unresolved").sum()
                    ),
                    "n_current_sample_unique": int(
                        proposed.loc[nonempty(proposed["current_sample_id"]), "current_sample_id"].nunique()
                    ),
                    "n_proposed_sample_unique": int(
                        proposed.loc[nonempty(proposed["proposed_sample_id"]), "proposed_sample_id"].nunique()
                    ),
                    "n_sample_values_changed": int(
                        proposed["current_sample_id"].ne(proposed["proposed_sample_id"]).sum()
                    ),
                    "n_current_library_blank": int((~nonempty(proposed["current_library_id"])).sum()),
                    "n_proposed_library_blank": int((~nonempty(proposed["proposed_library_id"])).sum()),
                    "n_current_library_unique": int(
                        proposed.loc[nonempty(proposed["current_library_id"]), "current_library_id"].nunique()
                    ),
                    "n_proposed_library_unique": int(
                        proposed.loc[nonempty(proposed["proposed_library_id"]), "proposed_library_id"].nunique()
                    ),
                    "n_library_values_changed": int(
                        proposed["current_library_id"].ne(proposed["proposed_library_id"]).sum()
                    ),
                    "n_current_donor_blank": int((~nonempty(proposed["current_donor_id"])).sum()),
                    "n_proposed_donor_blank": int((~nonempty(proposed["proposed_donor_id"])).sum()),
                    "n_donor_values_changed": int(
                        proposed["current_donor_id"].ne(proposed["proposed_donor_id"]).sum()
                    ),
                    "n_tcr_library_join_id": int(nonempty(proposed["tcr_library_join_id"]).sum()),
                }
            )

    source_table = pd.DataFrame(source_rows)
    structure_table = pd.DataFrame(
        structure_rows,
        columns=[
            "source_gse_id",
            "field",
            "direction",
            "value",
            "n_mapped_values",
            "n_cells",
            "mapped_values",
        ],
    )
    source_table.to_csv(output / "sample_identity_by_source.csv", index=False)
    structure_table.to_csv(output / "sample_identity_mapping_structure.csv", index=False)

    atlas_after = file_state(atlas)
    source_files_unchanged = all(
        source_states[source] == file_state(path.resolve())
        for source, path in SOURCE_IDENTITY_FILES.items()
    )
    checks = {
        "atlas_cell_count": n_cells == EXPECTED_ATLAS_CELLS,
        "atlas_unchanged": atlas_before == atlas_after,
        "source_identity_files_unchanged": source_files_unchanged,
        "sidecar_sha256": sidecar_sha == EXPECTED_SIDECAR_SHA256,
        "sidecar_row_count": len(sidecar) == EXPECTED_SIDECAR_ROWS,
        "sidecar_source_count": len(sidecar_sources) == EXPECTED_SIDECAR_SOURCES,
        "sidecar_identity_keys_unique": not sidecar.duplicated(
            ["source_gse_id", "source_obs_name"]
        ).any(),
        "all_affected_atlas_rows_have_tcr_library_identity": (
            affected_without_tcr_identity == 0
        ),
        "all_source_recovery_rows_match": source_recovery_unmatched == 0,
        "proposed_sample_ids_nonblank": int(source_table["n_proposed_sample_blank"].sum()) == 0,
        "proposed_library_ids_nonblank": int(source_table["n_proposed_library_blank"].sum()) == 0,
        "gse125527_recovers_30_samples": int(
            source_table.loc[source_table["source_gse_id"].eq("GSE125527"), "n_proposed_sample_unique"].iloc[0]
        ) == 30,
        "gse254249_recovers_92_samples": int(
            source_table.loc[source_table["source_gse_id"].eq("GSE254249"), "n_proposed_sample_unique"].iloc[0]
        ) == 92,
    }
    manifest = {
        "status": "PASS_SAMPLE_IDENTITY_READY" if all(checks.values()) else "FAIL_SAMPLE_IDENTITY",
        "atlas": str(atlas),
        "tcr_chain_calls_applied": False,
        "tcr_join_library_is_biological_sample": False,
        "n_cells": n_cells,
        "n_sources": int(len(source_table)),
        "n_sidecar_identity_sources": len(sidecar_sources),
        "n_affected_atlas_rows": affected_atlas_rows,
        "n_exact_tcr_library_identity_matches": exact_tcr_identity_matches,
        "n_source_recovery_rows": source_recovery_rows,
        "n_current_sample_blank": int(source_table["n_current_sample_blank"].sum()),
        "n_proposed_sample_blank": int(source_table["n_proposed_sample_blank"].sum()),
        "n_proposed_sample_unresolved": int(source_table["n_proposed_sample_unresolved"].sum()),
        "n_current_library_blank": int(source_table["n_current_library_blank"].sum()),
        "n_proposed_library_blank": int(source_table["n_proposed_library_blank"].sum()),
        "n_current_donor_blank": int(source_table["n_current_donor_blank"].sum()),
        "n_proposed_donor_blank": int(source_table["n_proposed_donor_blank"].sum()),
        "n_sample_values_changed": int(source_table["n_sample_values_changed"].sum()),
        "n_library_values_changed": int(source_table["n_library_values_changed"].sum()),
        "n_donor_values_changed": int(source_table["n_donor_values_changed"].sum()),
        "n_mapping_structure_rows": int(len(structure_table)),
        "checks": checks,
    }
    (output / "sample_identity_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    if not all(checks.values()):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
