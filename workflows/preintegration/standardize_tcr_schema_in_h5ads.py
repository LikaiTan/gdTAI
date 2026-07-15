#!/usr/bin/env python3
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


import sys
from pathlib import Path
from typing import Dict, Iterable, List

import anndata as ad
import pandas as pd

ad.settings.allow_write_nullable_strings = True


ROOT_SCANPY = _TNK_PROJECT_ROOT / "analysis_26GSE_V4/scanpy_projects"
ROOT_METADATA = _TNK_PROJECT_ROOT / "downloads/per_gse_h5ad_with_metadata"
OUT_REPORT = _TNK_PROJECT_ROOT / "analysis_26GSE_V4/reports/tcr_schema_standardization_report.csv"
SKIP_FILES = {"GSE221776_with_metadata.h5ad"}


CHAIN_SYNONYMS: Dict[str, List[str]] = {
    "TRA_cdr3": ["TRA_cdr3", "TCR_TRA_cdr3", "TCR_TRA_cdr3_x", "TCR_TRA_cdr3_y", "cdr3a"],
    "TRA_cdr3_nt": ["TRA_cdr3_nt", "TCR_TRA_cdr3_nt", "cdr3a_nt"],
    "TRA_v": ["TRA_v", "TCR_TRA_v", "TRA_v_gene"],
    "TRA_j": ["TRA_j", "TCR_TRA_j", "TRA_j_gene"],
    "TRA_clone_id": ["TRA_clone_id", "TCR_TRA_clone_id", "raw_alpha_clonotype_id"],
    "TRB_cdr3": ["TRB_cdr3", "TCR_TRB_cdr3", "TCR_TRB_cdr3_x", "TCR_TRB_cdr3_y", "cdr3b"],
    "TRB_cdr3_nt": ["TRB_cdr3_nt", "TCR_TRB_cdr3_nt", "cdr3b_nt", "tcr_TRB_cdr3_nt"],
    "TRB_v": ["TRB_v", "TCR_TRB_v", "TRB_v_gene", "tcr_TRB_v"],
    "TRB_d": ["TRB_d", "tcr_TRB_d"],
    "TRB_j": ["TRB_j", "TCR_TRB_j", "TRB_j_gene", "tcr_TRB_j"],
    "TRB_clone_id": ["TRB_clone_id", "TCR_TRB_clone_id", "raw_beta_clonotype_id", "tcr_TRB_clone_id"],
    "TRB_umis": ["TRB_umis", "tcr_TRB_umis"],
    "TRB_reads": ["TRB_reads", "tcr_TRB_reads"],
    "TRG_cdr3": ["TRG_cdr3", "cdr3g"],
    "TRG_cdr3_nt": ["TRG_cdr3_nt", "cdr3g_nt"],
    "TRG_clone_id": ["TRG_clone_id", "raw_gamma_clonotype_id"],
    "TRD_cdr3": ["TRD_cdr3", "cdr3d"],
    "TRD_cdr3_nt": ["TRD_cdr3_nt", "cdr3d_nt"],
    "TRD_clone_id": ["TRD_clone_id", "raw_delta_clonotype_id"],
}

TCRSEQ_ALTS = ["TCRseq", "withTCR", "has_tcr", "tcr_vdj_flag"]


def relevant_tcr_cols(cols: Iterable[str]) -> List[str]:
    out = []
    for c in cols:
        lc = c.lower()
        if (
            c == "TCRseq"
            or c.startswith("TRA")
            or c.startswith("TRB")
            or c.startswith("TRG")
            or c.startswith("TRD")
            or c.startswith("TCR_")
            or lc.startswith("tcr_")
            or "cdr3" in lc
            or "clon" in lc
            or lc in {"withtcr", "has_tcr", "tcr_vdj_flag"}
        ):
            out.append(c)
    return out


def empty_to_na(s: pd.Series) -> pd.Series:
    x = s.astype("string")
    x = x.str.strip()
    return x.replace({"": pd.NA, "nan": pd.NA, "None": pd.NA, "none": pd.NA, "NA": pd.NA, "na": pd.NA})


def first_nonempty(df: pd.DataFrame, cols: List[str]) -> pd.Series:
    tmp = pd.DataFrame({c: empty_to_na(df[c]) for c in cols})
    return tmp.bfill(axis=1).iloc[:, 0]


def normalize_tcrseq_value_series(s: pd.Series) -> pd.Series:
    x = s.astype("string").str.strip().str.lower()
    out = pd.Series(pd.NA, index=s.index, dtype="string")
    yes_vals = {"yes", "true", "1", "with tcr/vdj", "with tcr", "tcr", "vdj"}
    no_vals = {"no", "false", "0", "without tcr/vdj", "without tcr", ""}
    out[x.isin(yes_vals)] = "yes"
    out[x.isin(no_vals)] = "no"
    return out


def series_equal_as_string(a: pd.Series, b: pd.Series) -> bool:
    ax = a.astype("string").fillna("<NA>")
    bx = b.astype("string").fillna("<NA>")
    return ax.equals(bx)


def standardize_one_h5ad(path: Path) -> Dict[str, object]:
    adata = ad.read_h5ad(path)
    before_obs = adata.obs.copy()
    obs = before_obs.copy()
    before_cols = set(before_obs.columns)
    before_tcr_cols = relevant_tcr_cols(before_cols)
    before_yes = 0
    if "TCRseq" in before_obs.columns:
        before_yes = int((before_obs["TCRseq"].astype(str).str.lower() == "yes").sum())

    drop_cols = set()
    added_standard_cols = []
    changed = False

    for std_col, synonyms in CHAIN_SYNONYMS.items():
        present = [c for c in synonyms if c in obs.columns]
        if not present:
            continue
        new_series = first_nonempty(obs, present)
        if std_col not in before_obs.columns or not series_equal_as_string(before_obs[std_col], new_series):
            changed = True
        obs[std_col] = new_series
        if std_col not in before_cols:
            added_standard_cols.append(std_col)
        for c in present:
            if c != std_col:
                drop_cols.add(c)

    # Standardize TCRseq from existing indicators but do not infer from sample-level metadata.
    if any(c in obs.columns for c in TCRSEQ_ALTS):
        tcrseq = pd.Series(pd.NA, index=obs.index, dtype="string")
        for c in TCRSEQ_ALTS:
            if c in obs.columns:
                cur = normalize_tcrseq_value_series(obs[c])
                tcrseq = tcrseq.fillna(cur)
                if c != "TCRseq":
                    drop_cols.add(c)

        chain_cols = [c for c in CHAIN_SYNONYMS if c in obs.columns]
        if chain_cols:
            has_any_chain = pd.Series(False, index=obs.index)
            for c in chain_cols:
                has_any_chain = has_any_chain | empty_to_na(obs[c]).notna()
            tcrseq.loc[has_any_chain] = "yes"
            tcrseq = tcrseq.fillna("no")

        if "TCRseq" not in before_obs.columns or not series_equal_as_string(before_obs.get("TCRseq", pd.Series(pd.NA, index=obs.index, dtype="string")), tcrseq):
            changed = True
        obs["TCRseq"] = tcrseq

    # Drop old, incompatible chain columns after coalescing.
    safe_drop = [c for c in drop_cols if c in obs.columns and c not in CHAIN_SYNONYMS]
    if safe_drop:
        changed = True
        obs = obs.drop(columns=safe_drop)

    # Normalize TCR-related textual columns to stable categorical strings.
    categorical_cols = [
        "TCRseq",
        "TRA_cdr3",
        "TRA_cdr3_nt",
        "TRA_v",
        "TRA_j",
        "TRA_clone_id",
        "TRB_cdr3",
        "TRB_cdr3_nt",
        "TRB_v",
        "TRB_d",
        "TRB_j",
        "TRB_clone_id",
        "TRB_umis",
        "TRB_reads",
        "TRG_cdr3",
        "TRG_cdr3_nt",
        "TRG_clone_id",
        "TRD_cdr3",
        "TRD_cdr3_nt",
        "TRD_clone_id",
        "clonotype",
        "clonotype_number",
    ]
    for c in categorical_cols:
        if c in obs.columns:
            if c == "TCRseq":
                s = normalize_tcrseq_value_series(obs[c]).fillna("no")
            else:
                s = empty_to_na(obs[c])
            obs[c] = s.astype("category")

    # Reorder TCR columns to keep the standardized schema compact.
    preferred_front = [
        "TCRseq",
        "TRA_cdr3",
        "TRA_cdr3_nt",
        "TRA_v",
        "TRA_j",
        "TRA_clone_id",
        "TRB_cdr3",
        "TRB_cdr3_nt",
        "TRB_v",
        "TRB_d",
        "TRB_j",
        "TRB_clone_id",
        "TRB_umis",
        "TRB_reads",
        "TRG_cdr3",
        "TRG_cdr3_nt",
        "TRG_clone_id",
        "TRD_cdr3",
        "TRD_cdr3_nt",
        "TRD_clone_id",
        "clonotype",
        "clonotype_number",
    ]
    present_front = [c for c in preferred_front if c in obs.columns]
    remain = [c for c in obs.columns if c not in present_front]
    if list(obs.columns) != present_front + remain:
        changed = True
    obs = obs[present_front + remain]

    if "_index" in obs.columns:
        changed = True
        obs = obs.rename(columns={"_index": "obs__index"})

    if not changed:
        return {
            "path": str(path),
            "before_tcr_cols_n": len(before_tcr_cols),
            "after_tcr_cols_n": len(before_tcr_cols),
            "before_TCRseq_yes": before_yes,
            "after_TCRseq_yes": before_yes,
            "added_standard_cols": "",
            "dropped_old_cols": "",
            "changed": "no",
        }

    adata.obs = obs
    tmp_path = path.with_name(path.name + ".tmp")
    if tmp_path.exists():
        tmp_path.unlink()
    adata.write_h5ad(tmp_path)
    tmp_path.replace(path)

    after_cols = set(obs.columns)
    after_tcr_cols = relevant_tcr_cols(after_cols)
    after_yes = 0
    if "TCRseq" in obs.columns:
        after_yes = int((obs["TCRseq"].astype(str).str.lower() == "yes").sum())

    return {
        "path": str(path),
        "before_tcr_cols_n": len(before_tcr_cols),
        "after_tcr_cols_n": len(after_tcr_cols),
        "before_TCRseq_yes": before_yes,
        "after_TCRseq_yes": after_yes,
        "added_standard_cols": "|".join(added_standard_cols),
        "dropped_old_cols": "|".join(sorted(safe_drop)),
        "changed": "yes",
    }


def iter_targets() -> List[Path]:
    out: List[Path] = []
    for gdir in sorted(ROOT_SCANPY.glob("GSE*")):
        h5s = sorted((gdir / "outputs").glob("GSE*_scanpy_processed.h5ad"))
        out.extend(h5s)
    out.extend(sorted(p for p in ROOT_METADATA.glob("*.h5ad") if p.name not in SKIP_FILES))
    return out


def resolve_targets(args: List[str]) -> List[Path]:
    if not args:
        return iter_targets()

    resolved: List[Path] = []
    for arg in args:
        p = Path(arg)
        if p.exists():
            resolved.append(p)
            continue

        scanpy_matches = sorted(ROOT_SCANPY.glob(f"{arg}/outputs/GSE*_scanpy_processed.h5ad"))
        metadata_matches = sorted(ROOT_METADATA.glob(f"{arg}*.h5ad"))
        matches = scanpy_matches + metadata_matches
        if not matches:
            raise SystemExit(f"No h5ad target matched: {arg}")
        resolved.extend(matches)
    return resolved


def main() -> None:
    rows = []
    for path in resolve_targets(sys.argv[1:]):
        try:
            rows.append(standardize_one_h5ad(path))
        except Exception as e:
            rows.append(
                {
                    "path": str(path),
                    "before_tcr_cols_n": "",
                    "after_tcr_cols_n": "",
                    "before_TCRseq_yes": "",
                    "after_TCRseq_yes": "",
                    "added_standard_cols": "",
                    "dropped_old_cols": "",
                    "error": str(e),
                }
            )
    pd.DataFrame(rows).to_csv(OUT_REPORT, index=False)
    print(f"wrote {OUT_REPORT}")


if __name__ == "__main__":
    main()
