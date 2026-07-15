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


import io
import re
import sys
import tarfile
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd


ROOT = _TNK_PROJECT_ROOT / "downloads"
OUT_DIR = Path("analysis_output_v3")

# GSEs verified to have raw TCR evidence but were not integrated in prior pass.
TARGET_GSES = [
    "GSE145926",
    "GSE171037",
    "GSE178882",
    "GSE190870",
    "GSE212217",
    "GSE228597",
    "GSE243905",
    "GSE254249",
]

FINAL_COLS = [
    "cell_barcode",
    "GSE",
    "sample",
    "TRA_cdr3",
    "TRA_v",
    "TRA_d",
    "TRA_j",
    "TRA_cdr3_nt",
    "TRA_clone_id",
    "TRA_umis",
    "TRA_reads",
    "TRB_cdr3",
    "TRB_v",
    "TRB_d",
    "TRB_j",
    "TRB_cdr3_nt",
    "TRB_clone_id",
    "TRB_umis",
    "TRB_reads",
]


def _clean(x: object) -> str:
    s = str(x).strip()
    if s.lower() in {"", "nan", "none", "na"}:
        return ""
    return s


def _norm_barcode(x: object) -> str:
    s = _clean(x)
    if not s:
        return ""
    s = s.split("_")[0]
    m = re.match(r"^([ACGTN]+-\d+)(-\d+)+$", s, flags=re.I)
    if m:
        return m.group(1)
    return s


def _first_nonempty(vals: pd.Series) -> str:
    for v in vals:
        s = _clean(v)
        if s:
            return s
    return ""


def _num_or_zero(vals: pd.Series) -> int:
    x = pd.to_numeric(vals, errors="coerce").fillna(0)
    if len(x) == 0:
        return 0
    return int(x.max())


def _sample_from_name(name: str) -> str:
    base = Path(name).name
    base = re.sub(r"\.csv(\.gz)?$", "", base, flags=re.I)
    base = re.sub(r"\.txt(\.gz)?$", "", base, flags=re.I)
    base = re.sub(r"_(all|filtered)_contig_annotations(_tcr)?$", "", base, flags=re.I)
    base = re.sub(r"_scTCR$", "", base, flags=re.I)
    base = re.sub(r"_TCR(_enriched)?$", "", base, flags=re.I)
    base = re.sub(r"_10X_VDJ\\.merge$", "", base, flags=re.I)
    return base


def _preferred_file_map(files: List[Path]) -> List[Path]:
    # De-duplicate same basename across extracted_final/extracted_filtered/extracted_correct.
    score = {
        "extracted_correct": 3,
        "extracted_final": 2,
        "extracted_filtered": 1,
    }
    chosen: Dict[str, Path] = {}
    chosen_score: Dict[str, int] = {}
    for p in files:
        key = p.name
        s = 0
        pp = str(p)
        for k, v in score.items():
            if k in pp:
                s = max(s, v)
        if key not in chosen or s > chosen_score[key]:
            chosen[key] = p
            chosen_score[key] = s
    return list(chosen.values())


def _contig_to_compatible(df: pd.DataFrame, gse: str, sample: str) -> pd.DataFrame:
    cols = {c.lower(): c for c in df.columns}
    if "barcode" not in cols or "chain" not in cols:
        return pd.DataFrame(columns=FINAL_COLS)

    bcol = cols["barcode"]
    ccol = cols["chain"]
    vcol = cols.get("v_gene", cols.get("v_call"))
    dcol = cols.get("d_gene", cols.get("d_call"))
    jcol = cols.get("j_gene", cols.get("j_call"))
    cdr3_col = cols.get("cdr3")
    cdr3nt_col = cols.get("cdr3_nt", cols.get("junction"))
    clon_col = cols.get("raw_clonotype_id", cols.get("clonotype_id", cols.get("clone_id")))
    umis_col = cols.get("umis")
    reads_col = cols.get("reads")

    work = pd.DataFrame(
        {
            "cell_barcode": df[bcol].map(_norm_barcode),
            "chain": df[ccol].astype(str).str.upper(),
            "v": df[vcol] if vcol else "",
            "d": df[dcol] if dcol else "",
            "j": df[jcol] if jcol else "",
            "cdr3": df[cdr3_col] if cdr3_col else "",
            "cdr3_nt": df[cdr3nt_col] if cdr3nt_col else "",
            "clone_id": df[clon_col] if clon_col else "",
            "umis": pd.to_numeric(df[umis_col], errors="coerce").fillna(0) if umis_col else 0,
            "reads": pd.to_numeric(df[reads_col], errors="coerce").fillna(0) if reads_col else 0,
        }
    )
    work = work[work["cell_barcode"] != ""]
    work = work[work["chain"].isin(["TRA", "TRB"])]
    if work.empty:
        return pd.DataFrame(columns=FINAL_COLS)

    out_rows = []
    for (bc, ch), sub in work.groupby(["cell_barcode", "chain"], dropna=False):
        out_rows.append(
            {
                "cell_barcode": bc,
                "sample": sample,
                "chain": ch,
                "cdr3": _first_nonempty(sub["cdr3"]),
                "v": _first_nonempty(sub["v"]),
                "d": _first_nonempty(sub["d"]),
                "j": _first_nonempty(sub["j"]),
                "cdr3_nt": _first_nonempty(sub["cdr3_nt"]),
                "clone_id": _first_nonempty(sub["clone_id"]),
                "umis": _num_or_zero(sub["umis"]),
                "reads": _num_or_zero(sub["reads"]),
            }
        )
    flat = pd.DataFrame(out_rows)

    tra = flat[flat["chain"] == "TRA"].rename(
        columns={
            "cdr3": "TRA_cdr3",
            "v": "TRA_v",
            "d": "TRA_d",
            "j": "TRA_j",
            "cdr3_nt": "TRA_cdr3_nt",
            "clone_id": "TRA_clone_id",
            "umis": "TRA_umis",
            "reads": "TRA_reads",
        }
    )
    trb = flat[flat["chain"] == "TRB"].rename(
        columns={
            "cdr3": "TRB_cdr3",
            "v": "TRB_v",
            "d": "TRB_d",
            "j": "TRB_j",
            "cdr3_nt": "TRB_cdr3_nt",
            "clone_id": "TRB_clone_id",
            "umis": "TRB_umis",
            "reads": "TRB_reads",
        }
    )

    merged = pd.merge(
        tra[
            [
                "cell_barcode",
                "sample",
                "TRA_cdr3",
                "TRA_v",
                "TRA_d",
                "TRA_j",
                "TRA_cdr3_nt",
                "TRA_clone_id",
                "TRA_umis",
                "TRA_reads",
            ]
        ],
        trb[
            [
                "cell_barcode",
                "sample",
                "TRB_cdr3",
                "TRB_v",
                "TRB_d",
                "TRB_j",
                "TRB_cdr3_nt",
                "TRB_clone_id",
                "TRB_umis",
                "TRB_reads",
            ]
        ],
        on=["cell_barcode", "sample"],
        how="outer",
    )
    merged["GSE"] = gse
    for c in FINAL_COLS:
        if c not in merged.columns:
            merged[c] = ""
    merged = merged[FINAL_COLS]
    return merged


def _process_regular_contig_file(path: Path, gse: str) -> pd.DataFrame:
    try:
        # Most files are csv.gz
        df = pd.read_csv(path)
    except Exception:
        try:
            # Fallback for rare delimiter cases
            df = pd.read_csv(path, sep="\t")
        except Exception:
            return pd.DataFrame(columns=FINAL_COLS)
    sample = _sample_from_name(path.name)
    return _contig_to_compatible(df, gse, sample)


def _process_tar_with_contig(path: Path, gse: str) -> pd.DataFrame:
    rows = []
    sample = _sample_from_name(path.name)
    try:
        with tarfile.open(path, "r:*") as tar:
            for m in tar.getmembers():
                if not m.isfile():
                    continue
                name = m.name.lower()
                if "contig_annotations" not in name:
                    continue
                if not (name.endswith(".csv") or name.endswith(".csv.gz")):
                    continue
                f = tar.extractfile(m)
                if f is None:
                    continue
                payload = f.read()
                try:
                    if name.endswith(".gz"):
                        bio = io.BytesIO(payload)
                        df = pd.read_csv(bio, compression="gzip")
                    else:
                        bio = io.BytesIO(payload)
                        df = pd.read_csv(bio)
                except Exception:
                    continue
                rows.append(_contig_to_compatible(df, gse, sample))
    except Exception:
        return pd.DataFrame(columns=FINAL_COLS)
    if not rows:
        return pd.DataFrame(columns=FINAL_COLS)
    out = pd.concat(rows, ignore_index=True)
    return out


def _process_special_190870(path: Path) -> pd.DataFrame:
    # This file already holds TRA/TRB columns in metadata style.
    try:
        df = pd.read_csv(path)
    except Exception:
        return pd.DataFrame(columns=FINAL_COLS)

    if "cell_id" not in df.columns:
        return pd.DataFrame(columns=FINAL_COLS)

    out = pd.DataFrame()
    out["cell_barcode"] = (
        df["cell_id"].astype(str).str.replace(r"^[^_]+_", "", regex=True).map(_norm_barcode)
    )
    out["GSE"] = "GSE190870"
    out["sample"] = df.get("orig.ident", "").astype(str)

    out["TRA_cdr3"] = df.get("TRA_1_cdr3", "")
    out["TRA_v"] = df.get("TRA_1_v_gene", "")
    out["TRA_d"] = df.get("TRA_1_d_gene", "")
    out["TRA_j"] = df.get("TRA_1_j_gene", "")
    out["TRA_cdr3_nt"] = df.get("TRA_1_cdr3_nt", "")
    out["TRA_clone_id"] = df.get("clonotype", "")
    out["TRA_umis"] = df.get("TRA_1_expr", 0)
    out["TRA_reads"] = df.get("TRA_1_expr", 0)

    out["TRB_cdr3"] = df.get("TRB_1_cdr3", "")
    out["TRB_v"] = df.get("TRB_1_v_gene", "")
    out["TRB_d"] = df.get("TRB_1_d_gene", "")
    out["TRB_j"] = df.get("TRB_1_j_gene", "")
    out["TRB_cdr3_nt"] = df.get("TRB_1_cdr3_nt", "")
    out["TRB_clone_id"] = df.get("clonotype", "")
    out["TRB_umis"] = df.get("TRB_1_expr", 0)
    out["TRB_reads"] = df.get("TRB_1_expr", 0)

    out = out[out["cell_barcode"] != ""]
    for c in FINAL_COLS:
        if c not in out.columns:
            out[c] = ""
    return out[FINAL_COLS]


def _process_special_254249(path: Path) -> pd.DataFrame:
    # Although suffix is txt.gz, content is comma-separated with standard 10x VDJ columns.
    try:
        df = pd.read_csv(path, sep=",")
    except Exception:
        return pd.DataFrame(columns=FINAL_COLS)

    # Some parses can collapse into one giant column; try split fallback.
    if len(df.columns) == 1 and "," in df.columns[0]:
        raw = pd.read_csv(path, header=None)
        header = raw.iloc[0, 0].split(",")
        body = raw.iloc[1:, 0].str.split(",", expand=True)
        body.columns = header
        df = body.reset_index(drop=True)

    return _contig_to_compatible(df, "GSE254249", "GSE254249")


def _process_special_228597(path: Path) -> pd.DataFrame:
    # Prefer provided combined h5ad with integrated TCR metadata.
    try:
        import anndata as ad
    except Exception:
        return pd.DataFrame(columns=FINAL_COLS)
    try:
        adata = ad.read_h5ad(path, backed="r")
        obs = adata.obs.copy()
    except Exception:
        return pd.DataFrame(columns=FINAL_COLS)

    out = pd.DataFrame()
    out["cell_barcode"] = pd.Series(obs.index.astype(str)).map(_norm_barcode).values
    out["GSE"] = "GSE228597"
    out["sample"] = obs.get("sample_id", "").astype(str).values

    out["TRA_cdr3"] = obs.get("TRA_cdr3", "").astype(str).values
    out["TRA_v"] = obs.get("TRA_v_gene", "").astype(str).values
    out["TRA_d"] = ""
    out["TRA_j"] = obs.get("TRA_j_gene", "").astype(str).values
    out["TRA_cdr3_nt"] = obs.get("TRA_cdr3_nt", "").astype(str).values
    out["TRA_clone_id"] = ""
    out["TRA_umis"] = 0
    out["TRA_reads"] = 0

    out["TRB_cdr3"] = obs.get("TRB_cdr3", "").astype(str).values
    out["TRB_v"] = obs.get("TRB_v_gene", "").astype(str).values
    out["TRB_d"] = ""
    out["TRB_j"] = obs.get("TRB_j_gene", "").astype(str).values
    out["TRB_cdr3_nt"] = obs.get("TRB_cdr3_nt", "").astype(str).values
    out["TRB_clone_id"] = ""
    out["TRB_umis"] = 0
    out["TRB_reads"] = 0

    out = out[out["cell_barcode"] != ""]
    out = out[FINAL_COLS]
    return out


def _dedupe(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df
    key_cols = ["cell_barcode", "GSE", "sample"]
    keep_cols = [c for c in FINAL_COLS if c not in key_cols]
    # collapse possible duplicates from repeated extracted_* directories.
    out = (
        df.groupby(key_cols, as_index=False)[keep_cols]
        .agg(_first_nonempty)
        .reset_index(drop=True)
    )
    return out[FINAL_COLS]


def integrate_one_gse(gse: str) -> pd.DataFrame:
    gdir = ROOT / gse
    if not gdir.exists():
        return pd.DataFrame(columns=FINAL_COLS)

    rows = []

    if gse == "GSE190870":
        p = gdir / "GSE167036_meta_tcr.csv.gz"
        if p.exists():
            rows.append(_process_special_190870(p))
    elif gse == "GSE228597":
        p = gdir / "suppl" / "GSE228597_combined_pbmc_data_with_tcr.h5ad"
        if p.exists():
            rows.append(_process_special_228597(p))
    elif gse == "GSE254249":
        p = gdir / "suppl" / "GSE254249_10X_VDJ.merge.txt.gz"
        if p.exists():
            rows.append(_process_special_254249(p))
    else:
        all_files: List[Path] = []
        for sub in [gdir / "suppl"]:
            if not sub.exists():
                continue
            if gse == "GSE243905":
                # Use per-sample TCR tarballs so sample identity is preserved for remapping.
                for extracted_name in [
                    "extracted_correct_GSE243905",
                    "extracted_final_GSE243905",
                    "extracted_filtered_GSE243905",
                ]:
                    extracted_dir = sub / extracted_name
                    if extracted_dir.exists():
                        all_files.extend(extracted_dir.rglob("*_TCR.tar.gz"))
                all_files.extend(sub.rglob("*_TCR.tar.gz"))
            else:
                all_files.extend(sub.rglob("*filtered_contig_annotations*.csv.gz"))
                all_files.extend(sub.rglob("*all_contig_annotations*_tcr.csv.gz"))
                all_files.extend(sub.rglob("*all_contig_annotations*.csv.gz"))
                all_files.extend(sub.rglob("*_TCR.tar.gz"))
                all_files.extend(sub.rglob("*tcr_aggr.tar.gz"))

        files = _preferred_file_map(sorted(set(all_files)))
        for fp in files:
            low = fp.name.lower()
            if low.endswith(".tar.gz"):
                rows.append(_process_tar_with_contig(fp, gse))
            else:
                rows.append(_process_regular_contig_file(fp, gse))

    if not rows:
        return pd.DataFrame(columns=FINAL_COLS)
    out = pd.concat(rows, ignore_index=True)
    out = _dedupe(out)
    return out


def rebuild_global_tcr() -> pd.DataFrame:
    files = sorted(OUT_DIR.glob("tcr_GSE*.csv"))
    frames = []
    for f in files:
        try:
            df = pd.read_csv(f, low_memory=False)
        except Exception:
            continue
        if "GSE" not in df.columns:
            gse = f.stem.replace("tcr_", "")
            df["GSE"] = gse
        for c in FINAL_COLS:
            if c not in df.columns:
                df[c] = ""
        frames.append(df[FINAL_COLS])
    if not frames:
        return pd.DataFrame(columns=FINAL_COLS)
    merged = pd.concat(frames, ignore_index=True)
    # Keep first row per key for speed on large combined table.
    merged = merged.drop_duplicates(subset=["cell_barcode", "GSE", "sample"], keep="first")
    return merged


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    target_gses = TARGET_GSES
    if len(sys.argv) > 1:
        target_gses = [x.strip() for x in sys.argv[1:] if x.strip()]

    report_rows = []
    for gse in target_gses:
        df = integrate_one_gse(gse)
        outp = OUT_DIR / f"tcr_{gse}.csv"
        df.to_csv(outp, index=False)
        report_rows.append({"GSE": gse, "tcr_rows": int(len(df)), "output": str(outp)})
        print(f"{gse}: wrote {len(df)} rows -> {outp}")

    all_df = rebuild_global_tcr()
    all_path = OUT_DIR / "tcr_all_data_complete.csv"
    all_df.to_csv(all_path, index=False)
    print(f"Rebuilt {all_path} with {len(all_df)} rows")

    report = pd.DataFrame(report_rows).sort_values("GSE")
    report_path = _TNK_PROJECT_ROOT / "analysis_26GSE_V4/reports/tcr_newly_integrated_gse_report.csv"
    report.to_csv(report_path, index=False)
    print(f"Wrote {report_path}")
    print(report.to_string(index=False))


if __name__ == "__main__":
    main()
