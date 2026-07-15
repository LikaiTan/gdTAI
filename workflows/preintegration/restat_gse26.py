
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

#!/usr/bin/env python3
import json
from pathlib import Path
import re
import pandas as pd
import anndata as ad

ROOT = _TNK_PROJECT_ROOT / "analysis_26GSE_V4/scanpy_projects"
OUT = _TNK_PROJECT_ROOT / "analysis_26GSE_V4/reports/gse_cells_genes_donors_restat.csv"

DONOR_CANDIDATES = [
    'donor', 'donor_id', 'patient', 'patient_id', 'patientid',
    'subject', 'subject_id', 'individual', 'sample_donor', 'orig.ident',
    'case_id', 'pid', 'donorid', 'sample'
]


def pick_donor_col(cols):
    low = {c.lower(): c for c in cols}
    # exact preferred
    for cand in DONOR_CANDIDATES:
        if cand in low:
            return low[cand]
    # fuzzy contains
    for c in cols:
        cl = c.lower()
        if any(k in cl for k in ['donor', 'patient', 'subject', 'individual']):
            return c
    return ''


def nunique_nonempty(s):
    if s is None:
        return None
    vals = pd.Series(s).astype(str).str.strip()
    vals = vals[(vals != '') & (vals.str.lower() != 'nan') & (vals.str.lower() != 'none')]
    return int(vals.nunique()) if len(vals) else None

rows = []
for gdir in sorted(ROOT.glob('GSE*')):
    gse = gdir.name
    cfgp = gdir / 'config' / 'config.json'
    h5p = gdir / 'outputs' / 'scanpy_processed.h5ad'

    cfg = {}
    if cfgp.exists():
        try:
            cfg = json.loads(cfgp.read_text())
        except Exception:
            cfg = {}

    has_metadata = 'YES' if cfg.get('metadata_file') else 'NO'
    metadata_file = cfg.get('metadata_file', '')

    if not h5p.exists():
        rows.append({
            'GSE': gse,
            'has_metadata': has_metadata,
            'status': 'no_h5ad',
            'cells': None,
            'genes': None,
            'donor_count_h5ad_matched': None,
            'donor_col_h5ad': '',
            'donor_count_source_metadata': None,
            'donor_col_source_metadata': '',
            'metadata_file': metadata_file,
        })
        continue

    try:
        adata = ad.read_h5ad(h5p)
        cells, genes = int(adata.n_obs), int(adata.n_vars)

        donor_col_h5ad = pick_donor_col(list(adata.obs.columns))
        donor_count_h5ad = nunique_nonempty(adata.obs[donor_col_h5ad]) if donor_col_h5ad else None

        donor_count_src = None
        donor_col_src = ''
        if metadata_file:
            mp = Path(metadata_file)
            if mp.exists():
                try:
                    mdf = pd.read_csv(mp, sep=None, engine='python', low_memory=False)
                    donor_col_src = pick_donor_col(list(mdf.columns))
                    donor_count_src = nunique_nonempty(mdf[donor_col_src]) if donor_col_src else None
                except Exception:
                    donor_col_src = ''
                    donor_count_src = None

        rows.append({
            'GSE': gse,
            'has_metadata': has_metadata,
            'status': 'ok',
            'cells': cells,
            'genes': genes,
            'donor_count_h5ad_matched': donor_count_h5ad,
            'donor_col_h5ad': donor_col_h5ad,
            'donor_count_source_metadata': donor_count_src,
            'donor_col_source_metadata': donor_col_src,
            'metadata_file': metadata_file,
        })
    except Exception as e:
        rows.append({
            'GSE': gse,
            'has_metadata': has_metadata,
            'status': f'error: {type(e).__name__}',
            'cells': None,
            'genes': None,
            'donor_count_h5ad_matched': None,
            'donor_col_h5ad': '',
            'donor_count_source_metadata': None,
            'donor_col_source_metadata': '',
            'metadata_file': metadata_file,
        })

out = pd.DataFrame(rows).sort_values('GSE')
OUT.parent.mkdir(parents=True, exist_ok=True)
out.to_csv(OUT, index=False)
print(f'Wrote {OUT} with {len(out)} rows')
print(out[['GSE','status','cells','genes','donor_count_h5ad_matched','donor_col_h5ad']].to_string(index=False))
