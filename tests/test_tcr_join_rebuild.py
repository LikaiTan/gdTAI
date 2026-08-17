from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "workflows/intake/rebuild_flagged_tcr_joins.py"
SPEC = importlib.util.spec_from_file_location("rebuild_flagged_tcr_joins", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_barcode_core_handles_prefixed_and_decorated_barcodes() -> None:
    expected = "AAACCTGAGCGTGAGT"
    assert MODULE.barcode_core("AAACCTGAGCGTGAGT-1-0") == expected
    assert MODULE.barcode_core("P53-post-P-CD8-AAACCTGAGCGTGAGT-1") == expected
    assert MODULE.barcode_core("Rd1_N014V01_AAACCTGAGCGTGAGT-1") == expected


def test_sample_key_rules_align_rna_and_raw_modalities() -> None:
    obs = pd.DataFrame({"sample": ["gsm5399838_c14"]})
    assert MODULE.rna_sample_keys("GSE178882", obs).iloc[0] == "c14"
    assert MODULE.raw_sample_key_from_name(
        "GSE178882", "GSM5399868_c14_filtered_contig_annotations.csv.gz"
    ) == "c14"

    obs = pd.DataFrame({"sample": ["gsm8035466_om-hl-007-c8"]})
    assert MODULE.rna_sample_keys("GSE254176", obs).iloc[0] == "om_hl_007_c8"
    assert MODULE.raw_sample_key_from_name(
        "GSE254176", "GSM8035485_OM-HL-007-C8-TCR_filtered_contig_annotations.csv.gz"
    ) == "om_hl_007_c8"
    assert MODULE.raw_sample_key_from_name(
        "GSE311112", "GSM9317194_cll_5_btk_clone_tcr_all_contig_annotations.csv.gz"
    ) == "5_btk_clone"


def test_best_contig_uses_umi_then_reads_and_keeps_all_chain_families() -> None:
    rows = pd.DataFrame(
        {
            "sample_id": ["s1", "s1", "s1"],
            "barcode_core": ["A" * 16] * 3,
            "chain": ["TRD", "TRD", "TRG"],
            "cdr3": ["low", "high", "gamma"],
            "v": ["TRDV1", "TRDV2", "TRGV9"],
            "d": ["", "", ""],
            "j": ["TRDJ1", "TRDJ1", "TRGJ1"],
            "cdr3_nt": ["a", "b", "c"],
            "clone_id": ["1", "2", "3"],
            "contig_id": ["a", "b", "c"],
            "umis": [1, 4, 2],
            "reads": [100, 20, 30],
            "umi_available": [True] * 3,
            "read_available": [True] * 3,
            "full_length": [True] * 3,
            "source_file": ["raw"] * 3,
        }
    )
    selected = MODULE.select_best_contigs(rows)
    assert set(selected["chain"]) == {"TRD", "TRG"}
    assert selected.loc[selected["chain"].eq("TRD"), "cdr3"].iloc[0] == "high"
    assert selected.loc[selected["chain"].eq("TRD"), "n_productive_contigs"].iloc[0] == 2


def test_wide_tables_keep_unavailable_umi_as_null() -> None:
    frame = MODULE.empty_contig_frame()
    assert "umis" in frame.columns
    sidecar = pd.DataFrame(index=["cell"])
    sidecar = MODULE.add_empty_chain_columns(sidecar)
    assert pd.isna(sidecar.loc["cell", "TRA_umis"])


def test_truthy_handles_an_empty_boolean_series() -> None:
    result = MODULE.truthy(pd.Series([], dtype=bool))
    assert result.empty


def test_unresolved_donor_sources_do_not_report_cross_donor_reuse() -> None:
    sidecar = pd.DataFrame(
        {
            "has_TRA_TRB_paired_rebuilt": [True, True],
            "join_sample_id": ["s1", "s2"],
            "join_donor_id": ["d1", "d2"],
            "TRA_cdr3": ["TRA", "TRA"],
            "TRB_cdr3": ["TRB", "TRB"],
        }
    )
    result = MODULE.receptor_reuse_metrics("GSE178882", sidecar)
    assert pd.isna(result["n_pairs_across_donors"])
    assert result["n_pairs_across_samples"] == 1
