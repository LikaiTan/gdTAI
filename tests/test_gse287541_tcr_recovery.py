from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "workflows/intake/recover_gse287541_tcr_from_sra.py"
SPEC = importlib.util.spec_from_file_location("recover_gse287541_tcr_from_sra", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_public_manifest_has_46_unique_tcr_libraries() -> None:
    manifest = MODULE.load_manifest(MODULE.DEFAULT_MANIFEST, None, None)
    assert len(manifest) == 46
    assert manifest["run_accession"].is_unique
    assert manifest["join_sample_id"].is_unique


def test_round_is_part_of_duplicate_participant_join_key() -> None:
    assert MODULE.join_sample_id("N131V01_Rd4a_TCR") == "rd4a_n131v01"
    assert MODULE.join_sample_id("N131V01_Rd5_TCR") == "rd5_n131v01"
