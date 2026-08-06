import importlib.util
from pathlib import Path
import sys
import unittest

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "workflows/intake/audit_extension_tcr_source_chains.py"
SPEC = importlib.util.spec_from_file_location("audit_extension_tcr_source_chains", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class ExtensionTcrSourceAuditTests(unittest.TestCase):
    def test_only_productive_high_confidence_cell_contigs_with_cdr3_are_eligible(self) -> None:
        frame = pd.DataFrame(
            {
                "chain": ["TRD", "TRD", "TRG", "TRA", "TRB"],
                "productive": [True, False, True, True, True],
                "is_cell": [True, True, True, True, True],
                "high_confidence": [True, True, False, True, True],
                "cdr3": ["CACD", "", "CALW", "CAVA", "CASB"],
            }
        )
        result = MODULE.summarize_contig_frame(frame)
        self.assertEqual(result["total"]["TRD"], 2)
        self.assertEqual(result["productive_flag"]["TRD"], 1)
        self.assertEqual(result["eligible_productive_cdr3"]["TRD"], 1)
        self.assertEqual(result["eligible_productive_cdr3"]["TRG"], 0)
        self.assertEqual(result["eligible_productive_cdr3"]["TRA"], 1)
        self.assertEqual(result["eligible_productive_cdr3"]["TRB"], 1)

    def test_locus_and_junction_aliases_are_supported(self) -> None:
        frame = pd.DataFrame(
            {
                "locus": ["TRG", "TRD"],
                "is_productive": ["TRUE", "TRUE"],
                "junction_aa": ["CALX", ""],
            }
        )
        result = MODULE.summarize_contig_frame(frame)
        self.assertEqual(result["eligible_productive_cdr3"]["TRG"], 1)
        self.assertEqual(result["eligible_productive_cdr3"]["TRD"], 0)


if __name__ == "__main__":
    unittest.main()
