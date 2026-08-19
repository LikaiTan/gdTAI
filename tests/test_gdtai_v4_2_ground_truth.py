from __future__ import annotations

import numpy as np
import pandas as pd

from workflows.gdtai.run_gdtai_v4_2_ground_truth import qualifying_pair, source_umi_status


def test_source_umi_status_distinguishes_available_unavailable_and_mixed() -> None:
    source = np.array(["A", "A", "B", "B", "C", "C", "D"], dtype=object)
    paired = np.array([1, 1, 1, 1, 1, 1, 0], dtype=bool)
    left = np.array([1, 1, 0, 0, 1, 0, 0], dtype=bool)
    right = np.array([1, 1, 0, 0, 1, 0, 0], dtype=bool)
    result = source_umi_status(source, paired, left, right).set_index("source_gse_id")
    assert result.loc["A", "umi_status"] == "AVAILABLE"
    assert result.loc["B", "umi_status"] == "UNAVAILABLE"
    assert result.loc["C", "umi_status"] == "MIXED"
    assert result.loc["D", "umi_status"] == "JOIN_UNRESOLVED"


def test_qualifying_pair_excludes_umi_one_but_accepts_true_unavailable() -> None:
    source = np.array(["A", "A", "B", "C", "C"], dtype=object)
    paired = np.ones(5, dtype=bool)
    left_available = np.array([1, 1, 0, 1, 0], dtype=bool)
    right_available = np.array([1, 1, 0, 1, 0], dtype=bool)
    left_umi = np.array([2, 1, np.nan, 3, np.nan])
    right_umi = np.array([2, 4, np.nan, 2, np.nan])
    status = pd.DataFrame(
        {"source_gse_id": ["A", "B", "C"], "umi_status": ["AVAILABLE", "UNAVAILABLE", "MIXED"]}
    )
    observed = qualifying_pair(
        source, paired, left_available, right_available, left_umi, right_umi, status, 2
    )
    assert observed.tolist() == [True, False, True, True, False]
