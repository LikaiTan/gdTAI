from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflows/gdtai"))

from gdtai_v4_gpu_core import (  # noqa: E402
    AtomicCandidateCheckpoint,
    GpuEnvironmentError,
    stage1_threshold_frontier,
    stage2_threshold_frontier,
    validate_gpu_environment,
    weighted_standardize,
)


class TestGpuCoreWithoutCuda(unittest.TestCase):
    def test_weighted_standardization(self) -> None:
        matrix = np.asarray([[0, 1], [2, 3], [4, 8]], dtype=np.float32)
        weight = np.asarray([1, 2, 1], dtype=np.float32)
        scaled, mean, scale = weighted_standardize(matrix, weight)
        normalized = weight / weight.sum()
        self.assertTrue(np.allclose((scaled * normalized[:, None]).sum(axis=0), 0, atol=1e-6))
        self.assertTrue(np.all(scale > 0))
        self.assertEqual(mean.shape, (2,))

    def test_environment_rejects_missing_contract_before_cuda_import(self) -> None:
        with self.assertRaises(GpuEnvironmentError):
            validate_gpu_environment({
                "cuda_visible_devices": "0", "cublas_workspace_config": ":4096:8",
                "pythonhashseed": "20260807", "mps_path_prefix": "/tmp/gdtai-v41-",
                "forbidden_mps_path": "/tmp/nvidia-mps",
            })

    def test_stage1_per_source_nk_guardrail(self) -> None:
        score = np.asarray([.99, .98, .97, .96, .99, .98, .97, .96, .8, .1, .7, .1])
        source = np.asarray(["a"] * 4 + ["b"] * 4 + ["nk_a"] * 2 + ["nk_b"] * 2, object)
        gdt = np.asarray([1, 0, 0, 0, 1, 0, 0, 0] + [0] * 4, bool)
        abt = np.asarray([0, 1, 1, 1, 0, 1, 1, 1] + [0] * 4, bool)
        nk = np.asarray([0] * 8 + [1] * 4, bool)
        frontier, decision = stage1_threshold_frontier(score, source, gdt, abt, nk, ["a", "b"], {
            "gdt_recall_per_source": 1.0, "abt_recall_per_source": 1.0,
            "maximum_nk_passage_per_source": .5,
        })
        self.assertTrue(decision.passed)
        self.assertEqual(len(frontier), len(np.unique(np.r_[score, 0, 1])))
        self.assertLessEqual(frontier.loc[frontier.eligible, "maximum_source_nk_passage"].max(), .5)

    def test_stage2_frontier_and_guardrail(self) -> None:
        score = np.asarray([.9, .8, .1, .05, .95, .85, .1, .02])
        labels = np.asarray([1, 1, 0, 0, 1, 1, 0, 0])
        sources = np.asarray(["a"] * 4 + ["b"] * 4, object)
        paired = labels == 0
        nk = np.asarray([.04, .03, .02, .01])
        frontier, decision = stage2_threshold_frontier(
            score, labels, sources, paired, nk, np.ones(8), np.ones(4), .5,
            np.zeros(8, bool), np.zeros(4, bool), "balanced",
            {"objective": "f1", "minimum_macro_recall": .8, "maximum_paired_abt_fpr": 0, "maximum_strict_nk_fpr": 0},
        )
        self.assertTrue(decision.passed)
        self.assertEqual(decision.macro_recall, 1.0)
        self.assertEqual(len(frontier), len(np.unique(np.r_[score, nk, 0, 1])))

    def test_checkpoint_refuses_contract_change(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            root = Path(value)
            checkpoint = AtomicCandidateCheckpoint(root, "candidate", {"split": "a", "seed": 1})
            expected = np.asarray([.1, .2, .3])
            checkpoint.save({"probability": expected}, {"candidate": "test"})
            loaded = checkpoint.load()
            self.assertTrue(np.array_equal(loaded["arrays"]["probability"], expected))
            changed = AtomicCandidateCheckpoint(root, "candidate", {"split": "b", "seed": 1})
            with self.assertRaises(RuntimeError):
                changed.load()


if __name__ == "__main__":
    unittest.main()
