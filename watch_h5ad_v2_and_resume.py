#!/usr/bin/env python3
"""Watch for h5ad_v2.csv and resume the remaining rescue-enabled pipeline work.

This monitor polls every 120 seconds. Once `h5ad_v2.csv` appears and parses
cleanly, it reruns:
1. Phase 0 audit against h5ad_v2.csv
2. Phase 1 coarse T/NK extraction using the updated canonical audit table

The script is intentionally narrow: it does not advance to Phase 1b or later.
It only refreshes Phase 0 and Phase 1 after the repaired registry is available.
"""

from __future__ import annotations

import json
import os
import signal
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


# Config
PROJECT_ROOT = Path(__file__).resolve().parent
PYTHON_BIN = Path("/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python")
REGISTRY_V2 = PROJECT_ROOT / "h5ad_v2.csv"
LOG_DIR = PROJECT_ROOT / "Integrated_dataset" / "logs"
WATCH_LOG = LOG_DIR / "watch_h5ad_v2_and_resume.log"
STATE_JSON = LOG_DIR / "watch_h5ad_v2_and_resume_state.json"
PID_FILE = LOG_DIR / "watch_h5ad_v2_and_resume.pid"
POLL_SECONDS = 120

PHASE0_SCRIPT = PROJECT_ROOT / "phase0_dataset_audit.py"
PHASE1_SCRIPT = PROJECT_ROOT / "phase1_extract_tnk_candidates.py"
REQUIRED_REGISTRY_COLUMNS = {"h5ad_path", "gse_id", "source_root"}


def utc_now() -> str:
    """Return a stable UTC timestamp string for logs and state."""
    return datetime.now(timezone.utc).isoformat()


def append_log(message: str) -> None:
    """Append one timestamped line to the watcher log and stdout."""
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    line = f"[{utc_now()}] {message}"
    with WATCH_LOG.open("a", encoding="utf-8") as handle:
        handle.write(line + "\n")
    print(line, flush=True)


def write_state(status: str, **extra: object) -> None:
    """Persist watcher state for later inspection."""
    payload = {"timestamp_utc": utc_now(), "status": status}
    payload.update(extra)
    STATE_JSON.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def pid_is_running(pid: int) -> bool:
    """Return True if the given PID appears to be alive."""
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def ensure_single_instance() -> None:
    """Refuse to start a second watcher if one is already alive."""
    if PID_FILE.exists():
        content = PID_FILE.read_text(encoding="utf-8").strip()
        if content:
            try:
                existing_pid = int(content)
            except ValueError:
                existing_pid = -1
            if existing_pid > 0 and pid_is_running(existing_pid):
                raise RuntimeError(f"Watcher already running with PID {existing_pid}")

    PID_FILE.write_text(str(os.getpid()), encoding="utf-8")


def cleanup_pid_file() -> None:
    """Remove the PID file when the watcher exits."""
    if PID_FILE.exists():
        content = PID_FILE.read_text(encoding="utf-8").strip()
        if content == str(os.getpid()):
            PID_FILE.unlink()


def install_signal_handlers() -> None:
    """Handle termination cleanly so the PID file is not left stale."""
    def handle_signal(signum, _frame) -> None:
        append_log(f"Received signal {signum}; stopping watcher.")
        write_state("stopped", signal=signum)
        cleanup_pid_file()
        sys.exit(0)

    signal.signal(signal.SIGTERM, handle_signal)
    signal.signal(signal.SIGINT, handle_signal)


def registry_is_ready(registry_path: Path) -> tuple[bool, str]:
    """Check whether the repaired registry exists and has the required columns."""
    if not registry_path.exists():
        return False, "registry_missing"
    if registry_path.stat().st_size == 0:
        return False, "registry_empty"

    try:
        registry = pd.read_csv(registry_path)
    except Exception as exc:
        return False, f"registry_unreadable:{type(exc).__name__}:{exc}"

    missing = REQUIRED_REGISTRY_COLUMNS - set(registry.columns)
    if missing:
        return False, f"registry_missing_columns:{sorted(missing)}"
    if registry.empty:
        return False, "registry_has_no_rows"
    return True, f"registry_rows:{len(registry)}"


def run_step(label: str, cmd: list[str]) -> None:
    """Run one pipeline step and stream output into the watcher log."""
    append_log(f"Starting {label}: {' '.join(cmd)}")
    with WATCH_LOG.open("a", encoding="utf-8") as handle:
        handle.write(f"\n[{utc_now()}] >>> {label}\n")
        handle.flush()
        result = subprocess.run(
            cmd,
            cwd=PROJECT_ROOT,
            stdout=handle,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )
    if result.returncode != 0:
        raise RuntimeError(f"{label} failed with exit code {result.returncode}")
    append_log(f"Completed {label}.")


def run_resume_pipeline() -> None:
    """Run the approved automatic continuation once h5ad_v2.csv appears."""
    run_step(
        "phase0_audit_h5ad_v2",
        [str(PYTHON_BIN), str(PHASE0_SCRIPT), "--registry", str(REGISTRY_V2)],
    )
    run_step(
        "phase1_rebuild_candidates_from_updated_audit",
        [str(PYTHON_BIN), str(PHASE1_SCRIPT)],
    )


def main() -> None:
    """Poll for h5ad_v2.csv and continue the pipeline when it becomes available."""
    ensure_single_instance()
    install_signal_handlers()
    append_log("Watcher started. Polling every 120 seconds for h5ad_v2.csv.")
    write_state("waiting", poll_seconds=POLL_SECONDS, registry=str(REGISTRY_V2))

    try:
        while True:
            ready, detail = registry_is_ready(REGISTRY_V2)
            if not ready:
                write_state("waiting", poll_seconds=POLL_SECONDS, registry=str(REGISTRY_V2), detail=detail)
                append_log(f"Registry not ready: {detail}")
                time.sleep(POLL_SECONDS)
                continue

            append_log(f"Detected ready registry: {detail}")
            write_state("running", registry=str(REGISTRY_V2), detail=detail)
            run_resume_pipeline()
            write_state("completed", registry=str(REGISTRY_V2), detail=detail)
            append_log("Automatic continuation finished successfully.")
            return
    except Exception as exc:
        append_log(f"Watcher failed: {type(exc).__name__}: {exc}")
        write_state("failed", error=f"{type(exc).__name__}: {exc}", registry=str(REGISTRY_V2))
        raise
    finally:
        cleanup_pid_file()


if __name__ == "__main__":
    main()
