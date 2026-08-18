#!/usr/bin/env python3
"""Recover GSE287541 10x TCR contigs and UMI/read support from public SRA runs.

The workflow processes a bounded number of libraries concurrently, writes only
compact Cell Ranger contig outputs to the project, and removes generated FASTQ,
SRA, BAM, and assembly scratch after each successful library. Cell Ranger's UI
is always disabled, so this workflow never binds a network port.
"""

from __future__ import annotations

# TNK_WORKFLOW_BOOTSTRAP
import sys as _tnk_sys
from pathlib import Path as _TnkPath

_TNK_PROJECT_ROOT = _TnkPath(__file__).resolve().parents[2]
for _tnk_path in (_TNK_PROJECT_ROOT, _TNK_PROJECT_ROOT / "src"):
    _tnk_value = str(_tnk_path)
    if _tnk_value not in _tnk_sys.path:
        _tnk_sys.path.insert(0, _tnk_value)

import argparse
import concurrent.futures
import gzip
import hashlib
import json
import re
import shutil
import subprocess
import time
from dataclasses import asdict, dataclass
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATASET_ROOT = PROJECT_ROOT / "data/datasets/GSE287541"
DEFAULT_MANIFEST = DATASET_ROOT / "raw/geo_sra/ena_run_manifest.tsv"
DEFAULT_OUTPUT = DATASET_ROOT / "interim/tcr_recovery/cellranger"
DEFAULT_REFERENCE = (
    PROJECT_ROOT
    / "data/references/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0"
)
DEFAULT_WORK = Path("/tmp/gse287541_tcr_recovery")
CELLRANGER = Path("/home/tanlikai/cellranger-9.0.1/bin/cellranger")
EXPECTED_READ_LENGTHS = (10, 10, 27)


@dataclass
class RecoveryResult:
    run_accession: str
    sample_title: str
    join_sample_id: str
    status: str
    elapsed_seconds: float
    n_filtered_contigs: int = 0
    n_productive_contigs: int = 0
    n_productive_cells: int = 0
    filtered_contig_sha256: str = ""
    error: str = ""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE)
    parser.add_argument("--work-root", type=Path, default=DEFAULT_WORK)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--cores-per-sample", type=int, default=30)
    parser.add_argument("--memory-gb-per-sample", type=int, default=80)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--run", nargs="+", help="Optional SRR accessions to recover.")
    parser.add_argument("--keep-work", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def normalize_key(value: object) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(value).strip().lower()).strip("_")


def join_sample_id(sample_title: str) -> str:
    match = re.fullmatch(r"(.+?)_(Rd[^_]+)_TCR", sample_title, flags=re.I)
    if match is None:
        raise ValueError(f"Unexpected GSE287541 TCR title: {sample_title}")
    participant_visit, round_id = match.groups()
    return normalize_key(f"{round_id}_{participant_visit}")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run_command(command: list[str], *, cwd: Path, log_handle) -> None:
    log_handle.write("$ " + " ".join(command) + "\n")
    log_handle.flush()
    subprocess.run(
        command,
        cwd=cwd,
        stdout=log_handle,
        stderr=subprocess.STDOUT,
        text=True,
        check=True,
    )


def fastq_read_length(path: Path) -> int:
    with path.open("rt", encoding="ascii") as handle:
        handle.readline()
        return len(handle.readline().strip())


def copy_as_gzip(source: Path, destination: Path) -> None:
    with source.open("rb") as src, destination.open("wb") as raw_destination:
        with gzip.GzipFile(
            filename="",
            mode="wb",
            compresslevel=6,
            fileobj=raw_destination,
            mtime=0,
        ) as dst:
            shutil.copyfileobj(src, dst, length=8 * 1024 * 1024)


def load_manifest(path: Path, runs: list[str] | None, limit: int | None) -> pd.DataFrame:
    manifest = pd.read_csv(path, sep="\t", dtype=str)
    manifest = manifest.loc[manifest["sample_title"].str.endswith("_TCR", na=False)].copy()
    manifest["join_sample_id"] = manifest["sample_title"].map(join_sample_id)
    if manifest["run_accession"].duplicated().any():
        raise ValueError("The run manifest contains duplicate TCR run accessions.")
    if manifest["join_sample_id"].duplicated().any():
        duplicates = manifest.loc[
            manifest["join_sample_id"].duplicated(False), "join_sample_id"
        ].tolist()
        raise ValueError(f"Multiple TCR runs collapse to one RNA join key: {duplicates}")
    if runs:
        missing = sorted(set(runs) - set(manifest["run_accession"]))
        if missing:
            raise ValueError(f"Requested runs are absent from the TCR manifest: {missing}")
        manifest = manifest.loc[manifest["run_accession"].isin(runs)]
    manifest = manifest.sort_values("run_accession").reset_index(drop=True)
    return manifest.head(limit) if limit else manifest


def existing_result(row: pd.Series, output_root: Path) -> RecoveryResult | None:
    output_dir = output_root / row["join_sample_id"]
    metadata_path = output_dir / "recovery_metadata.json"
    contig_path = output_dir / "filtered_contig_annotations.csv.gz"
    if not metadata_path.exists() or not contig_path.exists():
        return None
    metadata = json.loads(metadata_path.read_text())
    if metadata.get("run_accession") != row["run_accession"]:
        raise ValueError(f"Existing output has the wrong SRR accession: {output_dir}")
    if sha256_file(contig_path) != metadata.get("filtered_contig_sha256"):
        raise ValueError(f"Existing output checksum failed: {contig_path}")
    return RecoveryResult(
        run_accession=row["run_accession"],
        sample_title=row["sample_title"],
        join_sample_id=row["join_sample_id"],
        status="skipped_validated",
        elapsed_seconds=0.0,
        n_filtered_contigs=int(metadata["n_filtered_contigs"]),
        n_productive_contigs=int(metadata["n_productive_contigs"]),
        n_productive_cells=int(metadata["n_productive_cells"]),
        filtered_contig_sha256=metadata["filtered_contig_sha256"],
    )


def recover_one(row: pd.Series, args: argparse.Namespace) -> RecoveryResult:
    start = time.monotonic()
    run = row["run_accession"]
    title = row["sample_title"]
    sample_key = row["join_sample_id"]
    if not args.overwrite:
        existing = existing_result(row, args.output_root)
        if existing is not None:
            return existing

    work_dir = args.work_root / run
    final_dir = args.output_root / sample_key
    staging_dir = args.output_root / f".{sample_key}.{run}.tmp"
    if work_dir.exists():
        shutil.rmtree(work_dir)
    if staging_dir.exists():
        shutil.rmtree(staging_dir)
    work_dir.mkdir(parents=True)
    staging_dir.mkdir(parents=True)
    log_path = staging_dir / "recovery.log"

    try:
        with log_path.open("w", encoding="utf-8") as log:
            run_command(["prefetch", run, "--output-directory", str(work_dir)], cwd=work_dir, log_handle=log)
            sra_path = work_dir / run / f"{run}.sra"
            if not sra_path.exists():
                raise FileNotFoundError(f"prefetch did not create {sra_path}")

            fastq_dir = work_dir / "fastq"
            fastq_dir.mkdir()
            run_command(
                [
                    "fasterq-dump",
                    str(sra_path),
                    "--include-technical",
                    "--split-files",
                    "--threads",
                    str(args.cores_per_sample),
                    "--outdir",
                    str(fastq_dir),
                ],
                cwd=work_dir,
                log_handle=log,
            )
            split_fastqs = [fastq_dir / f"{run}_{number}.fastq" for number in range(1, 5)]
            if not all(path.exists() for path in split_fastqs):
                raise FileNotFoundError(f"Expected four split FASTQs for {run}")
            lengths = tuple(fastq_read_length(path) for path in split_fastqs)
            if lengths[:3] != EXPECTED_READ_LENGTHS or lengths[3] < 75:
                raise ValueError(f"Unexpected 10x VDJ read structure for {run}: {lengths}")

            roles = ("I1", "I2", "R1", "R2")
            renamed: list[Path] = []
            for source, role in zip(split_fastqs, roles):
                destination = fastq_dir / f"{run}_S1_L001_{role}_001.fastq"
                source.rename(destination)
                renamed.append(destination)
            run_command(
                ["pigz", "-p", str(min(args.cores_per_sample, 12)), *map(str, renamed)],
                cwd=work_dir,
                log_handle=log,
            )

            run_command(
                [
                    str(CELLRANGER),
                    "vdj",
                    "--id=vdj",
                    f"--reference={args.reference}",
                    f"--fastqs={fastq_dir}",
                    f"--sample={run}",
                    "--chain=TR",
                    f"--localcores={args.cores_per_sample}",
                    f"--localmem={args.memory_gb_per_sample}",
                    "--disable-ui",
                ],
                cwd=work_dir,
                log_handle=log,
            )

        outs = work_dir / "vdj/outs"
        filtered = pd.read_csv(outs / "filtered_contig_annotations.csv", low_memory=False)
        productive = filtered.loc[
            filtered["productive"].astype(str).str.lower().eq("true")
            & filtered["chain"].isin(["TRA", "TRB", "TRG", "TRD"])
        ]
        copy_as_gzip(
            outs / "filtered_contig_annotations.csv",
            staging_dir / "filtered_contig_annotations.csv.gz",
        )
        copy_as_gzip(
            outs / "all_contig_annotations.csv",
            staging_dir / "all_contig_annotations.csv.gz",
        )
        shutil.copy2(outs / "metrics_summary.csv", staging_dir / "metrics_summary.csv")
        checksum = sha256_file(staging_dir / "filtered_contig_annotations.csv.gz")
        metadata = {
            "run_accession": run,
            "sample_title": title,
            "join_sample_id": sample_key,
            "read_lengths": list(lengths),
            "n_filtered_contigs": len(filtered),
            "n_productive_contigs": len(productive),
            "n_productive_cells": productive["barcode"].nunique(),
            "filtered_contig_sha256": checksum,
            "cellranger": str(CELLRANGER),
            "cellranger_version": subprocess.check_output(
                [str(CELLRANGER), "--version"], text=True
            ).strip(),
            "cellranger_reference": str(args.reference),
            "cellranger_ui_disabled": True,
            "source_manifest": str(args.manifest),
        }
        (staging_dir / "recovery_metadata.json").write_text(
            json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
        )
        if final_dir.exists():
            shutil.rmtree(final_dir)
        staging_dir.replace(final_dir)
        if not args.keep_work:
            shutil.rmtree(work_dir)
        return RecoveryResult(
            run_accession=run,
            sample_title=title,
            join_sample_id=sample_key,
            status="recovered",
            elapsed_seconds=time.monotonic() - start,
            n_filtered_contigs=len(filtered),
            n_productive_contigs=len(productive),
            n_productive_cells=productive["barcode"].nunique(),
            filtered_contig_sha256=checksum,
        )
    except Exception as error:
        return RecoveryResult(
            run_accession=run,
            sample_title=title,
            join_sample_id=sample_key,
            status="failed",
            elapsed_seconds=time.monotonic() - start,
            error=f"{type(error).__name__}: {error}",
        )


def main() -> None:
    args = parse_args()
    if args.workers < 1 or args.cores_per_sample < 1 or args.memory_gb_per_sample < 16:
        raise ValueError("Invalid worker, core, or memory request.")
    if args.workers * args.cores_per_sample > 80:
        raise ValueError("Requested recovery exceeds the 80-core project envelope.")
    if args.workers * args.memory_gb_per_sample > 400:
        raise ValueError("Requested recovery exceeds the 400-GB project envelope.")
    for required in (args.manifest, args.reference, CELLRANGER):
        if not required.exists():
            raise FileNotFoundError(required)
    args.output_root.mkdir(parents=True, exist_ok=True)
    args.work_root.mkdir(parents=True, exist_ok=True)
    manifest = load_manifest(args.manifest, args.run, args.limit)
    print(f"Recovering {len(manifest)} GSE287541 TCR libraries with {args.workers} workers")

    results: list[RecoveryResult] = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(recover_one, row, args): row["run_accession"]
            for _, row in manifest.iterrows()
        }
        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            results.append(result)
            print(
                f"{result.run_accession} {result.join_sample_id}: {result.status} "
                f"({result.elapsed_seconds / 60:.1f} min) {result.error}",
                flush=True,
            )

    result_frame = pd.DataFrame(asdict(result) for result in results).sort_values("run_accession")
    result_path = args.output_root.parent / "recovery_run_summary.tsv"
    result_frame.to_csv(result_path, sep="\t", index=False)
    failed = result_frame["status"].eq("failed")
    full_manifest = load_manifest(args.manifest, None, None)
    validated_outputs = sum(
        existing_result(row, args.output_root) is not None
        for _, row in full_manifest.iterrows()
    )
    completion = {
        "expected_tcr_libraries": 46,
        "libraries_in_this_run": len(manifest),
        "libraries_with_validated_output": validated_outputs,
        "complete": bool(validated_outputs == 46),
        "cellranger_ui_disabled": True,
        "run_summary": str(result_path),
    }
    (args.output_root.parent / "recovery_completion.json").write_text(
        json.dumps(completion, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(completion, indent=2))
    if failed.any():
        raise SystemExit(1)


if __name__ == "__main__":
    main()
