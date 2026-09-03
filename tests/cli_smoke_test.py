#!/usr/bin/env python3
"""Dependency-free smoke tests for the standalone LinearRaccess CLI."""

from __future__ import annotations

import argparse
import math
import subprocess
import tempfile
from pathlib import Path


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, text=True, capture_output=True, check=False)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--with-raccess", action="store_true")
    args = parser.parse_args()
    executable = args.executable.resolve()
    require(executable.is_file(), f"missing executable: {executable}")

    with tempfile.TemporaryDirectory(prefix="linearraccess_cli_") as raw_tmp:
        tmp = Path(raw_tmp)
        valid = tmp / "valid.fa"
        valid.write_text(">test\nGCGAAACGC\n", encoding="ascii")

        output = tmp / "probabilities.txt"
        completed = run(
            [
                str(executable),
                f"-seqfile={valid}",
                f"-outfile={output}",
                "-access_len=1,3",
                "-beam=0",
                "-c_multi=9",
                "-c_hairpin=9",
                "-no-fast-logsumexp",
                "-no-normalize",
                "-probabilities",
                "-source-revision=test-revision",
            ]
        )
        require(completed.returncode == 0, completed.stdout + completed.stderr)
        lines = output.read_text(encoding="ascii").splitlines()
        expected_metadata = {
            "# linearraccess_metadata=1",
            "# software_version=0.1.0",
            "# source_revision=test-revision",
            "# engine=lincapr",
            "# requested_energy_model=turner2004",
            "# effective_energy_model=turner2004",
            "# beam=0",
            "# c_multi=9",
            "# c_hairpin=9",
            "# logsumexp=exact",
            "# normalize_profiles=false",
            "# length_factor=0",
            "# output_mode=probability",
            "# value=probability",
        }
        missing = expected_metadata.difference(lines)
        require(not missing, f"missing metadata: {sorted(missing)}")
        values = [
            float(item.split(",", 1)[1])
            for line in lines
            if line and not line.startswith((">", "#"))
            for item in line.split()[1].split(";")
            if item
        ]
        require(len(values) == 16, f"unexpected probability count: {len(values)}")
        require(all(math.isfinite(value) for value in values), "non-finite probability")

        duplicate = run(
            [
                str(executable),
                f"-seqfile={valid}",
                f"-outfile={tmp / 'duplicate.txt'}",
                "-access_len=1,1",
            ]
        )
        require(duplicate.returncode != 0, "duplicate access lengths were accepted")

        optional_backend = run(
            [
                str(executable),
                f"-seqfile={valid}",
                f"-outfile={tmp / 'raccess.txt'}",
                "-access_len=1",
                "-engine=raccess",
            ]
        )
        if args.with_raccess:
            require(
                optional_backend.returncode == 0,
                optional_backend.stdout + optional_backend.stderr,
            )
        else:
            require(
                optional_backend.returncode != 0,
                "unavailable Raccess backend was accepted",
            )

        malformed = tmp / "malformed.fa"
        malformed.write_text("GCGAAACGC\n", encoding="ascii")
        malformed_run = run(
            [
                str(executable),
                f"-seqfile={malformed}",
                f"-outfile={tmp / 'malformed.txt'}",
                "-access_len=1",
            ]
        )
        require(malformed_run.returncode != 0, "malformed FASTA was accepted")

    print("cli_smoke_test\tpass")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
