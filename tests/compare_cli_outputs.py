#!/usr/bin/env python3
"""Check that two LinRacc builds produce byte-identical output.

Intended for performance-only changes: run a reference build (e.g. the
previous release) and a candidate build over a grid of CLI settings and
require identical output files.  Not part of ctest, since it needs a second
binary and a FASTA file.

  tests/compare_cli_outputs.py --reference old/LinRacc --candidate build/LinRacc \
      --fasta seqs.fa [--with-raccess] [--jobs 8]
"""

from __future__ import annotations

import argparse
import filecmp
import itertools
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records, name, seq = [], None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(seq)))
            name, seq = line[1:].split()[0], []
        elif line.strip():
            seq.append(line.strip())
    if name is not None:
        records.append((name, "".join(seq)))
    return records


def settings(with_raccess: bool, max_exact_len: int):
    access_lens = ["1,3,7,10", "5", "20,2,9"]
    beams = [0, 5, 100, 500]
    caps = ["30", "6", "len"]
    modes = [["-probabilities"], ["-byloop"], ["-probabilities", "-no-fast-logsumexp"]]
    energies = [["-energy=turner2004"], ["-energy=turner1999"]]
    if with_raccess:
        energies.append(["-engine=raccess"])
    for lens, beam, cap, mode, energy in itertools.product(access_lens, beams, caps, modes, energies):
        yield {"lens": lens, "beam": beam, "cap": cap, "flags": mode + energy,
               "max_len": max_exact_len if beam == 0 or cap == "len" else None}


def run(binary: Path, fasta: Path, out: Path, cfg: dict, seq_len: int) -> None:
    cap = str(seq_len) if cfg["cap"] == "len" else cfg["cap"]
    cmd = [str(binary), f"-seqfile={fasta}", f"-outfile={out}", f"-access_len={cfg['lens']}",
           f"-beam={cfg['beam']}", f"-c_hairpin={cap}", f"-c_multi={cap}", *cfg["flags"]]
    subprocess.run(cmd, check=True, capture_output=True)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--with-raccess", action="store_true")
    parser.add_argument("--max-exact-len", type=int, default=400,
                        help="skip beam=0 / uncapped runs on longer sequences")
    parser.add_argument("--jobs", type=int, default=8)
    args = parser.parse_args()

    records = read_fasta(args.fasta)
    with tempfile.TemporaryDirectory(prefix="linracc_cmp_") as raw:
        tmp = Path(raw)
        tasks = []
        for i, (name, seq) in enumerate(records):
            fa = tmp / f"s{i}.fa"
            fa.write_text(f">{name}\n{seq}\n")
            for j, cfg in enumerate(settings(args.with_raccess, args.max_exact_len)):
                if cfg["max_len"] is not None and len(seq) > cfg["max_len"]:
                    continue
                tasks.append((i, j, fa, cfg, len(seq)))

        def check(task):
            i, j, fa, cfg, n = task
            ref, cand = tmp / f"r{i}_{j}.out", tmp / f"c{i}_{j}.out"
            run(args.reference, fa, ref, cfg, n)
            run(args.candidate, fa, cand, cfg, n)
            same = filecmp.cmp(ref, cand, shallow=False)
            ref.unlink()
            cand.unlink()
            return same, records[i][0], cfg

        with ThreadPoolExecutor(max_workers=args.jobs) as ex:
            results = list(ex.map(check, tasks))

    failures = [(name, cfg) for same, name, cfg in results if not same]
    for name, cfg in failures[:20]:
        print(f"DIFF {name} {cfg}")
    print(f"compared {len(results)} runs, {len(failures)} differ")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
