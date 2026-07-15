#!/usr/bin/env python3
"""Extract FASTA records by accession while preserving source headers and sequences."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path


ACCESSION_TOKEN = re.compile(r"^(?:[^|]*\|)?([^|\s]+)")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="source FASTA")
    parser.add_argument("accessions", type=Path, help="one accession per line")
    parser.add_argument("output", type=Path, help="subset FASTA")
    return parser.parse_args()


def requested_accessions(path: Path) -> list[str]:
    accessions = [
        line.strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    if not accessions:
        raise ValueError(f"accession list is empty: {path}")
    duplicates = sorted({value for value in accessions if accessions.count(value) > 1})
    if duplicates:
        raise ValueError(f"duplicate accession(s): {', '.join(duplicates)}")
    return accessions


def records(text: str):
    header: str | None = None
    sequence: list[str] = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(sequence)
            header = line
            sequence = []
        elif header is None:
            raise ValueError("sequence data appears before the first FASTA header")
        else:
            sequence.append(line)
    if header is not None:
        yield header, "".join(sequence)


def accession(header: str) -> str:
    match = ACCESSION_TOKEN.match(header[1:])
    if match is None:
        raise ValueError(f"cannot parse accession from header: {header}")
    return match.group(1)


def main() -> int:
    args = parse_args()
    try:
        requested = requested_accessions(args.accessions)
        wanted = set(requested)
        selected = {
            accession(header): (header, sequence)
            for header, sequence in records(args.input.read_text(encoding="utf-8"))
            if accession(header) in wanted
        }
        missing = [value for value in requested if value not in selected]
        if missing:
            raise ValueError(f"accession(s) not found: {', '.join(missing)}")

        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("w", encoding="utf-8", newline="\n") as output:
            for value in requested:
                header, sequence = selected[value]
                output.write(f"{header}\n")
                for offset in range(0, len(sequence), 60):
                    output.write(f"{sequence[offset:offset + 60]}\n")
    except (OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    print(f"wrote {len(requested)} records to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
