#!/usr/bin/env python3
"""Remove FASTA records that contain xQuest's reserved pseudo-residue letters."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


RESERVED = frozenset("XUBJ")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    return parser.parse_args()


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
            sequence.append(line.upper())
    if header is not None:
        yield header, "".join(sequence)


def main() -> int:
    args = parse_args()
    try:
        source = args.input.read_text(encoding="utf-8")
        kept: list[tuple[str, str]] = []
        dropped: list[tuple[str, str]] = []
        for header, sequence in records(source):
            found = "".join(sorted(set(sequence) & RESERVED))
            if found:
                dropped.append((header, found))
            else:
                kept.append((header, sequence))
        if not kept:
            raise ValueError("no xQuest-compatible FASTA records remain")
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("w", encoding="utf-8", newline="\n") as output:
            for header, sequence in kept:
                output.write(f"{header}\n")
                for offset in range(0, len(sequence), 60):
                    output.write(f"{sequence[offset:offset + 60]}\n")
    except (OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    print(
        f"wrote {len(kept)} records to {args.output}; "
        f"dropped {len(dropped)} record(s) containing {''.join(sorted(RESERVED))}"
    )
    for header, found in dropped:
        print(f"dropped ({found}): {header}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
