#!/usr/bin/env python3
"""Build an xQuest-safe FASTA from proteins named in a GPx supplement workbook.

The default keeps the union of glycan-bearing Protein1 accessions and their
reported Protein2 partners. ``--glycoproteins-only`` keeps only Protein1, which
is smaller but cannot recover interactions with partner-only proteins.
"""

from __future__ import annotations

import argparse
import re
import sys
import zipfile
from pathlib import Path
from xml.etree import ElementTree


RESERVED = frozenset("XUBJ")
CONTAMINANT_HEADER = re.compile(
    r"(?:^>|\|)(?:rev(?:erse)?_|decoy_|con__|contaminant)", re.IGNORECASE
)
ACCESSION_TOKEN = re.compile(r"^(?:[^|]*\|)?([^|\s]+)")
SHEET_NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("supplement", type=Path, help="supplementary .xlsx workbook")
    parser.add_argument("output", type=Path, help="focused FASTA output")
    parser.add_argument(
        "--fasta",
        action="append",
        type=Path,
        required=True,
        help="source FASTA (repeat for additional databases)",
    )
    parser.add_argument(
        "--glycoproteins-only",
        action="store_true",
        help="keep Protein1 glycoproteins but omit Protein2-only partners",
    )
    parser.add_argument(
        "--exclude-accessions",
        type=Path,
        help="optional file with one explicitly excluded accession per line",
    )
    return parser.parse_args()


def qname(local: str) -> str:
    return f"{{{SHEET_NS}}}{local}"


def column_name(reference: str) -> str:
    match = re.match(r"([A-Z]+)", reference)
    if match is None:
        raise ValueError(f"invalid worksheet cell reference: {reference}")
    return match.group(1)


def shared_strings(archive: zipfile.ZipFile) -> list[str]:
    try:
        root = ElementTree.fromstring(archive.read("xl/sharedStrings.xml"))
    except KeyError:
        return []
    return [
        "".join(text.text or "" for text in item.iter(qname("t")))
        for item in root.findall(qname("si"))
    ]


def cell_value(cell: ElementTree.Element, strings: list[str]) -> str:
    cell_type = cell.attrib.get("t")
    value = cell.find(qname("v"))
    if cell_type == "inlineStr":
        return "".join(text.text or "" for text in cell.iter(qname("t")))
    if value is None or value.text is None:
        return ""
    if cell_type == "s":
        return strings[int(value.text)]
    return value.text


def reported_accessions(path: Path, glycoproteins_only: bool) -> tuple[list[str], int]:
    with zipfile.ZipFile(path) as archive:
        strings = shared_strings(archive)
        sheet = ElementTree.fromstring(archive.read("xl/worksheets/sheet1.xml"))

    rows: list[dict[str, str]] = []
    for row in sheet.iter(qname("row")):
        values = {
            column_name(cell.attrib["r"]): cell_value(cell, strings).strip()
            for cell in row.findall(qname("c"))
        }
        rows.append(values)
    if not rows or rows[0].get("E") != "Protein1 (ID)":
        raise ValueError("first worksheet does not contain Protein1 (ID) in column E")
    if not glycoproteins_only and rows[0].get("K") != "Protein2 (ID)":
        raise ValueError("first worksheet does not contain Protein2 (ID) in column K")

    protein1 = {row.get("E", "") for row in rows[1:] if row.get("E")}
    ordered: list[str] = []
    seen: set[str] = set()
    columns = ("E",) if glycoproteins_only else ("E", "K")
    for row in rows[1:]:
        for column in columns:
            value = row.get(column, "")
            if value and value not in seen:
                seen.add(value)
                ordered.append(value)
    return ordered, len(protein1)


def fasta_records(path: Path):
    header: str | None = None
    sequence: list[str] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(sequence).upper()
            header = line
            sequence = []
        elif header is None:
            raise ValueError(f"sequence before first FASTA header in {path}")
        else:
            sequence.append(line)
    if header is not None:
        yield header, "".join(sequence).upper()


def accession(header: str) -> str:
    match = ACCESSION_TOKEN.match(header[1:])
    if match is None:
        raise ValueError(f"cannot parse accession from header: {header}")
    return match.group(1)


def load_sources(paths: list[Path]) -> dict[str, tuple[str, str]]:
    records: dict[str, tuple[str, str]] = {}
    for path in paths:
        for header, sequence in fasta_records(path):
            key = accession(header)
            if key in records and records[key][1] != sequence:
                raise ValueError(f"conflicting sequences for accession {key}")
            records.setdefault(key, (header, sequence))
    return records


def exclusions(path: Path | None) -> set[str]:
    if path is None:
        return set()
    return {
        line.strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    }


def main() -> int:
    args = parse_args()
    try:
        wanted, glycoprotein_count = reported_accessions(
            args.supplement, args.glycoproteins_only
        )
        source = load_sources(args.fasta)
        excluded = exclusions(args.exclude_accessions)
        missing = [value for value in wanted if value not in source]
        if missing:
            raise ValueError(f"accession(s) missing from source FASTA: {', '.join(missing)}")

        kept: list[tuple[str, str, str]] = []
        dropped: list[tuple[str, str]] = []
        for value in wanted:
            header, sequence = source[value]
            if value in excluded:
                dropped.append((value, "explicit exclusion"))
                continue
            if CONTAMINANT_HEADER.search(header):
                dropped.append((value, "decoy/contaminant header"))
                continue
            reserved = "".join(sorted(set(sequence) & RESERVED))
            if reserved:
                dropped.append((value, f"xQuest-reserved residue(s) {reserved}"))
                continue
            kept.append((value, header, sequence))

        if not kept:
            raise ValueError("no compatible reported proteins remain")
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("w", encoding="utf-8", newline="\n") as output:
            for _, header, sequence in kept:
                output.write(f"{header}\n")
                for offset in range(0, len(sequence), 60):
                    output.write(f"{sequence[offset:offset + 60]}\n")
    except (OSError, ValueError, zipfile.BadZipFile) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    print(
        f"wrote {len(kept)} proteins to {args.output}: "
        f"{glycoprotein_count} reported glycoproteins plus "
        f"{len(kept) - glycoprotein_count} partner-only proteins"
    )
    print(f"discarded {len(dropped)} incompatible/explicit contaminant record(s)")
    for value, reason in dropped:
        print(f"dropped {value}: {reason}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
