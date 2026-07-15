#!/usr/bin/env python3
"""Convert centroided MGF spectra to the narrow mzXML form GlycoQuest consumes.

The converter is intentionally dependency-free and streams both input and output. It is
primarily for public repositories that publish MGF but not mzXML. Use ProteoWizard's
``msconvert`` when vendor RAW data and full acquisition metadata are available.
"""

from __future__ import annotations

import argparse
import base64
import re
import struct
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Iterator, TextIO
from xml.sax.saxutils import escape


SCAN_RE = re.compile(r"(?:\bscan=|\.)(\d+)(?:\D|$)", re.IGNORECASE)
CHARGE_RE = re.compile(r"(\d+)")


@dataclass
class Spectrum:
    title: str = ""
    scan: int | None = None
    rt_seconds: float | None = None
    precursor_mz: float | None = None
    charge: int | None = None
    peaks: list[tuple[float, float]] = field(default_factory=list)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="source MGF")
    parser.add_argument("output", type=Path, help="destination mzXML")
    parser.add_argument(
        "--scan-list",
        type=Path,
        help="optional text file containing one retained scan number per line",
    )
    parser.add_argument(
        "--require-ion",
        action="append",
        type=float,
        default=[],
        metavar="MZ",
        help="retain spectra matching any supplied fragment m/z (repeatable)",
    )
    parser.add_argument(
        "--ion-tolerance-ppm",
        type=float,
        default=20.0,
        help="fragment tolerance used by --require-ion (default: 20 ppm)",
    )
    parser.add_argument(
        "--max-scans",
        type=int,
        default=0,
        help="stop after writing this many retained scans (0 = unlimited)",
    )
    return parser.parse_args()


def load_scan_list(path: Path | None) -> set[int] | None:
    if path is None:
        return None
    scans: set[int] = set()
    for line_number, raw in enumerate(path.read_text().splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        try:
            scans.add(int(line))
        except ValueError as exc:
            raise ValueError(f"invalid scan number on {path}:{line_number}: {line}") from exc
    if not scans:
        raise ValueError(f"scan list is empty: {path}")
    return scans


def spectra(lines: Iterable[str]) -> Iterator[Spectrum]:
    current: Spectrum | None = None
    fallback_scan = 0
    for raw in lines:
        line = raw.strip()
        if line == "BEGIN IONS":
            if current is not None:
                raise ValueError("nested BEGIN IONS block")
            current = Spectrum()
            continue
        if line == "END IONS":
            if current is None:
                raise ValueError("END IONS without BEGIN IONS")
            fallback_scan += 1
            if current.scan is None:
                current.scan = fallback_scan
            yield current
            current = None
            continue
        if current is None or not line:
            continue

        if "=" in line:
            key, value = line.split("=", 1)
            key = key.upper()
            value = value.strip()
            if key == "TITLE":
                current.title = value
                native = re.search(r"\bscan=(\d+)", value, re.IGNORECASE)
                dotted = re.search(r"\.(\d+)\.\d+\.\d+(?:\s|$)", value)
                match = native or dotted or SCAN_RE.search(value)
                if match:
                    current.scan = int(match.group(1))
            elif key == "SCANS":
                current.scan = int(value.split()[0])
            elif key == "RTINSECONDS":
                current.rt_seconds = float(value.split()[0])
            elif key == "PEPMASS":
                current.precursor_mz = float(value.split()[0])
            elif key == "CHARGE":
                match = CHARGE_RE.search(value)
                current.charge = int(match.group(1)) if match else None
            continue

        fields = line.split()
        if len(fields) >= 2:
            try:
                current.peaks.append((float(fields[0]), float(fields[1])))
            except ValueError as exc:
                raise ValueError(f"invalid peak line: {line}") from exc

    if current is not None:
        raise ValueError("unterminated BEGIN IONS block")


def matches_required_ion(
    spectrum: Spectrum, targets: list[float], tolerance_ppm: float
) -> bool:
    if not targets:
        return True
    for peak_mz, _ in spectrum.peaks:
        for target in targets:
            if abs(peak_mz - target) <= target * tolerance_ppm * 1e-6:
                return True
    return False


def encode_peaks(peaks: list[tuple[float, float]]) -> str:
    packed = bytearray()
    for mz, intensity in peaks:
        packed.extend(struct.pack(">ff", mz, intensity))
    return base64.b64encode(packed).decode("ascii")


def write_scan(output: TextIO, spectrum: Spectrum) -> None:
    if spectrum.scan is None:
        raise ValueError("spectrum has no scan number")
    if spectrum.rt_seconds is None:
        raise ValueError(f"scan {spectrum.scan} has no RTINSECONDS")
    if spectrum.precursor_mz is None:
        raise ValueError(f"scan {spectrum.scan} has no PEPMASS")
    charge = (
        f' precursorCharge="{spectrum.charge}"' if spectrum.charge is not None else ""
    )
    title = escape(spectrum.title, {'"': "&quot;"})
    output.write(
        f'<scan num="{spectrum.scan}" msLevel="2" '
        f'retentionTime="PT{spectrum.rt_seconds:.6f}S" '
        f'peaksCount="{len(spectrum.peaks)}" filterLine="{title}">\n'
        f'<precursorMz{charge}>{spectrum.precursor_mz:.9f}</precursorMz>\n'
        '<peaks precision="32" byteOrder="network" contentType="m/z-int">'
        f"{encode_peaks(spectrum.peaks)}</peaks>\n"
        "</scan>\n"
    )


def convert(args: argparse.Namespace) -> int:
    if args.max_scans < 0:
        raise ValueError("--max-scans cannot be negative")
    if args.ion_tolerance_ppm <= 0:
        raise ValueError("--ion-tolerance-ppm must be positive")
    wanted_scans = load_scan_list(args.scan_list)
    args.output.parent.mkdir(parents=True, exist_ok=True)

    written = 0
    with args.input.open(encoding="utf-8", errors="strict") as source, args.output.open(
        "w", encoding="utf-8", newline="\n"
    ) as output:
        output.write(
            '<?xml version="1.0" encoding="UTF-8"?>\n'
            '<mzXML xmlns="http://sashimi.sourceforge.net/schema_revision/mzXML_3.1" '
            'version="3.1">\n<msRun>\n'
        )
        for spectrum in spectra(source):
            if wanted_scans is not None and spectrum.scan not in wanted_scans:
                continue
            if not matches_required_ion(
                spectrum, args.require_ion, args.ion_tolerance_ppm
            ):
                continue
            write_scan(output, spectrum)
            written += 1
            if args.max_scans and written >= args.max_scans:
                break
        output.write("</msRun>\n</mzXML>\n")

    if written == 0:
        args.output.unlink(missing_ok=True)
        raise ValueError("no spectra passed the requested filters")
    return written


def main() -> int:
    args = parse_args()
    try:
        count = convert(args)
    except (OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    print(f"wrote {count} MS2 scans to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
