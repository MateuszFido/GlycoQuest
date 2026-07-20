// Copyright (c) ETH Zurich, Mateusz Fido

export interface SpectrumSummaryInput {
  scan: number;
  precursorMz: number | null | undefined;
  charge: number | null | undefined;
  scanTimeMin: number | null | undefined;
}

/** Scan metadata only — retention/scan time never goes on the plot axes. */
export function formatSpectrumSummary(input: SpectrumSummaryInput): string {
  const parts = [`Scan ${input.scan}`];
  if (input.precursorMz != null && Number.isFinite(input.precursorMz) && input.precursorMz > 0) {
    parts.push(`precursor ${input.precursorMz.toFixed(4)} m/z`);
  }
  if (input.charge != null && Number.isFinite(input.charge) && input.charge > 0) {
    parts.push(`${input.charge}+`);
  }
  if (input.scanTimeMin != null && Number.isFinite(input.scanTimeMin)) {
    parts.push(`scan_time ${input.scanTimeMin.toFixed(2)} min`);
  }
  return parts.join(' | ');
}
