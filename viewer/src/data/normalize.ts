// Copyright (c) ETH Zurich, Mateusz Fido

import type {
  ViewerBundle,
  ViewerCrosslink,
  ViewerFiltering,
  ViewerIsotopePair,
  ViewerMeta,
  ViewerSpectrum,
} from '../types';

type RawBundle = Omit<ViewerBundle, 'filtering'> & {
  filtering?: Record<string, ViewerFiltering>;
  mirror_fragments?: unknown;
  fragments?: unknown;
  meta: Partial<ViewerMeta> & {
    generated_at?: string;
    generated_at_iso?: string;
    generated_at_unix?: number | null;
  };
  spectra: Record<string, Partial<ViewerSpectrum> & { mz: number[]; intensity: number[] }>;
  isotope_pairs?: Record<string, Partial<ViewerIsotopePair>>;
  crosslinks: Array<Partial<ViewerCrosslink> & Pick<ViewerCrosslink, 'id' | 'protein1' | 'protein2'>>;
};

export function normalizeViewerBundle(raw: unknown): ViewerBundle {
  const bundle = raw as RawBundle;
  const meta = normalizeMeta(bundle.meta);
  const spectra = normalizeSpectra(bundle.spectra ?? {});
  const rawCrosslinks = (bundle.crosslinks ?? []).map((xl) => normalizeCrosslink(xl, spectra));
  const crosslinks = uniquifyCrosslinkIds(rawCrosslinks);

  return {
    viewer_schema_version: 4,
    meta,
    proteins: bundle.proteins ?? [],
    crosslinks,
    qc: bundle.qc,
    spectra,
    isotope_pairs: normalizeIsotopePairs(bundle.isotope_pairs ?? {}),
    filtering: normalizeFiltering(bundle.filtering ?? {}, crosslinks),
  };
}

export function proteinPairKey(left: string, right: string): string {
  if (!left) return right;
  if (!right) return left;
  return left <= right ? `${left}|${right}` : `${right}|${left}`;
}

function uniquifyCrosslinkIds(crosslinks: ViewerCrosslink[]): ViewerCrosslink[] {
  const seen = new Map<string, number>();
  return crosslinks.map((xl, index) => {
    const baseId = xl.id || `xl_${index}`;
    const count = seen.get(baseId) ?? 0;
    seen.set(baseId, count + 1);
    if (count === 0 && xl.id) return xl;
    return {
      ...xl,
      id: `${baseId}__row_${index}`,
      source_id: baseId,
    };
  });
}

function normalizeMeta(meta: RawBundle['meta']): ViewerMeta {
  const unix = normalizeUnixTimestamp(meta.generated_at_unix, meta.generated_at);
  const iso = meta.generated_at_iso ?? (unix == null ? meta.generated_at ?? '' : unixToIso(unix));
  return {
    project: meta.project ?? '',
    input_label: meta.input_label ?? '',
    crosslinker: meta.crosslinker ?? '',
    crosslinker_mw: finiteOrNull(meta.crosslinker_mw),
    xlink_sites: meta.xlink_sites ?? '',
    glycan_library: meta.glycan_library ?? '',
    xquest_version: meta.xquest_version ?? null,
    generated_at: meta.generated_at ?? iso,
    generated_at_iso: iso,
    generated_at_unix: unix,
    total_hits: meta.total_hits ?? 0,
    passing_hits: meta.passing_hits ?? 0,
  };
}

function normalizeSpectra(rawSpectra: RawBundle['spectra']): Record<string, ViewerSpectrum> {
  const spectra: Record<string, ViewerSpectrum> = {};
  for (const [scan, spectrum] of Object.entries(rawSpectra)) {
    spectra[scan] = {
      mz: spectrum.mz ?? [],
      intensity: spectrum.intensity ?? [],
      retention_time_min: finiteOrNull(spectrum.retention_time_min),
      precursor_mz: spectrum.precursor_mz ?? 0,
      charge: spectrum.charge ?? 0,
    };
  }
  return spectra;
}

function normalizeIsotopePairs(
  rawPairs: Record<string, Partial<ViewerIsotopePair>>,
): Record<string, ViewerIsotopePair> {
  const pairs: Record<string, ViewerIsotopePair> = {};
  for (const [scan, pair] of Object.entries(rawPairs)) {
    const lightScan = finiteNumber(pair.light_scan) ?? 0;
    const heavyScan = finiteNumber(pair.heavy_scan) ?? 0;
    pairs[scan] = {
      id: pair.id ?? `${lightScan}:${heavyScan}`,
      source_artifact: pair.source_artifact ?? 'isotope_pairs.tsv',
      source_row: finiteNumber(pair.source_row) ?? null,
      light_file: pair.light_file ?? null,
      heavy_file: pair.heavy_file ?? null,
      light_scan: lightScan,
      heavy_scan: heavyScan,
      rt_light_min: finiteOrNull(pair.rt_light_min),
      rt_heavy_min: finiteOrNull(pair.rt_heavy_min),
      mz_light: finiteNumber(pair.mz_light) ?? 0,
      mz_heavy: finiteNumber(pair.mz_heavy) ?? 0,
      light_charge: finiteNumber(pair.light_charge) ?? 0,
      heavy_charge: finiteNumber(pair.heavy_charge) ?? 0,
    };
  }
  return pairs;
}

function normalizeCrosslink(
  xl: RawBundle['crosslinks'][number],
  spectra: Record<string, ViewerSpectrum>,
): ViewerCrosslink {
  const scanKey = xl.scan == null ? null : String(xl.scan);
  const spectrumRt = scanKey == null ? null : spectra[scanKey]?.retention_time_min ?? null;
  const linkType = normalizeLinkType(xl.link_type, xl.pep_seq2, xl.protein2);
  return {
    id: xl.id,
    source_id: xl.source_id,
    link_type: linkType,
    protein1: xl.protein1,
    pep_pos1: xl.pep_pos1 ?? null,
    pep_seq1: xl.pep_seq1 ?? '',
    link_pos1: xl.link_pos1 ?? null,
    abs_pos1: xl.abs_pos1 ?? null,
    protein2: xl.protein2,
    pep_pos2: xl.pep_pos2 ?? null,
    pep_seq2: xl.pep_seq2 ?? '',
    link_pos2: xl.link_pos2 ?? null,
    abs_pos2: xl.abs_pos2 ?? null,
    score: xl.score ?? 0,
    soft_score: xl.soft_score ?? 0,
    scan: xl.scan ?? null,
    retention_time_min: finiteOrNull(xl.retention_time_min) ?? spectrumRt,
    source_file: xl.source_file ?? null,
    charge: xl.charge ?? 0,
    precursor_mz: xl.precursor_mz ?? 0,
    precursor_error_ppm: xl.precursor_error_ppm ?? 0,
    xlinker_mass: finiteOrNull(xl.xlinker_mass),
    topology: xl.topology ?? '',
    protein_pair_key: xl.protein_pair_key ?? proteinPairKey(xl.protein1, xl.protein2),
    glycan_name: xl.glycan_name ?? null,
    glycan_composition: xl.glycan_composition ?? null,
    glyco_residue: xl.glyco_residue ?? null,
    glyco_peptide: xl.glyco_peptide ?? null,
    glyco_sites: (xl.glyco_sites ?? []).map((site) => ({
      peptide: finiteNumber(site.peptide) ?? 0,
      peptide_position: finiteNumber(site.peptide_position) ?? 0,
      residue: site.residue ?? '',
      sequon_present: site.sequon_present ?? null,
      plausible: site.plausible ?? false,
    })),
    diagnostic_ions: (xl.diagnostic_ions ?? []).map(normalizeDiagnosticIon),
    loss_label: xl.loss_label ?? null,
    postfilter_status: xl.postfilter_status ?? 'pass',
    mapped: xl.mapped ?? false,
  };
}

function normalizeDiagnosticIon(ion: Partial<ViewerCrosslink['diagnostic_ions'][number]>): ViewerCrosslink['diagnostic_ions'][number] {
  return {
    family: ion.family ?? '',
    expected_mz: finiteNumber(ion.expected_mz) ?? 0,
    observed_mz: finiteNumber(ion.observed_mz) ?? 0,
    loss_label: ion.loss_label ?? '',
    peak_index: finiteNumber(ion.peak_index) ?? 0,
    intensity: finiteNumber(ion.intensity) ?? 0,
    error_ppm: finiteNumber(ion.error_ppm) ?? 0,
  };
}

function normalizeLinkType(
  raw: string | undefined,
  seq2: string | undefined,
  protein2: string | undefined,
): string {
  const value = raw?.trim().toLowerCase();
  if (!value) return !seq2 && !protein2 ? 'monolink' : 'crosslink';
  if (value === 'mono') return 'monolink';
  if (value === 'xlink') return 'crosslink';
  return value;
}

function normalizeFiltering(
  rawFiltering: Record<string, ViewerFiltering>,
  crosslinks: ViewerCrosslink[],
): Record<string, ViewerFiltering> {
  const out: Record<string, ViewerFiltering> = {};
  for (const xl of crosslinks) {
    const filtering = rawFiltering[xl.id] ?? (xl.source_id ? rawFiltering[xl.source_id] : undefined);
    if (filtering) out[xl.id] = filtering;
  }
  return out;
}

function normalizeUnixTimestamp(value: number | null | undefined, legacy: string | undefined): number | null {
  if (typeof value === 'number' && Number.isFinite(value)) return value;
  if (!legacy || !/^\d+$/.test(legacy)) return null;
  return Number(legacy);
}

function unixToIso(seconds: number): string {
  return new Date(seconds * 1000).toISOString();
}

function finiteOrNull(value: number | null | undefined): number | null {
  return typeof value === 'number' && Number.isFinite(value) ? value : null;
}

function finiteNumber(value: number | null | undefined): number | null {
  return typeof value === 'number' && Number.isFinite(value) ? value : null;
}
