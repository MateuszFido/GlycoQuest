import type {
  ViewerBundle,
  ViewerCrosslink,
  ViewerFragments,
  ViewerMeta,
  ViewerMirrorFragments,
  ViewerSpectrum,
} from '../types';

type RawBundle = Omit<ViewerBundle, 'mirror_fragments'> & {
  mirror_fragments?: Record<string, ViewerMirrorFragments>;
  fragments?: Record<string, ViewerFragments>;
  meta: Partial<ViewerMeta> & {
    generated_at?: string;
    generated_at_iso?: string;
    generated_at_unix?: number | null;
  };
  spectra: Record<string, Partial<ViewerSpectrum> & { mz: number[]; intensity: number[] }>;
  crosslinks: Array<Partial<ViewerCrosslink> & Pick<ViewerCrosslink, 'id' | 'protein1' | 'protein2'>>;
};

const LEGACY_MATCH_TOLERANCE_DA = 0.5;

export function normalizeViewerBundle(raw: unknown): ViewerBundle {
  const bundle = raw as RawBundle;
  const meta = normalizeMeta(bundle.meta);
  const spectra = normalizeSpectra(bundle.spectra ?? {});
  const crosslinks = (bundle.crosslinks ?? []).map((xl) => normalizeCrosslink(xl, spectra));
  const fragments = bundle.fragments ?? {};

  return {
    viewer_schema_version: 2,
    meta,
    proteins: bundle.proteins ?? [],
    crosslinks,
    qc: bundle.qc,
    spectra,
    fragments,
    mirror_fragments:
      bundle.mirror_fragments ?? normalizeLegacyMirrorFragments(crosslinks, spectra, fragments),
  };
}

export function proteinPairKey(left: string, right: string): string {
  return left <= right ? `${left}|${right}` : `${right}|${left}`;
}

function normalizeMeta(meta: RawBundle['meta']): ViewerMeta {
  const unix = normalizeUnixTimestamp(meta.generated_at_unix, meta.generated_at);
  const iso = meta.generated_at_iso ?? (unix == null ? meta.generated_at ?? '' : unixToIso(unix));
  return {
    project: meta.project ?? '',
    input_label: meta.input_label ?? '',
    crosslinker: meta.crosslinker ?? '',
    xlink_sites: meta.xlink_sites ?? '',
    glycan_library: meta.glycan_library ?? '',
    resume: Boolean(meta.resume),
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

function normalizeCrosslink(
  xl: RawBundle['crosslinks'][number],
  spectra: Record<string, ViewerSpectrum>,
): ViewerCrosslink {
  const scanKey = xl.scan == null ? null : String(xl.scan);
  const spectrumRt = scanKey == null ? null : spectra[scanKey]?.retention_time_min ?? null;
  return {
    id: xl.id,
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
    topology: xl.topology ?? '',
    protein_pair_key: xl.protein_pair_key ?? proteinPairKey(xl.protein1, xl.protein2),
    glycan_name: xl.glycan_name ?? null,
    glycan_composition: xl.glycan_composition ?? null,
    glyco_residue: xl.glyco_residue ?? null,
    glyco_peptide: xl.glyco_peptide ?? null,
    loss_label: xl.loss_label ?? null,
    postfilter_status: xl.postfilter_status ?? 'pass',
    mapped: xl.mapped ?? false,
  };
}

function normalizeLegacyMirrorFragments(
  crosslinks: ViewerCrosslink[],
  spectra: Record<string, ViewerSpectrum>,
  fragments: Record<string, ViewerFragments>,
): Record<string, ViewerMirrorFragments> {
  const out: Record<string, ViewerMirrorFragments> = {};
  for (const xl of crosslinks) {
    const legacy = fragments[xl.id];
    const spectrum = xl.scan == null ? undefined : spectra[String(xl.scan)];
    if (!legacy) continue;

    const matchedObserved = legacy.matched_indices.filter((idx) => spectrum?.mz[idx] != null);
    const matchedTheoretical = matchedObserved.map((idx) =>
      nearestTheoreticalIndex(spectrum!.mz[idx], legacy.theoretical_mz),
    );

    out[xl.id] = {
      theoretical_mz: legacy.theoretical_mz,
      theoretical_intensity: legacy.theoretical_mz.map(() => 1),
      experimental_mz: matchedObserved.map((idx) => spectrum!.mz[idx]),
      experimental_intensity: matchedObserved.map((idx) => spectrum!.intensity[idx]),
      ion_types: legacy.labels.map((label) => ionType(label)),
      labels: legacy.labels,
      matched_indices_experimental: matchedObserved,
      matched_indices_theoretical: matchedTheoretical,
      annotation_source: 'glycoquest_approx',
    };
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

function nearestTheoreticalIndex(observedMz: number, theoreticalMz: number[]): number {
  let bestIndex = -1;
  let bestDelta = Number.POSITIVE_INFINITY;
  theoreticalMz.forEach((mz, index) => {
    const delta = Math.abs(mz - observedMz);
    if (delta <= LEGACY_MATCH_TOLERANCE_DA && delta < bestDelta) {
      bestIndex = index;
      bestDelta = delta;
    }
  });
  return bestIndex;
}

function ionType(label: string): string {
  if (label.includes('_b')) return 'b';
  if (label.includes('_y')) return 'y';
  return 'unknown';
}
