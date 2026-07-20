// Copyright (c) ETH Zurich, Mateusz Fido

import type { ViewerCrosslink, ViewerFiltering, ViewerSpectrum } from '../types';

const EXACT_MZ_TOLERANCE_DA = 0.0002;

export interface SpectrumPeakAnnotation {
  kind: 'xquest' | 'diagnostic' | 'monolink';
  marker: 'diamond';
  id: string;
  observedMz: number;
  theoreticalMz: number | null;
  massErrorDa: number | null;
  massErrorPpm: number | null;
  intensity: number;
  relativeIntensity: number;
  peakIndex: number;
  label: string;
  detailTitle: string;
  detailLines: string[];
  detailRows: PeakDetailRow[];
  sequence: string | null;
  crosslinkedSequence: string | null;
  fragmentSequence: string | null;
  protein: string | null;
}

export interface PeakDetailRow {
  label: string;
  value: string;
}

export interface SpectrumPeakMarker {
  id: string;
  kind: 'xquest' | 'diagnostic' | 'monolink' | 'mixed';
  marker: 'diamond';
  observedMz: number;
  intensity: number;
  relativeIntensity: number;
  peakIndex: number;
  label: string;
  detailTitle: string;
  detailRows: PeakDetailRow[];
  annotations: SpectrumPeakAnnotation[];
}

export interface SpectrumPeakMarkerContext {
  crosslink?: ViewerCrosslink;
  crosslinker?: string;
  crosslinkerMw?: number | null;
  xlinkSites?: string;
}

export interface SpectrumPeakLabel {
  mz: number;
  intensity: number;
  peakIndex: number;
}

export interface BuildPeakAnnotationOptions {
  remapMatchedFragments?: boolean;
}

export function buildPeakAnnotations(
  spectrum: ViewerSpectrum,
  filtering: ViewerFiltering | undefined,
  crosslink: ViewerCrosslink,
  _options: BuildPeakAnnotationOptions = {},
): SpectrumPeakAnnotation[] {
  const annotations: SpectrumPeakAnnotation[] = [];
  const maxIntensity = Math.max(...spectrum.intensity.filter((value) => Number.isFinite(value)), 1);

  for (const ion of filtering?.xquest_search.matched_ions ?? []) {
    const peakIndex = validPeakIndex(spectrum, ion.peak_index, ion.observed_mz);
    if (peakIndex == null) continue;
    const observedMz = spectrum.mz[peakIndex];
    const intensity = spectrum.intensity[peakIndex];
    if (!Number.isFinite(observedMz) || !Number.isFinite(intensity)) continue;
    const label = ion.label || `${ion.ion_type}${ion.position ?? ''}`;
    const theoreticalMz = ion.theoretical_mz;
      const fragment = fragmentContext(label, crosslink);
      annotations.push({
        kind: isMonolink(crosslink) ? 'monolink' : 'xquest',
        marker: 'diamond' as const,
        id: `xquest:${peakIndex}:${label}:${theoreticalMz}`,
        observedMz,
        theoreticalMz,
        massErrorDa: ion.error_da ?? observedMz - theoreticalMz,
        massErrorPpm:
          ion.error_ppm ??
          (theoreticalMz === 0 ? null : ((observedMz - theoreticalMz) / theoreticalMz) * 1_000_000),
        intensity,
        relativeIntensity: intensity / maxIntensity,
        peakIndex,
        label,
        detailTitle: label,
        detailLines: fragmentDetailLines(label, observedMz, intensity, theoreticalMz, fragment),
        detailRows: fragmentDetailRows(label, observedMz, intensity, theoreticalMz, fragment),
        sequence: fragment.peptide || null,
        crosslinkedSequence: fragment.crosslinkedSequence,
        fragmentSequence: fragment.fragmentSequence,
        protein: fragment.protein || null,
      });
  }

  for (const ion of filtering?.diagnostic_prefilter.matched_ions ?? []) {
    const peakIndex = validPeakIndex(spectrum, ion.peak_index, ion.observed_mz);
    if (peakIndex == null) continue;
    const observedMz = spectrum.mz[peakIndex];
    const intensity = spectrum.intensity[peakIndex];
    if (!Number.isFinite(observedMz) || !Number.isFinite(intensity)) continue;
    const label = `${ion.family} ${ion.observed_mz.toFixed(4)}`;
    annotations.push({
      kind: 'diagnostic',
      marker: 'diamond',
        id: `diagnostic:${peakIndex}:${ion.family}:${ion.observed_mz}`,
      observedMz,
      theoreticalMz: ion.expected_mz,
      massErrorDa: observedMz - ion.expected_mz,
      massErrorPpm: ion.error_ppm,
      intensity,
      relativeIntensity: intensity / maxIntensity,
      peakIndex,
      label,
      detailTitle: `${ion.family} diagnostic ion`,
      detailLines: [
        `Observed ${observedMz.toFixed(4)} m/z | intensity ${formatIntensity(intensity)}`,
        `Expected ${ion.expected_mz.toFixed(4)} m/z${ion.loss_label ? ` | ${ion.loss_label}` : ''}`,
        `Glycan ${crosslink.glycan_composition ?? crosslink.glycan_name ?? 'unknown'}`,
      ],
      detailRows: [
        { label: 'm/z', value: observedMz.toFixed(4) },
        { label: 'Intensity', value: formatIntensity(intensity) },
        { label: 'Expected', value: ion.expected_mz.toFixed(4) },
        { label: 'Error', value: formatError(observedMz - ion.expected_mz, ion.expected_mz) },
        { label: 'Glycan', value: crosslink.glycan_composition ?? crosslink.glycan_name ?? 'unknown' },
      ],
      sequence: null,
      crosslinkedSequence: null,
      fragmentSequence: null,
      protein: null,
    });
  }

  return annotations.sort((a, b) => a.observedMz - b.observedMz);
}

export function buildPeakMarkers(
  annotations: SpectrumPeakAnnotation[],
  context?: SpectrumPeakMarkerContext,
): SpectrumPeakMarker[] {
  const grouped = new Map<number, SpectrumPeakAnnotation[]>();
  for (const annotation of annotations) {
    const group = grouped.get(annotation.peakIndex) ?? [];
    group.push(annotation);
    grouped.set(annotation.peakIndex, group);
  }

  return Array.from(grouped.values())
    .map((group) => {
      group.sort((a, b) => kindRank(a.kind) - kindRank(b.kind) || a.label.localeCompare(b.label));
      const first = group[0];
      const kinds = new Set(group.map((annotation) => annotation.kind));
      const fragmentLabels = group
        .filter((annotation) => annotation.kind === 'xquest' || annotation.kind === 'monolink')
        .map((annotation) => annotation.label);
      const diagnosticLabels = group
        .filter((annotation) => annotation.kind === 'diagnostic')
        .map((annotation) => annotation.label);
      const kind: SpectrumPeakMarker['kind'] = kinds.size > 1 ? 'mixed' : first.kind;
      return {
        id: `peak:${first.peakIndex}`,
        kind,
        marker: 'diamond' as const,
        observedMz: first.observedMz,
        intensity: first.intensity,
        relativeIntensity: first.relativeIntensity,
        peakIndex: first.peakIndex,
        label: markerLabel(fragmentLabels, diagnosticLabels),
        detailTitle: markerTitle(fragmentLabels, diagnosticLabels),
        detailRows: markerDetailRows(group, context),
        annotations: group,
      };
    })
    .sort((a, b) => a.observedMz - b.observedMz);
}

export function pickTopPeakLabels(
  mz: number[],
  intensity: number[],
  xMin: number,
  xMax: number,
  limit = 8,
): SpectrumPeakLabel[] {
  const labels: SpectrumPeakLabel[] = [];
  for (let i = 0; i < mz.length; i++) {
    const mzValue = mz[i];
    const intensityValue = intensity[i];
    if (
      !Number.isFinite(mzValue) ||
      !Number.isFinite(intensityValue) ||
      mzValue < xMin ||
      mzValue > xMax ||
      intensityValue <= 0
    ) {
      continue;
    }
    labels.push({ mz: mzValue, intensity: intensityValue, peakIndex: i });
  }
  return labels.sort((a, b) => b.intensity - a.intensity).slice(0, Math.max(0, limit));
}

function fragmentDetailLines(
  label: string,
  observedMz: number,
  intensity: number,
  theoreticalMz: number | null,
  fragment: FragmentContext,
): string[] {
  const lines = [
    `Observed ${observedMz.toFixed(4)} m/z | intensity ${formatIntensity(intensity)}`,
    `P${fragment.arm} ${fragment.peptide || '?'} crosslink @${fragment.linkPos ?? '?'}`,
    fragment.protein,
  ];
  if (theoreticalMz != null) lines.splice(1, 0, `Theoretical ${theoreticalMz.toFixed(4)} m/z`);
  if (fragment.fragmentSequence) lines.splice(2, 0, `Fragment ${label} ${fragment.fragmentSequence}`);
  return lines;
}

function fragmentDetailRows(
  label: string,
  observedMz: number,
  intensity: number,
  theoreticalMz: number | null,
  fragment: FragmentContext,
): PeakDetailRow[] {
  const rows: PeakDetailRow[] = [
    { label: 'm/z', value: observedMz.toFixed(4) },
    { label: 'Intensity', value: formatIntensity(intensity) },
    { label: 'Fragment', value: fragment.fragmentSequence ? `${label} ${fragment.fragmentSequence}` : label },
    { label: 'Sequence', value: fragment.crosslinkedSequence ?? fragment.peptide ?? '?' },
  ];
  if (theoreticalMz != null) {
    rows.splice(1, 0, { label: 'Theoretical', value: theoreticalMz.toFixed(4) });
    rows.splice(2, 0, { label: 'Error', value: formatError(observedMz - theoreticalMz, theoreticalMz) });
  }
  if (fragment.protein) rows.push({ label: 'Protein', value: fragment.protein });
  return rows;
}

interface FragmentContext {
  arm: 1 | 2;
  ionType: string | null;
  ordinal: number | null;
  peptide: string;
  linkPos: number | null;
  protein: string;
  fragmentSequence: string | null;
  crosslinkedSequence: string | null;
}

function fragmentContext(label: string, crosslink: ViewerCrosslink): FragmentContext {
  const parsed = /^p([12])_([A-Za-z]+)(\d+)$/.exec(label);
  const arm = parsed?.[1] === '2' || (!parsed && label.startsWith('p2_')) ? 2 : 1;
  const peptide = arm === 1 ? crosslink.pep_seq1 : crosslink.pep_seq2;
  const linkPos = arm === 1 ? crosslink.link_pos1 : crosslink.link_pos2;
  const protein = arm === 1 ? crosslink.protein1 : crosslink.protein2;
  const ionType = parsed?.[2]?.toLowerCase() ?? null;
  const ordinal = parsed?.[3] ? Number(parsed[3]) : null;
  return {
    arm,
    ionType,
    ordinal,
    peptide,
    linkPos,
    protein,
    fragmentSequence: fragmentSequence(peptide, ionType, ordinal),
    crosslinkedSequence: markCrosslink(peptide, linkPos),
  };
}

function fragmentSequence(
  peptide: string,
  ionType: string | null,
  ordinal: number | null,
): string | null {
  if (!peptide || ordinal == null || !Number.isFinite(ordinal) || ordinal <= 0) return null;
  const bounded = Math.min(peptide.length, Math.floor(ordinal));
  if (ionType === 'b') return peptide.slice(0, bounded);
  if (ionType === 'y') return peptide.slice(peptide.length - bounded);
  return null;
}

function markCrosslink(peptide: string, linkPos: number | null): string | null {
  if (!peptide) return null;
  if (linkPos == null || linkPos < 1 || linkPos > peptide.length) return peptide;
  const index = linkPos - 1;
  return `${peptide.slice(0, index)}[${peptide[index]}]${peptide.slice(index + 1)}`;
}

function markerDetailRows(
  group: SpectrumPeakAnnotation[],
  context?: SpectrumPeakMarkerContext,
): PeakDetailRow[] {
  const first = group[0];
  const rows: PeakDetailRow[] = [
    { label: 'm/z', value: first.observedMz.toFixed(4) },
    { label: 'Intensity', value: formatIntensity(first.intensity) },
  ];
  rows.push(...crosslinkContextRows(context));
  const fragmentLabels = group
    .filter((annotation) => annotation.kind === 'xquest' || annotation.kind === 'monolink')
    .map((annotation) => annotation.label);
  if (fragmentLabels.length > 0) {
    rows.push({ label: 'Fragments', value: fragmentLabels.join(', ') });
  }
  const sequence = group.find((annotation) => annotation.crosslinkedSequence)?.crosslinkedSequence;
  if (sequence) rows.push({ label: 'Sequence', value: sequence });
  const fragmentSequenceText = group
    .filter((annotation) => annotation.fragmentSequence)
    .map((annotation) => `${annotation.label} ${annotation.fragmentSequence}`)
    .join(', ');
  if (fragmentSequenceText) rows.push({ label: 'Ions', value: fragmentSequenceText });
  const diagnosticLabels = group
    .filter((annotation) => annotation.kind === 'diagnostic')
    .map((annotation) => annotation.detailTitle);
  if (diagnosticLabels.length > 0) {
    rows.push({ label: 'Diagnostic', value: diagnosticLabels.join(', ') });
  }
  const error = group.find((annotation) => annotation.massErrorDa != null);
  if (error?.massErrorDa != null && error.theoreticalMz != null) {
    rows.push({ label: 'Error', value: formatError(error.massErrorDa, error.theoreticalMz) });
  }
  const protein = group.find((annotation) => annotation.protein)?.protein;
  if (protein) rows.push({ label: 'Protein', value: protein });
  return rows;
}

function crosslinkContextRows(context?: SpectrumPeakMarkerContext): PeakDetailRow[] {
  const rows: PeakDetailRow[] = [];
  const crosslink = context?.crosslink;
  const linker = [context?.crosslinker, context?.xlinkSites].filter(Boolean).join(' ');
  if (linker) rows.push({ label: 'Crosslinker', value: linker });
  if (!crosslink) return rows;
  const monolink = isMonolink(crosslink);
  rows.push({ label: 'Link type', value: monolink ? 'Monolink' : 'Crosslink' });
  const linkerMass = crosslink.xlinker_mass ?? context?.crosslinkerMw ?? null;
  if (linkerMass != null) rows.push({ label: 'MW', value: `${linkerMass.toFixed(4)} Da` });
  if (context?.xlinkSites) rows.push({ label: 'Chemistry', value: context.xlinkSites });

  rows.push({
    label: 'P1',
    value: `${shortProtein(crosslink.protein1)} ${markCrosslink(crosslink.pep_seq1, crosslink.link_pos1) ?? '?'}`,
  });
  if (!monolink) {
    rows.push({
      label: 'P2',
      value: `${shortProtein(crosslink.protein2)} ${markCrosslink(crosslink.pep_seq2, crosslink.link_pos2) ?? '?'}`,
    });
  }
  if (crosslink.abs_pos1 != null || crosslink.abs_pos2 != null) {
    rows.push({
      label: 'Sites',
      value: monolink
        ? `${crosslink.abs_pos1 ?? '?'}`
        : `${crosslink.abs_pos1 ?? '?'} <-> ${crosslink.abs_pos2 ?? '?'}`,
    });
  }
  const glycan = crosslink.glycan_composition ?? crosslink.glycan_name;
  if (glycan) {
    rows.push({
      label: 'Glycan',
      value: crosslink.loss_label ? `${glycan} ${crosslink.loss_label}` : glycan,
    });
  }
  return rows;
}

function shortProtein(id: string): string {
  const parts = id.split('|');
  return parts.length >= 3 ? parts[2] : id;
}

function markerLabel(fragmentLabels: string[], diagnosticLabels: string[]): string {
  const labels = [...fragmentLabels, ...diagnosticLabels];
  if (labels.length === 0) return 'matched peak';
  if (labels.length === 1) return labels[0];
  return `${labels[0]} +${labels.length - 1}`;
}

function markerTitle(fragmentLabels: string[], diagnosticLabels: string[]): string {
  const total = fragmentLabels.length + diagnosticLabels.length;
  if (total === 1) return fragmentLabels[0] ?? diagnosticLabels[0];
  const parts: string[] = [];
  if (fragmentLabels.length > 0) parts.push(`${fragmentLabels.length} matched fragment${fragmentLabels.length === 1 ? '' : 's'}`);
  if (diagnosticLabels.length > 0) parts.push(`${diagnosticLabels.length} diagnostic ion${diagnosticLabels.length === 1 ? '' : 's'}`);
  return parts.join(' + ');
}

function kindRank(kind: SpectrumPeakAnnotation['kind']): number {
  if (kind === 'diagnostic') return 0;
  if (kind === 'monolink') return 1;
  return 2;
}

function isMonolink(xl: ViewerCrosslink): boolean {
  return xl.link_type === 'monolink' || !xl.protein2;
}

function validPeakIndex(
  spectrum: ViewerSpectrum,
  peakIndex: number | null,
  observedMz: number,
): number | null {
  if (peakIndex == null || peakIndex < 0 || peakIndex >= spectrum.mz.length) return null;
  const mz = spectrum.mz[peakIndex];
  if (!Number.isFinite(mz) || Math.abs(mz - observedMz) > EXACT_MZ_TOLERANCE_DA) return null;
  return peakIndex;
}

function formatIntensity(value: number): string {
  if (Math.abs(value) >= 1000) return value.toFixed(0);
  return value.toFixed(2);
}

function formatError(deltaDa: number, theoreticalMz: number): string {
  const ppm = theoreticalMz === 0 ? null : (deltaDa / theoreticalMz) * 1_000_000;
  return ppm == null ? `${deltaDa.toFixed(4)} Da` : `${deltaDa.toFixed(4)} Da (${ppm.toFixed(1)} ppm)`;
}
