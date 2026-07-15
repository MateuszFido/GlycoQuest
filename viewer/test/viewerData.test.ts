import { normalizeViewerBundle } from '../src/data/normalize';
import { parseGlycanComposition, renderGlycanSvg } from '../src/glycan/snfg';
import {
  buildFullSequenceLayout,
  displayProteinLabel,
  groupCrosslinksByEndpoint,
  residuePointForPosition,
  residueSegmentsForRange,
  stableProteinIdsForSelection,
} from '../src/panels/pairMapLayout';
import { buildPeakAnnotations, buildPeakMarkers, pickTopPeakLabels } from '../src/panels/spectrumPeaks';
import { formatSpectrumSummary } from '../src/panels/spectrumSummary';
import { SelectionStore } from '../src/store/selection';
import type { ViewerBundle, ViewerCrosslink, ViewerFiltering } from '../src/types';

const v1Bundle = {
  viewer_schema_version: 1,
  meta: {
    project: 'demo',
    input_label: 'run.mzXML',
    crosslinker: 'dss',
    xlink_sites: 'K:K',
    glycan_library: 'nglyc',
    resume: false,
    generated_at: '60',
    total_hits: 3,
    passing_hits: 3,
  },
  proteins: [
    { id: 'P1', display_name: 'P1', sequence: 'ACDEFGHIK' },
    { id: 'P2', display_name: 'P2', sequence: 'LMNPQRSTV' },
    { id: 'P3', display_name: 'P3', sequence: 'WYACDEFGH' },
  ],
  crosslinks: [
    crosslink('x1', 'P1', 'P2', 2, 5, 7),
    crosslink('x2', 'P2', 'P1', 4, 7, 8),
    crosslink('x3', 'P1', 'P3', 6, 8, 9),
  ],
  qc: {
    funnel: [],
    outcomes: [],
    glycan_top: [],
    site_dist: [],
    score_hist: { bins: 1, min: 0, max: 1, counts: [1], n: 1 },
    ppm_hist: { bins: 1, min: 0, max: 1, counts: [1], n: 1 },
  },
  spectra: {
    '7': {
      mz: [100.1, 200.2],
      intensity: [10, 50],
      precursor_mz: 600.2,
      charge: 3,
    },
  },
  isotope_pairs: {
    '7': {
      id: '7:9',
      light_scan: 7,
      heavy_scan: 9,
      rt_light_min: 20.1,
      rt_heavy_min: 20.4,
      mz_light: 600.2,
      mz_heavy: 604.2,
      light_charge: 3,
      heavy_charge: 3,
    },
  },
  fragments: {
    x1: {
      theoretical_mz: [100.2, 300.3],
      labels: ['p1_b1', 'p2_y2'],
      matched_indices: [0],
    },
  },
};

const bundle = normalizeViewerBundle(v1Bundle) as ViewerBundle;

test('normalizes legacy bundles into schema v4 viewer data without approximate annotations', () => {
  assertEqual(bundle.viewer_schema_version, 4);
  assertEqual(bundle.meta.generated_at_iso, '1970-01-01T00:01:00.000Z');
  assertEqual(bundle.meta.generated_at_unix, 60);
  assertEqual(bundle.spectra['7'].retention_time_min, null);
  assertEqual(bundle.isotope_pairs['7'].heavy_scan, 9);
  assertEqual(bundle.isotope_pairs['7'].rt_light_min, 20.1);
  assertEqual(bundle.crosslinks[0].retention_time_min, null);
  assertEqual(bundle.crosslinks[0].protein_pair_key, 'P1|P2');
  assertArrayEqual(bundle.crosslinks[0].glyco_sites, []);
  assertEqual(Object.keys(bundle.filtering).length, 0);
});

test('does not convert legacy fragments into approximate mirror fragments', () => {
  assertEqual(Object.keys((bundle as any).mirror_fragments ?? {}).length, 0);
});

test('selection focuses the sequence map on the selected protein pair', () => {
  const store = new SelectionStore(bundle);
  store.selectCrosslink('x1');

  assertArrayEqual(store.focusedProteinIds, ['P1', 'P2']);
  assertArrayEqual(store.selectedPairCrosslinks.map((xl) => xl.id), ['x1', 'x2']);
  assertEqual(store.selectedPairCrosslinks.every((xl) => xl.protein_pair_key === 'P1|P2'), true);
});

test('normalizes monolinks as single-protein selections', () => {
  const raw = {
    ...structuredClone(v1Bundle),
    meta: {
      ...structuredClone(v1Bundle.meta),
      crosslinker_mw: 138.0680796,
    },
    crosslinks: [
      {
        ...crosslink('m1', 'P1', '', 2, 0, 7),
        link_type: 'monolink',
        protein2: '',
        pep_seq2: '',
        pep_pos2: null,
        link_pos2: null,
        abs_pos2: null,
        protein_pair_key: 'P1',
        xlinker_mass: 156.07864,
      },
    ],
  };

  const normalized = normalizeViewerBundle(raw) as ViewerBundle;
  const store = new SelectionStore(normalized);
  store.selectCrosslink('m1');

  assertEqual(normalized.meta.crosslinker_mw, 138.0680796);
  assertEqual(normalized.crosslinks[0].link_type, 'monolink');
  assertEqual(normalized.crosslinks[0].protein_pair_key, 'P1');
  assertEqual(normalized.crosslinks[0].xlinker_mass, 156.07864);
  assertArrayEqual(store.focusedProteinIds, ['P1']);
  assertArrayEqual(store.selectedPairCrosslinks.map((xl) => xl.id), ['m1']);
});

test('preserves every glycosylation site in multi-glycan results', () => {
  const raw = {
    ...structuredClone(v1Bundle),
    crosslinks: [
      {
        ...crosslink('multi1', 'P1', 'P2', 2, 5, 7),
        glycan_composition: 'HexNAc(2)Hex(5)',
        glyco_residue: 'N',
        glyco_peptide: 1,
        glyco_sites: [
          { peptide: 1, peptide_position: 2, residue: 'N', sequon_present: true, plausible: true },
          { peptide: 2, peptide_position: 3, residue: 'N', sequon_present: false, plausible: false },
        ],
      },
    ],
  };

  const normalized = normalizeViewerBundle(raw) as ViewerBundle;
  assertEqual(normalized.crosslinks[0].glyco_sites.length, 2);
  assertEqual(normalized.crosslinks[0].glyco_sites[0].peptide_position, 2);
  assertEqual(normalized.crosslinks[0].glyco_sites[1].plausible, false);
});

test('spectrum summary is Scan | precursor | charge | scan_time, never an axis label', () => {
  assertEqual(
    formatSpectrumSummary({
      scan: 6715,
      precursorMz: 938.4606,
      charge: 4,
      scanTimeMin: 19.71,
    }),
    'Scan 6715 | precursor 938.4606 m/z | 4+ | scan_time 19.71 min',
  );
  assertEqual(
    formatSpectrumSummary({
      scan: 100,
      precursorMz: 500.1,
      charge: 2,
      scanTimeMin: null,
    }),
    'Scan 100 | precursor 500.1000 m/z | 2+',
  );
});

test('spectrum peak annotations use Filtering diagnostic and xQuest rows only', () => {
  const xl = asViewerCrosslink(crosslink('x1', 'P1', 'P2', 2, 5, 7));
  xl.pep_seq1 = 'ACDK';
  xl.link_pos1 = 3;
  xl.glycan_composition = 'HexNAc(1)';

  const peaks = buildPeakAnnotations(
    {
      mz: [100.1, 204.0868, 500.5],
      intensity: [10, 80, 20],
      retention_time_min: null,
      precursor_mz: 600.2,
      charge: 3,
    },
    filteringFor('x1'),
    xl,
  );

  assertEqual(peaks.length, 2);
  assertEqual(peaks[0].kind, 'xquest');
  assertEqual(peaks[0].label, 'p1_b1');
  assertEqual(peaks[0].detailLines.some((line) => line.includes('P1 ACDK crosslink @3')), true);
  assertEqual(peaks[1].kind, 'diagnostic');
  assertEqual(peaks[1].label, 'HexNAc 204.0868');
  assertEqual(peaks[1].detailLines.some((line) => line.includes('HexNAc(1)')), true);
});

test('spectrum peak annotations do not remap Filtering rows onto isotope partner scans', () => {
  const xl = asViewerCrosslink(crosslink('x1', 'P1', 'P2', 2, 5, 7));
  xl.pep_seq1 = 'ACDK';
  xl.link_pos1 = 3;
  xl.glycan_composition = 'HexNAc(1)';

  const peaks = (buildPeakAnnotations as any)(
    {
      mz: [50.0, 100.23, 204.087],
      intensity: [500, 40, 80],
      retention_time_min: null,
      precursor_mz: 604.2,
      charge: 3,
    },
    filteringFor('x1'),
    xl,
    { remapMatchedFragments: true },
  );

  assertEqual(peaks.length, 0);
});

test('spectrum clickable diamonds represent every matched peak with grouped card rows', () => {
  const xl = asViewerCrosslink(crosslink('x1', 'P1', 'P2', 2, 5, 7));
  xl.pep_seq1 = 'AKVFKD';
  xl.link_pos1 = 5;

  const spectrum = {
    mz: [100.1, 200.2, 300.3, 400.4],
    intensity: [10, 50, 75, 20],
    retention_time_min: null,
    precursor_mz: 600.2,
    charge: 3,
  };
  const annotations = buildPeakAnnotations(
    spectrum,
    filteringFor('x1', {
      xquestRows: [
        xquestIon('p1_b2', 0, 100.11, 100.1, 10),
        xquestIon('p1_b5', 1, 200.21, 200.2, 50),
        xquestIon('p1_y2', 1, 200.22, 200.2, 50),
        xquestIon('p1_b6', 2, 300.31, 300.3, 75),
      ],
      diagnosticRows: [],
    }),
    xl,
  );

  xl.glycan_composition = 'HexNAc(2)';
  xl.loss_label = '-H2O';
  const markers = buildPeakMarkers(annotations, {
    crosslink: xl,
    crosslinker: 'DSS',
    xlinkSites: 'K:K',
  });

  assertEqual(markers.length, 3);
  assertArrayEqual(
    markers.map((marker) => marker.observedMz),
    [100.1, 200.2, 300.3],
  );
  assertEqual(markers[1].annotations.length, 2);
  assertEqual(markers[1].detailRows.some((row) => row.label === 'm/z' && row.value === '200.2000'), true);
  assertEqual(
    markers[1].detailRows.some(
      (row) => row.label === 'Sequence' && row.value.includes('AKVF[K]D'),
    ),
    true,
  );
  assertEqual(
    markers[1].detailRows.some((row) => row.label === 'Crosslinker' && row.value === 'DSS K:K'),
    true,
  );
  assertEqual(
    markers[1].detailRows.some((row) => row.label === 'Glycan' && row.value === 'HexNAc(2) -H2O'),
    true,
  );
});

test('monolink spectrum markers use red cards with crosslinker mass and chemistry', () => {
  const xl = asViewerCrosslink({
    ...crosslink('m1', 'P1', '', 2, 0, 7),
    link_type: 'monolink',
    pep_seq1: 'AKVFKD',
    link_pos1: 5,
    protein2: '',
    pep_seq2: '',
    pep_pos2: null,
    link_pos2: null,
    abs_pos2: null,
    protein_pair_key: 'P1',
    xlinker_mass: 156.07864,
  });

  const annotations = buildPeakAnnotations(
    {
      mz: [100.1],
      intensity: [10],
      retention_time_min: null,
      precursor_mz: 600.2,
      charge: 3,
    },
    filteringFor('m1', {
      xquestRows: [xquestIon('p1_b2', 0, 100.11, 100.1, 10)],
      diagnosticRows: [],
    }),
    xl,
  );
  const markers = buildPeakMarkers(annotations, {
    crosslink: xl,
    crosslinker: 'DSS',
    xlinkSites: 'K:K',
    crosslinkerMw: 138.0680796,
  });

  assertEqual(annotations[0].kind, 'monolink');
  assertEqual(markers[0].kind, 'monolink');
  assertEqual(
    markers[0].detailRows.some((row) => row.label === 'Link type' && row.value === 'Monolink'),
    true,
  );
  assertEqual(
    markers[0].detailRows.some((row) => row.label === 'MW' && row.value === '156.0786 Da'),
    true,
  );
  assertEqual(
    markers[0].detailRows.some((row) => row.label === 'Chemistry' && row.value === 'K:K'),
    true,
  );
  assertEqual(markers[0].detailRows.some((row) => row.label === 'P2'), false);
});

test('top peak labels use visible m/z range and intensity order', () => {
  const labels = pickTopPeakLabels(
    [50, 100, 150, 200, 250],
    [0.8, 0.2, 1.0, 0.6, 0.9],
    75,
    225,
    2,
  );

  assertArrayEqual(
    labels.map((label) => label.mz),
    [150, 200],
  );
});

test('normalizes duplicate legacy crosslink ids into unique selectable ids', () => {
  const raw = {
    ...structuredClone(v1Bundle),
    crosslinks: [
      crosslink('dup', 'P1', 'P2', 2, 5, 7),
      crosslink('dup', 'P1', 'P3', 4, 6, 7),
    ],
    fragments: {
      dup: {
        theoretical_mz: [100.2],
        labels: ['p1_b1'],
        matched_indices: [0],
      },
    },
  };

  const normalized = normalizeViewerBundle(raw) as ViewerBundle;

  assertArrayEqual(
    normalized.crosslinks.map((xl) => xl.id),
    ['dup', 'dup__row_1'],
  );
  assertEqual(normalized.crosslinks[1].source_id, 'dup');
  assertEqual(Object.keys((normalized as any).mirror_fragments ?? {}).length, 0);

  const store = new SelectionStore(normalized);
  store.selectCrosslink('dup__row_1');
  assertEqual(store.selectedCrosslink?.protein2, 'P3');
});

test('pair map full sequence layout maps absolute residues to wrapped row coordinates', () => {
  const layout = buildFullSequenceLayout(12, {
    left: 100,
    top: 20,
    residuesPerRow: 5,
    residuePitch: 10,
    rowHeight: 30,
    rowTrackHeight: 20,
  });

  const first = residuePointForPosition(layout, 1);
  const sixth = residuePointForPosition(layout, 6);

  assertEqual(first?.x, 105);
  assertEqual(first?.y, 30);
  assertEqual(first?.row, 0);
  assertEqual(first?.column, 0);
  assertEqual(sixth?.x, 105);
  assertEqual(sixth?.y, 60);
  assertEqual(sixth?.row, 1);
  assertEqual(sixth?.column, 0);
  assertEqual(residuePointForPosition(layout, 0), null);

  const segments = residueSegmentsForRange(layout, 4, 8);
  assertArrayEqual(
    segments.map((segment) => `${segment.start}-${segment.end}@${segment.x},${segment.y},${segment.width}`),
    ['4-5@130,20,20', '6-8@100,50,30'],
  );
});

test('pair map preserves lane order when selecting another crosslink in the same pair', () => {
  const first = asViewerCrosslink(crosslink('stable-a', 'P1', 'P2', 2, 5, 1));
  const reversed = asViewerCrosslink(crosslink('stable-b', 'P2', 'P1', 5, 2, 2));

  const initial = stableProteinIdsForSelection(first, null, []);
  const next = stableProteinIdsForSelection(reversed, initial.pairKey, initial.ids);

  assertArrayEqual(initial.ids, ['P1', 'P2']);
  assertArrayEqual(next.ids, ['P1', 'P2']);
  assertEqual(next.pairKey, 'P1|P2');
});

test('pair map uses compact labels without losing full protein identity', () => {
  const label = displayProteinLabel('sp|P02768|ALBU_HUMAN');

  assertEqual(label.short, 'ALBU_HUMAN');
  assertEqual(label.full, 'sp|P02768|ALBU_HUMAN');
});

test('pair map groups exact endpoint stacks for disambiguation', () => {
  const links = [
    asViewerCrosslink(crosslink('a', 'P1', 'P2', 2, 5, 1)),
    asViewerCrosslink(crosslink('b', 'P1', 'P2', 2, 5, 2)),
    asViewerCrosslink(crosslink('c', 'P1', 'P2', 3, 5, 3)),
  ];

  const groups = groupCrosslinksByEndpoint(links, ['P1', 'P2']);

  assertEqual(groups.length, 2);
  assertArrayEqual(
    groups[0].crosslinks.map((xl) => xl.id),
    ['a', 'b'],
  );
  assertArrayEqual(
    groups[1].crosslinks.map((xl) => xl.id),
    ['c'],
  );
});

test('parses glycan compositions into ordered SNFG components', () => {
  const parsed = parseGlycanComposition('HexNAc(3)Hex(6)NeuAc(1)');

  assertArrayEqual(
    parsed.map((part) => `${part.name}:${part.count}`),
    ['HexNAc:3', 'Hex:6', 'NeuAc:1'],
  );
});

test('renders glycan composition as inline svg markup', () => {
  const svg = renderGlycanSvg('HexNAc(2)Fuc(1)', { size: 18 });

  assertEqual(svg.includes('<svg'), true);
  assertEqual(svg.includes('HexNAc(2)Fuc(1)'), true);
  assertEqual(svg.includes('#0072BC'), true);
  assertEqual(svg.includes('#ED1C24'), true);
});

function crosslink(
  id: string,
  protein1: string,
  protein2: string,
  absPos1: number,
  absPos2: number,
  scan: number,
) {
  return {
    id,
    link_type: 'crosslink',
    protein1,
    pep_pos1: 1,
    pep_seq1: 'ACD',
    link_pos1: 1,
    abs_pos1: absPos1,
    protein2,
    pep_pos2: 1,
    pep_seq2: 'LMN',
    link_pos2: 1,
    abs_pos2: absPos2,
    score: 10,
    soft_score: 11,
    scan,
    charge: 3,
    precursor_mz: 600.2,
    precursor_error_ppm: 1.1,
    topology: protein1 === protein2 ? 'intra' : 'inter',
    glycan_name: null,
    glycan_composition: null,
    glyco_residue: null,
    glyco_peptide: null,
    glyco_sites: [],
    diagnostic_ions: [],
    loss_label: null,
    postfilter_status: 'pass',
    mapped: true,
  };
}

function asViewerCrosslink(xl: any): ViewerCrosslink {
  const left = xl.protein1;
  const right = xl.protein2;
  return {
    ...xl,
    retention_time_min: xl.retention_time_min ?? null,
    source_file: xl.source_file ?? null,
    protein_pair_key: xl.protein_pair_key ?? (!left ? right : !right ? left : left <= right ? `${left}|${right}` : `${right}|${left}`),
    xlinker_mass: xl.xlinker_mass ?? null,
  };
}

function filteringFor(
  crosslinkId: string,
  options: {
    diagnosticRows?: ViewerFiltering['diagnostic_prefilter']['matched_ions'];
    xquestRows?: ViewerFiltering['xquest_search']['matched_ions'];
  } = {},
): ViewerFiltering {
  return {
    input_scan: {
      status: 'available',
      source_file: 'run.mzXML',
      source_artifact: 'spectra/run.mzXML',
      source_row: null,
      scan: 7,
      retention_time_min: 20.1,
      precursor_mz: 600.2,
      charge: 3,
      peak_count: 3,
    },
    diagnostic_prefilter: {
      status: 'matched',
      source_artifact: 'filtered_spectra.tsv',
      source_row: 2,
      matched_family_count: 1,
      matched_families: ['HexNAc'],
      matched_ions:
        options.diagnosticRows ?? [
          {
            family: 'HexNAc',
            expected_mz: 204.0867,
            observed_mz: 204.0868,
            loss_label: '',
            peak_index: 1,
            intensity: 80,
            error_ppm: 0.49,
          },
        ],
    },
    isotope_pair: null,
    glycan_pruning: {
      status: 'retained',
      source_artifact: 'glycan_pruning.tsv',
      source_row: 3,
      selected_glycan: 'HexNAc(1)',
      selected_composition: 'HexNAc(1)',
      retained_count_for_scan: 1,
      required_families: ['HexNAc'],
    },
    xquest_search: {
      status: 'matched',
      source_artifact: 'jobs/HexNAc_1_/results/xquest.xml',
      source_row: null,
      xquest_version: 'xquest 2.1.7',
      score: 10,
      rank: 1,
      xlinkions_matched: '1/10',
      backboneions_matched: '3/20',
      num_matched_ions_alpha: null,
      num_matched_ions_beta: null,
      num_matched_common_ions_alpha: null,
      num_matched_common_ions_beta: null,
      num_matched_xlink_ions_alpha: null,
      num_matched_xlink_ions_beta: null,
      matched_ions: options.xquestRows ?? [xquestIon('p1_b1', 0, 100.2, 100.1, 10)],
      unavailable_reason: null,
    },
    postfilter: {
      status: 'pass',
      source_artifact: 'results/glycoquest_xquest.csv',
      source_row: 2,
      hard_status: 'pass',
      rules: [
        {
          name: 'diagnostic_positive',
          status: 'pass',
          value: 'true',
          threshold: 'required',
        },
      ],
    },
  };
}

function xquestIon(
  label: string,
  peakIndex: number,
  theoreticalMz: number,
  observedMz: number,
  intensity: number,
) {
  return {
    label,
    ion_type: label.replace(/^p[12]_/, ''),
    peptide: label.startsWith('p2_') ? 'beta' : 'alpha',
    position: label.replace(/\D/g, '') || null,
    theoretical_mz: theoreticalMz,
    observed_mz: observedMz,
    error_da: observedMz - theoreticalMz,
    error_ppm: ((observedMz - theoreticalMz) / theoreticalMz) * 1_000_000,
    intensity,
    peak_index: peakIndex,
  };
}

function test(name: string, fn: () => void): void {
  try {
    fn();
  } catch (error) {
    console.error(`not ok - ${name}`);
    throw error;
  }
  console.log(`ok - ${name}`);
}

function assertEqual<T>(actual: T, expected: T): void {
  if (actual !== expected) {
    throw new Error(`expected ${String(expected)}, got ${String(actual)}`);
  }
}

function assertArrayEqual<T>(actual: T[], expected: T[]): void {
  const same =
    actual.length === expected.length && actual.every((value, index) => value === expected[index]);
  if (!same) {
    throw new Error(`expected ${JSON.stringify(expected)}, got ${JSON.stringify(actual)}`);
  }
}
