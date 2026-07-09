import { normalizeViewerBundle } from '../src/data/normalize';
import { SelectionStore } from '../src/store/selection';
import type { ViewerBundle } from '../src/types';

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
  fragments: {
    x1: {
      theoretical_mz: [100.2, 300.3],
      labels: ['p1_b1', 'p2_y2'],
      matched_indices: [0],
    },
  },
};

const bundle = normalizeViewerBundle(v1Bundle) as ViewerBundle;

test('normalizes v1 bundles into schema v2 viewer data', () => {
  assertEqual(bundle.viewer_schema_version, 2);
  assertEqual(bundle.meta.generated_at_iso, '1970-01-01T00:01:00.000Z');
  assertEqual(bundle.meta.generated_at_unix, 60);
  assertEqual(bundle.spectra['7'].retention_time_min, null);
  assertEqual(bundle.crosslinks[0].retention_time_min, null);
  assertEqual(bundle.crosslinks[0].protein_pair_key, 'P1|P2');
});

test('converts legacy fragments into approximate mirror fragments', () => {
  const fragments = bundle.mirror_fragments.x1;
  assertEqual(fragments.annotation_source, 'glycoquest_approx');
  assertArrayEqual(fragments.experimental_mz, [100.1]);
  assertArrayEqual(fragments.experimental_intensity, [10]);
  assertArrayEqual(fragments.theoretical_mz, [100.2, 300.3]);
  assertArrayEqual(fragments.theoretical_intensity, [1, 1]);
  assertArrayEqual(fragments.matched_indices_experimental, [0]);
  assertArrayEqual(fragments.matched_indices_theoretical, [0]);
});

test('selection focuses the sequence map on the selected protein pair', () => {
  const store = new SelectionStore(bundle);
  store.selectCrosslink('x1');

  assertArrayEqual(store.focusedProteinIds, ['P1', 'P2']);
  assertArrayEqual(store.selectedPairCrosslinks.map((xl) => xl.id), ['x1', 'x2']);
  assertEqual(store.selectedPairCrosslinks.every((xl) => xl.protein_pair_key === 'P1|P2'), true);
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
    loss_label: null,
    postfilter_status: 'pass',
    mapped: true,
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
