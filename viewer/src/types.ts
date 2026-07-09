export interface ViewerBundle {
  viewer_schema_version: number;
  meta: ViewerMeta;
  proteins: ViewerProtein[];
  crosslinks: ViewerCrosslink[];
  qc: ViewerQc;
  spectra: Record<string, ViewerSpectrum>;
  fragments?: Record<string, ViewerFragments>;
  mirror_fragments: Record<string, ViewerMirrorFragments>;
}

export interface ViewerMeta {
  project: string;
  input_label: string;
  crosslinker: string;
  xlink_sites: string;
  glycan_library: string;
  resume: boolean;
  generated_at: string;
  generated_at_iso: string;
  generated_at_unix: number | null;
  total_hits: number;
  passing_hits: number;
}

export interface ViewerProtein {
  id: string;
  display_name: string;
  sequence: string;
}

export interface ViewerCrosslink {
  id: string;
  protein1: string;
  pep_pos1: number | null;
  pep_seq1: string;
  link_pos1: number | null;
  abs_pos1: number | null;
  protein2: string;
  pep_pos2: number | null;
  pep_seq2: string;
  link_pos2: number | null;
  abs_pos2: number | null;
  score: number;
  soft_score: number;
  scan: number | null;
  retention_time_min: number | null;
  source_file: string | null;
  charge: number;
  precursor_mz: number;
  precursor_error_ppm: number;
  topology: string;
  protein_pair_key: string;
  glycan_name: string | null;
  glycan_composition: string | null;
  glyco_residue: string | null;
  glyco_peptide: number | null;
  loss_label: string | null;
  postfilter_status: string;
  mapped: boolean;
}

export interface ViewerQc {
  funnel: NamedCount[];
  outcomes: NamedCount[];
  glycan_top: NamedCount[];
  site_dist: NamedCount[];
  score_hist: Histogram;
  ppm_hist: Histogram;
}

export interface NamedCount {
  label: string;
  count: number;
}

export interface Histogram {
  bins: number;
  min: number;
  max: number;
  counts: number[];
  n: number;
}

export interface ViewerSpectrum {
  mz: number[];
  intensity: number[];
  retention_time_min: number | null;
  precursor_mz: number;
  charge: number;
}

export interface ViewerFragments {
  theoretical_mz: number[];
  labels: string[];
  matched_indices: number[];
}

export interface ViewerMirrorFragments {
  theoretical_mz: number[];
  theoretical_intensity: number[];
  experimental_mz: number[];
  experimental_intensity: number[];
  ion_types: string[];
  labels: string[];
  matched_indices_experimental: number[];
  matched_indices_theoretical: number[];
  annotation_source: 'glycoquest_approx' | 'xquest' | 'external';
}

export type ViewerListener = () => void;

export interface ViewerFilters {
  showFailed: boolean;
  proteinId: string | null;
  minScore: number;
}
