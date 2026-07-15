export interface ViewerBundle {
  viewer_schema_version: number;
  meta: ViewerMeta;
  proteins: ViewerProtein[];
  crosslinks: ViewerCrosslink[];
  qc: ViewerQc;
  spectra: Record<string, ViewerSpectrum>;
  isotope_pairs: Record<string, ViewerIsotopePair>;
  filtering: Record<string, ViewerFiltering>;
}

export interface ViewerMeta {
  project: string;
  input_label: string;
  crosslinker: string;
  crosslinker_mw: number | null;
  xlink_sites: string;
  glycan_library: string;
  xquest_version: string | null;
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
  /** Original id before frontend uniquification of legacy duplicate rows. */
  source_id?: string;
  link_type: string;
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
  xlinker_mass: number | null;
  topology: string;
  protein_pair_key: string;
  glycan_name: string | null;
  glycan_composition: string | null;
  glyco_residue: string | null;
  glyco_peptide: number | null;
  glyco_sites: ViewerGlycoSite[];
  diagnostic_ions: ViewerDiagnosticIon[];
  loss_label: string | null;
  postfilter_status: string;
  mapped: boolean;
}

export interface ViewerGlycoSite {
  peptide: number;
  peptide_position: number;
  residue: string;
  sequon_present: boolean | null;
  plausible: boolean;
}

export interface ViewerDiagnosticIon {
  family: string;
  expected_mz: number;
  observed_mz: number;
  loss_label: string;
  peak_index: number;
  intensity: number;
  error_ppm: number;
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

export interface ViewerIsotopePair {
  id: string;
  source_artifact: string;
  source_row: number | null;
  light_file: string | null;
  heavy_file: string | null;
  light_scan: number;
  heavy_scan: number;
  rt_light_min: number | null;
  rt_heavy_min: number | null;
  mz_light: number;
  mz_heavy: number;
  light_charge: number;
  heavy_charge: number;
}

export interface ViewerFiltering {
  input_scan: ViewerFilteringInputScan;
  diagnostic_prefilter: ViewerFilteringDiagnosticPrefilter;
  isotope_pair: ViewerFilteringIsotopePair | null;
  glycan_pruning: ViewerFilteringGlycanPruning;
  xquest_search: ViewerFilteringXquestSearch;
  postfilter: ViewerFilteringPostfilter;
}

export interface ViewerFilteringInputScan {
  status: string;
  source_file: string | null;
  source_artifact: string;
  source_row: number | null;
  scan: number | null;
  retention_time_min: number | null;
  precursor_mz: number;
  charge: number;
  peak_count: number;
}

export interface ViewerFilteringDiagnosticPrefilter {
  status: string;
  source_artifact: string;
  source_row: number | null;
  matched_family_count: number;
  matched_families: string[];
  matched_ions: ViewerDiagnosticIon[];
}

export interface ViewerFilteringIsotopePair {
  status: string;
  source_artifact: string;
  source_row: number | null;
  light_scan: number;
  heavy_scan: number;
  rt_light_min: number | null;
  rt_heavy_min: number | null;
  mz_light: number;
  mz_heavy: number;
  light_charge: number;
  heavy_charge: number;
}

export interface ViewerFilteringGlycanPruning {
  status: string;
  source_artifact: string;
  source_row: number | null;
  selected_glycan: string | null;
  selected_composition: string | null;
  retained_count_for_scan: number;
  required_families: string[];
}

export interface ViewerFilteringXquestSearch {
  status: string;
  source_artifact: string;
  source_row: number | null;
  xquest_version: string | null;
  score: number;
  rank: number;
  xlinkions_matched: string | null;
  backboneions_matched: string | null;
  num_matched_ions_alpha: number | null;
  num_matched_ions_beta: number | null;
  num_matched_common_ions_alpha: number | null;
  num_matched_common_ions_beta: number | null;
  num_matched_xlink_ions_alpha: number | null;
  num_matched_xlink_ions_beta: number | null;
  matched_ions: ViewerXquestMatchedIon[];
  unavailable_reason: string | null;
}

export interface ViewerXquestMatchedIon {
  label: string;
  ion_type: string;
  peptide: string | null;
  position: string | null;
  theoretical_mz: number;
  observed_mz: number;
  error_da: number | null;
  error_ppm: number | null;
  intensity: number | null;
  peak_index: number | null;
}

export interface ViewerFilteringPostfilter {
  status: string;
  source_artifact: string;
  source_row: number | null;
  hard_status: string;
  rules: ViewerFilteringRule[];
}

export interface ViewerFilteringRule {
  name: string;
  status: string;
  value: string;
  threshold: string;
}

export type ViewerListener = () => void;

export interface ViewerFilters {
  showFailed: boolean;
  proteinId: string | null;
  minScore: number;
}
