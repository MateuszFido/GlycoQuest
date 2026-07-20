// Copyright (c) ETH Zurich, Mateusz Fido

import type { SelectionStore } from '../store/selection';
import type { ViewerFiltering } from '../types';

export function renderFilteringPanel(container: HTMLElement, store: SelectionStore): () => void {
  const panel = document.createElement('section');
  panel.className = 'gq-panel gq-panel--filtering';
  panel.innerHTML = '<h2>Filtering</h2>';
  const body = document.createElement('div');
  body.className = 'gq-panel-body';
  panel.appendChild(body);
  container.appendChild(panel);

  const render = () => {
    const xl = store.selectedCrosslink;
    if (!xl) {
      body.innerHTML = '<p class="gq-empty">Select a crosslink to inspect Filtering.</p>';
      return;
    }
    const filtering = store.bundle.filtering[xl.id] ?? (xl.source_id ? store.bundle.filtering[xl.source_id] : undefined);
    if (!filtering) {
      body.innerHTML = '<p class="gq-empty">No Filtering record was bundled for this crosslink.</p>';
      return;
    }
    body.replaceChildren(renderFilteringTree(filtering));
  };

  const unsub = store.subscribe(render);
  render();
  return () => {
    unsub();
    panel.remove();
  };
}

function renderFilteringTree(filtering: ViewerFiltering): HTMLElement {
  const tree = document.createElement('div');
  tree.className = 'gq-filtering-tree';
  tree.append(
    step('Input scan', filtering.input_scan.status, [
      row('Artifact', filtering.input_scan.source_artifact),
      row('Source file', filtering.input_scan.source_file ?? 'unavailable'),
      row('Scan', filtering.input_scan.scan),
      row('RT', formatNullable(filtering.input_scan.retention_time_min, ' min')),
      row('Precursor', formatNullable(filtering.input_scan.precursor_mz, ' m/z')),
      row('Charge', filtering.input_scan.charge ? `${filtering.input_scan.charge}+` : 'unavailable'),
      row('Peaks', filtering.input_scan.peak_count),
    ]),
    step('Diagnostic prefilter', filtering.diagnostic_prefilter.status, [
      row('Artifact', sourceLabel(filtering.diagnostic_prefilter.source_artifact, filtering.diagnostic_prefilter.source_row)),
      row('Families', filtering.diagnostic_prefilter.matched_families.join(', ') || 'none'),
      row('Matched ions', filtering.diagnostic_prefilter.matched_ions.length),
      ...filtering.diagnostic_prefilter.matched_ions.slice(0, 12).map((ion) =>
        row(
          `${ion.family} ${ion.observed_mz.toFixed(4)}`,
          `peak ${ion.peak_index}, ${formatNumber(ion.intensity)} intensity, ${ion.error_ppm.toFixed(2)} ppm`,
        ),
      ),
    ]),
    step('Isotope evidence', filtering.isotope_pair?.status ?? 'unavailable', isotopeRows(filtering)),
    step('Glycan pruning', filtering.glycan_pruning.status, [
      row('Artifact', sourceLabel(filtering.glycan_pruning.source_artifact, filtering.glycan_pruning.source_row)),
      row('Selected', filtering.glycan_pruning.selected_composition ?? filtering.glycan_pruning.selected_glycan ?? 'unavailable'),
      row('Retained for scan', filtering.glycan_pruning.retained_count_for_scan),
      row('Required families', filtering.glycan_pruning.required_families.join(', ') || 'none'),
    ]),
    step('xQuest search', filtering.xquest_search.status, [
      row('Artifact', sourceLabel(filtering.xquest_search.source_artifact, filtering.xquest_search.source_row)),
      row('Version', filtering.xquest_search.xquest_version ?? 'unavailable'),
      row('Rank', filtering.xquest_search.rank),
      row('Score', filtering.xquest_search.score.toFixed(2)),
      row('xlink ions', filtering.xquest_search.xlinkions_matched ?? 'unavailable'),
      row('backbone ions', filtering.xquest_search.backboneions_matched ?? 'unavailable'),
      row('Exact rows', filtering.xquest_search.matched_ions.length),
      ...(filtering.xquest_search.unavailable_reason ? [row('Reason', filtering.xquest_search.unavailable_reason)] : []),
      ...filtering.xquest_search.matched_ions.slice(0, 12).map((ion) =>
        row(
          ion.label,
          `${ion.observed_mz.toFixed(4)} m/z, ${ion.error_ppm == null ? 'error unavailable' : `${ion.error_ppm.toFixed(2)} ppm`}`,
        ),
      ),
    ]),
    step('Postfilter', filtering.postfilter.status, [
      row('Artifact', sourceLabel(filtering.postfilter.source_artifact, filtering.postfilter.source_row)),
      row('Hard status', filtering.postfilter.hard_status),
      ...filtering.postfilter.rules.map((rule) =>
        row(rule.name, `${rule.status}: ${rule.value} (${rule.threshold})`),
      ),
    ]),
  );
  return tree;
}

function isotopeRows(filtering: ViewerFiltering): Array<[string, string]> {
  const pair = filtering.isotope_pair;
  if (!pair) return [row('State', 'unavailable')];
  return [
    row('Artifact', sourceLabel(pair.source_artifact, pair.source_row)),
    row('Light scan', `${pair.light_scan} | ${pair.mz_light.toFixed(4)} m/z | ${pair.light_charge}+`),
    row('Heavy scan', `${pair.heavy_scan} | ${pair.mz_heavy.toFixed(4)} m/z | ${pair.heavy_charge}+`),
    row('RT delta', `${Math.abs((pair.rt_heavy_min ?? 0) - (pair.rt_light_min ?? 0)).toFixed(3)} min`),
  ];
}

function step(title: string, status: string, rows: Array<[string, string]>): HTMLElement {
  const item = document.createElement('details');
  item.className = `gq-filter-step gq-filter-step--${status}`;
  item.open = true;
  const summary = document.createElement('summary');
  summary.innerHTML = `<span>${escapeHtml(title)}</span><strong>${escapeHtml(status)}</strong>`;
  item.appendChild(summary);
  const table = document.createElement('dl');
  table.className = 'gq-filter-rows';
  for (const [label, value] of rows) {
    const dt = document.createElement('dt');
    dt.textContent = label;
    const dd = document.createElement('dd');
    dd.textContent = value;
    table.append(dt, dd);
  }
  item.appendChild(table);
  return item;
}

function row(label: string, value: unknown): [string, string] {
  return [label, value == null || value === '' ? 'unavailable' : String(value)];
}

function sourceLabel(artifact: string, sourceRow: number | null): string {
  return sourceRow == null ? artifact : `${artifact} row ${sourceRow}`;
}

function formatNullable(value: number | null, suffix: string): string {
  return value == null || !Number.isFinite(value) ? 'unavailable' : `${value.toFixed(4)}${suffix}`;
}

function formatNumber(value: number): string {
  return Math.abs(value) >= 1000 ? value.toFixed(0) : value.toFixed(2);
}

function escapeHtml(s: string): string {
  return s.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
}
