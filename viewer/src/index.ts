import type { ViewerBundle } from './types';
import { SelectionStore } from './store/selection';
import { renderQcPanel } from './panels/QcPanel';
import { renderSequencePanel } from './panels/SequencePanel';
import { renderNetworkPanel } from './panels/NetworkPanel';
import { renderSpectrumPanel } from './panels/SpectrumPanel';
import { renderHitsTable } from './panels/HitsTable';

export type { ViewerBundle } from './types';

export interface MountOptions {
  bundle?: ViewerBundle;
  bundleUrl?: string;
}

/**
 * Mount the GlycoQuest crosslink viewer into `container`.
 * Loads `viewer.json` from `bundleUrl` when `bundle` is not provided.
 */
export async function mountViewer(container: HTMLElement, options: MountOptions = {}): Promise<() => void> {
  let bundle: ViewerBundle;
  if (options.bundle) {
    bundle = options.bundle;
  } else {
    const url = options.bundleUrl ?? './viewer.json';
    const response = await fetch(url);
    if (!response.ok) throw new Error(`Failed to load viewer.json: ${response.status}`);
    bundle = (await response.json()) as ViewerBundle;
  }

  const store = new SelectionStore(bundle);
  const cleanups: Array<() => void> = [];

  container.innerHTML = '';
  container.id = 'app';

  const header = buildHeader(bundle);
  const toolbar = buildToolbar(store, bundle);
  const main = document.createElement('main');
  main.className = 'gq-main';
  const footer = buildFooter();

  container.append(header, toolbar, main, footer);

  cleanups.push(renderQcPanel(main, store));
  cleanups.push(renderNetworkPanel(main, store));
  cleanups.push(renderHitsTable(main, store));
  cleanups.push(renderSequencePanel(main, store));
  cleanups.push(renderSpectrumPanel(main, store));

  if (bundle.crosslinks.length > 0) {
    const first = store.visibleCrosslinks[0] ?? bundle.crosslinks[0];
    store.selectCrosslink(first.id);
  }

  return () => {
    for (const fn of cleanups) fn();
    container.innerHTML = '';
  };
}

function buildHeader(bundle: ViewerBundle): HTMLElement {
  const header = document.createElement('header');
  header.className = 'gq-header';
  const { meta } = bundle;
  header.innerHTML = `
    <h1>GlycoQuest · ${esc(meta.project)}</h1>
    <div class="gq-meta">
      <span>Input: ${esc(meta.input_label)}</span>
      <span>Crosslinker: ${esc(meta.crosslinker)} (${esc(meta.xlink_sites)})</span>
      <span>Glycans: ${esc(meta.glycan_library)}</span>
      <span>${meta.passing_hits} passing / ${meta.total_hits} total</span>
      ${meta.resume ? '<span>Resume mode</span>' : ''}
    </div>`;
  return header;
}

function buildToolbar(store: SelectionStore, bundle: ViewerBundle): HTMLElement {
  const bar = document.createElement('div');
  bar.className = 'gq-toolbar';

  const proteinSelect = document.createElement('label');
  proteinSelect.textContent = 'Protein ';
  const select = document.createElement('select');
  const allOpt = document.createElement('option');
  allOpt.value = '';
  allOpt.textContent = 'All proteins';
  select.appendChild(allOpt);
  for (const p of bundle.proteins) {
    const opt = document.createElement('option');
    opt.value = p.id;
    opt.textContent = p.display_name;
    select.appendChild(opt);
  }
  select.value = store.selectedProteinId ?? '';
  select.addEventListener('change', () => {
    store.selectProtein(select.value || null);
  });
  proteinSelect.appendChild(select);
  bar.appendChild(proteinSelect);

  const failedLabel = document.createElement('label');
  const failedCb = document.createElement('input');
  failedCb.type = 'checkbox';
  failedCb.addEventListener('change', () => store.setShowFailed(failedCb.checked));
  failedLabel.append(failedCb, ' Show failed hits');
  bar.appendChild(failedLabel);

  const scoreLabel = document.createElement('label');
  scoreLabel.textContent = 'Min score ';
  const scoreInput = document.createElement('input');
  scoreInput.type = 'number';
  scoreInput.min = '0';
  scoreInput.step = '0.5';
  scoreInput.value = '0';
  scoreInput.addEventListener('change', () => store.setMinScore(parseFloat(scoreInput.value) || 0));
  scoreLabel.appendChild(scoreInput);
  bar.appendChild(scoreLabel);

  store.subscribe(() => {
    select.value = store.filters.proteinId ?? '';
  });

  return bar;
}

function buildFooter(): HTMLElement {
  const footer = document.createElement('footer');
  footer.className = 'gq-footer';
  footer.innerHTML =
    'GlycoQuest crosslink viewer (MIT). Network layout inspired by xiNET (Apache-2.0). ' +
    'Spectra use approximate b/y annotation.';
  return footer;
}

function esc(s: string): string {
  return s.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
}
