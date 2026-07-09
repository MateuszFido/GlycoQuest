import type { ViewerBundle } from './types';
import { normalizeViewerBundle } from './data/normalize';
import { SelectionStore } from './store/selection';
import { renderQcPanel } from './panels/QcPanel';
import { renderNetworkPanel } from './panels/NetworkPanel';
import { renderSpectrumPanel } from './panels/SpectrumPanel';
import { renderHitsTable } from './panels/HitsTable';

export type { ViewerBundle, ViewerCrosslink } from './types';

export interface MountOptions {
  bundle?: unknown;
  bundleUrl?: string;
  theme?: ViewerTheme;
  onSelectCrosslink?: (crosslink: ViewerBundle['crosslinks'][number] | null) => void;
}

export type ViewerTheme = 'default' | 'lcmspector' | Partial<Record<ViewerThemeToken, string>>;

export type ViewerThemeToken =
  | '--gq-surface'
  | '--gq-panel'
  | '--gq-panel-muted'
  | '--gq-border'
  | '--gq-text'
  | '--gq-muted'
  | '--gq-accent'
  | '--gq-experimental'
  | '--gq-theoretical'
  | '--gq-glycan'
  | '--gq-crosslink'
  | '--gq-warning'
  | '--gq-fail';

/**
 * Mount the GlycoQuest crosslink viewer into `container`.
 * Loads `viewer.json` from `bundleUrl` when `bundle` is not provided.
 */
export async function mountViewer(container: HTMLElement, options: MountOptions = {}): Promise<() => void> {
  let rawBundle: unknown;
  if (options.bundle) {
    rawBundle = options.bundle;
  } else {
    const url = options.bundleUrl ?? './viewer.json';
    const response = await fetch(url);
    if (!response.ok) throw new Error(`Failed to load viewer.json: ${response.status}`);
    rawBundle = await response.json();
  }

  const bundle = normalizeViewerBundle(rawBundle);
  const store = new SelectionStore(bundle);
  const cleanups: Array<() => void> = [];

  container.innerHTML = '';

  const root = document.createElement('div');
  root.className = 'gq-viewer';
  applyTheme(root, options.theme);

  const header = buildHeader(bundle);
  const toolbar = buildToolbar(store, bundle);
  const main = document.createElement('main');
  main.className = 'gq-main';
  const left = document.createElement('section');
  left.className = 'gq-column gq-column--left';
  const right = document.createElement('section');
  right.className = 'gq-column gq-column--right';
  main.append(left, right);
  const footer = buildFooter();

  root.append(header, toolbar, main, footer);
  container.appendChild(root);

  if (options.onSelectCrosslink) {
    cleanups.push(
      store.subscribe(() => {
        options.onSelectCrosslink?.(store.selectedCrosslink);
      }),
    );
  }

  cleanups.push(renderHitsTable(left, store));
  cleanups.push(renderQcPanel(left, store));
  cleanups.push(renderNetworkPanel(right, store));
  cleanups.push(renderSpectrumPanel(right, store));

  if (bundle.crosslinks.length > 0) {
    const first = store.visibleCrosslinks[0] ?? bundle.crosslinks[0];
    store.selectCrosslink(first.id);
  }

  return () => {
    for (const fn of cleanups) fn();
    container.innerHTML = '';
  };
}

export { normalizeViewerBundle } from './data/normalize';

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
      ${meta.generated_at_iso ? `<span>Generated: ${esc(formatDate(meta.generated_at_iso))}</span>` : ''}
      ${meta.resume ? '<span>Resume mode</span>' : ''}
    </div>`;
  return header;
}

function buildToolbar(store: SelectionStore, bundle: ViewerBundle): HTMLElement {
  const bar = document.createElement('div');
  bar.className = 'gq-toolbar';

  const proteinSelect = document.createElement('label');
  proteinSelect.className = 'gq-field';
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
  failedLabel.className = 'gq-check';
  const failedCb = document.createElement('input');
  failedCb.type = 'checkbox';
  failedCb.checked = store.filters.showFailed;
  failedCb.addEventListener('change', () => store.setShowFailed(failedCb.checked));
  failedLabel.append(failedCb, ' Show failed hits');
  bar.appendChild(failedLabel);

  const scoreLabel = document.createElement('label');
  scoreLabel.className = 'gq-field';
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
    'GlycoQuest viewer (MIT). Approximate annotations are labeled when exact fragment evidence is unavailable.';
  return footer;
}

function applyTheme(root: HTMLElement, theme: ViewerTheme | undefined): void {
  if (!theme || theme === 'default') return;
  if (theme === 'lcmspector') {
    root.classList.add('gq-viewer--lcmspector');
    return;
  }
  for (const [key, value] of Object.entries(theme)) {
    if (value) root.style.setProperty(key, value);
  }
}

function formatDate(iso: string): string {
  const date = new Date(iso);
  if (Number.isNaN(date.getTime())) return iso;
  return date.toLocaleString(undefined, {
    year: 'numeric',
    month: 'short',
    day: '2-digit',
    hour: '2-digit',
    minute: '2-digit',
  });
}

function esc(s: string): string {
  return s.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
}
