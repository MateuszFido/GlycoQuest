// Copyright (c) ETH Zurich, Mateusz Fido

import type { SelectionStore } from '../store/selection';
import type { NamedCount, Histogram } from '../types';

export function renderQcPanel(container: HTMLElement, store: SelectionStore): () => void {
  const panel = document.createElement('section');
  panel.className = 'gq-panel gq-panel--qc';
  panel.innerHTML = '<h2>Quality control</h2>';
  const body = document.createElement('div');
  body.className = 'gq-panel-body';
  panel.appendChild(body);
  container.appendChild(panel);

  const render = () => {
    const { qc } = store.bundle;
    const grid = document.createElement('div');
    grid.className = 'gq-chart-grid';

    grid.appendChild(barChart('Prefilter funnel', qc.funnel));
    grid.appendChild(barChart('Post-filter outcome', qc.outcomes));
    grid.appendChild(barChart('Glycans (top 15)', qc.glycan_top));
    grid.appendChild(barChart('Glycosylation sites', qc.site_dist));
    grid.appendChild(histChart('xQuest score', qc.score_hist));
    grid.appendChild(histChart('Precursor error (ppm)', qc.ppm_hist));

    body.replaceChildren(grid);
  };

  const unsub = store.subscribe(render);
  render();
  return () => {
    unsub();
    panel.remove();
  };
}

function barChart(title: string, data: NamedCount[]): HTMLElement {
  const wrap = document.createElement('div');
  wrap.className = 'gq-chart-card';
  const h = document.createElement('h3');
  h.textContent = title;
  h.className = 'gq-chart-title';
  wrap.appendChild(h);

  if (data.length === 0) {
    const empty = document.createElement('p');
    empty.className = 'gq-empty';
    empty.textContent = 'No data';
    wrap.appendChild(empty);
    return wrap;
  }

  const max = Math.max(...data.map((d) => d.count), 1);
  for (const row of data) {
    const el = document.createElement('div');
    el.className = 'gq-bar-row';
    const pct = (row.count / max) * 100;
    el.innerHTML = `
      <span class="gq-bar-label" title="${escapeHtml(row.label)}">${escapeHtml(row.label)}</span>
      <div class="gq-bar-track"><div class="gq-bar-fill" style="width:${pct}%"></div></div>
      <span class="gq-bar-val">${formatCount(row.count)}</span>`;
    wrap.appendChild(el);
  }
  return wrap;
}

function histChart(title: string, hist: Histogram): HTMLElement {
  const wrap = document.createElement('div');
  wrap.className = 'gq-chart-card';
  const h = document.createElement('h3');
  h.textContent = title;
  h.className = 'gq-chart-title';
  wrap.appendChild(h);

  if (hist.n === 0) {
    const empty = document.createElement('p');
    empty.className = 'gq-empty';
    empty.textContent = 'No data';
    wrap.appendChild(empty);
    return wrap;
  }

  const max = Math.max(...hist.counts, 1);
  const bars = document.createElement('div');
  bars.className = 'gq-hist';
  for (const c of hist.counts) {
    const bar = document.createElement('div');
    bar.className = 'gq-hist-bar';
    bar.style.height = `${(c / max) * 100}%`;
    bar.title = String(c);
    bars.appendChild(bar);
  }
  wrap.appendChild(bars);

  const axis = document.createElement('div');
  axis.className = 'gq-hist-axis';
  axis.innerHTML = `<span>${hist.min.toFixed(1)}</span><span>n=${hist.n}</span><span>${hist.max.toFixed(1)}</span>`;
  wrap.appendChild(axis);
  return wrap;
}

function formatCount(n: number): string {
  return Number.isInteger(n) ? String(n) : n.toFixed(1);
}

function escapeHtml(s: string): string {
  return s.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
}
