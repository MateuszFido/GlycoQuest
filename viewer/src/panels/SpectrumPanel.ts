import uPlot from 'uplot';
import 'uplot/dist/uPlot.min.css';

import type { SelectionStore } from '../store/selection';
import type { ViewerMirrorFragments, ViewerSpectrum } from '../types';

export function renderSpectrumPanel(container: HTMLElement, store: SelectionStore): () => void {
  const panel = document.createElement('section');
  panel.className = 'gq-panel gq-panel--spectrum';
  panel.innerHTML = '<h2>MS/MS Mirror Plot</h2>';
  const body = document.createElement('div');
  body.className = 'gq-panel-body';
  const summary = document.createElement('div');
  summary.className = 'gq-spectrum-summary';
  const plotHost = document.createElement('div');
  plotHost.className = 'gq-spectrum-plot';
  const note = document.createElement('p');
  note.className = 'gq-spectrum-note';
  body.append(summary, plotHost, note);
  panel.appendChild(body);
  container.appendChild(panel);

  let plotInstance: uPlot | null = null;
  let resizeObserver: ResizeObserver | null = null;

  const destroyPlot = () => {
    resizeObserver?.disconnect();
    resizeObserver = null;
    plotInstance?.destroy();
    plotInstance = null;
  };

  const render = () => {
    destroyPlot();
    const xl = store.selectedCrosslink;
    if (!xl || xl.scan == null) {
      summary.textContent = '';
      plotHost.innerHTML = '<p class="gq-empty">Select a crosslink with a scan number to view MS/MS.</p>';
      note.textContent = '';
      return;
    }

    const spectrum = store.bundle.spectra[String(xl.scan)];
    const fragments = store.bundle.mirror_fragments[xl.id];
    if (!spectrum || spectrum.mz.length === 0) {
      summary.textContent = formatSpectrumSummary(xl.scan, xl.retention_time_min, null);
      plotHost.innerHTML = `<p class="gq-empty">No peak list bundled for scan ${xl.scan}. Reduced mzXML may be missing from spectra/.</p>`;
      note.textContent = '';
      return;
    }

    summary.replaceChildren(
      document.createTextNode(formatSpectrumSummary(xl.scan, xl.retention_time_min, spectrum)),
      resetButton(() => resetZoom(plotInstance, spectrum, fragments)),
    );
    plotHost.innerHTML = '';
    plotInstance = createMirrorPlot(plotHost, spectrum, fragments);
    resizeObserver = new ResizeObserver(() => {
      if (!plotInstance) return;
      plotInstance.setSize({ width: plotWidth(plotHost), height: 300 });
    });
    resizeObserver.observe(plotHost);
    note.textContent = fragments
      ? fragmentNote(fragments)
      : 'Raw MS/MS only. No theoretical mirror annotations were bundled for this crosslink.';
  };

  const unsub = store.subscribe(render);
  render();
  return () => {
    unsub();
    destroyPlot();
    panel.remove();
  };
}

function createMirrorPlot(
  host: HTMLElement,
  spectrum: ViewerSpectrum,
  fragments: ViewerMirrorFragments | undefined,
): uPlot {
  const experimental = normalizePositive(spectrum.intensity);
  const theoretical = fragments ? normalizePositive(fragments.theoretical_intensity) : [];
  const theoreticalMz = fragments?.theoretical_mz ?? [];
  const allMz = uniqueSorted([...spectrum.mz, ...theoreticalMz]);
  const xMin = allMz[0] ?? 0;
  const xMax = allMz[allMz.length - 1] ?? 2000;

  const opts: uPlot.Options = {
    width: plotWidth(host),
    height: 300,
    scales: {
      x: { min: xMin * 0.98, max: xMax * 1.02 },
      y: { min: -1.15, max: 1.15 },
    },
    cursor: {
      drag: { x: true, y: false, setScale: true },
    },
    axes: [
      { stroke: cssVar(host, '--gq-muted'), grid: { show: true, stroke: cssVar(host, '--gq-grid') } },
      {
        stroke: cssVar(host, '--gq-muted'),
        grid: { show: true, stroke: cssVar(host, '--gq-grid') },
        values: (_, vals) => vals.map((v) => Math.abs(Number(v)).toFixed(1)),
      },
    ],
    series: [{}, { label: 'Mirror intensity', stroke: 'transparent', points: { show: false } }],
    hooks: {
      draw: [
        (u) => {
          drawMirror(u, spectrum.mz, experimental, theoreticalMz, theoretical, fragments);
        },
      ],
    },
  };

  return new uPlot(opts, [allMz, allMz.map(() => null)] as uPlot.AlignedData, host);
}

function drawMirror(
  u: uPlot,
  experimentalMz: number[],
  experimental: number[],
  theoreticalMz: number[],
  theoretical: number[],
  fragments: ViewerMirrorFragments | undefined,
): void {
  const ctx = u.ctx;
  const zeroY = u.valToPos(0, 'y', true);
  const topY = u.bbox.top;
  const bottomY = u.bbox.top + u.bbox.height;

  ctx.save();
  ctx.strokeStyle = cssVar(u.root, '--gq-border');
  ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.moveTo(u.bbox.left, zeroY);
  ctx.lineTo(u.bbox.left + u.bbox.width, zeroY);
  ctx.stroke();
  ctx.restore();

  drawSticks(u, experimentalMz, experimental, 0, cssVar(u.root, '--gq-experimental'));
  drawSticks(u, theoreticalMz, theoretical.map((value) => -value), 0, cssVar(u.root, '--gq-theoretical'));

  if (fragments) {
    drawMatchedPeaks(u, experimentalMz, experimental, fragments);
    drawLabels(u, theoreticalMz, theoretical, fragments, bottomY, topY);
  }
}

function drawSticks(
  u: uPlot,
  mz: number[],
  intensity: number[],
  baseline: number,
  color: string,
): void {
  const ctx = u.ctx;
  const zeroY = u.valToPos(baseline, 'y', true);
  ctx.save();
  ctx.strokeStyle = color;
  ctx.lineWidth = 1.2;
  for (let i = 0; i < mz.length; i++) {
    const x = u.valToPos(mz[i], 'x', true);
    const y = u.valToPos(intensity[i], 'y', true);
    ctx.beginPath();
    ctx.moveTo(x, zeroY);
    ctx.lineTo(x, y);
    ctx.stroke();
  }
  ctx.restore();
}

function drawMatchedPeaks(
  u: uPlot,
  experimentalMz: number[],
  experimental: number[],
  fragments: ViewerMirrorFragments,
): void {
  const ctx = u.ctx;
  ctx.save();
  ctx.fillStyle = cssVar(u.root, '--gq-crosslink');
  for (const idx of fragments.matched_indices_experimental) {
    if (idx >= experimentalMz.length) continue;
    const x = u.valToPos(experimentalMz[idx], 'x', true);
    const y = u.valToPos(experimental[idx], 'y', true);
    ctx.beginPath();
    ctx.arc(x, y, 3.5, 0, Math.PI * 2);
    ctx.fill();
  }
  ctx.restore();
}

function drawLabels(
  u: uPlot,
  theoreticalMz: number[],
  theoretical: number[],
  fragments: ViewerMirrorFragments,
  bottomY: number,
  topY: number,
): void {
  const ctx = u.ctx;
  ctx.save();
  ctx.fillStyle = cssVar(u.root, '--gq-muted');
  ctx.font = '11px system-ui, -apple-system, Segoe UI, sans-serif';
  ctx.textAlign = 'center';
  for (let i = 0; i < fragments.matched_indices_theoretical.length; i++) {
    const theoIdx = fragments.matched_indices_theoretical[i];
    if (theoIdx < 0 || theoIdx >= theoreticalMz.length) continue;
    const x = u.valToPos(theoreticalMz[theoIdx], 'x', true);
    const y = Math.min(bottomY - 5, Math.max(topY + 12, u.valToPos(-theoretical[theoIdx], 'y', true) + 14));
    ctx.fillText(fragments.labels[theoIdx] ?? '', x, y);
  }
  ctx.restore();
}

function resetZoom(
  plot: uPlot | null,
  spectrum: ViewerSpectrum,
  fragments: ViewerMirrorFragments | undefined,
): void {
  if (!plot) return;
  const allMz = uniqueSorted([...spectrum.mz, ...(fragments?.theoretical_mz ?? [])]);
  const xMin = allMz[0] ?? 0;
  const xMax = allMz[allMz.length - 1] ?? 2000;
  plot.setScale('x', { min: xMin * 0.98, max: xMax * 1.02 });
  plot.setScale('y', { min: -1.15, max: 1.15 });
}

function resetButton(onClick: () => void): HTMLButtonElement {
  const button = document.createElement('button');
  button.type = 'button';
  button.className = 'gq-icon-button';
  button.title = 'Reset zoom';
  button.textContent = 'Reset';
  button.addEventListener('click', onClick);
  return button;
}

function normalizePositive(values: number[]): number[] {
  const max = Math.max(...values.filter((value) => Number.isFinite(value)), 1);
  return values.map((value) => (Number.isFinite(value) ? value / max : 0));
}

function uniqueSorted(values: number[]): number[] {
  return Array.from(new Set(values.filter((value) => Number.isFinite(value)))).sort((a, b) => a - b);
}

function plotWidth(host: HTMLElement): number {
  return Math.max(420, Math.floor(host.clientWidth || 640));
}

function cssVar(el: Element, name: string): string {
  return getComputedStyle(el).getPropertyValue(name).trim();
}

function fragmentNote(fragments: ViewerMirrorFragments): string {
  const matched = fragments.matched_indices_experimental.length;
  const theoretical = fragments.theoretical_mz.length;
  const source =
    fragments.annotation_source === 'glycoquest_approx'
      ? 'Approximate GlycoQuest b/y annotations'
      : `${fragments.annotation_source} annotations`;
  return `${source}: ${theoretical} theoretical ions, ${matched} matched experimental peak(s).`;
}

function formatSpectrumSummary(
  scan: number,
  rt: number | null,
  spectrum: ViewerSpectrum | null,
): string {
  const parts = [`Scan ${scan}`];
  if (rt != null) parts.push(`RT ${rt.toFixed(2)} min`);
  if (spectrum) {
    parts.push(`precursor ${spectrum.precursor_mz.toFixed(4)} m/z`);
    if (spectrum.charge > 0) parts.push(`${spectrum.charge}+`);
  }
  return parts.join(' | ');
}
