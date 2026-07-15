import uPlot from 'uplot';
import 'uplot/dist/uPlot.min.css';

import type { SelectionStore } from '../store/selection';
import type { ViewerCrosslink, ViewerFiltering, ViewerIsotopePair, ViewerSpectrum } from '../types';
import {
  buildPeakAnnotations,
  buildPeakMarkers,
  pickTopPeakLabels,
  type SpectrumPeakAnnotation,
  type SpectrumPeakMarker,
} from './spectrumPeaks';
import { formatSpectrumSummary } from './spectrumSummary';

export { formatSpectrumSummary } from './spectrumSummary';

export function renderSpectrumPanel(container: HTMLElement, store: SelectionStore): () => void {
  const panel = document.createElement('section');
  panel.className = 'gq-panel gq-panel--spectrum';
  panel.innerHTML = '<h2>MS/MS Plot</h2>';
  const body = document.createElement('div');
  body.className = 'gq-panel-body';
  const summary = document.createElement('div');
  summary.className = 'gq-spectrum-summary';
  const plotHost = document.createElement('div');
  plotHost.className = 'gq-spectrum-plot';
  const plotSurface = document.createElement('div');
  plotSurface.className = 'gq-spectrum-surface';
  const peakDetail = document.createElement('div');
  peakDetail.className = 'gq-peak-detail';
  plotHost.append(plotSurface, peakDetail);
  const note = document.createElement('p');
  note.className = 'gq-spectrum-note';
  body.append(summary, plotHost, note);
  panel.appendChild(body);
  container.appendChild(panel);

  let plotInstance: uPlot | null = null;
  let resizeObserver: ResizeObserver | null = null;
  let selectedCrosslinkId: string | null = null;
  let selectedSpectrumScan: number | null = null;

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
      selectedCrosslinkId = null;
      selectedSpectrumScan = null;
      summary.textContent = '';
      plotSurface.innerHTML = '<p class="gq-empty">Select a crosslink with a scan number to view MS/MS.</p>';
      peakDetail.replaceChildren();
      note.textContent = '';
      return;
    }

    if (selectedCrosslinkId !== xl.id) {
      selectedCrosslinkId = xl.id;
      selectedSpectrumScan = xl.scan;
    }

    const primaryScan = xl.scan;
    const isotopePair = store.bundle.isotope_pairs[String(primaryScan)];
    let activeScan = selectedSpectrumScan ?? primaryScan;
    if (!store.bundle.spectra[String(activeScan)]) {
      activeScan = primaryScan;
      selectedSpectrumScan = primaryScan;
    }

    const spectrum = store.bundle.spectra[String(activeScan)];
    const filtering = store.bundle.filtering[xl.id] ?? (xl.source_id ? store.bundle.filtering[xl.source_id] : undefined);
    const isPrimaryScan = activeScan === primaryScan;
    const scanTimeMin = spectrum?.retention_time_min ?? xl.retention_time_min;

    if (!spectrum || spectrum.mz.length === 0) {
      summary.textContent = formatSpectrumSummary({
        scan: activeScan,
        precursorMz: xl.precursor_mz,
        charge: xl.charge,
        scanTimeMin,
      });
      plotSurface.innerHTML = `<p class="gq-empty">No peak list bundled for scan ${activeScan}. Reduced mzXML may be missing from spectra/.</p>`;
      peakDetail.replaceChildren();
      note.textContent = '';
      return;
    }

    const summaryText = document.createElement('span');
    summaryText.className = 'gq-spectrum-summary-text';
    summaryText.textContent = formatSpectrumSummary({
      scan: activeScan,
      precursorMz: spectrum.precursor_mz || xl.precursor_mz,
      charge: spectrum.charge || xl.charge,
      scanTimeMin,
    });
    const summaryControls = document.createElement('div');
    summaryControls.className = 'gq-spectrum-controls';
    if (isotopePair) {
      summaryControls.append(
        isotopePairSelector(isotopePair, activeScan, store.bundle.spectra, (scan) => {
          selectedSpectrumScan = scan;
          render();
        }),
      );
    }
    summaryControls.append(resetButton(() => resetZoom(plotInstance, spectrum, filtering)));
    summary.replaceChildren(
      summaryText,
      summaryControls,
    );
    plotSurface.innerHTML = '';
    peakDetail.replaceChildren();
    const annotations = buildPeakAnnotations(spectrum, filtering, xl, {
      remapMatchedFragments: !isPrimaryScan,
    });
    const markers = buildPeakMarkers(annotations, {
      crosslink: xl,
      crosslinker: store.bundle.meta.crosslinker,
      crosslinkerMw: store.bundle.meta.crosslinker_mw,
      xlinkSites: store.bundle.meta.xlink_sites,
    });
    plotInstance = createMirrorPlot(plotSurface, spectrum, filtering, markers, (peak) => {
      if (peak) peakDetail.replaceChildren(buildPeakDetail(peak));
      else peakDetail.replaceChildren();
    });
    resizeObserver = new ResizeObserver(() => {
      if (!plotInstance) return;
      plotInstance.setSize({ width: plotWidth(plotSurface), height: 300 });
    });
    resizeObserver.observe(plotSurface);
    note.textContent =
      isPrimaryScan
        ? filteringNote(filtering, annotations)
        : isotopePartnerNote(xl, isotopePair, activeScan, filtering, annotations);
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
  filtering: ViewerFiltering | undefined,
  markers: SpectrumPeakMarker[],
  onPeakSelect: (peak: SpectrumPeakMarker | null) => void,
): uPlot {
  const experimental = normalizePositive(spectrum.intensity);
  const theoreticalMz = xquestTheoreticalMz(filtering);
  const theoretical = normalizePositive(xquestTheoreticalIntensity(filtering));
  const allMz = uniqueSorted([...spectrum.mz, ...theoreticalMz]);
  const xMin = allMz[0] ?? 0;
  const xMax = allMz[allMz.length - 1] ?? 2000;
  const muted = cssVar(host, '--gq-muted');
  const grid = cssVar(host, '--gq-grid');
  let syncPeakMarkers = () => {};

  const opts: uPlot.Options = {
    width: plotWidth(host),
    height: 300,
    legend: { show: true },
    scales: {
      x: { time: false, min: xMin * 0.98, max: xMax * 1.02 },
      y: { min: -1.15, max: 1.15 },
    },
    cursor: {
      show: true,
      x: true,
      y: false,
      points: { show: false },
      drag: { x: true, y: true, uni: 50, setScale: false },
      bind: {
        dblclick: () => null,
      },
    },
    axes: [
      {
        label: 'm/z',
        stroke: muted,
        grid: { show: true, stroke: grid },
        values: (_u: uPlot, vals: number[]) => vals.map((v) => Number(v).toFixed(0)),
      },
      {
        label: 'Rel. intensity',
        stroke: muted,
        grid: { show: true, stroke: grid },
        values: (_u: uPlot, vals: number[]) => vals.map((v) => Math.abs(Number(v)).toFixed(1)),
      },
    ],
    series: [
      {},
      {
        label: 'Experimental',
        stroke: cssVar(host, '--gq-experimental'),
        width: 0,
        points: { show: false },
        paths: () => null,
      },
      {
        label: 'Theoretical',
        stroke: cssVar(host, '--gq-theoretical'),
        width: 0,
        points: { show: false },
        paths: () => null,
      },
    ],
    hooks: {
      draw: [
        (u: uPlot) => {
          drawMirror(u, spectrum.mz, experimental, theoreticalMz, theoretical);
          syncPeakMarkers();
        },
      ],
      setSelect: [
        (u: uPlot) => {
          const { left, width, height } = u.select;
          if (width > 5 || height > 5) {
            const min = u.posToVal(left, 'x');
            const max = u.posToVal(left + width, 'x');
            if (Number.isFinite(min) && Number.isFinite(max) && max > min) {
              u.setScale('x', { min, max });
            }
          }
          u.setSelect({ left: 0, top: 0, width: 0, height: 0 }, false);
        },
      ],
      ready: [
        (u: uPlot) => {
          u.over.addEventListener('dblclick', (event) => {
            event.preventDefault();
            resetZoom(u, spectrum, filtering);
          });
        },
      ],
    },
  };

  const plot = new uPlot(
    opts,
    [
      allMz,
      alignSeries(allMz, spectrum.mz, experimental),
      alignSeries(allMz, theoreticalMz, theoretical),
    ] as uPlot.AlignedData,
    host,
  );
  syncPeakMarkers = installPeakMarkers(plot, markers, onPeakSelect);
  syncPeakMarkers();
  return plot;
}

function drawMirror(
  u: uPlot,
  experimentalMz: number[],
  experimental: number[],
  theoreticalMz: number[],
  theoretical: number[],
): void {
  const ctx = u.ctx;
  const zeroY = u.valToPos(0, 'y', true);
  const topY = u.bbox.top;
  const left = u.bbox.left;
  const width = u.bbox.width;

  ctx.save();
  ctx.beginPath();
  ctx.rect(left, topY, width, u.bbox.height);
  ctx.clip();

  ctx.strokeStyle = cssVar(u.root, '--gq-border') || '#d7dee8';
  ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.moveTo(left, zeroY);
  ctx.lineTo(left + width, zeroY);
  ctx.stroke();

  const experimentalColor = cssVar(u.root, '--gq-experimental') || '#0284c7';
  const theoreticalColor = cssVar(u.root, '--gq-theoretical') || '#be123c';

  drawSticks(u, experimentalMz, experimental, 0, experimentalColor);
  drawSticks(u, theoreticalMz, theoretical.map((value) => -value), 0, theoreticalColor);
  drawPeakLabels(u, experimentalMz, experimental);

  ctx.restore();
}

function drawPeakLabels(u: uPlot, mz: number[], intensity: number[]): void {
  const xMin = u.scales.x.min ?? Number.NEGATIVE_INFINITY;
  const xMax = u.scales.x.max ?? Number.POSITIVE_INFINITY;
  const labels = pickTopPeakLabels(mz, intensity, xMin, xMax, 8);
  if (labels.length === 0) return;

  const ctx = u.ctx;
  const left = u.bbox.left;
  const right = u.bbox.left + u.bbox.width;
  const top = u.bbox.top;
  ctx.save();
  ctx.font = '11px system-ui';
  ctx.lineJoin = 'round';
  for (const label of labels) {
    const x = u.valToPos(label.mz, 'x', true);
    const y = u.valToPos(label.intensity, 'y', true) - 14;
    if (!Number.isFinite(x) || !Number.isFinite(y) || x < left || x > right) continue;
    const clampedX = Math.max(left + 26, Math.min(right - 26, x));
    const clampedY = Math.max(top + 12, y);
    ctx.textAlign = 'center';
    ctx.strokeStyle = '#ffffff';
    ctx.lineWidth = 3;
    ctx.strokeText(label.mz.toFixed(4), clampedX, clampedY);
    ctx.fillStyle = cssVar(u.root, '--gq-text') || '#1f2937';
    ctx.fillText(label.mz.toFixed(4), clampedX, clampedY);
  }
  ctx.restore();
}

function drawSticks(
  u: uPlot,
  mz: number[],
  intensity: number[],
  baseline: number,
  color: string,
): void {
  if (mz.length === 0 || intensity.length === 0) return;
  const ctx = u.ctx;
  const zeroY = u.valToPos(baseline, 'y', true);
  const xMin = u.scales.x.min ?? Number.NEGATIVE_INFINITY;
  const xMax = u.scales.x.max ?? Number.POSITIVE_INFINITY;
  ctx.save();
  ctx.strokeStyle = color;
  ctx.lineWidth = 1.2;
  for (let i = 0; i < mz.length; i++) {
    const mzValue = mz[i];
    if (!Number.isFinite(mzValue) || mzValue < xMin || mzValue > xMax) continue;
    const intensityValue = intensity[i];
    if (!Number.isFinite(intensityValue) || intensityValue === 0) continue;
    const x = u.valToPos(mzValue, 'x', true);
    const y = u.valToPos(intensityValue, 'y', true);
    ctx.beginPath();
    ctx.moveTo(x, zeroY);
    ctx.lineTo(x, y);
    ctx.stroke();
  }
  ctx.restore();
}

function installPeakMarkers(
  plot: uPlot,
  markers: SpectrumPeakMarker[],
  onPeakSelect: (peak: SpectrumPeakMarker | null) => void,
): () => void {
  const layer = document.createElement('div');
  layer.className = 'gq-peak-marker-layer';
  plot.over.appendChild(layer);
  let selectedMarkerId: string | null = null;

  const selectMarker = (marker: SpectrumPeakMarker | null) => {
    selectedMarkerId = marker?.id ?? null;
    onPeakSelect(marker);
    sync();
  };

  plot.over.addEventListener('click', (event) => {
    if ((event.target as Element | null)?.closest?.('.gq-peak-marker')) return;
    selectMarker(null);
  });

  function sync(): void {
    renderPeakMarkers(layer, plot, markers, selectedMarkerId, selectMarker);
  }

  return sync;
}

function renderPeakMarkers(
  layer: HTMLDivElement,
  plot: uPlot,
  markers: SpectrumPeakMarker[],
  selectedMarkerId: string | null,
  onSelect: (marker: SpectrumPeakMarker) => void,
): void {
  const xMin = plot.scales.x.min ?? Number.NEGATIVE_INFINITY;
  const xMax = plot.scales.x.max ?? Number.POSITIVE_INFINITY;
  const fragment = document.createDocumentFragment();
  for (const marker of markers) {
    if (marker.observedMz < xMin || marker.observedMz > xMax) continue;
    const left = plot.valToPos(marker.observedMz, 'x');
    const top = plot.valToPos(marker.relativeIntensity, 'y');
    if (!Number.isFinite(left) || !Number.isFinite(top)) continue;
    if (
      left < -8 ||
      top < -8 ||
      left > plot.over.clientWidth + 8 ||
      top > plot.over.clientHeight + 8
    ) {
      continue;
    }

    const button = document.createElement('button');
    button.type = 'button';
    button.className = [
      'gq-peak-marker',
      `gq-peak-marker--${marker.kind}`,
      marker.id === selectedMarkerId ? 'gq-peak-marker--selected' : '',
    ]
      .filter(Boolean)
      .join(' ');
    button.dataset.peakId = marker.id;
    button.style.left = `${left}px`;
    button.style.top = `${top}px`;
    button.title = `${marker.label} | ${marker.observedMz.toFixed(4)} m/z | intensity ${formatPeakIntensity(marker.intensity)}`;
    button.setAttribute(
      'aria-label',
      `Peak ${marker.label} at ${marker.observedMz.toFixed(4)} m/z, intensity ${formatPeakIntensity(marker.intensity)}`,
    );
    button.addEventListener('pointerdown', (event) => {
      event.preventDefault();
      event.stopPropagation();
    });
    button.addEventListener('mousedown', (event) => {
      event.preventDefault();
      event.stopPropagation();
    });
    button.addEventListener('click', (event) => {
      event.preventDefault();
      event.stopPropagation();
      onSelect(marker);
    });
    fragment.appendChild(button);
  }
  layer.replaceChildren(fragment);
}

function buildPeakDetail(peak: SpectrumPeakMarker): HTMLElement {
  const detail = document.createElement('div');
  detail.className = `gq-peak-card gq-peak-card--${peak.kind}`;
  const title = document.createElement('div');
  title.className = 'gq-peak-card-title';
  title.textContent = peak.detailTitle;
  detail.appendChild(title);
  for (const rowData of peak.detailRows) {
    const row = document.createElement('div');
    row.className = 'gq-peak-card-row';
    const label = document.createElement('span');
    label.textContent = rowData.label;
    const value = document.createElement('strong');
    value.textContent = rowData.value;
    row.append(label, value);
    detail.appendChild(row);
  }
  return detail;
}

function resetZoom(
  plot: uPlot | null,
  spectrum: ViewerSpectrum,
  filtering: ViewerFiltering | undefined,
): void {
  if (!plot) return;
  const allMz = uniqueSorted([...spectrum.mz, ...xquestTheoreticalMz(filtering)]);
  const xMin = allMz[0] ?? 0;
  const xMax = allMz[allMz.length - 1] ?? 2000;
  plot.setScale('x', { min: xMin * 0.98, max: xMax * 1.02 });
  plot.setScale('y', { min: -1.15, max: 1.15 });
}

function resetButton(onClick: () => void): HTMLButtonElement {
  const button = document.createElement('button');
  button.type = 'button';
  button.className = 'gq-icon-button';
  button.title = 'Reset m/z zoom';
  button.textContent = 'Reset';
  button.addEventListener('click', onClick);
  return button;
}

function isotopePairSelector(
  pair: ViewerIsotopePair,
  selectedScan: number,
  spectra: Record<string, ViewerSpectrum>,
  onChange: (scan: number) => void,
): HTMLSelectElement {
  const select = document.createElement('select');
  select.className = 'gq-spectrum-select';
  select.title = 'Show the matched isotope-pair scan';
  for (const role of ['light', 'heavy'] as const) {
    const scan = role === 'light' ? pair.light_scan : pair.heavy_scan;
    const option = document.createElement('option');
    option.value = String(scan);
    option.disabled = !spectra[String(scan)];
    option.textContent = isotopeOptionLabel(pair, role);
    select.appendChild(option);
  }
  select.value = String(selectedScan);
  select.addEventListener('change', () => {
    const scan = Number(select.value);
    if (Number.isFinite(scan)) onChange(scan);
  });
  return select;
}

function isotopeOptionLabel(pair: ViewerIsotopePair, role: 'light' | 'heavy'): string {
  const scan = role === 'light' ? pair.light_scan : pair.heavy_scan;
  const mz = role === 'light' ? pair.mz_light : pair.mz_heavy;
  const charge = role === 'light' ? pair.light_charge : pair.heavy_charge;
  const rt = role === 'light' ? pair.rt_light_min : pair.rt_heavy_min;
  const rtLabel = rt == null ? '' : ` | RT ${rt.toFixed(2)} min`;
  return `${role === 'light' ? 'Light' : 'Heavy'} scan ${scan} | ${mz.toFixed(4)} m/z | ${charge}+${rtLabel}`;
}

function isotopePartnerNote(
  crosslink: ViewerCrosslink,
  pair: ViewerIsotopePair | undefined,
  activeScan: number,
  filtering: ViewerFiltering | undefined,
  annotations: SpectrumPeakAnnotation[],
): string {
  const role = pair ? isotopeRole(pair, activeScan) : null;
  const roleLabel = role ? `${role} ` : '';
  const xquestCount = annotations.filter((annotation) => isFragmentLike(annotation)).length;
  const diagnosticCount = annotations.filter((annotation) => annotation.kind === 'diagnostic').length;
  const exactRows = filtering?.xquest_search.matched_ions.length ?? 0;
  return `Raw ${roleLabel}DSS isotope partner scan for selected crosslink ${crosslink.scan}. ${xquestCount} exact xQuest marker(s), ${diagnosticCount} diagnostic marker(s) from Filtering. ${exactRows} xQuest row(s) are recorded for the selected crosslink scan.`;
}

function isotopeRole(pair: ViewerIsotopePair, scan: number): 'light' | 'heavy' | null {
  if (pair.light_scan === scan) return 'light';
  if (pair.heavy_scan === scan) return 'heavy';
  return null;
}

function normalizePositive(values: number[]): number[] {
  const max = Math.max(...values.filter((value) => Number.isFinite(value)), 1);
  return values.map((value) => (Number.isFinite(value) ? value / max : 0));
}

function xquestTheoreticalMz(filtering: ViewerFiltering | undefined): number[] {
  return (filtering?.xquest_search.matched_ions ?? [])
    .map((ion) => ion.theoretical_mz)
    .filter((value) => Number.isFinite(value));
}

function xquestTheoreticalIntensity(filtering: ViewerFiltering | undefined): number[] {
  return (filtering?.xquest_search.matched_ions ?? [])
    .map((ion) => ion.intensity ?? 1)
    .filter((value) => Number.isFinite(value));
}

function alignSeries(xValues: number[], mz: number[], values: number[]): Array<number | null> {
  const byMz = new Map<number, number>();
  for (let i = 0; i < mz.length; i++) {
    byMz.set(mz[i], values[i]);
  }
  return xValues.map((mzValue) => byMz.get(mzValue) ?? null);
}

function uniqueSorted(values: number[]): number[] {
  return Array.from(new Set(values.filter((value) => Number.isFinite(value)))).sort((a, b) => a - b);
}

function formatPeakIntensity(value: number): string {
  if (Math.abs(value) >= 1000) return value.toFixed(0);
  return value.toFixed(2);
}

function plotWidth(host: HTMLElement): number {
  return Math.max(420, Math.floor(host.clientWidth || 640));
}

function cssVar(el: Element, name: string): string {
  return getComputedStyle(el).getPropertyValue(name).trim();
}

function filteringNote(
  filtering: ViewerFiltering | undefined,
  annotations: SpectrumPeakAnnotation[],
): string {
  const matched = annotations.filter((peak) => isFragmentLike(peak)).length;
  const diagnostic = annotations.filter((peak) => peak.kind === 'diagnostic').length;
  const xquestRows = filtering?.xquest_search.matched_ions.length ?? 0;
  if (!filtering) return 'Raw MS/MS only. No Filtering record was bundled for this crosslink.';
  return `Filtering: ${xquestRows} exact xQuest row(s), ${matched} plotted xQuest marker(s), ${diagnostic} plotted diagnostic marker(s).`;
}

function isFragmentLike(annotation: SpectrumPeakAnnotation): boolean {
  return annotation.kind === 'xquest' || annotation.kind === 'monolink';
}
