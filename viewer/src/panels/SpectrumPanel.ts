import uPlot from 'uplot';
import 'uplot/dist/uPlot.min.css';

import type { SelectionStore } from '../store/selection';

let plotInstance: uPlot | null = null;

export function renderSpectrumPanel(container: HTMLElement, store: SelectionStore): () => void {
  const panel = document.createElement('section');
  panel.className = 'gq-panel';
  panel.innerHTML = '<h2>MS/MS spectrum</h2>';
  const body = document.createElement('div');
  body.className = 'gq-panel-body';
  const plotHost = document.createElement('div');
  plotHost.className = 'gq-spectrum-plot';
  const note = document.createElement('p');
  note.className = 'gq-spectrum-note';
  note.textContent = 'Approximate b/y mirror annotation (not identical to xQuest).';
  body.appendChild(plotHost);
  body.appendChild(note);
  panel.appendChild(body);
  container.appendChild(panel);

  const render = () => {
    destroyPlot();
    const xl = store.selectedCrosslink;
    if (!xl || xl.scan == null) {
      plotHost.innerHTML = '<p class="gq-empty">Select a crosslink with a scan number to view MS/MS.</p>';
      return;
    }

    const spectrum = store.bundle.spectra[String(xl.scan)];
    const fragments = store.bundle.fragments[xl.id];
    if (!spectrum || spectrum.mz.length === 0) {
      plotHost.innerHTML = `<p class="gq-empty">No peak list bundled for scan ${xl.scan}. Reduced mzXML may be missing from spectra/.</p>`;
      return;
    }

    plotHost.innerHTML = '';
    plotInstance = createMirrorPlot(plotHost, spectrum.mz, spectrum.intensity, fragments);
  };

  const unsub = store.subscribe(render);
  render();
  return () => {
    unsub();
    destroyPlot();
    panel.remove();
  };
}

function destroyPlot(): void {
  if (plotInstance) {
    plotInstance.destroy();
    plotInstance = null;
  }
}

function createMirrorPlot(
  host: HTMLElement,
  mz: number[],
  intensity: number[],
  fragments: { theoretical_mz: number[]; labels: string[]; matched_indices: number[] } | undefined,
): uPlot {
  const maxObs = Math.max(...intensity, 1);
  const normObs = intensity.map((i) => (i / maxObs) * 100);

  const theoMz = fragments?.theoretical_mz ?? [];
  const theoInt = theoMz.map(() => 50);

  const allMz = [...mz, ...theoMz].sort((a, b) => a - b);
  const xMin = allMz[0] ?? 0;
  const xMax = allMz[allMz.length - 1] ?? 2000;

  const width = Math.max(host.clientWidth || 480, 400);
  const height = 220;

  const data: uPlot.AlignedData = [mz, normObs];

  const opts: uPlot.Options = {
    width,
    height,
    scales: {
      x: { min: xMin * 0.98, max: xMax * 1.02 },
      y: { min: -10, max: 110 },
    },
    axes: [
      { stroke: '#94a3b8', grid: { show: true, stroke: '#334155' } },
      { stroke: '#94a3b8', grid: { show: true, stroke: '#334155' } },
    ],
    series: [
      {},
      {
        label: 'Experimental',
        stroke: '#38bdf8',
        width: 0,
        points: { show: false },
      },
    ],
    hooks: {
      draw: [
        (u) => {
          const ctx = u.ctx;
          const plotLeft = u.bbox.left;
          const plotTop = u.bbox.top;
          const plotHeight = u.bbox.height;
          const baseline = plotTop + plotHeight;

          ctx.save();
          ctx.strokeStyle = '#38bdf8';
          ctx.fillStyle = 'rgba(56, 189, 248, 0.35)';
          ctx.lineWidth = 1;
          const barW = Math.max(1, plotLeft > 0 ? (u.bbox.width / mz.length) * 0.6 : 2);
          for (let i = 0; i < mz.length; i++) {
            const x = u.valToPos(mz[i], 'x', true);
            const y = u.valToPos(normObs[i], 'y', true);
            ctx.beginPath();
            ctx.moveTo(x, baseline);
            ctx.lineTo(x, y);
            ctx.stroke();
            ctx.fillRect(x - barW / 2, y, barW, baseline - y);
          }
          ctx.restore();

          if (!fragments || theoMz.length === 0) return;
          const mirrorBase = plotTop + plotHeight * 0.15;
          ctx.save();
          ctx.strokeStyle = '#fbbf24';
          ctx.lineWidth = 1;
          for (const tmz of theoMz) {
            const x = u.valToPos(tmz, 'x', true);
            ctx.beginPath();
            ctx.moveTo(x, mirrorBase);
            ctx.lineTo(x, mirrorBase + 30);
            ctx.stroke();
          }
          ctx.restore();

          ctx.save();
          ctx.fillStyle = '#22c55e';
          for (const idx of fragments.matched_indices) {
            if (idx >= mz.length) continue;
            const x = u.valToPos(mz[idx], 'x', true);
            const y = u.valToPos(normObs[idx], 'y', true);
            ctx.beginPath();
            ctx.arc(x, y, 3, 0, Math.PI * 2);
            ctx.fill();
          }
          ctx.restore();
        },
      ],
    },
  };

  const plot = new uPlot(opts, data, host);

  if (theoMz.length > 0) {
    const legend = document.createElement('div');
    legend.style.cssText = 'font-size:11px;color:var(--muted);margin-top:4px';
    legend.textContent = `${theoMz.length} theoretical ions · ${fragments?.matched_indices.length ?? 0} matched peaks`;
    host.appendChild(legend);
  }

  return plot;
}
