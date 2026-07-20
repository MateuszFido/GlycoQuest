// Copyright (c) ETH Zurich, Mateusz Fido

import { renderGlycanSvg } from '../glycan/snfg';
import type { SelectionStore } from '../store/selection';

export function renderHitsTable(container: HTMLElement, store: SelectionStore): () => void {
  const panel = document.createElement('section');
  panel.className = 'gq-panel gq-panel--hits';
  panel.innerHTML = '<h2>Links</h2>';
  const body = document.createElement('div');
  body.className = 'gq-panel-body';
  panel.appendChild(body);
  container.appendChild(panel);

  const render = () => {
    const rows = store.visibleCrosslinks;
    if (rows.length === 0) {
      body.innerHTML = '<p class="gq-empty">No crosslinks match current filters.</p>';
      return;
    }

    const table = document.createElement('table');
    table.className = 'gq-table';
    table.innerHTML = `<thead><tr>
      <th>Scan</th><th>RT</th><th>Type</th><th>P1</th><th>P2</th><th>Glycan</th><th>Score</th><th>Status</th>
    </tr></thead>`;
    const tbody = document.createElement('tbody');

    for (const xl of rows) {
      const tr = document.createElement('tr');
      tr.dataset.crosslinkId = xl.id;
      if (xl.id === store.selectedCrosslinkId) tr.classList.add('selected');
      const statusCls = xl.postfilter_status === 'pass' ? 'pass' : 'fail';

      const scanTd = document.createElement('td');
      scanTd.textContent = xl.scan == null ? '' : String(xl.scan);
      const rtTd = document.createElement('td');
      rtTd.textContent = formatRt(xl.retention_time_min);
      const typeTd = document.createElement('td');
      typeTd.textContent = xl.link_type === 'monolink' ? 'monolink' : 'crosslink';
      const p1Td = document.createElement('td');
      p1Td.textContent = xl.protein1;
      const p2Td = document.createElement('td');
      p2Td.textContent = xl.protein2 || '-';
      const glycanTd = document.createElement('td');
      if (xl.glycan_composition) {
        const button = document.createElement('button');
        button.type = 'button';
        button.className = 'gq-glycan-chip';
        button.innerHTML = renderGlycanSvg(xl.glycan_composition, { size: 14 });
        button.title = xl.glycan_composition;
        button.addEventListener('click', (event) => {
          event.stopPropagation();
          store.selectCrosslink(xl.id);
          store.selectGlycan(xl.glycan_composition);
        });
        glycanTd.appendChild(button);
      } else {
        glycanTd.textContent = '-';
      }
      const scoreTd = document.createElement('td');
      scoreTd.textContent = xl.score.toFixed(2);
      const statusTd = document.createElement('td');
      statusTd.className = statusCls;
      statusTd.textContent = xl.postfilter_status;

      tr.append(scanTd, rtTd, typeTd, p1Td, p2Td, glycanTd, scoreTd, statusTd);
      tr.addEventListener('click', () => store.selectCrosslink(xl.id));
      tbody.appendChild(tr);
    }
    table.appendChild(tbody);
    body.replaceChildren(table);
  };

  const unsub = store.subscribe(render);
  render();
  return () => {
    unsub();
    panel.remove();
  };
}

function formatRt(value: number | null): string {
  return value == null ? '-' : `${value.toFixed(2)} min`;
}
