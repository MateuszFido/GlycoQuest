import type { SelectionStore } from '../store/selection';

export function renderHitsTable(container: HTMLElement, store: SelectionStore): () => void {
  const panel = document.createElement('section');
  panel.className = 'gq-panel';
  panel.innerHTML = '<h2>Crosslinks</h2>';
  const body = document.createElement('div');
  body.className = 'gq-panel-body';
  body.style.maxHeight = '320px';
  body.style.overflow = 'auto';
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
      <th>Scan</th><th>P1</th><th>P2</th><th>Glycan</th><th>Score</th><th>Status</th>
    </tr></thead>`;
    const tbody = document.createElement('tbody');

    for (const xl of rows) {
      const tr = document.createElement('tr');
      if (xl.id === store.selectedCrosslinkId) tr.classList.add('selected');
      const statusCls = xl.postfilter_status === 'pass' ? 'pass' : 'fail';
      tr.innerHTML = `
        <td>${xl.scan ?? ''}</td>
        <td>${esc(xl.protein1)}</td>
        <td>${esc(xl.protein2)}</td>
        <td>${esc(xl.glycan_composition ?? '-')}</td>
        <td>${xl.score.toFixed(2)}</td>
        <td class="${statusCls}">${xl.postfilter_status}</td>`;
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

function esc(s: string): string {
  return s.replace(/&/g, '&amp;').replace(/</g, '&lt;');
}
