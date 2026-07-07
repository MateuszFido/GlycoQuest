import type { SelectionStore } from '../store/selection';
import type { ViewerCrosslink, ViewerProtein } from '../types';

const BAR_HEIGHT = 14;
const BAR_GAP = 36;
const MARGIN = { top: 24, left: 120, right: 24 };

interface LayoutProtein {
  protein: ViewerProtein;
  y: number;
  barWidth: number;
}

/**
 * Clean-room SVG crosslink network (xiNET-inspired protein bars + residue-resolution edges).
 */
export function renderNetworkPanel(container: HTMLElement, store: SelectionStore): () => void {
  const panel = document.createElement('section');
  panel.className = 'gq-panel gq-panel--network';
  panel.innerHTML = '<h2>Crosslink network</h2>';
  const body = document.createElement('div');
  body.className = 'gq-panel-body';
  panel.appendChild(body);
  container.appendChild(panel);

  const render = () => {
    const crosslinks = store.visibleCrosslinks.filter((xl) => xl.mapped);
    const proteinIds = new Set<string>();
    for (const xl of crosslinks) {
      proteinIds.add(xl.protein1);
      proteinIds.add(xl.protein2);
    }
    const proteins = store.bundle.proteins.filter((p) => proteinIds.has(p.id));
    if (proteins.length === 0) {
      body.innerHTML =
        '<p class="gq-empty">No mapped crosslinks to display. Check FASTA protein IDs match xQuest output.</p>';
      return;
    }

    const maxLen = Math.max(...proteins.map((p) => p.sequence.length), 1);
    const barAreaWidth = 600;
    const height = MARGIN.top + proteins.length * BAR_GAP + 40;

    const svg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
    svg.setAttribute('class', 'gq-network-svg');
    svg.setAttribute('viewBox', `0 0 ${MARGIN.left + barAreaWidth + MARGIN.right} ${height}`);

    const layouts: LayoutProtein[] = proteins.map((protein, i) => ({
      protein,
      y: MARGIN.top + i * BAR_GAP,
      barWidth: Math.max(40, (protein.sequence.length / maxLen) * barAreaWidth),
    }));

    const layoutById = new Map(layouts.map((l) => [l.protein.id, l]));

    for (const layout of layouts) {
      const { protein, y, barWidth } = layout;
      const bar = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
      bar.setAttribute('class', 'protein-bar');
      bar.setAttribute('x', String(MARGIN.left));
      bar.setAttribute('y', String(y));
      bar.setAttribute('width', String(barWidth));
      bar.setAttribute('height', String(BAR_HEIGHT));
      bar.setAttribute('rx', '3');
      svg.appendChild(bar);

      const label = document.createElementNS('http://www.w3.org/2000/svg', 'text');
      label.setAttribute('class', 'protein-label');
      label.setAttribute('x', '8');
      label.setAttribute('y', String(y + BAR_HEIGHT - 2));
      label.textContent = truncateId(protein.display_name, 14);
      svg.appendChild(label);

      const lenLabel = document.createElementNS('http://www.w3.org/2000/svg', 'text');
      lenLabel.setAttribute('class', 'protein-label');
      lenLabel.setAttribute('x', String(MARGIN.left + barWidth + 6));
      lenLabel.setAttribute('y', String(y + BAR_HEIGHT - 2));
      lenLabel.setAttribute('fill', '#94a3b8');
      lenLabel.textContent = `${protein.sequence.length} aa`;
      svg.appendChild(lenLabel);
    }

    for (const xl of crosslinks) {
      drawEdge(svg, xl, layoutById, store.selectedCrosslinkId === xl.id, () =>
        store.selectCrosslink(xl.id),
      );
    }

    body.replaceChildren(svg);
  };

  const unsub = store.subscribe(render);
  render();
  return () => {
    unsub();
    panel.remove();
  };
}

function drawEdge(
  svg: SVGSVGElement,
  xl: ViewerCrosslink,
  layouts: Map<string, LayoutProtein>,
  selected: boolean,
  onClick: () => void,
): void {
  const l1 = layouts.get(xl.protein1);
  const l2 = layouts.get(xl.protein2);
  if (!l1 || !l2 || !xl.abs_pos1 || !xl.abs_pos2) return;

  const x1 = residueX(l1, xl.abs_pos1);
  const y1 = l1.y + BAR_HEIGHT / 2;
  const x2 = residueX(l2, xl.abs_pos2);
  const y2 = l2.y + BAR_HEIGHT / 2;

  const path = document.createElementNS('http://www.w3.org/2000/svg', 'path');
  const midY = (y1 + y2) / 2;
  const d =
    y1 === y2
      ? `M ${x1} ${y1} L ${x2} ${y2}`
      : `M ${x1} ${y1} C ${x1} ${midY}, ${x2} ${midY}, ${x2} ${y2}`;
  path.setAttribute('d', d);
  let cls = 'xlink-edge';
  if (xl.glycan_composition) cls += ' xlink-edge--glycan';
  if (selected) cls += ' xlink-edge--selected';
  if (xl.postfilter_status !== 'pass') cls += ' xlink-edge--failed';
  path.setAttribute('class', cls);
  path.addEventListener('click', onClick);
  const title = document.createElementNS('http://www.w3.org/2000/svg', 'title');
  title.textContent = edgeTitle(xl);
  path.appendChild(title);
  svg.appendChild(path);

  for (const [x, y] of [
    [x1, y1],
    [x2, y2],
  ] as const) {
    const dot = document.createElementNS('http://www.w3.org/2000/svg', 'circle');
    dot.setAttribute('class', selected ? 'xlink-site xlink-site--selected' : 'xlink-site');
    dot.setAttribute('cx', String(x));
    dot.setAttribute('cy', String(y));
    dot.setAttribute('r', '4');
    dot.addEventListener('click', onClick);
    svg.appendChild(dot);
  }
}

function residueX(layout: LayoutProtein, absPos: number): number {
  const len = layout.protein.sequence.length;
  const frac = len > 1 ? (absPos - 1) / (len - 1) : 0.5;
  return MARGIN.left + frac * layout.barWidth;
}

function edgeTitle(xl: ViewerCrosslink): string {
  const parts = [
    `${xl.protein1}:${xl.abs_pos1} ↔ ${xl.protein2}:${xl.abs_pos2}`,
    `score ${xl.score.toFixed(2)}`,
  ];
  if (xl.glycan_composition) parts.push(xl.glycan_composition);
  return parts.join(' · ');
}

function truncateId(id: string, max: number): string {
  return id.length > max ? id.slice(0, max - 1) + '…' : id;
}
