import type { SelectionStore } from '../store/selection';
import type { ViewerCrosslink, ViewerProtein } from '../types';

const SVG_NS = 'http://www.w3.org/2000/svg';
const LANE_LEFT = 136;
const LANE_RIGHT = 40;
const LANE_WIDTH = 760;
const LANE_HEIGHT = 12;
const TOP_LANE_Y = 92;
const BOTTOM_LANE_Y = 204;

interface Lane {
  protein: ViewerProtein;
  y: number;
}

interface Site {
  absPos: number | null;
  pepPos: number | null;
  peptide: string;
}

export function renderNetworkPanel(container: HTMLElement, store: SelectionStore): () => void {
  const panel = document.createElement('section');
  panel.className = 'gq-panel gq-panel--pair-map';
  panel.innerHTML = '<h2>Sequence Pair Map</h2>';
  const body = document.createElement('div');
  body.className = 'gq-panel-body';
  panel.appendChild(body);
  container.appendChild(panel);

  const render = () => {
    const selected = store.selectedCrosslink;
    if (!selected) {
      body.innerHTML = '<p class="gq-empty">Select a crosslink to view sequence context.</p>';
      return;
    }

    const proteins = focusedProteins(store, selected);
    if (proteins.length === 0 || !selected.mapped) {
      body.innerHTML =
        '<p class="gq-empty">No mapped sequence coordinates for the selected crosslink.</p>';
      return;
    }

    const pairCrosslinks = store.selectedPairCrosslinks.filter((xl) => xl.mapped);
    const svg = selected.protein1 === selected.protein2
      ? renderIntraproteinMap(proteins[0], pairCrosslinks, selected, store)
      : renderInterproteinMap(proteins, pairCrosslinks, selected, store);
    const detail = buildDetail(selected, pairCrosslinks.length);
    body.replaceChildren(svg, detail);
  };

  const unsub = store.subscribe(render);
  render();
  return () => {
    unsub();
    panel.remove();
  };
}

function renderInterproteinMap(
  proteins: ViewerProtein[],
  crosslinks: ViewerCrosslink[],
  selected: ViewerCrosslink,
  store: SelectionStore,
): SVGSVGElement {
  const lanes: Lane[] = [
    { protein: proteins[0], y: TOP_LANE_Y },
    { protein: proteins[1], y: BOTTOM_LANE_Y },
  ];
  const svg = baseSvg(BOTTOM_LANE_Y + 70);
  lanes.forEach((lane) => drawLane(svg, lane, crosslinks, selected));

  drawInterproteinLinks(svg, lanes, crosslinks, selected, store);
  return svg;
}

function renderIntraproteinMap(
  protein: ViewerProtein,
  crosslinks: ViewerCrosslink[],
  selected: ViewerCrosslink,
  store: SelectionStore,
): SVGSVGElement {
  const lane = { protein, y: BOTTOM_LANE_Y - 28 };
  const svg = baseSvg(BOTTOM_LANE_Y + 76);
  drawLane(svg, lane, crosslinks, selected);
  drawIntraproteinArcs(svg, lane, crosslinks, selected, store);
  return svg;
}

function baseSvg(height: number): SVGSVGElement {
  const width = LANE_LEFT + LANE_WIDTH + LANE_RIGHT;
  const svg = document.createElementNS(SVG_NS, 'svg');
  svg.setAttribute('class', 'gq-pair-svg');
  svg.setAttribute('viewBox', `0 0 ${width} ${height}`);
  svg.setAttribute('width', String(width));
  svg.setAttribute('height', String(height));
  svg.setAttribute('role', 'img');
  return svg;
}

function drawLane(
  svg: SVGSVGElement,
  lane: Lane,
  crosslinks: ViewerCrosslink[],
  selected: ViewerCrosslink,
): void {
  const y = lane.y;
  appendText(svg, lane.protein.display_name, 14, y + 9, 'gq-pair-label');
  appendText(svg, `${lane.protein.sequence.length} aa`, LANE_LEFT + LANE_WIDTH + 8, y + 9, 'gq-pair-length');

  for (const xl of crosslinks) {
    for (const coverage of peptideCoverage(lane.protein.id, xl)) {
      const start = residueX(lane.protein, coverage.start);
      const end = residueX(lane.protein, coverage.end);
      appendRect(svg, start, y - 8, Math.max(2, end - start), LANE_HEIGHT + 16, 'gq-peptide-band');
    }
  }

  appendRect(svg, LANE_LEFT, y, LANE_WIDTH, LANE_HEIGHT, 'gq-sequence-lane');

  for (const pos of xlinkPositions(lane.protein.id, crosslinks)) {
    const dot = appendCircle(svg, residueX(lane.protein, pos), y + LANE_HEIGHT / 2, 4, 'gq-site-dot');
    dot.appendChild(title(`Crosslink residue ${lane.protein.id}:${pos}`));
  }

  for (const marker of glycanMarkers(lane.protein.id, crosslinks)) {
    const x = residueX(lane.protein, marker.position);
    const diamond = document.createElementNS(SVG_NS, 'path');
    diamond.setAttribute('class', 'gq-glycan-marker');
    diamond.setAttribute('d', `M ${x} ${y - 16} l 6 6 l -6 6 l -6 -6 Z`);
    diamond.appendChild(title(marker.label));
    svg.appendChild(diamond);
  }

  const selectedPositions = selectedSitesForProtein(selected, lane.protein.id);
  for (const pos of selectedPositions) {
    appendCircle(svg, residueX(lane.protein, pos), y + LANE_HEIGHT / 2, 6, 'gq-site-dot gq-site-dot--selected');
  }
}

function drawInterproteinLinks(
  svg: SVGSVGElement,
  lanes: Lane[],
  crosslinks: ViewerCrosslink[],
  selected: ViewerCrosslink,
  store: SelectionStore,
): void {
  const draw = (xl: ViewerCrosslink) => {
    const topSite = siteForProtein(xl, lanes[0].protein.id);
    const bottomSite = siteForProtein(xl, lanes[1].protein.id);
    if (!topSite.absPos || !bottomSite.absPos) return;
    const x1 = residueX(lanes[0].protein, topSite.absPos);
    const y1 = lanes[0].y + LANE_HEIGHT / 2;
    const x2 = residueX(lanes[1].protein, bottomSite.absPos);
    const y2 = lanes[1].y + LANE_HEIGHT / 2;
    const line = document.createElementNS(SVG_NS, 'line');
    line.setAttribute('x1', String(x1));
    line.setAttribute('y1', String(y1));
    line.setAttribute('x2', String(x2));
    line.setAttribute('y2', String(y2));
    line.setAttribute('class', linkClass(xl, selected, 'gq-xlink-line'));
    line.addEventListener('click', () => store.selectCrosslink(xl.id));
    line.appendChild(title(linkTitle(xl)));
    svg.appendChild(line);
  };

  crosslinks.filter((xl) => xl.id !== selected.id).forEach(draw);
  draw(selected);
}

function drawIntraproteinArcs(
  svg: SVGSVGElement,
  lane: Lane,
  crosslinks: ViewerCrosslink[],
  selected: ViewerCrosslink,
  store: SelectionStore,
): void {
  const draw = (xl: ViewerCrosslink, index: number) => {
    const s1 = siteForProtein(xl, lane.protein.id, 1);
    const s2 = siteForProtein(xl, lane.protein.id, 2);
    if (!s1.absPos || !s2.absPos) return;
    const x1 = residueX(lane.protein, s1.absPos);
    const x2 = residueX(lane.protein, s2.absPos);
    const y = lane.y + LANE_HEIGHT / 2;
    const distance = Math.abs(x2 - x1);
    const arcHeight = Math.min(118, Math.max(36, distance * 0.34 + index * 9));
    const path = document.createElementNS(SVG_NS, 'path');
    path.setAttribute(
      'd',
      `M ${x1} ${y} C ${x1} ${y - arcHeight}, ${x2} ${y - arcHeight}, ${x2} ${y}`,
    );
    path.setAttribute('class', linkClass(xl, selected, 'gq-xlink-arc'));
    path.addEventListener('click', () => store.selectCrosslink(xl.id));
    path.appendChild(title(linkTitle(xl)));
    svg.appendChild(path);
  };

  crosslinks.filter((xl) => xl.id !== selected.id).forEach(draw);
  draw(selected, crosslinks.length);
}

function focusedProteins(store: SelectionStore, selected: ViewerCrosslink): ViewerProtein[] {
  const ids = selected.protein1 === selected.protein2
    ? [selected.protein1]
    : [selected.protein1, selected.protein2];
  return ids
    .slice(0, 2)
    .map((id) => store.bundle.proteins.find((protein) => protein.id === id))
    .filter((protein): protein is ViewerProtein => Boolean(protein));
}

function siteForProtein(xl: ViewerCrosslink, proteinId: string, preferredArm?: 1 | 2): Site {
  if (preferredArm !== 2 && xl.protein1 === proteinId) {
    return { absPos: xl.abs_pos1, pepPos: xl.pep_pos1, peptide: xl.pep_seq1 };
  }
  if (preferredArm !== 1 && xl.protein2 === proteinId) {
    return { absPos: xl.abs_pos2, pepPos: xl.pep_pos2, peptide: xl.pep_seq2 };
  }
  return { absPos: null, pepPos: null, peptide: '' };
}

function selectedSitesForProtein(xl: ViewerCrosslink, proteinId: string): number[] {
  const positions: number[] = [];
  if (xl.protein1 === proteinId && xl.abs_pos1 != null) positions.push(xl.abs_pos1);
  if (xl.protein2 === proteinId && xl.abs_pos2 != null) positions.push(xl.abs_pos2);
  return positions;
}

function peptideCoverage(proteinId: string, xl: ViewerCrosslink): Array<{ start: number; end: number }> {
  const sites = [siteForProtein(xl, proteinId, 1), siteForProtein(xl, proteinId, 2)];
  return sites.flatMap((site) => {
    if (!site.pepPos || !site.peptide) return [];
    return [{ start: site.pepPos, end: site.pepPos + site.peptide.length - 1 }];
  });
}

function xlinkPositions(proteinId: string, crosslinks: ViewerCrosslink[]): number[] {
  return Array.from(
    new Set(
      crosslinks.flatMap((xl) =>
        selectedSitesForProtein(xl, proteinId).filter((pos) => Number.isFinite(pos)),
      ),
    ),
  );
}

function glycanMarkers(
  proteinId: string,
  crosslinks: ViewerCrosslink[],
): Array<{ position: number; label: string }> {
  const markers: Array<{ position: number; label: string }> = [];
  for (const xl of crosslinks) {
    if (!xl.glyco_residue || !xl.glyco_peptide) continue;
    const site = xl.glyco_peptide === 1 ? siteForProtein(xl, proteinId, 1) : siteForProtein(xl, proteinId, 2);
    if (!site.pepPos || !site.peptide) continue;
    for (let i = 0; i < site.peptide.length; i++) {
      if (site.peptide[i].toUpperCase() === xl.glyco_residue.toUpperCase()) {
        markers.push({
          position: site.pepPos + i,
          label: `Glycan ${xl.glycan_composition ?? xl.glycan_name ?? 'site'}`,
        });
      }
    }
  }
  return markers;
}

function residueX(protein: ViewerProtein, absPos: number): number {
  const len = Math.max(protein.sequence.length, 1);
  const frac = len > 1 ? (absPos - 1) / (len - 1) : 0.5;
  return LANE_LEFT + Math.max(0, Math.min(1, frac)) * LANE_WIDTH;
}

function linkClass(xl: ViewerCrosslink, selected: ViewerCrosslink, base: string): string {
  const classes = [base];
  if (xl.glycan_composition) classes.push('gq-xlink--glycan');
  if (xl.postfilter_status !== 'pass') classes.push('gq-xlink--failed');
  if (xl.id === selected.id) classes.push('gq-xlink--selected');
  return classes.join(' ');
}

function linkTitle(xl: ViewerCrosslink): string {
  const parts = [
    `${xl.protein1}:${xl.abs_pos1 ?? '?'} <-> ${xl.protein2}:${xl.abs_pos2 ?? '?'}`,
    `score ${xl.score.toFixed(2)}`,
  ];
  if (xl.glycan_composition) parts.push(xl.glycan_composition);
  return parts.join(' | ');
}

function buildDetail(selected: ViewerCrosslink, pairCount: number): HTMLElement {
  const detail = document.createElement('div');
  detail.className = 'gq-detail';
  detail.textContent = [
    `${selected.protein1}:${selected.abs_pos1 ?? '?'} <-> ${selected.protein2}:${selected.abs_pos2 ?? '?'}`,
    `score ${selected.score.toFixed(2)}`,
    `scan ${selected.scan ?? '?'}`,
    selected.retention_time_min == null ? null : `RT ${selected.retention_time_min.toFixed(2)} min`,
    selected.glycan_composition ? `glycan ${selected.glycan_composition}` : null,
    `${pairCount} visible link(s) for this pair`,
  ]
    .filter((part): part is string => Boolean(part))
    .join(' | ');
  return detail;
}

function appendRect(
  svg: SVGSVGElement,
  x: number,
  y: number,
  width: number,
  height: number,
  className: string,
): SVGRectElement {
  const rect = document.createElementNS(SVG_NS, 'rect');
  rect.setAttribute('class', className);
  rect.setAttribute('x', String(x));
  rect.setAttribute('y', String(y));
  rect.setAttribute('width', String(width));
  rect.setAttribute('height', String(height));
  rect.setAttribute('rx', '4');
  svg.appendChild(rect);
  return rect;
}

function appendCircle(
  svg: SVGSVGElement,
  cx: number,
  cy: number,
  r: number,
  className: string,
): SVGCircleElement {
  const circle = document.createElementNS(SVG_NS, 'circle');
  circle.setAttribute('class', className);
  circle.setAttribute('cx', String(cx));
  circle.setAttribute('cy', String(cy));
  circle.setAttribute('r', String(r));
  svg.appendChild(circle);
  return circle;
}

function appendText(
  svg: SVGSVGElement,
  value: string,
  x: number,
  y: number,
  className: string,
): SVGTextElement {
  const text = document.createElementNS(SVG_NS, 'text');
  text.setAttribute('class', className);
  text.setAttribute('x', String(x));
  text.setAttribute('y', String(y));
  text.textContent = value;
  svg.appendChild(text);
  return text;
}

function title(value: string): SVGTitleElement {
  const titleEl = document.createElementNS(SVG_NS, 'title');
  titleEl.textContent = value;
  return titleEl;
}
