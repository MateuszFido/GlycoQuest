// Copyright (c) ETH Zurich, Mateusz Fido

import { renderGlycanSvg } from '../glycan/snfg';
import type { SelectionStore } from '../store/selection';
import type { ViewerCrosslink, ViewerProtein } from '../types';
import {
  buildFullSequenceLayout,
  displayProteinLabel,
  groupCrosslinksByEndpoint,
  residuePointForPosition,
  residueSegmentsForRange,
  selectedSitesForProtein,
  stableProteinIdsForSelection,
  type EdgeGroup,
  type FullSequenceLayout,
} from './pairMapLayout';

const SVG_NS = 'http://www.w3.org/2000/svg';
const LABEL_WIDTH = 168;
const LANE_LEFT = LABEL_WIDTH;
const LANE_RIGHT = 48;
const RESIDUE_PITCH = 14;
const RESIDUES_PER_ROW = 72;
const LANE_WIDTH = RESIDUES_PER_ROW * RESIDUE_PITCH;
const ROW_TRACK_HEIGHT = 22;
const ROW_HEIGHT = 30;
const LETTER_BASELINE_OFFSET = 16;
const MAP_TOP = 24;
const LANE_HEADER_HEIGHT = 42;
const LANE_GAP = 56;

interface Lane {
  protein: ViewerProtein;
  top: number;
  layout: FullSequenceLayout;
}

interface Site {
  absPos: number | null;
  pepPos: number | null;
  linkPos: number | null;
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

  let pendingChoiceGroup: EdgeGroup | null = null;
  let focusedPairKey: string | null = null;
  let focusedProteinIds: string[] = [];

  const render = () => {
    const selected = store.selectedCrosslink;
    if (!selected) {
      body.innerHTML = '<p class="gq-empty">Select a crosslink to view sequence context.</p>';
      return;
    }

    const focus = stableProteinIdsForSelection(selected, focusedPairKey, focusedProteinIds);
    focusedPairKey = focus.pairKey;
    focusedProteinIds = focus.ids;
    const proteins = proteinsById(store, focus.ids);
    if (proteins.length === 0 || proteins.length !== focus.ids.length || !selected.mapped) {
      body.innerHTML =
        '<p class="gq-empty">No mapped sequence coordinates for the selected crosslink.</p>';
      return;
    }

    const pairCrosslinks = store.selectedPairCrosslinks.filter((xl) => xl.mapped);
    const svg =
      focus.ids.length === 1
        ? renderIntraproteinMap(proteins[0], pairCrosslinks, selected, store, (group) => {
            pendingChoiceGroup = group;
            render();
          })
        : renderInterproteinMap(proteins, pairCrosslinks, selected, store, (group) => {
            pendingChoiceGroup = group;
            render();
          });

    if (
      pendingChoiceGroup &&
      !pendingChoiceGroup.crosslinks.some((xl) =>
        pairCrosslinks.some((visible) => visible.id === xl.id),
      )
    ) {
      pendingChoiceGroup = null;
    }

    const picker = pendingChoiceGroup ? buildEdgePicker(pendingChoiceGroup, store, () => {
      pendingChoiceGroup = null;
      render();
    }) : null;
    const detail = buildDetail(selected, pairCrosslinks.length, store);
    const glycanDetail = store.selectedGlycanComposition
      ? buildGlycanDetail(store.selectedGlycanComposition)
      : null;
    body.replaceChildren(
      ...([svg, picker, detail, glycanDetail].filter(Boolean) as Array<HTMLElement | SVGSVGElement>),
    );
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
  onStackedClick: (group: EdgeGroup) => void,
): SVGSVGElement {
  const topLane = buildLane(proteins[0], MAP_TOP);
  const bottomLane = buildLane(proteins[1], topLane.top + laneHeight(topLane) + LANE_GAP);
  const lanes: Lane[] = [topLane, bottomLane];
  const svg = baseSvg(bottomLane.top + laneHeight(bottomLane) + 28);
  lanes.forEach((lane) => drawLane(svg, lane, crosslinks, selected));
  drawInterproteinLinks(svg, lanes, crosslinks, selected, store, onStackedClick);
  lanes.forEach((lane) => drawLaneFeatures(svg, lane, crosslinks, selected, store));
  return svg;
}

function renderIntraproteinMap(
  protein: ViewerProtein,
  crosslinks: ViewerCrosslink[],
  selected: ViewerCrosslink,
  store: SelectionStore,
  onStackedClick: (group: EdgeGroup) => void,
): SVGSVGElement {
  const lane = buildLane(protein, MAP_TOP);
  const svg = baseSvg(lane.top + laneHeight(lane) + 28);
  drawLane(svg, lane, crosslinks, selected);
  drawIntraproteinArcs(svg, lane, crosslinks, selected, store, onStackedClick);
  drawLaneFeatures(svg, lane, crosslinks, selected, store);
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

function buildLane(protein: ViewerProtein, top: number): Lane {
  return {
    protein,
    top,
    layout: buildFullSequenceLayout(protein.sequence.length, {
      left: LANE_LEFT,
      top: top + LANE_HEADER_HEIGHT,
      residuesPerRow: RESIDUES_PER_ROW,
      residuePitch: RESIDUE_PITCH,
      rowHeight: ROW_HEIGHT,
      rowTrackHeight: ROW_TRACK_HEIGHT,
    }),
  };
}

function laneHeight(lane: Lane): number {
  return LANE_HEADER_HEIGHT + lane.layout.height;
}

function drawLane(
  svg: SVGSVGElement,
  lane: Lane,
  crosslinks: ViewerCrosslink[],
  selected: ViewerCrosslink,
): void {
  const label = displayProteinLabel(lane.protein.display_name);
  const labelText = appendText(svg, label.short, 14, lane.top + 15, 'gq-pair-label');
  labelText.appendChild(title(label.full));
  appendText(svg, `${lane.protein.sequence.length} aa`, 14, lane.top + 33, 'gq-pair-length');

  drawSequenceRows(svg, lane);
  drawPeptideBands(
    svg,
    lane,
    crosslinks.filter((xl) => xl.id !== selected.id),
    'gq-peptide-band',
  );
  drawPeptideBands(svg, lane, [selected], 'gq-peptide-band gq-peptide-band--selected');
  drawResidueLetters(svg, lane, selected, crosslinks);
}

function drawSequenceRows(svg: SVGSVGElement, lane: Lane): void {
  for (let row = 0; row < lane.layout.rowCount; row++) {
    const rowStart = row * lane.layout.residuesPerRow + 1;
    const rowEnd = Math.min(lane.protein.sequence.length, rowStart + lane.layout.residuesPerRow - 1);
    const rowLength = Math.max(1, rowEnd - rowStart + 1);
    const y = lane.layout.top + row * lane.layout.rowHeight;
    const rowIndex = appendText(svg, String(rowStart), LANE_LEFT - 10, y + LETTER_BASELINE_OFFSET, 'gq-residue-row-index');
    rowIndex.setAttribute('text-anchor', 'end');
    appendRect(
      svg,
      lane.layout.left,
      y,
      rowLength * lane.layout.residuePitch,
      lane.layout.rowTrackHeight,
      'gq-residue-track',
    );
  }
}

function drawPeptideBands(
  svg: SVGSVGElement,
  lane: Lane,
  crosslinks: ViewerCrosslink[],
  className: string,
): void {
  for (const xl of crosslinks) {
    for (const coverage of peptideCoverage(lane.protein.id, xl)) {
      for (const segment of residueSegmentsForRange(lane.layout, coverage.start, coverage.end)) {
        appendRect(
          svg,
          segment.x,
          segment.y,
          Math.max(2, segment.width),
          lane.layout.rowTrackHeight,
          className,
        );
      }
    }
  }
}

function drawLaneFeatures(
  svg: SVGSVGElement,
  lane: Lane,
  crosslinks: ViewerCrosslink[],
  selected: ViewerCrosslink,
  store: SelectionStore,
): void {
  for (const pos of xlinkPositions(lane.protein.id, crosslinks)) {
    const point = residuePointForPosition(lane.layout, pos);
    if (!point) continue;
    const dot = appendCircle(svg, point.x, point.y, 4, 'gq-site-dot');
    dot.appendChild(title(`Crosslink residue ${lane.protein.id}:${pos}`));
  }

  for (const marker of monolinkMarkers(lane.protein.id, crosslinks)) {
    const point = residuePointForPosition(lane.layout, marker.position);
    if (!point) continue;
    const diamond = document.createElementNS(SVG_NS, 'path');
    diamond.setAttribute('class', 'gq-monolink-marker');
    diamond.setAttribute(
      'd',
      `M ${point.x} ${point.y - 7} l 7 7 l -7 7 l -7 -7 Z`,
    );
    diamond.setAttribute('tabindex', '0');
    diamond.setAttribute('role', 'button');
    diamond.appendChild(title(marker.label));
    const select = () => {
      const target = marker.crosslinks.find((xl) => xl.id === selected.id) ?? marker.crosslinks[0];
      store.selectCrosslink(target.id);
    };
    diamond.addEventListener('click', (event) => {
      event.stopPropagation();
      select();
    });
    diamond.addEventListener('keydown', (event) => {
      if (event.key === 'Enter' || event.key === ' ') {
        event.preventDefault();
        select();
      }
    });
    svg.appendChild(diamond);
  }

  for (const marker of glycanMarkers(lane.protein.id, crosslinks)) {
    const point = residuePointForPosition(lane.layout, marker.position);
    if (!point) continue;
    const diamond = document.createElementNS(SVG_NS, 'path');
    diamond.setAttribute('class', 'gq-glycan-marker');
    diamond.setAttribute(
      'd',
      `M ${point.x} ${point.y - 16} l 6 6 l -6 6 l -6 -6 Z`,
    );
    diamond.setAttribute('tabindex', '0');
    diamond.setAttribute('role', 'button');
    diamond.appendChild(title(marker.label));
    diamond.addEventListener('click', (event) => {
      event.stopPropagation();
      store.selectGlycan(marker.composition);
    });
    diamond.addEventListener('keydown', (event) => {
      if (event.key === 'Enter' || event.key === ' ') {
        event.preventDefault();
        store.selectGlycan(marker.composition);
      }
    });
    svg.appendChild(diamond);
  }

  const selectedPositions = selectedSitesForProtein(selected, lane.protein.id);
  for (const pos of selectedPositions) {
    const point = residuePointForPosition(lane.layout, pos);
    if (!point) continue;
    appendCircle(svg, point.x, point.y, 6, 'gq-site-dot gq-site-dot--selected');
  }
}

function drawResidueLetters(
  svg: SVGSVGElement,
  lane: Lane,
  selected: ViewerCrosslink,
  crosslinks: ViewerCrosslink[],
): void {
  const xlinkSet = new Set(xlinkPositions(lane.protein.id, crosslinks));
  const monolinkSet = new Set(monolinkMarkers(lane.protein.id, crosslinks).map((marker) => marker.position));
  const selectedSet = new Set(selectedSitesForProtein(selected, lane.protein.id));
  for (let i = 0; i < lane.protein.sequence.length; i++) {
    const pos = i + 1;
    const point = residuePointForPosition(lane.layout, pos);
    if (!point) continue;
    const classes = ['gq-residue-letter'];
    if (selectedSet.has(pos)) classes.push('gq-residue-letter--selected');
    else if (monolinkSet.has(pos)) classes.push('gq-residue-letter--monolink');
    else if (xlinkSet.has(pos)) classes.push('gq-residue-letter--xlink');
    const text = appendText(
      svg,
      lane.protein.sequence[i],
      point.x,
      point.trackY + LETTER_BASELINE_OFFSET,
      classes.join(' '),
    );
    text.setAttribute('text-anchor', 'middle');
    text.appendChild(title(`${lane.protein.id}:${pos}`));
  }
}

function drawInterproteinLinks(
  svg: SVGSVGElement,
  lanes: Lane[],
  crosslinks: ViewerCrosslink[],
  selected: ViewerCrosslink,
  store: SelectionStore,
  onStackedClick: (group: EdgeGroup) => void,
): void {
  const groups = groupCrosslinksByEndpoint(crosslinks, [
    lanes[0].protein.id,
    lanes[1].protein.id,
  ]);
  const ordered = [
    ...groups.filter((group) => !group.crosslinks.some((xl) => xl.id === selected.id)),
    ...groups.filter((group) => group.crosslinks.some((xl) => xl.id === selected.id)),
  ];

  for (const group of ordered) {
    const representative =
      group.crosslinks.find((xl) => xl.id === selected.id) ?? group.crosslinks[0];
    const topSite = siteForProtein(representative, lanes[0].protein.id);
    const bottomSite = siteForProtein(representative, lanes[1].protein.id);
    if (!topSite.absPos || !bottomSite.absPos) continue;
    const topPoint = residuePointForPosition(lanes[0].layout, topSite.absPos);
    const bottomPoint = residuePointForPosition(lanes[1].layout, bottomSite.absPos);
    if (!topPoint || !bottomPoint) continue;
    const line = document.createElementNS(SVG_NS, 'line');
    line.setAttribute('x1', String(topPoint.x));
    line.setAttribute('y1', String(topPoint.y));
    line.setAttribute('x2', String(bottomPoint.x));
    line.setAttribute('y2', String(bottomPoint.y));
    line.setAttribute('class', linkClass(representative, selected, 'gq-xlink-line'));
    line.appendChild(title(linkTitle(representative, group.crosslinks.length)));
    attachEdgeHandlers(svg, line, group, store, onStackedClick);
  }
}

function drawIntraproteinArcs(
  svg: SVGSVGElement,
  lane: Lane,
  crosslinks: ViewerCrosslink[],
  selected: ViewerCrosslink,
  store: SelectionStore,
  onStackedClick: (group: EdgeGroup) => void,
): void {
  const groups = groupCrosslinksByEndpoint(crosslinks, [lane.protein.id]);
  const ordered = [
    ...groups.filter((group) => !group.crosslinks.some((xl) => xl.id === selected.id)),
    ...groups.filter((group) => group.crosslinks.some((xl) => xl.id === selected.id)),
  ];

  ordered.forEach((group, index) => {
    const representative =
      group.crosslinks.find((xl) => xl.id === selected.id) ?? group.crosslinks[0];
    const s1 = siteForProtein(representative, lane.protein.id, 1);
    const s2 = siteForProtein(representative, lane.protein.id, 2);
    if (!s1.absPos || !s2.absPos) return;
    const p1 = residuePointForPosition(lane.layout, s1.absPos);
    const p2 = residuePointForPosition(lane.layout, s2.absPos);
    if (!p1 || !p2) return;
    const distance = Math.hypot(p2.x - p1.x, p2.y - p1.y);
    const arcHeight = Math.min(86, Math.max(28, distance * 0.16 + index * 8));
    const controlY =
      p1.row === p2.row
        ? Math.max(lane.layout.top - 12, p1.y - arcHeight)
        : (p1.y + p2.y) / 2;
    const path = document.createElementNS(SVG_NS, 'path');
    path.setAttribute(
      'd',
      `M ${p1.x} ${p1.y} C ${p1.x} ${controlY}, ${p2.x} ${controlY}, ${p2.x} ${p2.y}`,
    );
    path.setAttribute('class', linkClass(representative, selected, 'gq-xlink-arc'));
    path.appendChild(title(linkTitle(representative, group.crosslinks.length)));
    attachEdgeHandlers(svg, path, group, store, onStackedClick);
  });
}

function attachEdgeHandlers(
  svg: SVGSVGElement,
  shape: SVGLineElement | SVGPathElement,
  group: EdgeGroup,
  store: SelectionStore,
  onStackedClick: (group: EdgeGroup) => void,
): void {
  const onClick = () => {
    if (group.crosslinks.length === 1) {
      store.selectCrosslink(group.crosslinks[0].id);
      return;
    }
    onStackedClick(group);
  };
  const hit = shape.cloneNode(false) as SVGLineElement | SVGPathElement;
  hit.setAttribute('class', 'gq-xlink-hit-target');
  hit.addEventListener('click', (event) => {
    event.stopPropagation();
    onClick();
  });
  shape.addEventListener('click', (event) => {
    event.stopPropagation();
    onClick();
  });
  svg.appendChild(hit);
  svg.appendChild(shape);
}

function buildEdgePicker(
  group: EdgeGroup,
  store: SelectionStore,
  onClose: () => void,
): HTMLElement {
  const picker = document.createElement('div');
  picker.className = 'gq-edge-picker';
  const heading = document.createElement('div');
  heading.className = 'gq-edge-picker-title';
  heading.textContent = `${group.crosslinks.length} crosslinks at this endpoint — choose one:`;
  picker.appendChild(heading);
  for (const xl of group.crosslinks) {
    const button = document.createElement('button');
    button.type = 'button';
    button.className = 'gq-edge-choice';
    if (xl.id === store.selectedCrosslinkId) button.classList.add('gq-edge-choice--selected');
    button.textContent = [
      `scan ${xl.scan ?? '?'}`,
      xl.glycan_composition ?? 'no glycan',
      `score ${xl.score.toFixed(2)}`,
      `${xl.protein1}:${xl.abs_pos1 ?? '?'} -> ${xl.protein2}:${xl.abs_pos2 ?? '?'}`,
    ].join(' | ');
    button.addEventListener('click', () => {
      store.selectCrosslink(xl.id);
      onClose();
    });
    picker.appendChild(button);
  }
  return picker;
}

function proteinsById(store: SelectionStore, ids: string[]): ViewerProtein[] {
  return ids
    .slice(0, 2)
    .map((id) => store.bundle.proteins.find((protein) => protein.id === id))
    .filter((protein): protein is ViewerProtein => Boolean(protein));
}

function siteForProtein(xl: ViewerCrosslink, proteinId: string, preferredArm?: 1 | 2): Site {
  if (preferredArm !== 2 && xl.protein1 === proteinId) {
    return {
      absPos: xl.abs_pos1,
      pepPos: xl.pep_pos1,
      linkPos: xl.link_pos1,
      peptide: xl.pep_seq1,
    };
  }
  if (preferredArm !== 1 && xl.protein2 === proteinId) {
    return {
      absPos: xl.abs_pos2,
      pepPos: xl.pep_pos2,
      linkPos: xl.link_pos2,
      peptide: xl.pep_seq2,
    };
  }
  return { absPos: null, pepPos: null, linkPos: null, peptide: '' };
}

function peptideCoverage(
  proteinId: string,
  xl: ViewerCrosslink,
): Array<{ start: number; end: number }> {
  const sites = [siteForProtein(xl, proteinId, 1), siteForProtein(xl, proteinId, 2)];
  return sites.flatMap((site) => {
    const start = peptideStart(site);
    if (!start || !site.peptide) return [];
    return [{ start, end: start + site.peptide.length - 1 }];
  });
}

function xlinkPositions(proteinId: string, crosslinks: ViewerCrosslink[]): number[] {
  return Array.from(
    new Set(
      crosslinks.filter((xl) => !isMonolink(xl)).flatMap((xl) =>
        selectedSitesForProtein(xl, proteinId).filter((pos) => Number.isFinite(pos)),
      ),
    ),
  );
}

function monolinkMarkers(
  proteinId: string,
  crosslinks: ViewerCrosslink[],
): Array<{ position: number; label: string; crosslinks: ViewerCrosslink[] }> {
  const groups = new Map<number, ViewerCrosslink[]>();
  for (const xl of crosslinks.filter(isMonolink)) {
    for (const pos of selectedSitesForProtein(xl, proteinId)) {
      const group = groups.get(pos) ?? [];
      group.push(xl);
      groups.set(pos, group);
    }
  }
  return Array.from(groups.entries()).map(([position, grouped]) => ({
    position,
    label: `Monolink ${proteinId}:${position}${grouped.length > 1 ? ` (${grouped.length} hits)` : ''}`,
    crosslinks: grouped,
  }));
}

function isMonolink(xl: ViewerCrosslink): boolean {
  return xl.link_type === 'monolink' || !xl.protein2;
}

function glycanMarkers(
  proteinId: string,
  crosslinks: ViewerCrosslink[],
): Array<{ position: number; label: string; composition: string }> {
  const markers = new Map<number, { position: number; compositions: string[] }>();
  const addMarker = (position: number, composition: string) => {
    const marker = markers.get(position) ?? { position, compositions: [] };
    if (!marker.compositions.includes(composition)) marker.compositions.push(composition);
    markers.set(position, marker);
  };

  for (const xl of crosslinks) {
    if (!xl.glycan_composition) continue;
    if (xl.glyco_sites.length > 0) {
      for (const glycoSite of xl.glyco_sites) {
        const site = siteForProtein(xl, proteinId, glycoSite.peptide === 1 ? 1 : 2);
        const start = peptideStart(site);
        if (start && glycoSite.peptide_position > 0) {
          addMarker(start + glycoSite.peptide_position - 1, xl.glycan_composition);
        }
      }
      continue;
    }
    if (xl.glyco_residue && xl.glyco_peptide) {
      const site =
        xl.glyco_peptide === 1
          ? siteForProtein(xl, proteinId, 1)
          : siteForProtein(xl, proteinId, 2);
      const start = peptideStart(site);
      if (!start || !site.peptide) continue;
      for (let i = 0; i < site.peptide.length; i++) {
        if (site.peptide[i].toUpperCase() === xl.glyco_residue.toUpperCase()) {
          addMarker(start + i, xl.glycan_composition);
        }
      }
      continue;
    }
    for (const pos of selectedSitesForProtein(xl, proteinId)) {
      addMarker(pos, xl.glycan_composition);
    }
  }
  return Array.from(markers.values()).map((marker) => ({
    position: marker.position,
    label: `Glycan ${marker.compositions.join(', ')}`,
    composition: marker.compositions[0],
  }));
}

function peptideStart(site: Site): number | null {
  if (site.pepPos != null) return site.pepPos;
  if (site.absPos != null && site.linkPos != null) return site.absPos - site.linkPos + 1;
  return null;
}

function linkClass(xl: ViewerCrosslink, selected: ViewerCrosslink, base: string): string {
  const classes = [base];
  if (xl.glycan_composition) classes.push('gq-xlink--glycan');
  if (xl.postfilter_status !== 'pass') classes.push('gq-xlink--failed');
  if (xl.id === selected.id) classes.push('gq-xlink--selected');
  return classes.join(' ');
}

function linkTitle(xl: ViewerCrosslink, stackCount: number): string {
  const parts = [
    `${xl.protein1}:${xl.abs_pos1 ?? '?'} <-> ${xl.protein2}:${xl.abs_pos2 ?? '?'}`,
    `score ${xl.score.toFixed(2)}`,
  ];
  if (xl.glycan_composition) parts.push(xl.glycan_composition);
  if (stackCount > 1) parts.push(`${stackCount} hits at endpoint`);
  return parts.join(' | ');
}

function buildDetail(
  selected: ViewerCrosslink,
  pairCount: number,
  store: SelectionStore,
): HTMLElement {
  const detail = document.createElement('div');
  detail.className = 'gq-detail';
  const meta = document.createElement('div');
  meta.textContent = [
    `${selected.protein1}:${selected.abs_pos1 ?? '?'} <-> ${selected.protein2}:${selected.abs_pos2 ?? '?'}`,
    `score ${selected.score.toFixed(2)}`,
    `scan ${selected.scan ?? '?'}`,
    selected.retention_time_min == null
      ? null
      : `scan_time ${selected.retention_time_min.toFixed(2)} min`,
    `${pairCount} visible link(s) for this pair`,
  ]
    .filter((part): part is string => Boolean(part))
    .join(' | ');
  detail.appendChild(meta);

  if (selected.glycan_composition) {
    const button = document.createElement('button');
    button.type = 'button';
    button.className = 'gq-glycan-chip';
    button.innerHTML = renderGlycanSvg(selected.glycan_composition, { size: 16 });
    button.appendChild(document.createTextNode(` ${selected.glycan_composition}`));
    button.addEventListener('click', () => store.selectGlycan(selected.glycan_composition));
    detail.appendChild(button);
  }
  return detail;
}

function buildGlycanDetail(composition: string): HTMLElement {
  const detail = document.createElement('div');
  detail.className = 'gq-glycan-detail';
  detail.innerHTML = renderGlycanSvg(composition, { size: 26 });
  const label = document.createElement('span');
  label.textContent = composition;
  detail.appendChild(label);
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
