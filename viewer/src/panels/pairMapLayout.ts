// Copyright (c) ETH Zurich, Mateusz Fido

import type { ViewerCrosslink, ViewerProtein } from '../types';

export interface ProteinLabel {
  short: string;
  full: string;
}

export interface EdgeGroup {
  key: string;
  crosslinks: ViewerCrosslink[];
}

export interface FullSequenceLayoutOptions {
  left: number;
  top: number;
  residuesPerRow: number;
  residuePitch: number;
  rowHeight: number;
  rowTrackHeight: number;
}

export interface FullSequenceLayout extends FullSequenceLayoutOptions {
  sequenceLength: number;
  rowCount: number;
  width: number;
  height: number;
}

export interface ResiduePoint {
  x: number;
  y: number;
  row: number;
  column: number;
  trackY: number;
}

export interface ResidueSegment {
  start: number;
  end: number;
  x: number;
  y: number;
  width: number;
  row: number;
}

export interface StableProteinSelection {
  pairKey: string;
  ids: string[];
}

export function displayProteinLabel(id: string): ProteinLabel {
  const parts = id.split('|');
  const short = parts.length >= 3 && parts[2] ? parts[2] : id;
  return { short, full: id };
}

export function groupCrosslinksByEndpoint(
  crosslinks: ViewerCrosslink[],
  laneProteinIds: string[],
): EdgeGroup[] {
  const groups = new Map<string, ViewerCrosslink[]>();
  for (const xl of crosslinks) {
    const key = endpointKey(xl, laneProteinIds);
    const group = groups.get(key) ?? [];
    group.push(xl);
    groups.set(key, group);
  }
  return Array.from(groups.entries()).map(([key, grouped]) => ({ key, crosslinks: grouped }));
}

export function endpointKey(xl: ViewerCrosslink, laneProteinIds: string[]): string {
  return laneProteinIds
    .map((proteinId) => {
      const positions = selectedSitesForProtein(xl, proteinId).join(',');
      return `${proteinId}:${positions || '?'}`;
    })
    .join('->');
}

export function selectedSitesForProtein(xl: ViewerCrosslink, proteinId: string): number[] {
  const positions: number[] = [];
  if (xl.protein1 === proteinId && xl.abs_pos1 != null) positions.push(xl.abs_pos1);
  if (xl.protein2 === proteinId && xl.abs_pos2 != null) positions.push(xl.abs_pos2);
  return positions;
}

export function buildFullSequenceLayout(
  sequenceLength: number,
  options: FullSequenceLayoutOptions,
): FullSequenceLayout {
  const residuesPerRow = Math.max(1, Math.floor(options.residuesPerRow));
  const length = Math.max(0, Math.floor(sequenceLength));
  const rowCount = Math.max(1, Math.ceil(length / residuesPerRow));
  const width = Math.min(Math.max(length, 1), residuesPerRow) * options.residuePitch;
  return {
    ...options,
    residuesPerRow,
    sequenceLength: length,
    rowCount,
    width,
    height: rowCount * options.rowHeight,
  };
}

export function residuePointForPosition(
  layout: FullSequenceLayout,
  absPos: number | null | undefined,
): ResiduePoint | null {
  if (absPos == null || !Number.isFinite(absPos)) return null;
  const pos = Math.floor(absPos);
  if (pos < 1 || pos > layout.sequenceLength) return null;
  const index = pos - 1;
  const row = Math.floor(index / layout.residuesPerRow);
  const column = index % layout.residuesPerRow;
  const trackY = layout.top + row * layout.rowHeight;
  return {
    x: layout.left + column * layout.residuePitch + layout.residuePitch / 2,
    y: trackY + layout.rowTrackHeight / 2,
    row,
    column,
    trackY,
  };
}

export function residueSegmentsForRange(
  layout: FullSequenceLayout,
  start: number | null | undefined,
  end: number | null | undefined,
): ResidueSegment[] {
  if (start == null || end == null || !Number.isFinite(start) || !Number.isFinite(end)) return [];
  const from = Math.max(1, Math.floor(Math.min(start, end)));
  const to = Math.min(layout.sequenceLength, Math.floor(Math.max(start, end)));
  if (from > to || layout.sequenceLength === 0) return [];

  const segments: ResidueSegment[] = [];
  let current = from;
  while (current <= to) {
    const row = Math.floor((current - 1) / layout.residuesPerRow);
    const rowEnd = Math.min(to, (row + 1) * layout.residuesPerRow);
    const startColumn = (current - 1) % layout.residuesPerRow;
    const endColumn = (rowEnd - 1) % layout.residuesPerRow;
    segments.push({
      start: current,
      end: rowEnd,
      x: layout.left + startColumn * layout.residuePitch,
      y: layout.top + row * layout.rowHeight,
      width: (endColumn - startColumn + 1) * layout.residuePitch,
      row,
    });
    current = rowEnd + 1;
  }
  return segments;
}

export function stableProteinIdsForSelection(
  selected: ViewerCrosslink,
  previousPairKey: string | null,
  previousIds: string[],
): StableProteinSelection {
  const pairKey = crosslinkPairKey(selected);
  const selectedIds = proteinIdsForCrosslink(selected);
  const previousCanBeReused =
    pairKey === previousPairKey &&
    previousIds.length === selectedIds.length &&
    previousIds.every((id) => selectedIds.includes(id)) &&
    selectedIds.every((id) => previousIds.includes(id));

  return {
    pairKey,
    ids: previousCanBeReused ? previousIds.slice() : selectedIds,
  };
}

function crosslinkPairKey(xl: ViewerCrosslink): string {
  if (xl.protein_pair_key) return xl.protein_pair_key;
  if (!xl.protein1) return xl.protein2;
  if (!xl.protein2) return xl.protein1;
  return xl.protein1 <= xl.protein2
    ? `${xl.protein1}|${xl.protein2}`
    : `${xl.protein2}|${xl.protein1}`;
}

function proteinIdsForCrosslink(xl: ViewerCrosslink): string[] {
  if (xl.link_type === 'monolink' || !xl.protein2) return [xl.protein1].filter(Boolean);
  return xl.protein1 === xl.protein2
    ? [xl.protein1]
    : [xl.protein1, xl.protein2].filter(Boolean);
}
