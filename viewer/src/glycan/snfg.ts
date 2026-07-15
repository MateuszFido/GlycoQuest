export interface GlycanComponent {
  name: string;
  count: number;
  shape: SnfgShape;
  color: string;
}

export type SnfgShape =
  | 'circle'
  | 'square'
  | 'diamond'
  | 'triangle'
  | 'star'
  | 'crossed_square'
  | 'divided_diamond'
  | 'flat_hexagon';

export interface RenderOptions {
  size?: number;
}

const COMPONENT_RE = /([A-Za-z]+)\((\d+)\)/g;

const SNFG: Record<string, { shape: SnfgShape; color: string }> = {
  Hex: { shape: 'circle', color: '#00A651' },
  Glc: { shape: 'circle', color: '#0072BC' },
  Man: { shape: 'circle', color: '#00A651' },
  Gal: { shape: 'circle', color: '#FFD400' },
  HexNAc: { shape: 'square', color: '#0072BC' },
  GlcNAc: { shape: 'square', color: '#0072BC' },
  GalNAc: { shape: 'square', color: '#FFD400' },
  Fuc: { shape: 'triangle', color: '#ED1C24' },
  dHex: { shape: 'triangle', color: '#ED1C24' },
  NeuAc: { shape: 'diamond', color: '#A54399' },
  NeuGc: { shape: 'diamond', color: '#8FCCE0' },
  Kdn: { shape: 'diamond', color: '#00A651' },
  Pent: { shape: 'star', color: '#F47920' },
  Xyl: { shape: 'star', color: '#F47920' },
};

export function parseGlycanComposition(composition: string): GlycanComponent[] {
  const components: GlycanComponent[] = [];
  for (const match of composition.matchAll(COMPONENT_RE)) {
    const name = match[1];
    const symbol = SNFG[name] ?? { shape: 'flat_hexagon' as const, color: '#64748b' };
    components.push({
      name,
      count: Number(match[2]),
      shape: symbol.shape,
      color: symbol.color,
    });
  }
  return components;
}

export function renderGlycanSvg(composition: string, options: RenderOptions = {}): string {
  const size = options.size ?? 18;
  const gap = 5;
  const components = parseGlycanComposition(composition);
  if (components.length === 0) {
    return `<span class="gq-glycan-text">${escapeHtml(composition)}</span>`;
  }
  const width = components.reduce(
    (sum, component) => sum + size + (component.count > 1 ? 9 : 0) + gap,
    0,
  );
  const symbols = components
    .map((component, index) => {
      const x = indexOffset(components, index, size, gap);
      return `${shapeMarkup(component.shape, x, 1, size, component.color)}${countMarkup(component, x, size)}`;
    })
    .join('');
  return `<svg class="gq-glycan-svg" viewBox="0 0 ${width} ${size + 4}" width="${width}" height="${size + 4}" role="img" aria-label="${escapeHtml(composition)}"><title>${escapeHtml(composition)}</title>${symbols}</svg>`;
}

function indexOffset(components: GlycanComponent[], index: number, size: number, gap: number): number {
  let x = 0;
  for (let i = 0; i < index; i++) {
    x += size + (components[i].count > 1 ? 9 : 0) + gap;
  }
  return x;
}

function shapeMarkup(shape: SnfgShape, x: number, y: number, size: number, color: string): string {
  const stroke = '#333333';
  const mid = x + size / 2;
  const bottom = y + size;
  if (shape === 'circle') {
    return `<circle cx="${mid}" cy="${y + size / 2}" r="${size / 2 - 1}" fill="${color}" stroke="${stroke}"/>`;
  }
  if (shape === 'square') {
    return `<rect x="${x + 1}" y="${y + 1}" width="${size - 2}" height="${size - 2}" fill="${color}" stroke="${stroke}"/>`;
  }
  if (shape === 'diamond') {
    return `<path d="M ${mid} ${y + 1} L ${x + size - 1} ${y + size / 2} L ${mid} ${bottom - 1} L ${x + 1} ${y + size / 2} Z" fill="${color}" stroke="${stroke}"/>`;
  }
  if (shape === 'triangle') {
    return `<path d="M ${mid} ${y + 1} L ${x + size - 1} ${bottom - 1} L ${x + 1} ${bottom - 1} Z" fill="${color}" stroke="${stroke}"/>`;
  }
  if (shape === 'star') {
    return starPath(mid, y + size / 2, size / 2 - 1, color, stroke);
  }
  if (shape === 'crossed_square') {
    return `<rect x="${x + 1}" y="${y + 1}" width="${size - 2}" height="${size - 2}" fill="${color}" stroke="${stroke}"/><line x1="${x + 1}" y1="${y + 1}" x2="${x + size - 1}" y2="${bottom - 1}" stroke="${stroke}"/><line x1="${x + size - 1}" y1="${y + 1}" x2="${x + 1}" y2="${bottom - 1}" stroke="${stroke}"/>`;
  }
  if (shape === 'divided_diamond') {
    return `<path d="M ${mid} ${y + 1} L ${x + size - 1} ${y + size / 2} L ${mid} ${y + size / 2} L ${x + 1} ${y + size / 2} Z" fill="${color}" stroke="${stroke}"/><path d="M ${x + 1} ${y + size / 2} L ${mid} ${bottom - 1} L ${x + size - 1} ${y + size / 2} Z" fill="#FFFFFF" stroke="${stroke}"/>`;
  }
  return `<path d="M ${x + 4} ${y + 1} H ${x + size - 4} L ${x + size - 1} ${y + size / 2} L ${x + size - 4} ${bottom - 1} H ${x + 4} L ${x + 1} ${y + size / 2} Z" fill="${color}" stroke="${stroke}"/>`;
}

function starPath(cx: number, cy: number, rOuter: number, color: string, stroke: string): string {
  const rInner = rOuter * 0.38;
  const points: string[] = [];
  for (let j = 0; j < 10; j++) {
    const angle = ((j * 36 - 90) * Math.PI) / 180;
    const r = j % 2 === 0 ? rOuter : rInner;
    points.push(`${cx + r * Math.cos(angle)},${cy + r * Math.sin(angle)}`);
  }
  return `<polygon points="${points.join(' ')}" fill="${color}" stroke="${stroke}"/>`;
}

function countMarkup(component: GlycanComponent, x: number, size: number): string {
  if (component.count <= 1) return '';
  return `<text x="${x + size - 1}" y="${size + 3}" class="gq-glycan-count">${component.count}</text>`;
}

function escapeHtml(value: string): string {
  return value
    .replace(/&/g, '&amp;')
    .replace(/</g, '&lt;')
    .replace(/>/g, '&gt;')
    .replace(/"/g, '&quot;');
}
