import type { SelectionStore } from '../store/selection';
import type { ViewerCrosslink } from '../types';

export function renderSequencePanel(container: HTMLElement, store: SelectionStore): () => void {
  const panel = document.createElement('section');
  panel.className = 'gq-panel';
  panel.innerHTML = '<h2>Sequence map</h2>';
  const body = document.createElement('div');
  body.className = 'gq-panel-body';
  panel.appendChild(body);
  container.appendChild(panel);

  const render = () => {
    const proteinId = store.selectedProteinId;
    const protein = store.bundle.proteins.find((p) => p.id === proteinId);
    if (!protein) {
      body.innerHTML = '<p class="gq-empty">Select a protein to view sequence context.</p>';
      return;
    }

    const crosslinks = store.visibleCrosslinks.filter(
      (xl) => xl.protein1 === proteinId || xl.protein2 === proteinId,
    );
    const selected = store.selectedCrosslink;

    const xlinkPositions = new Set<number>();
    const glycoPositions = new Map<number, string>();
    for (const xl of crosslinks) {
      collectPositions(xl, proteinId!, xlinkPositions, glycoPositions);
    }

    let selectedPositions = new Set<number>();
    if (selected) {
      selectedPositions = new Set(
        [selected.abs_pos1, selected.abs_pos2].filter(
          (p): p is number => p != null && isOnProtein(selected, proteinId!, p),
        ),
      );
    }

    const scroll = document.createElement('div');
    scroll.className = 'gq-seq-scroll';
    const groups = document.createElement('div');
    groups.className = 'gq-seq-groups';

    const seq = protein.sequence;
    for (let i = 0; i < seq.length; i += 10) {
      const group = document.createElement('div');
      group.className = 'gq-seq-group';
      const idx = document.createElement('div');
      idx.className = 'gq-seq-idx';
      idx.textContent = String(i + 1);
      group.appendChild(idx);

      const residues = document.createElement('div');
      residues.className = 'gq-seq-residues';
      for (let j = i; j < Math.min(i + 10, seq.length); j++) {
        const pos = j + 1;
        const res = document.createElement('span');
        res.className = 'gq-residue';
        res.textContent = seq[j];
        if (xlinkPositions.has(pos)) res.classList.add('gq-residue--xlink');
        else if (crosslinks.some((xl) => peptideCovers(xl, proteinId!, pos))) {
          res.classList.add('gq-residue--covered');
        }
        if (glycoPositions.has(pos)) res.classList.add('gq-residue--glyco');
        if (selectedPositions.has(pos)) res.classList.add('gq-residue--selected');
        res.title = glycoPositions.get(pos) ?? `Residue ${pos}`;
        residues.appendChild(res);
      }
      group.appendChild(residues);
      groups.appendChild(group);
    }
    scroll.appendChild(groups);

    const detail = document.createElement('div');
    detail.className = 'gq-detail';
    if (selected && (selected.protein1 === proteinId || selected.protein2 === proteinId)) {
      detail.innerHTML = formatCrosslinkDetail(selected, proteinId!);
    } else {
      detail.textContent = `${crosslinks.length} crosslink(s) on ${protein.display_name} · ${seq.length} residues`;
    }

    body.replaceChildren(scroll, detail);
  };

  const unsub = store.subscribe(render);
  render();
  return () => {
    unsub();
    panel.remove();
  };
}

function collectPositions(
  xl: ViewerCrosslink,
  proteinId: string,
  xlinks: Set<number>,
  glycos: Map<number, string>,
): void {
  if (xl.protein1 === proteinId && xl.abs_pos1) xlinks.add(xl.abs_pos1);
  if (xl.protein2 === proteinId && xl.abs_pos2) xlinks.add(xl.abs_pos2);

  if (xl.glyco_residue && xl.glyco_peptide) {
    const onP1 = xl.glyco_peptide === 1 && xl.protein1 === proteinId;
    const onP2 = xl.glyco_peptide === 2 && xl.protein2 === proteinId;
    if (onP1 || onP2) {
      const pep = onP1 ? xl.pep_seq1 : xl.pep_seq2;
      const pepPos = onP1 ? xl.pep_pos1 : xl.pep_pos2;
      if (pepPos && pep) {
        for (let i = 0; i < pep.length; i++) {
          if (pep[i].toUpperCase() === xl.glyco_residue.toUpperCase()) {
            glycos.set(pepPos + i, `Glycan: ${xl.glycan_composition ?? xl.glycan_name ?? '?'}`);
          }
        }
      }
    }
  }
}

function peptideCovers(xl: ViewerCrosslink, proteinId: string, pos: number): boolean {
  if (xl.protein1 === proteinId && xl.pep_pos1 && xl.pep_seq1) {
    const end = xl.pep_pos1 + xl.pep_seq1.length - 1;
    if (pos >= xl.pep_pos1 && pos <= end) return true;
  }
  if (xl.protein2 === proteinId && xl.pep_pos2 && xl.pep_seq2) {
    const end = xl.pep_pos2 + xl.pep_seq2.length - 1;
    if (pos >= xl.pep_pos2 && pos <= end) return true;
  }
  return false;
}

function isOnProtein(xl: ViewerCrosslink, proteinId: string, pos: number): boolean {
  if (xl.protein1 === proteinId && xl.abs_pos1 === pos) return true;
  if (xl.protein2 === proteinId && xl.abs_pos2 === pos) return true;
  return false;
}

function formatCrosslinkDetail(xl: ViewerCrosslink, proteinId: string): string {
  const parts = [
    `Score ${xl.score.toFixed(2)} · scan ${xl.scan ?? '?'}`,
    `P1 ${xl.protein1} ${xl.pep_seq1} @${xl.abs_pos1 ?? '?'}`,
    `P2 ${xl.protein2} ${xl.pep_seq2} @${xl.abs_pos2 ?? '?'}`,
  ];
  if (xl.glycan_composition) parts.push(`Glycan ${xl.glycan_composition}`);
  if (xl.glyco_residue) parts.push(`Site pep${xl.glyco_peptide}:${xl.glyco_residue}`);
  parts.push(`Viewing ${proteinId}`);
  return parts.join(' · ');
}
