import { proteinPairKey } from '../data/normalize';
import type { ViewerBundle, ViewerCrosslink, ViewerFilters, ViewerListener } from '../types';

export class SelectionStore {
  bundle: ViewerBundle;
  selectedCrosslinkId: string | null = null;
  selectedProteinId: string | null = null;
  filters: ViewerFilters = { showFailed: false, proteinId: null, minScore: 0 };

  private listeners = new Set<ViewerListener>();

  constructor(bundle: ViewerBundle) {
    this.bundle = bundle;
  }

  subscribe(listener: ViewerListener): () => void {
    this.listeners.add(listener);
    return () => this.listeners.delete(listener);
  }

  private notify(): void {
    for (const l of this.listeners) l();
  }

  get visibleCrosslinks(): ViewerCrosslink[] {
    return this.bundle.crosslinks.filter((xl) => {
      if (!this.filters.showFailed && xl.postfilter_status !== 'pass') return false;
      if (xl.score < this.filters.minScore) return false;
      if (this.filters.proteinId) {
        return xl.protein1 === this.filters.proteinId || xl.protein2 === this.filters.proteinId;
      }
      return true;
    });
  }

  get selectedCrosslink(): ViewerCrosslink | null {
    if (!this.selectedCrosslinkId) return null;
    return this.bundle.crosslinks.find((xl) => xl.id === this.selectedCrosslinkId) ?? null;
  }

  get focusedProteinIds(): string[] {
    const selected = this.selectedCrosslink;
    if (!selected) return [];
    return selected.protein1 === selected.protein2
      ? [selected.protein1]
      : [selected.protein1, selected.protein2];
  }

  get selectedPairCrosslinks(): ViewerCrosslink[] {
    const selected = this.selectedCrosslink;
    if (!selected) return [];
    const selectedKey = selected.protein_pair_key ?? proteinPairKey(selected.protein1, selected.protein2);
    return this.visibleCrosslinks.filter((xl) => {
      const key = xl.protein_pair_key ?? proteinPairKey(xl.protein1, xl.protein2);
      return key === selectedKey;
    });
  }

  selectCrosslink(id: string | null): void {
    this.selectedCrosslinkId = id;
    if (id) {
      const xl = this.bundle.crosslinks.find((x) => x.id === id);
      if (xl) {
        this.selectedProteinId = xl.protein1;
      }
    }
    this.notify();
  }

  selectProtein(id: string | null): void {
    this.selectedProteinId = id;
    this.filters.proteinId = id;
    this.notify();
  }

  setShowFailed(show: boolean): void {
    this.filters.showFailed = show;
    this.notify();
  }

  setMinScore(score: number): void {
    this.filters.minScore = score;
    this.notify();
  }
}
