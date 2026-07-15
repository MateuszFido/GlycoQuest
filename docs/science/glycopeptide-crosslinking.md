# Glycopeptide crosslinking


## The biological question

Glycoproteins carry oligosaccharides on specific amino acid side chains. In **crosslinking mass spectrometry (XL-MS)**, a bifunctional reagent captures molecules that are close in space. **GlycoQuest** targets cases where:

- One peptide in the pair carries an **intact glycan**
- The other peptide is typically unglycosylated
- A **chemical crosslinker** connects either the peptide backbones or a glycan handle to a peptide amine

This reveals **spatial proximity** between a glycosylated region and another peptide — useful for mapping glycan-mediated cell-surface interactions as well as glycoprotein architecture.

For isotopically labelled crosslinkers, heavy/light isotope prefiltering applies to the **crosslinker**, not the glycan. 

## N-linked vs O-linked glycans

| Type | Attachment | Sequon | Default library |
|------|------------|--------|-----------------|
| N-linked | Asparagine (Asn) | N-X-S/T (X ≠ Pro) | `nglyc309` |
| O-linked | Serine / Threonine | No strict sequon | `oglyc78` |

GlycoQuest's score for N-glycan matches is influenced by sequon presence; O-glycans rely on diagnostic ions and database context.

## Why prefilter before xQuest?

Full crosslink searches are computationally expensive. Glycopeptide spectra often show characteristic **oxonium ions** (e.g. HexNAc at ~204 m/z). Requiring diagnostic evidence:

- Reduces false searches on non-glycan spectra
- Prunes the glycan library per spectrum
- Focuses xQuest on plausible glycopeptide–peptide candidates

## Typical experimental setup (DSS)

1. Proteins crosslinked with DSS (light/heavy mixture)
2. Enzymatic digest (trypsin)
3. LC–MS/MS
4. mzXML → GlycoQuest → xQuest

See [Crosslinkers](../configuration/crosslinkers.md) for other options.

## Glycan-to-peptide setup (SiaNAz + NHS linker)

1. Feed cells ManNAz so their N-glycans contain azido-sialic acid (SiaNAz)
2. Click a cyclooctyne-bearing NHS linker onto SiaNAz
3. Let the NHS ester capture a nearby primary amine, normally Lys
4. Digest, enrich, and acquire LC–MS/MS
5. Encode the physical glycan-SiaNAz-linker-Lys bridge as xQuest `X:K`

`X` is xQuest's first variable-modification pseudo-residue. GlycoQuest assigns
the N-glycan to that slot, so `X:K` forces the glycan and the crosslink onto the
same Asn. The physical crosslink still runs through SiaNAz rather than directly
through the Asn side chain.

## Published reference datasets

Datasets that directly report glycopeptide–peptide crosslinks (GPx) do exist:

| Study | GPx experiment | Public data |
|-------|----------------|-------------|
| Xie et al., *Chemical Science* (2021), [doi:10.1039/D1SC00814E](https://doi.org/10.1039/D1SC00814E) | NHS–cyclooctyne workflow for sialic-acid-mediated cell-surface interactions in PNT2 cells; Supplementary Table S1 contains 494 GPx rows. | Public MGF, mzML, FASTA, mzIdentML, and mzTab: [MassIVE MSV000087442](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?accession=MSV000087442), [doi:10.25345/C5VV5S](https://doi.org/10.25345/C5VV5S) |
| Chen et al., *Analytical Chemistry* (2025), [doi:10.1021/acs.analchem.4c04134](https://doi.org/10.1021/acs.analchem.4c04134) | Enrichable disulfide and photocleavable linkers on PNT2 cells; 12,779 SSBXL and 4,186 PCBXL GPx are reported. | Deposited as [MassIVE MSV000093174](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?accession=MSV000093174), [doi:10.25345/C5416T929](https://doi.org/10.25345/C5416T929), but anonymous access redirected to a private task when checked on 2026-07-15. The article and synthetic supporting information are public. |

### Paper-derived search masses

All three presets assume the glycan library contains ordinary NeuAc. Their
`xlinkermw` therefore includes the SiaNAz-for-NeuAc difference as well as the
retained linker fragment:

| Preset | Sample state searched | Relative bridge formula | Monoisotopic mass | Sites |
|--------|-----------------------|-------------------------|------------------:|-------|
| `nhs-cyclooctyne` | Intact product from Xie et al. | C10H11N3O2 | 205.085126607 Da | `X:K` |
| `ssbxl` | Chen et al. after TCEP cleavage and iodoacetamide alkylation | C28H27N7O5S | 573.179438173 Da | `X:K` |
| `pcbxl` | Chen et al. after photocleavage | C25H24N6O3 | 456.190988659 Da | `X:K` |

For SSBXL, bridge + NeuAc residue + proton is 865.2821, independently matching
the m/z 865.28 linker signature required by the paper's MeroX search. The SSBXL
mass is specifically for the post-TCEP/**IAA** material, not the intact enrichment
reagent. Likewise, PCBXL uses the post-photorelease bridge shown in the synthetic
scheme.

### Executed MSV000087442 proof of concept

The repository includes reproducible commands and small reference files in
`examples/MSV000087442/`. Four spectra
were reconstructed by matching high-scoring Supplementary Table S1 rows to the
in-situ MGF using precursor mass (within 10 ppm) plus both m/z 204.0867 and
366.1395. GlycoQuest retained all four scans, all three xQuest jobs completed,
and four hits were consolidated and passed the GlycoQuest post-filter. The
generated xQuest definitions use `AArequired X:K`, ensuring that the glycan
pseudo-residue is also the crosslink site.

This validates the raw-data conversion, linker mass/site definition, prefilter,
xQuest launch, and result extraction on real GPx data. It does **not** establish
search-engine concordance: the paper used MeroX with glycan-specific offset ions,
and the two xQuest assignments differ from Supplementary Table S1. The paper also
does not place scan identifiers beside those table rows, so the checked-in scan
selection is explicitly a mass/fragment-ion reconstruction.

## Related

- [Diagnostic ions](diagnostic-ions.md)
- [Workflow overview](../workflow/overview.md)
- [Interpreting hits](../results/interpreting-hits.md)
