# GlycoQuest V1 Design

## 1. Scope and Non-Goals

GlycoQuest V1 is a standalone Rust CLI package that prepares and runs xQuest searches for crosslinked glycopeptide-peptide spectra (DSS by default; DMTMM, published glycan-to-peptide linkers, and custom chemistries via settings). It accepts xQuest-compatible MS/MS input, a FASTA database, and a bundled or explicit glycan library. It filters spectra for glycan diagnostic ions, annotates neutral losses, optionally filters for DSS light/heavy isotopic pairs when `crosslinker.label=light-heavy`, prunes per-spectrum glycan candidates, writes inspectable xQuest job folders, optionally runs xQuest, and consolidates xQuest XML results with GlycoQuest-specific annotations and post-filters.

V1 does not include LCMSpector web or desktop integration, raw vendor conversion, automatic glycan database discovery, a database, a service, or a web UI. It must fail clearly for raw vendor formats and require explicit conversion to xQuest-compatible mzXML. V1 must not claim FDR/q-value support until an xQuest/xProphet path is explicitly validated for this glycan-filtered search mode.

## 2. RNxQuest/xQuest Findings

RNxQuest patterns to reuse:

- Cargo package layout with a CLI binary and bundled templates, reusing RNxQuest's visible package structure without its original packaging mechanics.
- Parameter-generation CLI that validates inputs, records arguments, creates job directories, and writes `.def` files.
- Transparent generated files and commands: RNxQuest writes xQuest definition files, shell scripts, and mass mapping CSVs rather than hiding xQuest.
- Mass/loss expansion pattern: RNxQuest generates neutral-loss mass variants and keeps a mapping table for result annotation.
- Result extraction pattern: RNxQuest flattens `spectrum_search` and `search_hit` XML attributes to CSV, then adds domain-specific cleanup and annotations.

xQuest patterns and limits observed locally:

- xQuest jobs are configured through `xquest.def` and `xmm.def`.
- xQuest defaults include DSS settings such as `AArequired K:K`, `xlinkermw 138.0680796`, `xlinktypes 1111`, fixed C carbamidomethylation, and optional `variable_mod`.
- `pQuest.pl` prepares per-mzXML folders and pseudoSH master maps; `runXquest.pl` orchestrates compare-peaks and xQuest execution.
- xQuest's default DSS path expects a light/heavy isotope-coded crosslinker distribution with `isotopeshift 12.075321`, `printisotopicscanpairs 1`, and `cp_isotopediff 12.075321`.
- The local light-only instructions set `isotopeshift 0`, `printlightonlypairs 1`, and `cp_isotopediff 0`, but that is not the V1 default.
- xQuest comments document N-terminus support through `ntermxlinkable 1` and amino acid code `Z`, but this needs integration testing.
- GlycoQuest should use the patched local runtime `glycoquest/xQuest/V2.1.6/xquest`. `V2.1.5` stays archival/original.
- xQuest represents variable modifications by rewriting residues to pseudo-residues. `X`, `U`, `B`, and `J` map to the first four `variable_mod` entries. For example, `variable_mod N,203.079373,M,15.994915` means `X` is glycosylated `N` and `U` is oxidized `M`.
- `V2.1.6` fixes the local multi-variable-mod path so pseudo-residue masses are populated for `X/U/B/J` and `nvariable_mod` is enforced as a maximum number of variable modifications per peptide.
- Since xQuest treats pseudo-residues as ordinary residues after expansion, GlycoQuest must parse `seq1` and `seq2` in xQuest output, map pseudo-residues back to their configured source residues and mass deltas, and enforce GlycoQuest-specific glycan counts in post-filtering.
- Input FASTA files must be rejected or sanitized if they contain literal `X`, `U`, `B`, or `J`, because those letters are reserved by xQuest variable-mod expansion.
- xQuest does not score glycan diagnostic ions, glycan-family evidence, MS1 charge plausibility, or glycosylation sequons. GlycoQuest owns those prefilters, annotations, and post-filters.

## 3. CLI Contract

Primary command:

```bash
glycoquest search \
  --msms ./mzxmls \
  --fasta proteins.fasta \
  --glycans glycans.tsv \
  --xquest-root ./glycoquest/xQuest/V2.1.6/xquest \
  --out glycoquest_out \
  --diagnostic-tolerance-ppm 10 \
  --neutral-loss-tolerance-da 0.05 \
  --crosslinker dss \
  --crosslinker-label light-heavy \
  --ms1-tolerance-ppm 10 \
  --ms2-tolerance-da 0.2 \
  --dry-run
```

Run mode omits `--dry-run`. Crosslinker presets: `dss` (default, light/heavy isotope prefilter required), `dmtmm` (zero-length Lys–Asp/Glu linkage), and the unlabeled `nhs-cyclooctyne`, `ssbxl`, and `pcbxl` glycan-to-peptide models. Advanced crosslinker fields live in `settings.ini` under `[crosslinker]`: `label` (`light-heavy`, `light-only`, `none`), `shift_da`, `xlinkermw`, and `xlink_sites`.

Exit codes: `0` success, `1` validation/config error, `2` no spectra pass filters, `3` xQuest execution failure, `4` result extraction failure.

Raw vendor files such as `.raw`, `.wiff`, `.d`, `.baf`, and `.tdf` fail with: `Unsupported raw vendor input: <path>. Convert explicitly to mzXML before running GlycoQuest V1.`

Real GPx reference data are available. Xie et al. (*Chemical Science*, 2021,
doi:10.1039/D1SC00814E) deposited the public `MSV000087442` dataset
(doi:10.25345/C5VV5S). Its `nhs-cyclooctyne` model is encoded as an unlabeled
`X:K` bridge of 205.085126607 Da relative to a NeuAc glycan library, and an
the full-spectrum search preset uses `data/MSV000087442/`,
`configs/msv000087442-full.ini`, and the bundled `msv000087442-sianaz` library.
Chen et al. (*Analytical Chemistry*, 2025, doi:10.1021/acs.analchem.4c04134)
reported `MSV000093174` (doi:10.25345/C5416T929); the paper-derived post-release
SSBXL and PCBXL models are bundled, although anonymous dataset access was private
when checked on 2026-07-15. Local real-data tests use explicit
`GLYCOQUEST_INTEGRATION_MZXML` and `GLYCOQUEST_INTEGRATION_FASTA` paths.

## 4. File Formats

Glycan library CSV/TSV schema:

```text
name,composition,monoisotopic_mass,diagnostic_ions,residue_targets
HexNAc1,HexNAc(1),203.079373,HexNAc@204.0867;HexNAc@186.0760[-H2O];HexNAc@168.0654[-2H2O],N
NeuAc1HexNAc1,NeuAc(1)HexNAc(1),494.174789,NeuAc@292.1027;NeuAc@274.0921[-H2O],N;S;T
```

Each `diagnostic_ions` entry is an expanded search target (`family@mz` with optional
`[-loss_label]`). Neutral-loss deltas are applied at load time from the global
loss table in `diagnostic_ion_catalog.txt` (M base templates × N losses).

Validation errors should be specific: missing required column, duplicate `name`, non-numeric mass/ion, non-positive mass, empty diagnostic ions, unreadable file, or unsupported delimiter. Empty diagnostic ions are invalid unless the user explicitly opts in, because V1 must not silently skip diagnostic filtering.

Bundled seed data should live under `glycan_databases/fragpipe/` as raw upstream
references. The starter inputs are FragPipe's Byonic-style
`Nglyc309_Byonic.glyc` and `Oglyc78_Byonic.glyc`, plus `glycan_residues.txt`
for residue masses and `diagnostic_ion_catalog.txt` for diagnostic-ion templates
and neutral-loss deltas. These `.glyc` files are composition-only lists such as
`HexNAc(2)Hex(5)Fuc(1)`. The GlycoQuest converter should parse residue counts,
sum monoisotopic masses from `glycan_residues.txt`, assign default residue
targets from the source list type (`N` for N-glycan, `S;T` for O-glycan unless
overridden), and expand per-glycan search ions from the catalog (M × N).
The converter should report parsed row counts rather than trusting counts
embedded in filenames.

GlycoQuest should include this FragPipe-to-GlycoQuest conversion as an
in-built converter used during the package build or data-install step, with an
explicit CLI path for regeneration. V1 should not require a checked-in generated
TSV to be hand-maintained when the raw upstream `.glyc`, residue, and diagnostic
catalog files are present.

Outputs:

- `plan.json`: normalized inputs, options, generated jobs, and xQuest commands.
- `filtered_spectra.tsv`: retained spectra with matched diagnostic ions, neutral losses, and candidate glycans.
- `isotope_pairs.tsv`: diagnostic-positive spectra that pass the DSS light/heavy isotope-pair prefilter when enabled (`label=light-heavy`); header-only when isotope prefilter is disabled.
- `rejected_spectra.tsv`: spectra rejected by diagnostic or DSS isotope-pair filters.
- `glycan_pruning.tsv`: spectrum-to-glycan candidate table.
- `prefiltered_mzxml/`: reduced mzXML files generated from the supplied full mzXML and used as xQuest inputs.
- `jobs/<job_id>/xquest.def`, `xmm.def`, `glycoquest_matchlist.txt`, symlinked mzXML, and `run.sh`.
- `results/glycoquest_xquest.csv`: flattened xQuest hits with GlycoQuest annotations.
  Each hit carries its glycan (`glycan_name`, `glycan_composition`, `glycan_mass`,
  `loss_label`), attachment residue and peptide (`glyco_residue`, `glyco_peptide`),
  originating scan and diagnostic evidence (`source_file`, `scan`,
  `matched_families`, `matched_ion_count`), and post-filter outcome (`hard_status`,
  `soft_score`, `postfilter_status`). Hits are de-duplicated across per-glycan jobs,
  keeping the highest-scoring copy, and sorted by soft score.
- `results/xiview.csv`: passing crosslinks in CLMS-CSV layout
  (`Protein1,PepPos1,PepSeq1,LinkPos1,AbsPos1,Protein2,PepPos2,PepSeq2,LinkPos2,AbsPos2,Score,…`).
  Peptides are resolved from xQuest pseudo-residues via the job manifest and located in
  the FASTA to compute 1-based peptide start (`PepPos*`) and absolute cross-link residue
  (`AbsPos* = PepPos + LinkPos - 1`) positions; glycan annotations are appended as trailing
  columns.
- `results/report.html`: self-contained QC/glycan report (inline CSS + SVG, no network or
  JavaScript) covering the prefilter funnel, post-filter outcome breakdown, xQuest score and
  precursor-error distributions, glycan and glycosylation-site distributions, and a hits table.
- `results/viewer/`: interactive crosslink viewer bundle (`viewer.json` schema version 1,
  `database.fasta`, static `index.html` + JS/CSS). Coordinated views: QC charts, residue-resolution
  network, per-protein sequence map with glycan sites, and approximate MS/MS mirror plots from
  reduced `spectra/` mzXML. Built from the MIT `viewer/` TypeScript package.

The glycan library input (`--glycans`) accepts either a bundled database id
(`nglyc309`, `oglyc78`) or a path to a user-provided CSV/TSV file matching the
schema above. A custom O-glycan library with `residue_targets` `S;T` drives the
xQuest `variable_mod` onto Ser/Thr rather than Asn. `--resume` re-consolidates
existing job XML without prefilter/manifest state and therefore produces reduced
annotation (glycan reconstructed from the job id; diagnostic/sequon soft features
skipped).

## 5. Workflow

Validate inputs and xQuest executables, parse the glycan library, reject raw vendor formats, and parse mzXML MS2 scans. For each scan, match diagnostic-ion families within `--diagnostic-tolerance-ppm`, default 10 ppm. Any intensity is enough for pruning. Keep spectra with at least one configured glycan diagnostic family, and record any matching neutral-loss variants as annotations rather than hard requirements.

Prune glycans per spectrum by diagnostic families. If a spectrum has only HexNAc evidence, keep all glycans containing HexNAc. If it has only NeuAc evidence, keep all glycans containing NeuAc. If it has HexNAc and NeuAc evidence together, keep glycans containing both families. This is spectrum-local pruning; different spectra may keep different subsets of the glycan library.

After glycan-evidence filtering, optionally run a file-side prefilter for DSS light/heavy isotopic pairing when `crosslinker.label=light-heavy`. The default DSS workflow requires a partner spectrum separated by the configured isotope shift (`shift_da`, default 12.075321 Da), within MS1 and retention-time tolerances. When `label=none` (DMTMM preset) or `label=light-only`, diagnostic-positive spectra pass through without an isotope partner requirement; `isotope_pairs.tsv` remains empty and xQuest jobs use single-scan matchlists with `isotopeshift 0`.

Unlabeled crosslinkers skip isotope pairing. DMTMM uses `xlinkermw=-18.0109` and `K:E,K:D`. The GPx presets use `X:K`, where `X` is xQuest's first variable-modification pseudo-residue and therefore the glycan-bearing Asn: NHS–cyclooctyne 205.085126607 Da, post-TCEP/IAA SSBXL 573.179438173 Da, and post-photorelease PCBXL 456.190988659 Da. The GPx values are relative to an ordinary NeuAc glycan library.

The full mzXML should be supplied to GlycoQuest for real-data runs. Creating
small, xQuest-sized mzXML inputs is the job of the diagnostic and DSS
prefiltering step, not a manual prerequisite. An explicitly configured real-data
mzXML should remain available for optional integration testing so the prefilter
is exercised against actual scan density, precursor charge metadata, and noise.

Generate xQuest job inputs for retained spectra. When isotope prefilter is enabled, matchlist rows use the paired layout (light scan, heavy scan, m/z pair). When disabled, rows list individual scans. Dry-run writes all files, `plan.json`, and commands. Run mode executes `run.sh` per job and logs to `logs/`. Result extraction flattens xQuest XML and applies GlycoQuest post-filters; isotope-pair evidence is required only for `light-heavy` mode.

Post-xQuest filtering has hard and soft requirements. Hard requirements are configured crosslink evidence through the xQuest hit, glycan mass residual compatibility, a configured glycan residue present in the matched peptide sequence, and diagnostic-ion evidence in the originating spectrum. Soft scoring features are precursor charge verified from the MS1 isotopic envelope, glycosylation sequon presence, mass error, and diagnostic-ion count. Precursor charge must be scientifically reasonable when present: tryptic peptide-peptide hits are commonly +2/+3, while glycopeptide-peptide and glycopeptide-glycopeptide crosslinks are expected to skew higher, often +4 to +6, with +2 to +7 treated as possible.

## 6. xQuest Parameter Strategy

Use xQuest DSS crosslinking rather than RNxQuest RNA monolink mode. Default settings are `xkinkerID DSS`, `crosslinkername DSS`, `xlinkermw 138.0680796`, `isotopeshift 12.075321`, `cp_isotopediff 12.075321`, `printisotopicscanpairs 1`, and `xlinktypes 0011` for intraprotein and interprotein peptide-peptide crosslinks. Default sites are `K:K`. If enabled, N-termini use `ntermxlinkable 1` plus `Z` pairs such as `K:K,K:Z,Z:Z`, pending integration verification.

The crosslinker and labeling mode stay configurable. V1 defaults to `--crosslinker dss` with `label=light-heavy`. `label=none` skips isotope prefilter for DMTMM and the GPx presets. `label=light-only` skips isotope prefilter but writes xQuest light-only pair settings.

GlycoQuest V1 should run against the patched `V2.1.6` runtime or an external xQuest installation that is explicitly validated to have equivalent variable-modification behavior. At startup, GlycoQuest should verify that the selected xQuest tree has the expected pseudo-residue variable-mod support before generating glycan jobs.

Fixed C carbamidomethylation remains a fixed modification (`C 57.02146`). Glycan masses are configurable variable modifications on glycan-bearing residues. DSS/DMTMM keep the glycan separate from the peptide-backbone crosslink; GPx uses `X:K` so xQuest crosslinks the glycan-modified Asn pseudo-residue itself to Lys, while the preset mass accounts for the SiaNAz/linker delta relative to NeuAc. Oxidation may be specified in the same job as a glycan, for example with glycan `N,<mass>` as the first variable modification and oxidation `M,15.994915` as the second.

The default job model is one glycan composition/loss variant per xQuest job. Glycan parent masses and water-loss variants should be generated as separate search variants and mapped back to the parent glycan composition in GlycoQuest output. Do not group different glycan masses into one xQuest job unless that grouping is separately validated.

Set `variable_mod` and `nvariable_mod` per job from the selected glycan and optional peptide modifications. The first pseudo-residue should be reserved for the glycan variant where possible, so output parsing can treat `X` as the glycan-bearing residue. `nvariable_mod` is an xQuest per-peptide maximum, not a GlycoQuest hit-level glycan-count rule, so GlycoQuest must post-filter result pairs for the intended V1 class. V1 defaults to peptide-glycopeptide crosslinks and should require exactly one configured glycan pseudo-residue across the matched peptide pair. Glycopeptide-glycopeptide crosslinks are a later requirement.

Resource control is part of the contract. Before execution, dry-run and run mode should compute the number of generated xQuest jobs, retained spectra, and exact spectrum-job assignments. If those exceed `--max-jobs`, `--max-pruned-spectra`, or `--max-total-job-spectrum-comparisons`, GlycoQuest should fail before launching xQuest unless the user raises the limits.

## 7. Testing Plan

Unit tests should cover glycan library parsing, FragPipe glycan conversion, validation errors, mzXML peak decoding, diagnostic-ion matching, neutral-loss annotation, DSS light/heavy isotope-pair matching, pruning, raw-format rejection, FASTA pseudo-residue rejection, xQuest variable-mod parameter generation, pseudo-residue output parsing, post-filter hard/soft rules, and generated `.def`/matchlist/command correctness. Fixture data should include tiny mzXML files with diagnostic-positive paired scans, diagnostic-positive unpaired scans, failing scans, a tiny FASTA, a glycan library with HexNAc and sialic-acid examples, and minimal xQuest XML for extraction tests. Integration tests should include dry-run golden files, an opt-in real mzXML/FASTA pair configured through environment variables, a patched-xQuest variable-mod regression test, a glycan-plus-oxidation job, and an optional xQuest-installed search that verifies regular light/heavy execution and XML extraction.

## 8. Implementation Plan

1. Scaffold Rust crate, CLI, Cargo version metadata, and template embedding.
2. Add input validation, xQuest executable discovery, patched-runtime validation, and FASTA pseudo-residue rejection.
3. Add glycan CSV/TSV parser, schema validation, and in-built FragPipe `.glyc` converter.
4. Add mzXML parser and diagnostic/neutral-loss matcher.
5. Add DSS light/heavy isotope-pair prefilter, spectrum/pruning outputs, and pruned mzXML writing.
6. Add glycan/loss variant expansion and resource-limit planning.
7. Add xQuest job generator for regular light/heavy matchlists and `.def` files.
8. Add pseudo-residue-aware variable-mod builder for glycan plus optional oxidation.
9. Add dry-run `plan.json` and `run.sh` generation.
10. Add run mode with subprocess logging.
11. Add XML extraction and GlycoQuest post-filters.
12. Add fixture tests and optional xQuest integration tests.

## 9. Risks and Open Questions

- xQuest variable-modification capacity is now understood as four pseudo-residue slots in the patched `V2.1.6` path, but the package must ensure users run that runtime or an equivalent external installation.
- The extracted patched xQuest tree must be packaged or distributed deliberately later. The first build can proceed against the local `V2.1.6` runtime, but the parent project should not depend on a locally edited untracked runtime by accident.
- FASTA collisions with literal `X/U/B/J` would corrupt xQuest's variable-mod semantics unless rejected or sanitized.
- One glycan/loss variant per xQuest job is scientifically clear but can create many jobs. Resource limits and dry-run summaries are required before execution.
- N-terminus crosslink support through `Z` must be verified with generated jobs.
- Glycan residue targets are chemistry assumptions and must stay configurable.
- DSS light/heavy isotope pairing is required only when `crosslinker.label=light-heavy`; tolerance defaults should be verified against local instrument behavior.
- False-discovery handling remains an open question; V1 should not report FDR/q-values until xProphet use is designed and validated for GlycoQuest output.
- mzXML parsing should be intentionally narrow in V1. mzML support should wait until xQuest input handling is locally verified.
