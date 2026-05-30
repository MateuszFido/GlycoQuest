# GlycoQuest Design Prompt

## Role
You are designing GlycoQuest, a Python CLI extension for xQuest.

GlycoQuest should function as a wrapper around xQuest in the same spirit as
RNxQuest, but for glycopeptide-peptide crosslinks. The requested output is a
technical design and implementation plan, not immediate code.

## Required Local Context
Before proposing the design, inspect the local reference projects:

- `RNxQuest/`: package layout, CLI entry points, template files, parameter
  generation, result extraction, and post-processing.
- `xQuest/`: available documentation, expected input formats, default files,
  and how xQuest jobs are configured and launched.

Summarize the relevant RNxQuest/xQuest patterns you will reuse. Do not invent
xQuest behavior or parameters that are not supported by the local files or
documentation; list any uncertainty explicitly.

## V1 Goal
Design a standalone Python CLI package that accepts xQuest-compatible MS/MS
input files and an explicit glycan library, filters spectra using glycan
diagnostic ions and neutral losses, prunes the glycan search space, generates
xQuest jobs, runs xQuest or produces an inspectable dry-run plan, and
consolidates xQuest results with GlycoQuest-specific post-filters.

The V1 tool should be usable from the command line when xQuest is installed.
Do not design LCMSpector web or desktop integration for V1, but keep the CLI
and file outputs process-friendly so LCMSpector can call it later.

## Runtime and Packaging Constraints
- Implement as a Python CLI package modeled after RNxQuest.
- xQuest is the required search engine dependency.
- Prefer Python standard library code unless a dependency is clearly justified.
- Do not require LCMSpector, Rust, TypeScript, a database, or a web service.
- Do not hide xQuest behind an opaque abstraction; users should be able to
  inspect generated xQuest files and commands.
- If raw vendor formats are supplied, fail with a clear message requiring
  explicit conversion to an xQuest-compatible format.

## Required Input Contracts
Define a concrete V1 CLI and file contract.

At minimum, include MS/MS input path or directory, FASTA/database path, glycan
library path, xQuest executable or installation path, output directory, mass
tolerances for diagnostic-ion filtering, and dry-run versus run behavior.

The glycan library should be an explicit CSV or TSV file. Propose columns such
as `name`, `composition`, `monoisotopic_mass`, `diagnostic_ions`, and
`neutral_losses`.

Define validation rules and error messages for malformed input.

## Scientific and Search Requirements
The primary crosslinker is DSS. It crosslinks primary amines, mainly lysines
and protein N-termini. The covalent crosslink is between peptide residues; the
glycan is not part of the crosslinking bond.

The design must account for intraprotein and interprotein crosslinks; K-K,
K-N-terminus, and N-terminus-N-terminus edge cases where supported by xQuest;
carbamidomethylation as fixed/common modification; oxidation as variable
modification; enough variable modification capacity for glycan masses and
oxidation; glycan fragments, partial compositions, and neutral water losses;
diagnostic ions such as HexNAc and sialic acid ions; and post-xQuest filtering
analogous to RNxQuest-style result cleanup.

If a chemical assumption is uncertain, call it out and propose how the V1
design should make it configurable rather than hard-coded.

## Required Design Output
Produce a structured design with these sections:

1. **Scope and non-goals**: what V1 will and will not do.
2. **RNxQuest/xQuest findings**: concrete local patterns to reuse.
3. **CLI contract**: commands, arguments, examples, and exit behavior.
4. **File formats**: glycan library schema, filtered-spectrum output, generated
   xQuest files, and consolidated result files.
5. **Workflow**: validation, spectrum parsing, diagnostic filtering, glycan
   pruning, xQuest job generation, xQuest execution, result extraction, and
   post-filtering.
6. **xQuest parameter strategy**: glycan masses, fixed/variable modifications,
   neutral losses, DSS settings, and job splitting.
7. **Testing plan**: unit tests, fixture data, integration tests with xQuest,
   and checks for generated command/file correctness.
8. **Implementation plan**: small ordered tasks suitable for commits.
9. **Risks and open questions**: especially xQuest limitations, chemistry
   assumptions, false-discovery handling, and input parsing constraints.

## Forbidden Output
- Do not write implementation code yet.
- Do not design a web UI.
- Do not introduce LCMSpector-specific APIs.
- Do not silently skip diagnostic-ion filtering.
- Do not claim FDR support unless the xQuest/xProphet path is explicitly
  designed and verified.
- Do not assume automatic glycan database discovery; V1 uses the explicit
  user-provided glycan library.
