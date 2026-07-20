# Installation

GlycoQuest is a Rust CLI using a local **xQuest V2.1.7** installation. You also need **mzXML** MS/MS input (convert vendor raw files separately) and a protein **FASTA** database.

## Download pre-built binaries

Pre-built `glycoquest` binaries are published as release artifacts for **Windows**, **macOS**, and **Linux**. This is the recommended install path.

1. Open the [GlycoQuest releases](https://git.proxiomics.com/matt/GlycoQuest/-/releases) page.
2. Download the archive for your platform from the latest release.
3. Extract it and add the `glycoquest` binary to your `PATH`, or invoke it by full path.

### Linux

Release binaries are built on **Ubuntu 22.04** (glibc 2.35) so they run on ETH Euler and other 22.04-based systems. Binaries from newer distros (e.g. Ubuntu 24.04 / glibc 2.39) will fail with `GLIBC_2.39 not found`.

```bash
tar xzf glycoquest-*-x86_64-unknown-linux-gnu.tar.gz
./glycoquest --help
```

### macOS

```bash
tar xzf glycoquest-*-aarch64-apple-darwin.tar.gz   # Apple Silicon (CI builds)
./glycoquest --help
```

Intel Macs: run the aarch64 binary under Rosetta, or build locally with `cargo build --release` (GitHub Actions no longer provides macOS x86_64 runners).

### Windows

Extract `glycoquest-*-x86_64-pc-windows-msvc.zip` and run from PowerShell or Command Prompt:

```powershell
.\glycoquest.exe --help
```

After installation, run `glycoquest --help` from any directory (if on `PATH`) or from the extract folder.

## Requirements

| Component | Purpose | Notes |
|-----------|---------|-------|
| **Rust toolchain** (edition 2024) | Build `glycoquest` | [rustup.rs](https://rustup.rs) |
| **xQuest V2.1.7+** | Crosslink search engine | Path passed as `--xquest-root`, available on [Gitlab](https://gitlab.ethz.ch/leitner_lab/xquest_xprophet) |
| **Perl + DB_File + XML::Parser** | xQuest indexing and mzXML/XML parsing | Package managers, CPAN, or conda. On ETH Euler: `scripts/bootstrap-euler-perl.sh` then `scripts/check-xquest-perl.pl` |
| **ProteoWizard msconvert** | Raw → mzXML (optional) | GlycoQuest does **not** convert vendor formats |

### xQuest layout

Point `--xquest-root` at the directory that contains `bin/xquest.pl` and bundled Perl libraries. 

```text
xQuest/V2.1.7/xquest/
├── bin/xquest.pl
├── 1209/lib/perl5/   # bundled Perl modules (PERL5LIB set in run.sh)
└── ...
```

The default in `settings.ini` at the repository root is `xquest_bin = bin/xquest.pl` (relative to `--xquest-root`).

## Build from source

If you'd like to build from source, 

### Requirements

| Component | Purpose | Notes |
|-----------|---------|-------|
| **Rust toolchain** (edition 2024) | Build `glycoquest` | [rustup.rs](https://rustup.rs) |

### Build GlycoQuest

```bash
git clone <repository-url> GlycoQuest
cd GlycoQuest
cargo build --release
./target/release/glycoquest --help
```

The binary is `./target/release/glycoquest`. If you add it to `PATH`, you can just call it as `glycoquest --options`, otherwise by full path `/your/path/here/target/release/glycoquest --options`.

### Build the viewer (optional)

GlycoQuest copies prebuilt viewer assets into `results/viewer/` after each run. Rebuild only when developing the viewer:

```bash
cd viewer
npm install
npm run build
```

## Input formats

| Format | Supported |
|--------|-----------|
| mzXML (MS/MS) | Yes |
| Directory of mzXML files | Yes |
| `.raw`, `.wiff`, `.d`, `.baf`, `.tdf` | **No** — convert with msconvert first |

Example msconvert command:

```bash
msconvert your-input-file.raw --mzXML --filter "peakPicking true 1-"
```

## Verify installation

GlycoQuest includes a dry-run option (without xQuest execution):

```bash
glycoquest tests/fixtures/mzxml/dss_pair.mzXML \
  --database data/example.fasta \
  --xquest-root xQuest/V2.1.7/xquest \
  --out /tmp/glycoquest_test \
  --dry-run
```

If you built from source and have not added the binary to `PATH`, use `./target/release/glycoquest` instead of `glycoquest`.

A successful run prints a **readiness** block with green `✓ PASS` for MS input, FASTA, glycan library, xQuest runtime, and output directory. See [First run](first-run.md) for the full workflow.

## Environment variables

| Variable | Effect |
|----------|--------|
| `GLYCOQUEST_KEEP_TMP=1` | Keep `tmp/jobs/` and logs after a successful run (debugging) |
| `GLYCOQUEST_GLYCAN_DATA_DIR` | Override path to bundled glycan `.glyc` databases |
| `GLYCOQUEST_XQUEST_ROOT` | Integration tests: xQuest root path |

## Next steps

- [First run](first-run.md) — dry-run, full search, open results
- [Settings reference](../configuration/settings-reference.md) — `settings.ini`
