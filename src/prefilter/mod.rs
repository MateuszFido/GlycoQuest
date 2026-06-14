//! Spectrum prefilter orchestration and TSV output writers.

mod diagnostic;
mod isotope;
mod prune;

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::cli::settings::Settings;
use crate::crosslinker::CrosslinkerProfile;
use crate::glyco::GlycanLibrary;
use crate::mzxml::{self, Ms2Scan};

pub use diagnostic::{DiagnosticMatch, MatchedIon};
pub use isotope::{IsotopePair, IsotopePairOutcome, ScanRef};
pub use prune::PrunedGlycan;

#[derive(Debug, Clone, PartialEq)]
pub struct FilteredSpectrum {
    pub source_file: PathBuf,
    pub scan_number: u32,
    pub retention_time_min: f64,
    pub precursor_mz: f64,
    pub precursor_charge: Option<u8>,
    pub matched_families: Vec<String>,
    pub matched_ions: Vec<MatchedIon>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct RejectedSpectrum {
    pub source_file: PathBuf,
    pub scan_number: u32,
    pub reason: RejectReason,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RejectReason {
    NoDiagnostic,
    NoIsotopePair,
}

impl RejectReason {
    fn as_str(self) -> &'static str {
        match self {
            Self::NoDiagnostic => "no_diagnostic",
            Self::NoIsotopePair => "no_isotope_pair",
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct GlycanPruningRow {
    pub source_file: PathBuf,
    pub scan_number: u32,
    pub glycan_name: String,
    pub composition: String,
    pub matched_families: Vec<String>,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct PrefilterStats {
    pub scans_total: usize,
    pub diagnostic_positive: usize,
    pub isotope_pairs: usize,
    pub filtered_scans: usize,
    pub rejected: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub struct PrefilterResult {
    pub filtered: Vec<FilteredSpectrum>,
    pub isotope_pairs: Vec<IsotopePair>,
    pub rejected: Vec<RejectedSpectrum>,
    pub pruning: Vec<GlycanPruningRow>,
    pub stats: PrefilterStats,
}

pub fn run_prefilter(
    files: &[PathBuf],
    library: &GlycanLibrary,
    settings: &Settings,
    crosslinker: &CrosslinkerProfile,
) -> Result<PrefilterResult, String> {
    let mut rejected = Vec::new();
    let mut diagnostic_positive = Vec::new();
    let mut stats = PrefilterStats::default();

    for file in files {
        let scans = mzxml::parse_scans(file)?;
        stats.scans_total += scans.len();

        for scan in scans {
            let diag = diagnostic::match_diagnostic_ions(
                &scan,
                library,
                settings.diagnostic_tolerance_ppm,
            );
            if !diag.passes {
                rejected.push(RejectedSpectrum {
                    source_file: file.clone(),
                    scan_number: scan.scan_number,
                    reason: RejectReason::NoDiagnostic,
                });
                continue;
            }

            stats.diagnostic_positive += 1;
            diagnostic_positive.push(ScanWithMatch {
                source_file: file.clone(),
                scan,
                diagnostic: diag,
            });
        }
    }

    let mut isotope_pairs = Vec::new();
    let mut filtered = Vec::new();
    let mut pruning = Vec::new();

    if crosslinker.requires_isotope_pair_prefilter() {
        apply_isotope_pair_filter(
            &diagnostic_positive,
            settings,
            &mut isotope_pairs,
            &mut filtered,
            &mut pruning,
            &mut rejected,
            &mut stats,
            library,
        )?;
    } else {
        pass_diagnostic_positive(
            &diagnostic_positive,
            &mut filtered,
            &mut pruning,
            library,
        )?;
        stats.filtered_scans = filtered.len();
    }

    stats.rejected = rejected.len();

    if settings.max_pruned_spectra > 0 && stats.filtered_scans > settings.max_pruned_spectra as usize
    {
        return Err(format!(
            "retained spectra ({}) exceed max_pruned_spectra ({})",
            stats.filtered_scans, settings.max_pruned_spectra
        ));
    }

    Ok(PrefilterResult {
        filtered,
        isotope_pairs,
        rejected,
        pruning,
        stats,
    })
}

fn apply_isotope_pair_filter(
    diagnostic_positive: &[ScanWithMatch],
    settings: &Settings,
    isotope_pairs: &mut Vec<IsotopePair>,
    filtered: &mut Vec<FilteredSpectrum>,
    pruning: &mut Vec<GlycanPruningRow>,
    rejected: &mut Vec<RejectedSpectrum>,
    stats: &mut PrefilterStats,
    library: &GlycanLibrary,
) -> Result<(), String> {
    let scan_refs: Vec<ScanRef> = diagnostic_positive
        .iter()
        .map(|entry| ScanRef {
            source_file: entry.source_file.clone(),
            scan: entry.scan.clone(),
        })
        .collect();

    for outcome in isotope::match_isotope_pairs(&scan_refs, settings) {
        match outcome {
            IsotopePairOutcome::Paired(pair) => {
                stats.isotope_pairs += 1;
                isotope_pairs.push(pair.clone());

                if let Some(light) =
                    find_entry(diagnostic_positive, &pair.light_file, pair.light_scan)
                {
                    append_filtered_and_pruning(filtered, pruning, light, library)?;
                }
                if let Some(heavy) =
                    find_entry(diagnostic_positive, &pair.heavy_file, pair.heavy_scan)
                {
                    append_filtered_and_pruning(filtered, pruning, heavy, library)?;
                }
            }
            IsotopePairOutcome::Unpaired(scan_ref) => {
                rejected.push(RejectedSpectrum {
                    source_file: scan_ref.source_file,
                    scan_number: scan_ref.scan.scan_number,
                    reason: RejectReason::NoIsotopePair,
                });
            }
        }
    }

    stats.filtered_scans = filtered.len();
    Ok(())
}

fn pass_diagnostic_positive(
    diagnostic_positive: &[ScanWithMatch],
    filtered: &mut Vec<FilteredSpectrum>,
    pruning: &mut Vec<GlycanPruningRow>,
    library: &GlycanLibrary,
) -> Result<(), String> {
    for entry in diagnostic_positive {
        append_filtered_and_pruning(filtered, pruning, entry, library)?;
    }
    Ok(())
}

struct ScanWithMatch {
    source_file: PathBuf,
    scan: Ms2Scan,
    diagnostic: DiagnosticMatch,
}

fn find_entry<'a>(
    entries: &'a [ScanWithMatch],
    file: &Path,
    scan_number: u32,
) -> Option<&'a ScanWithMatch> {
    entries.iter().find(|e| e.source_file == file && e.scan.scan_number == scan_number)
}

fn append_filtered_and_pruning(
    filtered: &mut Vec<FilteredSpectrum>,
    pruning: &mut Vec<GlycanPruningRow>,
    entry: &ScanWithMatch,
    library: &GlycanLibrary,
) -> Result<(), String> {
    if filtered.iter().any(|f| {
        f.source_file == entry.source_file && f.scan_number == entry.scan.scan_number
    }) {
        return Ok(());
    }

    filtered.push(FilteredSpectrum {
        source_file: entry.source_file.clone(),
        scan_number: entry.scan.scan_number,
        retention_time_min: entry.scan.retention_time_min,
        precursor_mz: entry.scan.precursor_mz,
        precursor_charge: entry.scan.precursor_charge,
        matched_families: entry.diagnostic.matched_families.clone(),
        matched_ions: entry.diagnostic.matched_ions.clone(),
    });

    let candidates = prune::prune_glycans(&entry.diagnostic.matched_families, library)?;
    for glycan in candidates {
        pruning.push(GlycanPruningRow {
            source_file: entry.source_file.clone(),
            scan_number: entry.scan.scan_number,
            glycan_name: glycan.name,
            composition: glycan.composition,
            matched_families: entry.diagnostic.matched_families.clone(),
        });
    }

    Ok(())
}

pub fn write_outputs(out_dir: &Path, result: &PrefilterResult) -> Result<(), String> {
    write_filtered_spectra(&out_dir.join("filtered_spectra.tsv"), &result.filtered)?;
    write_isotope_pairs(&out_dir.join("isotope_pairs.tsv"), &result.isotope_pairs)?;
    write_rejected_spectra(&out_dir.join("rejected_spectra.tsv"), &result.rejected)?;
    write_glycan_pruning(&out_dir.join("glycan_pruning.tsv"), &result.pruning)?;
    Ok(())
}

fn write_io(result: std::io::Result<()>) -> Result<(), String> {
    result.map_err(|err| err.to_string())
}

fn write_filtered_spectra(path: &Path, rows: &[FilteredSpectrum]) -> Result<(), String> {
    let mut w = tsv_writer(path)?;
    write_io(writeln!(
        w,
        "source_file\tscan\trt_min\tprecursor_mz\tcharge\tmatched_families\tmatched_ions"
    ))?;
    for row in rows {
        write_io(writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            row.source_file.display(),
            row.scan_number,
            row.retention_time_min,
            row.precursor_mz,
            row.precursor_charge
                .map(|c| c.to_string())
                .unwrap_or_else(|| ".".into()),
            row.matched_families.join(";"),
            format_matched_ions(&row.matched_ions),
        ))?;
    }
    Ok(())
}

fn write_isotope_pairs(path: &Path, rows: &[IsotopePair]) -> Result<(), String> {
    let mut w = tsv_writer(path)?;
    write_io(writeln!(
        w,
        "light_file\tlight_scan\theavy_file\theavy_scan\trt_light\trt_heavy\tmz_light\tmz_heavy\tcharge"
    ))?;
    for row in rows {
        write_io(writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            row.light_file.display(),
            row.light_scan,
            row.heavy_file.display(),
            row.heavy_scan,
            row.rt_light_min,
            row.rt_heavy_min,
            row.mz_light,
            row.mz_heavy,
            row.charge,
        ))?;
    }
    Ok(())
}

fn write_rejected_spectra(path: &Path, rows: &[RejectedSpectrum]) -> Result<(), String> {
    let mut w = tsv_writer(path)?;
    write_io(writeln!(w, "source_file\tscan\treason"))?;
    for row in rows {
        write_io(writeln!(
            w,
            "{}\t{}\t{}",
            row.source_file.display(),
            row.scan_number,
            row.reason.as_str(),
        ))?;
    }
    Ok(())
}

fn write_glycan_pruning(path: &Path, rows: &[GlycanPruningRow]) -> Result<(), String> {
    let mut w = tsv_writer(path)?;
    write_io(writeln!(
        w,
        "source_file\tscan\tglycan_name\tcomposition\tmatched_families"
    ))?;
    for row in rows {
        write_io(writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}",
            row.source_file.display(),
            row.scan_number,
            row.glycan_name,
            row.composition,
            row.matched_families.join(";"),
        ))?;
    }
    Ok(())
}

fn format_matched_ions(ions: &[MatchedIon]) -> String {
    ions.iter()
        .map(|ion| {
            if ion.loss_label.is_empty() {
                format!("{}@{:.4}", ion.family, ion.observed_mz)
            } else {
                format!("{}@{:.4}[{}]", ion.family, ion.observed_mz, ion.loss_label)
            }
        })
        .collect::<Vec<_>>()
        .join(";")
}

fn tsv_writer(path: &Path) -> Result<BufWriter<File>, String> {
    let file = File::create(path).map_err(|err| {
        format!("cannot write output file {}: {err}", path.display())
    })?;
    Ok(BufWriter::new(file))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::crosslinker::CrosslinkerProfile;
    use crate::glyco::load_glycan_database;
    use std::fs;

    fn dss_profile() -> CrosslinkerProfile {
        CrosslinkerProfile::resolve(&Settings::defaults(), Some("dss")).unwrap()
    }

    fn unlabeled_profile() -> CrosslinkerProfile {
        CrosslinkerProfile::resolve(&Settings::defaults(), Some("dmtmm")).unwrap()
    }

    fn fixture(name: &str) -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/mzxml")
            .join(name)
    }

    #[test]
    fn end_to_end_dss_pair_fixture() {
        let library = load_glycan_database("nglyc309").unwrap();
        let settings = Settings::defaults();
        let files = vec![fixture("dss_pair.mzXML")];
        let result = run_prefilter(&files, &library, &settings, &dss_profile()).unwrap();
        assert_eq!(result.stats.diagnostic_positive, 2);
        assert_eq!(result.stats.isotope_pairs, 1);
        assert_eq!(result.filtered.len(), 2);
        assert!(!result.pruning.is_empty());
    }

    #[test]
    fn no_diagnostic_fixture_rejects_all() {
        let library = load_glycan_database("nglyc309").unwrap();
        let settings = Settings::defaults();
        let files = vec![fixture("no_diagnostic.mzXML")];
        let result = run_prefilter(&files, &library, &settings, &dss_profile()).unwrap();
        assert!(result.filtered.is_empty());
        assert_eq!(result.rejected.len(), 1);
        assert!(matches!(
            result.rejected[0].reason,
            RejectReason::NoDiagnostic
        ));
    }

    #[test]
    fn unlabeled_crosslinker_passes_single_diagnostic_scan() {
        let library = load_glycan_database("nglyc309").unwrap();
        let settings = Settings::defaults();
        let files = vec![fixture("hexnac_positive.mzXML")];
        let result =
            run_prefilter(&files, &library, &settings, &unlabeled_profile()).unwrap();
        assert_eq!(result.stats.diagnostic_positive, 1);
        assert_eq!(result.stats.isotope_pairs, 0);
        assert_eq!(result.filtered.len(), 1);
        assert!(result.isotope_pairs.is_empty());
    }

    #[test]
    fn writes_all_tsv_outputs() {
        let library = load_glycan_database("nglyc309").unwrap();
        let settings = Settings::defaults();
        let files = vec![fixture("dss_pair.mzXML")];
        let result = run_prefilter(&files, &library, &settings, &dss_profile()).unwrap();

        let out = std::env::temp_dir().join(format!(
            "glycoquest_prefilter_out_{}",
            std::process::id()
        ));
        fs::create_dir_all(&out).unwrap();
        write_outputs(&out, &result).unwrap();

        for name in [
            "filtered_spectra.tsv",
            "isotope_pairs.tsv",
            "rejected_spectra.tsv",
            "glycan_pruning.tsv",
        ] {
            assert!(out.join(name).is_file(), "missing {name}");
        }
        let _ = fs::remove_dir_all(out);
    }
}
