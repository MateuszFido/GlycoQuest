//! Load advanced options from `settings.ini`.

use std::path::{Path, PathBuf};

use configparser::ini::Ini;

/// Advanced GlycoQuest options read from `settings.ini`.
#[derive(Debug, Clone, PartialEq)]
pub struct Settings {
    pub xquest_bin: Option<PathBuf>,
    pub diagnostic_tolerance_ppm: f64,
    pub neutral_loss_tolerance_da: f64,
    pub ms1_tolerance_ppm: f64,
    pub ms2_tolerance_da: f64,
    pub isotope_pair_ms1_tolerance_ppm: f64,
    pub isotope_pair_rt_tolerance_min: f64,
    pub crosslinker_name: String,
    pub crosslinker_label: String,
    pub crosslinker_shift_da: f64,
    pub xlink_sites: String,
    pub nterm_xlinkable: bool,
    pub fixed_carbamidomethyl_cys: bool,
    pub variable_oxidation: bool,
    pub glycan_targets: String,
    pub max_jobs: u32,
    pub max_pruned_spectra: u32,
    pub max_total_job_spectrum_comparisons: u64,
    pub min_score: f64,
    pub max_precursor_error_ppm: f64,
}

impl Settings {
    pub fn defaults() -> Self {
        Self {
            xquest_bin: None,
            diagnostic_tolerance_ppm: 10.0,
            neutral_loss_tolerance_da: 0.05,
            ms1_tolerance_ppm: 10.0,
            ms2_tolerance_da: 0.2,
            isotope_pair_ms1_tolerance_ppm: 10.0,
            isotope_pair_rt_tolerance_min: 2.0,
            crosslinker_name: "dss".into(),
            crosslinker_label: "light-heavy".into(),
            crosslinker_shift_da: 12.075321,
            xlink_sites: "K:K".into(),
            nterm_xlinkable: false,
            fixed_carbamidomethyl_cys: true,
            variable_oxidation: false,
            glycan_targets: "N,S,T".into(),
            max_jobs: 0,
            max_pruned_spectra: 0,
            max_total_job_spectrum_comparisons: 0,
            min_score: 0.0,
            max_precursor_error_ppm: 20.0,
        }
    }

    /// Load settings from `path`, falling back to built-in defaults for missing keys.
    pub fn load_from_file(path: impl AsRef<Path>) -> Result<Self, String> {
        let path = path.as_ref();
        let path_str = path
            .to_str()
            .ok_or_else(|| format!("settings path is not valid UTF-8: {}", path.display()))?;

        let mut ini = Ini::new();
        ini.load(path_str)
            .map_err(|e| format!("failed to read settings file {}: {e}", path.display()))?;

        Ok(Self::from_ini(&ini))
    }

    fn from_ini(ini: &Ini) -> Self {
        let mut s = Self::defaults();

        s.xquest_bin = non_empty_path(ini.get("xquest", "xquest_bin"));

        s.diagnostic_tolerance_ppm =
            get_f64(ini, "tolerances", "diagnostic_tolerance_ppm", s.diagnostic_tolerance_ppm);
        s.neutral_loss_tolerance_da =
            get_f64(ini, "tolerances", "neutral_loss_tolerance_da", s.neutral_loss_tolerance_da);
        s.ms1_tolerance_ppm = get_f64(ini, "tolerances", "ms1_tolerance_ppm", s.ms1_tolerance_ppm);
        s.ms2_tolerance_da = get_f64(ini, "tolerances", "ms2_tolerance_da", s.ms2_tolerance_da);
        s.isotope_pair_ms1_tolerance_ppm = get_f64(
            ini,
            "tolerances",
            "isotope_pair_ms1_tolerance_ppm",
            s.isotope_pair_ms1_tolerance_ppm,
        );
        s.isotope_pair_rt_tolerance_min = get_f64(
            ini,
            "tolerances",
            "isotope_pair_rt_tolerance_min",
            s.isotope_pair_rt_tolerance_min,
        );

        if let Some(v) = non_empty_string(ini.get("crosslinker", "name")) {
            s.crosslinker_name = v;
        }
        if let Some(v) = non_empty_string(ini.get("crosslinker", "label")) {
            s.crosslinker_label = v;
        }
        s.crosslinker_shift_da = get_f64(ini, "crosslinker", "shift_da", s.crosslinker_shift_da);
        if let Some(v) = non_empty_string(ini.get("crosslinker", "xlink_sites")) {
            s.xlink_sites = v;
        }
        s.nterm_xlinkable = get_bool(ini, "crosslinker", "nterm_xlinkable", s.nterm_xlinkable);

        s.fixed_carbamidomethyl_cys =
            get_bool(ini, "modifications", "fixed_carbamidomethyl_cys", s.fixed_carbamidomethyl_cys);
        s.variable_oxidation =
            get_bool(ini, "modifications", "variable_oxidation", s.variable_oxidation);

        if let Some(v) = non_empty_string(ini.get("glycan", "targets")) {
            s.glycan_targets = v;
        }

        s.max_jobs = get_u32(ini, "limits", "max_jobs", s.max_jobs);
        s.max_pruned_spectra = get_u32(ini, "limits", "max_pruned_spectra", s.max_pruned_spectra);
        s.max_total_job_spectrum_comparisons = get_u64(
            ini,
            "limits",
            "max_total_job_spectrum_comparisons",
            s.max_total_job_spectrum_comparisons,
        );
        s.min_score = get_f64(ini, "limits", "min_score", s.min_score);
        s.max_precursor_error_ppm =
            get_f64(ini, "limits", "max_precursor_error_ppm", s.max_precursor_error_ppm);

        s
    }
}

pub fn default_settings_path() -> PathBuf {
    PathBuf::from("settings.ini")
}

fn non_empty_path(value: Option<String>) -> Option<PathBuf> {
    non_empty_string(value).map(PathBuf::from)
}

fn non_empty_string(value: Option<String>) -> Option<String> {
    value.filter(|s| !s.trim().is_empty())
}

fn get_f64(ini: &Ini, section: &str, key: &str, default: f64) -> f64 {
    ini.get(section, key)
        .and_then(|v| v.parse().ok())
        .unwrap_or(default)
}

fn get_u32(ini: &Ini, section: &str, key: &str, default: u32) -> u32 {
    ini.get(section, key)
        .and_then(|v| v.parse().ok())
        .unwrap_or(default)
}

fn get_u64(ini: &Ini, section: &str, key: &str, default: u64) -> u64 {
    ini.get(section, key)
        .and_then(|v| v.parse().ok())
        .unwrap_or(default)
}

fn get_bool(ini: &Ini, section: &str, key: &str, default: bool) -> bool {
    match ini.get(section, key).as_deref().map(str::trim) {
        Some("1" | "true" | "yes" | "on") => true,
        Some("0" | "false" | "no" | "off") => false,
        Some(_) => default,
        None => default,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn loads_crosslinker_from_ini_string() {
        let mut ini = Ini::new();
        ini.read(
            r#"
[crosslinker]
name = dss
label = light-heavy
"#
            .into(),
        )
        .unwrap();
        let settings = Settings::from_ini(&ini);
        assert_eq!(settings.crosslinker_name, "dss");
        assert_eq!(settings.crosslinker_label, "light-heavy");
    }
}
