//! Crosslinker profiles and bundled chemistry presets.

mod preset;

use crate::cli::settings::Settings;

pub use preset::{CrosslinkerLabel, CrosslinkerPreset};

/// Resolved crosslinker configuration used by prefilter, job generation, and post-filters.
#[derive(Debug, Clone, PartialEq)]
pub struct CrosslinkerProfile {
    pub name: String,
    pub label: CrosslinkerLabel,
    pub shift_da: f64,
    pub xlinkermw: f64,
    pub xlink_sites: String,
    pub nterm_xlinkable: bool,
}

impl CrosslinkerProfile {
    /// Build a profile from settings, applying a bundled preset when `cli_name` matches.
    pub fn resolve(settings: &Settings, cli_name: Option<&str>) -> Result<Self, String> {
        let preset_name = cli_name
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .or_else(|| {
                if settings.crosslinker_name.trim().eq_ignore_ascii_case("dss")
                    || settings
                        .crosslinker_name
                        .trim()
                        .eq_ignore_ascii_case("dmtmm")
                {
                    Some(settings.crosslinker_name.as_str())
                } else {
                    None
                }
            });

        let mut profile = if let Some(name) = preset_name {
            if let Some(preset) = CrosslinkerPreset::by_name(name) {
                preset.to_profile()
            } else {
                profile_from_settings(settings)
            }
        } else {
            profile_from_settings(settings)
        };

        if let Some(name) = cli_name.filter(|s| !s.trim().is_empty()) {
            profile.name = name.trim().to_string();
        }

        profile.validate()?;
        Ok(profile)
    }

    pub fn requires_isotope_pair_prefilter(&self) -> bool {
        self.label == CrosslinkerLabel::LightHeavy
    }

    pub fn validate(&self) -> Result<(), String> {
        if self.label == CrosslinkerLabel::LightHeavy && self.shift_da <= 0.0 {
            return Err(format!(
                "crosslinker label light-heavy requires shift_da > 0 (got {})",
                self.shift_da
            ));
        }
        Ok(())
    }
}

fn profile_from_settings(settings: &Settings) -> CrosslinkerProfile {
    CrosslinkerProfile {
        name: settings.crosslinker_name.clone(),
        label: CrosslinkerLabel::parse(&settings.crosslinker_label),
        shift_da: settings.crosslinker_shift_da,
        xlinkermw: settings.xlinkermw,
        xlink_sites: settings.xlink_sites.clone(),
        nterm_xlinkable: settings.nterm_xlinkable,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dmtmm_preset_disables_isotope_prefilter() {
        let settings = Settings::defaults();
        let profile = CrosslinkerProfile::resolve(&settings, Some("dmtmm")).unwrap();
        assert_eq!(profile.name, "dmtmm");
        assert_eq!(profile.label, CrosslinkerLabel::None);
        assert!(!profile.requires_isotope_pair_prefilter());
        assert!((profile.xlinkermw - (-18.0109)).abs() < 1e-4);
        assert_eq!(profile.xlink_sites, "K:E,K:D");
    }

    #[test]
    fn dss_preset_requires_isotope_prefilter() {
        let settings = Settings::defaults();
        let profile = CrosslinkerProfile::resolve(&settings, Some("dss")).unwrap();
        assert!(profile.requires_isotope_pair_prefilter());
        assert!((profile.xlinkermw - 138.0680796).abs() < 1e-4);
    }

    #[test]
    fn rejects_light_heavy_with_zero_shift() {
        let mut settings = Settings::defaults();
        settings.crosslinker_name = "custom".into();
        settings.crosslinker_label = "light-heavy".into();
        settings.crosslinker_shift_da = 0.0;
        let err = CrosslinkerProfile::resolve(&settings, None).unwrap_err();
        assert!(err.contains("shift_da"));
    }
}
