// Copyright (c) ETH Zurich, Mateusz Fido

//! Crosslinker preset construction and verification.

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
                if CrosslinkerPreset::by_name(&settings.crosslinker_name).is_some() {
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

/// Construct profiles from settings.ini file.
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
    // test for the pseudoresidues from example dataset
    fn glycan_protein_link_glycan_to_lysine() {
        let settings = Settings::defaults();
        for (name, mass) in [
            ("nhs-cyclooctyne", 205.085126607),
            ("ssbxl", 573.179438173),
            ("pcbxl", 456.190988659),
        ] {
            let profile = CrosslinkerProfile::resolve(&settings, Some(name)).unwrap();
            assert_eq!(profile.label, CrosslinkerLabel::None);
            assert_eq!(profile.xlink_sites, "X:K");
            assert!((profile.xlinkermw - mass).abs() < 1e-9);
            assert!(!profile.requires_isotope_pair_prefilter());
        }
    }

    #[test]
    fn ssbxl_bridge_signature_ion() {
        const NEUAC_RESIDUE: f64 = 291.095416527;
        const PROTON: f64 = 1.007276466621;
        let signature_mz = 573.179438173 + NEUAC_RESIDUE + PROTON;
        assert!((signature_mz - 865.282131166).abs() < 1e-9);
    }

    #[test]
    fn settings_name_selects_glycan_protein_preset() {
        let mut settings = Settings::defaults();
        settings.crosslinker_name = "ssbxl".into();
        let profile = CrosslinkerProfile::resolve(&settings, None).unwrap();
        assert!((profile.xlinkermw - 573.179438173).abs() < 1e-9);
        assert_eq!(profile.xlink_sites, "X:K");
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
