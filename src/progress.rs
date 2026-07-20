// Copyright (c) ETH Zurich, Mateusz Fido

//! Terminal-aware progress rendering shared by the Rust and xQuest phases.

use std::io::{self, IsTerminal};
use std::time::Duration;

use indicatif::{ProgressBar, ProgressDrawTarget, ProgressState, ProgressStyle};

use crate::cli::ProgressMode;

const REFRESH_HZ: u8 = 10;

#[derive(Debug, Clone)]
pub(crate) struct ProgressReporter {
    enabled: bool,
    force: bool,
}

impl ProgressReporter {
    pub(crate) fn new(mode: ProgressMode) -> Self {
        let enabled = match mode {
            ProgressMode::Auto => terminal_supports_progress(),
            ProgressMode::Always => true,
            ProgressMode::Never => false,
        };
        Self {
            enabled,
            force: mode == ProgressMode::Always,
        }
    }

    pub(crate) fn spinner(&self, step: usize, steps: usize, label: &str) -> PhaseProgress {
        let bar = if self.force {
            ProgressBar::with_draw_target(
                None,
                ProgressDrawTarget::term_like_with_hz(
                    Box::new(console::Term::buffered_stderr()),
                    REFRESH_HZ,
                ),
            )
        } else if self.enabled {
            ProgressBar::with_draw_target(None, ProgressDrawTarget::stderr_with_hz(REFRESH_HZ))
        } else {
            ProgressBar::hidden()
        };
        bar.set_style(spinner_style());
        bar.set_prefix(format!("[{step}/{steps}] {label}"));
        if self.enabled {
            bar.enable_steady_tick(Duration::from_millis(100));
        }
        PhaseProgress {
            bar,
            enabled: self.enabled,
        }
    }

    pub(crate) fn determinate(
        &self,
        step: usize,
        steps: usize,
        label: &str,
        length: u64,
    ) -> PhaseProgress {
        let phase = self.spinner(step, steps, label);
        phase.make_determinate(length);
        phase
    }
}

fn terminal_supports_progress() -> bool {
    io::stderr().is_terminal()
        && std::env::var("TERM").is_ok_and(|term| !term.is_empty() && term != "dumb")
}

#[derive(Debug, Clone)]
pub(crate) struct PhaseProgress {
    bar: ProgressBar,
    enabled: bool,
}

impl PhaseProgress {
    pub(crate) fn enabled(&self) -> bool {
        self.enabled
    }

    pub(crate) fn make_determinate(&self, length: u64) {
        self.bar.set_position(0);
        self.bar.set_length(length.max(1));
        self.bar.set_style(determinate_style());
        self.bar.reset_elapsed();
        self.bar.reset_eta();
    }

    pub(crate) fn set_position(&self, position: u64) {
        self.bar.set_position(position);
    }

    pub(crate) fn inc(&self, delta: u64) {
        self.bar.inc(delta);
    }

    pub(crate) fn set_message(&self, message: impl Into<String>) {
        self.bar.set_message(message.into());
    }

    pub(crate) fn println(&self, message: impl AsRef<str>) {
        if self.enabled {
            let _ = self.bar.println(message.as_ref());
        } else {
            eprintln!("{}", message.as_ref());
        }
    }

    pub(crate) fn finish(&self, message: impl Into<String>) {
        self.bar.set_position(self.bar.length().unwrap_or(0));
        self.bar.finish_with_message(message.into());
    }

    pub(crate) fn abandon(&self, message: impl Into<String>) {
        self.bar.abandon_with_message(message.into());
    }

    #[cfg(test)]
    pub(crate) fn test_snapshot(&self) -> (u64, String) {
        (self.bar.position(), self.bar.message().to_string())
    }
}

fn spinner_style() -> ProgressStyle {
    ProgressStyle::with_template("{prefix:.bold} {spinner:.green} {elapsed_precise} {wide_msg}")
        .expect("valid spinner progress template")
        .tick_strings(&["⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏"])
}

fn determinate_style() -> ProgressStyle {
    ProgressStyle::with_template(
        "{prefix:.bold} [{bar:20.cyan/blue}] {pos}/{len} · {elapsed_precise} · ETA {phase_eta}\n  {spinner:.green} {rate} · {wide_msg}",
    )
    .expect("valid determinate progress template")
    .with_key(
        "rate",
        |state: &ProgressState, writer: &mut dyn std::fmt::Write| {
            let rate = state.per_sec();
            if rate >= 1_000_000.0 {
                let _ = write!(writer, "{:.1}M/s", rate / 1_000_000.0);
            } else if rate >= 1_000.0 {
                let _ = write!(writer, "{:.1}k/s", rate / 1_000.0);
            } else if rate >= 10.0 {
                let _ = write!(writer, "{rate:.0}/s");
            } else {
                let _ = write!(writer, "{rate:.1}/s");
            }
        },
    )
    .with_key(
        "phase_eta",
        |state: &ProgressState, writer: &mut dyn std::fmt::Write| {
            if state.pos() == 0 {
                let _ = writer.write_str("—");
            } else {
                let seconds = state.eta().as_secs();
                if seconds >= 3600 {
                    let _ = write!(writer, "{}h {:02}m", seconds / 3600, (seconds % 3600) / 60);
                } else if seconds >= 60 {
                    let _ = write!(writer, "{}m {:02}s", seconds / 60, seconds % 60);
                } else {
                    let _ = write!(writer, "{seconds}s");
                }
            }
        },
    )
    .progress_chars("█▓░")
    .tick_strings(&["⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏"])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn explicit_progress_modes_override_terminal_detection() {
        assert!(ProgressReporter::new(ProgressMode::Always).enabled);
        assert!(!ProgressReporter::new(ProgressMode::Never).enabled);
    }

    #[test]
    fn changing_determinate_work_units_resets_position() {
        let progress = ProgressReporter::new(ProgressMode::Never);
        let phase = progress.determinate(2, 4, "Preparing xQuest jobs", 4);
        phase.inc(4);
        phase.make_determinate(3);
        assert_eq!(phase.test_snapshot().0, 0);
    }
}
