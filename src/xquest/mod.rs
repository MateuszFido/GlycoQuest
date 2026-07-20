// Copyright (c) ETH Zurich, Mateusz Fido

//! xQuest runtime discovery, job generation, and definitions.

mod defs;
mod generate;
mod matchlist;
mod perl_deps;
mod run;
mod runtime;

pub use defs::{JobDefs, write_job_defs};
pub use generate::{GeneratedJob, generate_jobs};
pub use matchlist::{MatchlistRow, write_matchlist};
pub use perl_deps::check_perl_search_deps;
pub(crate) use run::execute_jobs_with_progress;
pub use run::{JobRunRecord, log_job_summary};
pub use runtime::{XQuestRuntime, resolve_runtime};
