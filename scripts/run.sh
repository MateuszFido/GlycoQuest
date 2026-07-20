#!/usr/bin/env bash
# Submit (or locally run) a GlycoQuest search on a Slurm cluster (e.g. ETH Euler).
#
# Usage:
#   scripts/run.sh [wrapper options] -- <glycoquest args...>
#   scripts/run.sh [wrapper options] <glycoquest args...>
#
# Examples:
#   scripts/run.sh data/run.mzXML --database proteins.fasta --dry-run
#   scripts/run.sh --cpus 16 --time 08:00:00 -- \
#     data/run.mzXML --database proteins.fasta --crosslinker dss --out "$SCRATCH/gq_out"
#   scripts/run.sh --local data/run.mzXML --database proteins.fasta --jobs 4
#
# Wrapper options (not passed to GlycoQuest):
#   --cpus N          CPUs / parallel xQuest jobs (default: 16)
#   --time T          Wall time (default: 04:00:00)
#   --mem-per-cpu M   Memory per CPU (default: 4G)
#   --tmp SIZE        Node-local scratch request (default: 50G)
#   --job-name NAME   Slurm job name (default: glycoquest)
#   --mail-type TYPE  Slurm mail events (default: END,FAIL; empty to disable)
#   --partition NAME  Slurm partition (optional)
#   --account NAME    Slurm account (optional)
#   --local           Run on this machine; do not sbatch
#   --print           Print the sbatch command and exit
#   -h, --help        Show this help
#
# Environment:
#   GLYCOQUEST        Path to the glycoquest binary (optional)
#   GLYCOQUEST_STACK  Euler module stack (default: stack/2024-06)
#   GLYCOQUEST_PERL   Perl module (default: perl/5.38.0)
#
# Under Slurm, --jobs defaults to $SLURM_CPUS_PER_TASK and --progress to never
# when those flags are omitted. --xquest-root defaults to <repo>/V2.1.7/xquest
# when omitted and that directory exists.

set -euo pipefail

SCRIPT_PATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

CPUS=16
TIME=04:00:00
MEM_PER_CPU=4G
TMP=50G
JOB_NAME=glycoquest
MAIL_TYPE=END,FAIL
PARTITION=""
ACCOUNT=""
LOCAL=0
PRINT_ONLY=0
LOG_DIR=jobs

usage() {
  sed -n '2,36p' "$SCRIPT_PATH" | sed -E 's/^# ?//'
}

args_contain() {
  local needle=$1
  shift
  local arg
  for arg in "$@"; do
    if [[ "$arg" == "$needle" || "$arg" == "$needle"=* ]]; then
      return 0
    fi
    # clap-style short attached value: -j8
    if [[ "$needle" == "-j" && "$arg" =~ ^-j[0-9] ]]; then
      return 0
    fi
  done
  return 1
}

resolve_glycoquest() {
  if [[ -n "${GLYCOQUEST:-}" ]]; then
    if [[ ! -x "$GLYCOQUEST" ]]; then
      echo "error: GLYCOQUEST is set but not executable: $GLYCOQUEST" >&2
      exit 1
    fi
    printf '%s\n' "$GLYCOQUEST"
    return
  fi

  local candidate
  for candidate in \
    "$REPO_ROOT/target/release/glycoquest" \
    "$REPO_ROOT/target/debug/glycoquest"
  do
    if [[ -x "$candidate" ]]; then
      printf '%s\n' "$candidate"
      return
    fi
  done

  if command -v glycoquest >/dev/null 2>&1; then
    command -v glycoquest
    return
  fi

  echo "error: could not find glycoquest binary." >&2
  echo "  Build with: cargo build --release" >&2
  echo "  Or set GLYCOQUEST=/path/to/glycoquest" >&2
  exit 1
}

# --- parse wrapper options ---------------------------------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)
      usage
      exit 0
      ;;
    --cpus)
      CPUS=$2
      shift 2
      ;;
    --time)
      TIME=$2
      shift 2
      ;;
    --mem-per-cpu)
      MEM_PER_CPU=$2
      shift 2
      ;;
    --tmp)
      TMP=$2
      shift 2
      ;;
    --job-name)
      JOB_NAME=$2
      shift 2
      ;;
    --mail-type)
      MAIL_TYPE=$2
      shift 2
      ;;
    --partition)
      PARTITION=$2
      shift 2
      ;;
    --account)
      ACCOUNT=$2
      shift 2
      ;;
    --local)
      LOCAL=1
      shift
      ;;
    --print)
      PRINT_ONLY=1
      shift
      ;;
    --)
      shift
      break
      ;;
    -*)
      # First GlycoQuest flag (e.g. --database before INPUT is unusual but allowed by clap after INPUT).
      # Treat unknown dashed tokens as start of GlycoQuest argv.
      break
      ;;
    *)
      # Positional INPUT or other GlycoQuest argument.
      break
      ;;
  esac
done

if [[ $# -eq 0 ]]; then
  echo "error: missing GlycoQuest arguments." >&2
  echo "Run: scripts/run.sh --help" >&2
  exit 1
fi

if ! args_contain --database "$@"; then
  echo "error: GlycoQuest requires --database <FASTA>." >&2
  exit 1
fi

# --- submitter mode ----------------------------------------------------------
if [[ -z "${SLURM_JOB_ID:-}" && "$LOCAL" -eq 0 ]]; then
  sbatch_args=(
    --job-name="$JOB_NAME"
    --ntasks=1
    --cpus-per-task="$CPUS"
    --time="$TIME"
    --mem-per-cpu="$MEM_PER_CPU"
    --tmp="$TMP"
    --output="$LOG_DIR/${JOB_NAME}.%j.out"
    --error="$LOG_DIR/${JOB_NAME}.%j.err"
  )
  if [[ -n "$MAIL_TYPE" ]]; then
    sbatch_args+=(--mail-type="$MAIL_TYPE")
  fi
  if [[ -n "$PARTITION" ]]; then
    sbatch_args+=(--partition="$PARTITION")
  fi
  if [[ -n "$ACCOUNT" ]]; then
    sbatch_args+=(--account="$ACCOUNT")
  fi

  if [[ "$PRINT_ONLY" -eq 1 ]]; then
    printf 'sbatch'
    printf ' %q' "${sbatch_args[@]}" "$SCRIPT_PATH" "$@"
    printf '\n'
    exit 0
  fi

  if ! command -v sbatch >/dev/null 2>&1; then
    echo "error: sbatch not found. Use --local to run without Slurm, or run on a cluster login node." >&2
    exit 1
  fi

  GLYCOQUEST_BIN="$(resolve_glycoquest)"
  mkdir -p "$LOG_DIR"

  echo "Submitting GlycoQuest job (${CPUS} CPUs, time=${TIME}, mem-per-cpu=${MEM_PER_CPU})"
  echo "  binary: $GLYCOQUEST_BIN"
  sbatch "${sbatch_args[@]}" "$SCRIPT_PATH" "$@"
  exit 0
fi

# --- worker / local mode -----------------------------------------------------
GLYCOQUEST_BIN="$(resolve_glycoquest)"
extra=()

if ! args_contain --jobs "$@" && ! args_contain -j "$@"; then
  if [[ -n "${SLURM_CPUS_PER_TASK:-}" ]]; then
    extra+=(--jobs "$SLURM_CPUS_PER_TASK")
  elif [[ "$LOCAL" -eq 1 ]]; then
    extra+=(--jobs "$CPUS")
  fi
fi

if ! args_contain --progress "$@"; then
  extra+=(--progress never)
fi

default_xquest="$REPO_ROOT/V2.1.7/xquest"
if ! args_contain --xquest-root "$@" && [[ -d "$default_xquest" ]]; then
  extra+=(--xquest-root "$default_xquest")
fi

if [[ -n "${SLURM_JOB_ID:-}" ]]; then
  STACK="${GLYCOQUEST_STACK:-stack/2024-06}"
  PERL_MOD="${GLYCOQUEST_PERL:-perl/5.38.0}"

  if command -v module >/dev/null 2>&1; then
    # shellcheck disable=SC1091
    module load "$STACK" 2>/dev/null || true
    module load "$PERL_MOD"
  fi

  if ! perl -MDB_File -e 1 2>/dev/null; then
    echo "error: Perl DB_File is not available (required by xQuest indexing)." >&2
    echo "  Tried module: ${PERL_MOD} (stack: ${STACK})" >&2
    echo "  Verify with: module load ${STACK}; module load ${PERL_MOD}; perl -MDB_File -e 1" >&2
    exit 1
  fi

  export OMP_NUM_THREADS=1
  echo "GlycoQuest Slurm job ${SLURM_JOB_ID} on $(hostname)"
  echo "  cpus=${SLURM_CPUS_PER_TASK:-?}  tmpdir=${TMPDIR:-n/a}"
fi

echo "+ $(printf '%q ' "$GLYCOQUEST_BIN" "${extra[@]}" "$@")"
exec "$GLYCOQUEST_BIN" "${extra[@]}" "$@"
