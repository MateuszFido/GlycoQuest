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
#   GLYCOQUEST              Path to the glycoquest binary (optional)
#   GLYCOQUEST_STACK        Euler module stack (default: stack/2024-04)
#   GLYCOQUEST_COMPILER     Compiler module for Perl (default: gcc/8.5.0)
#   GLYCOQUEST_PERL         Perl module (default: perl/5.38.0)
#   GLYCOQUEST_BERKELEY_DB  Berkeley DB module: auto (default), a full module
#                           name, or empty to skip
#   GLYCOQUEST_PERL5        local::lib / cpanm prefix with DB_File
#                           (default: $HOME/perl5 if it exists)
#
# Under Slurm, --jobs defaults to $SLURM_CPUS_PER_TASK and --progress to never
# when those flags are omitted. --xquest-root defaults to <repo>/V2.1.7/xquest
# when omitted and that directory exists.
#
# Euler notes:
#   - perl/5.38.0 needs stack/2024-04 + gcc/8.5.0 (`module spider perl/5.38.0`).
#   - berkeley-db is hierarchy-specific (hashed names). Default "auto" tries to
#     load one if visible; if not, continues (OK when DB_File already has RPATH).
#   - Spack Perl has no DB_File; install once against module perl + berkeley-db:
#       module load stack/2024-04 gcc/8.5.0 perl/5.38.0 eth_proxy
#       cpanm --local-lib=$HOME/perl5 DB_File
#   - DB_File.so needs libdb-18.1.so at runtime. If the berkeley-db module is not
#     visible, this script globs the stack's Spack tree and prepends the lib dir
#     to LD_LIBRARY_PATH.
#   - xQuest job scripts overwrite PERL5LIB, so this wrapper exposes a local
#     install via PERL5OPT=-I... (not PERL5LIB alone).

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
  sed -n '2,52p' "$SCRIPT_PATH" | sed -E 's/^# ?//'
}

# Print berkeley-db modules visible in the *current* Lmod hierarchy (stdout).
# Prefer Spack-hashed names; the bare berkeley-db/X.Y.Z often fails to load.
available_berkeley_db_modules() {
  local all hashed bare
  all=$(
    {
      module -t avail berkeley-db 2>&1 || true
      module --show_hidden -t avail berkeley-db 2>&1 || true
      module -t avail 2>&1 || true
    } | sed -E 's/\([^)]*\)//g; s/^[[:space:]]+//; s/[[:space:]]+$//' \
      | grep -E '^berkeley-db/' \
      | grep -v '/$' \
      | sort -u
  )
  [[ -n "$all" ]] || return 0
  hashed=$(echo "$all" | grep -E 'berkeley-db/[0-9.]+-[A-Za-z0-9]+$' || true)
  bare=$(echo "$all" | grep -vE 'berkeley-db/[0-9.]+-[A-Za-z0-9]+$' || true)
  [[ -n "$hashed" ]] && echo "$hashed"
  [[ -n "$bare" ]] && echo "$bare"
}

# Load a berkeley-db module. Arg: auto | empty(skip) | explicit module name.
# Sets BERKELEY_DB to the module that was loaded (or empty if skipped).
# In "auto" mode, failure is soft: DB_File may already be linked with an RPATH.
load_berkeley_db() {
  local requested=$1
  BERKELEY_DB=""

  if [[ -z "$requested" ]]; then
    return 0
  fi

  local candidates=()
  local soft=0
  if [[ "$requested" == "auto" ]]; then
    soft=1
  else
    candidates+=("$requested")
  fi

  local mod
  while IFS= read -r mod; do
    [[ -n "$mod" ]] || continue
    local seen=0
    local c
    for c in "${candidates[@]+"${candidates[@]}"}"; do
      if [[ "$c" == "$mod" ]]; then
        seen=1
        break
      fi
    done
    if [[ "$seen" -eq 0 ]]; then
      candidates+=("$mod")
    fi
  done < <(available_berkeley_db_modules)

  # Last-resort names seen on Euler for common toolchains (may 404 — that's OK).
  if [[ "$soft" -eq 1 ]]; then
    for mod in \
      berkeley-db/18.1.40-c3kwxwi \
      berkeley-db/18.1.40-v4v4ea3 \
      berkeley-db/18.1.40-cbtnnjh \
      berkeley-db/18.1.40
    do
      local seen=0
      local c
      for c in "${candidates[@]+"${candidates[@]}"}"; do
        [[ "$c" == "$mod" ]] && seen=1 && break
      done
      [[ "$seen" -eq 0 ]] && candidates+=("$mod")
    done
  fi

  if [[ ${#candidates[@]} -eq 0 ]]; then
    if [[ "$soft" -eq 1 ]]; then
      echo "warning: no berkeley-db module visible after ${STACK} ${COMPILER}; continuing." >&2
      echo "  If DB_File later fails to load libdb, set GLYCOQUEST_BERKELEY_DB from:" >&2
      echo "    module spider berkeley-db" >&2
      return 0
    fi
    echo "error: no berkeley-db module visible after ${STACK} ${COMPILER}." >&2
    echo "  Run: module spider berkeley-db" >&2
    return 1
  fi

  local err=""
  for mod in "${candidates[@]}"; do
    if err=$(module load "$mod" 2>&1); then
      BERKELEY_DB=$mod
      return 0
    fi
  done

  if [[ "$soft" -eq 1 ]]; then
    echo "warning: could not load berkeley-db after ${STACK} ${COMPILER}; continuing." >&2
    echo "  Tried: ${candidates[*]}" >&2
    echo "  If perl -MDB_File fails on libdb, set GLYCOQUEST_BERKELEY_DB explicitly." >&2
    return 0
  fi

  echo "error: could not load berkeley-db module '${requested}' after ${STACK} ${COMPILER}." >&2
  echo "  Last module error:" >&2
  echo "$err" >&2
  echo "  Run: module spider berkeley-db" >&2
  return 1
}

# If DB_File.so was built against Spack berkeley-db, ensure libdb-*.so is findable.
# Sets LD_LIBRARY_PATH when a matching lib is found under the active stack.
ensure_libdb_runtime_path() {
  # Already resolvable?
  if perl -MDB_File -e 1 2>/dev/null; then
    return 0
  fi

  local stack_name=${STACK#stack/}
  local stack_root="/cluster/software/stacks/${stack_name}"
  local lib=""
  local candidate

  # Prefer the compiler directory that matches the loaded gcc module.
  local gcc_ver=${COMPILER#gcc/}
  shopt -s nullglob
  local globs=(
    "${stack_root}/spack/opt/spack/"*"/gcc-${gcc_ver}/berkeley-db-"*/lib/libdb-18.1.so
    "${stack_root}/spack/opt/spack/"*"/gcc-${gcc_ver}/berkeley-db-"*/lib/libdb.so
    "${stack_root}/spack/opt/spack/"*"/berkeley-db-"*/lib/libdb-18.1.so
  )
  for candidate in "${globs[@]}"; do
    if [[ -f "$candidate" ]]; then
      lib=$candidate
      break
    fi
  done
  shopt -u nullglob

  if [[ -z "$lib" ]]; then
    return 1
  fi

  local libdir
  libdir=$(dirname "$lib")
  case ":${LD_LIBRARY_PATH:-}:" in
    *":${libdir}:"*) ;;
    *)
      export LD_LIBRARY_PATH="${libdir}${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
      echo "note: added berkeley-db libs to LD_LIBRARY_PATH: ${libdir}" >&2
      ;;
  esac
  return 0
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
  # perl/5.38.0 on Euler requires stack/2024-04 + gcc/8.5.0 (module spider).
  # berkeley-db hashes differ per toolchain; default "auto" picks from module avail.
  STACK="${GLYCOQUEST_STACK:-stack/2024-04}"
  # Single-dash defaults: empty GLYCOQUEST_COMPILER= / GLYCOQUEST_BERKELEY_DB= skips.
  COMPILER="${GLYCOQUEST_COMPILER-gcc/8.5.0}"
  PERL_MOD="${GLYCOQUEST_PERL:-perl/5.38.0}"
  BERKELEY_DB_REQUEST="${GLYCOQUEST_BERKELEY_DB-auto}"
  BERKELEY_DB=""

  if command -v module >/dev/null 2>&1; then
    module load "$STACK"
    if [[ -n "$COMPILER" ]]; then
      module load "$COMPILER"
    fi
    # Load berkeley-db before Perl so LD_LIBRARY_PATH is set for DB_File.so.
    load_berkeley_db "$BERKELEY_DB_REQUEST"
    module load "$PERL_MOD"
  else
    echo "warning: module command not found; using whatever perl is on PATH" >&2
  fi

  # Spack Perl has no DB_File; a cpanm install under $HOME/perl5 is the usual fix.
  # xQuest run.sh overwrites PERL5LIB, so expose the local tree via PERL5OPT.
  PERL5_ROOT="${GLYCOQUEST_PERL5-}"
  if [[ -z "$PERL5_ROOT" && -d "${HOME}/perl5/lib/perl5" ]]; then
    PERL5_ROOT="${HOME}/perl5"
  fi
  if [[ -n "$PERL5_ROOT" && -d "${PERL5_ROOT}/lib/perl5" ]]; then
    perl5_lib="${PERL5_ROOT}/lib/perl5"
    case ":${PERL5LIB:-}:" in
      *":${perl5_lib}:"*) ;;
      *) export PERL5LIB="${perl5_lib}${PERL5LIB:+:$PERL5LIB}" ;;
    esac
    case " ${PERL5OPT:-} " in
      *" -I${perl5_lib} "*|*" -I${perl5_lib}"*) ;;
      *) export PERL5OPT="-I${perl5_lib}${PERL5OPT:+ ${PERL5OPT}}" ;;
    esac
  fi

  # Module may be invisible in this hierarchy; still need libdb for the XS .so.
  ensure_libdb_runtime_path || true

  if ! perl -MDB_File -e 1 2>/dev/null; then
    db_err=$(perl -MDB_File -e 1 2>&1 || true)
    echo "error: Perl DB_File is not available (required by xQuest indexing)." >&2
    echo "  Loaded: ${STACK} ${COMPILER} ${PERL_MOD} berkeley-db=${BERKELEY_DB:-none}" >&2
    echo "  PERL5OPT=${PERL5OPT:-<unset>}  LD_LIBRARY_PATH=${LD_LIBRARY_PATH:-<unset>}" >&2
    echo "  perl said: ${db_err}" >&2
    echo "  If the module is installed but libdb is missing:" >&2
    echo "    ls /cluster/software/stacks/${STACK#stack/}/spack/opt/spack/*/gcc-*/berkeley-db-*/lib/libdb-18.1.so" >&2
    echo "    export LD_LIBRARY_PATH=<that-lib-dir>:\$LD_LIBRARY_PATH" >&2
    echo "  Or install once (with berkeley-db / libdb available), then resubmit:" >&2
    echo "    module load ${STACK} ${COMPILER} ${PERL_MOD} eth_proxy" >&2
    echo "    cpanm --local-lib=\$HOME/perl5 DB_File" >&2
    exit 1
  fi

  export OMP_NUM_THREADS=1
  echo "GlycoQuest Slurm job ${SLURM_JOB_ID} on $(hostname)"
  echo "  modules: ${STACK} ${COMPILER} ${BERKELEY_DB:-none} ${PERL_MOD}"
  if [[ -n "${PERL5_ROOT:-}" ]]; then
    echo "  perl5 local: ${PERL5_ROOT} (PERL5OPT=${PERL5OPT:-})"
  fi
  if [[ -n "${LD_LIBRARY_PATH:-}" ]]; then
    echo "  LD_LIBRARY_PATH includes berkeley-db (len=${#LD_LIBRARY_PATH})"
  fi
  echo "  cpus=${SLURM_CPUS_PER_TASK:-?}  tmpdir=${TMPDIR:-n/a}"
fi

echo "+ $(printf '%q ' "$GLYCOQUEST_BIN" "${extra[@]}" "$@")"
exec "$GLYCOQUEST_BIN" "${extra[@]}" "$@"
