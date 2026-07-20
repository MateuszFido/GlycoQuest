#!/usr/bin/env bash
# One-shot install of XS Perl modules required by GlycoQuest/xQuest on ETH Euler.
# Run on a login node (needs eth_proxy for CPAN).
#
#   scripts/bootstrap-euler-perl.sh
#   scripts/bootstrap-euler-perl.sh && scripts/check-xquest-perl.pl
#
# Installs into $HOME/perl5 (override with GLYCOQUEST_PERL5):
#   DB_File, XML::Parser, and GD/GD::Graph when the bundled 1209/lib64 GD.so
#   cannot load against the loaded module perl.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
STACK="${GLYCOQUEST_STACK:-stack/2024-04}"
COMPILER="${GLYCOQUEST_COMPILER:-gcc/8.5.0}"
PERL_MOD="${GLYCOQUEST_PERL:-perl/5.38.0}"
PERL5_ROOT="${GLYCOQUEST_PERL5:-$HOME/perl5}"
XQUEST_ROOT="${GLYCOQUEST_XQUEST_ROOT:-$REPO_ROOT/V2.1.7/xquest}"

if ! command -v module >/dev/null 2>&1; then
  echo "error: Lmod 'module' command not found (are you on Euler?)" >&2
  exit 1
fi

module load "$STACK"
[[ -n "$COMPILER" ]] && module load "$COMPILER"
module load eth_proxy 2>/dev/null || true
if module spider berkeley-db >/dev/null 2>&1; then
  module load "${GLYCOQUEST_BERKELEY_DB:-berkeley-db}" 2>/dev/null \
    || echo "note: could not module-load berkeley-db; continuing" >&2
fi
# Soft-load gd/libgd when visible (needed to build/run GD.pm).
for gdmod in gd libgd; do
  module load "$gdmod" 2>/dev/null && break || true
done
module load "$PERL_MOD"

if ! command -v cpanm >/dev/null 2>&1; then
  echo "error: cpanm not found after loading ${PERL_MOD}." >&2
  echo "  Install App::cpanminus once, or use: curl -L https://cpanmin.us | perl - App::cpanminus" >&2
  exit 1
fi

export PERL5OPT="-I${PERL5_ROOT}/lib/perl5${PERL5OPT:+ ${PERL5OPT}}"
# Prefer bundled 1209 paths (including lib64 where legacy GD lives).
export PERL5LIB="${XQUEST_ROOT}/1209/lib64/perl5:${XQUEST_ROOT}/1209/lib/perl5:${XQUEST_ROOT}/1209/share/perl5:${XQUEST_ROOT}/modules${PERL5LIB:+:$PERL5LIB}"

stack_name=${STACK#stack/}
gcc_ver=${COMPILER#gcc/}
shopt -s nullglob
for so in \
  "/cluster/software/stacks/${stack_name}/spack/opt/spack/"*"/gcc-${gcc_ver}/berkeley-db-"*/lib/libdb-18.1.so \
  "/cluster/software/stacks/${stack_name}/spack/opt/spack/"*"/gcc-${gcc_ver}/libexpat-"*/lib/libexpat.so* \
  "/cluster/software/stacks/${stack_name}/spack/opt/spack/"*"/gcc-${gcc_ver}/libgd-"*/lib/libgd.so* \
  "/cluster/software/stacks/${stack_name}/spack/opt/spack/"*"/gcc-${gcc_ver}/gd-"*/lib/libgd.so*
do
  [[ -f "$so" ]] || continue
  libdir=$(dirname "$so")
  case ":${LD_LIBRARY_PATH:-}:" in
    *":${libdir}:"*) ;;
    *) export LD_LIBRARY_PATH="${libdir}${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}" ;;
  esac
done
shopt -u nullglob

pkgs=(DB_File XML::Parser)
# Rebuild GD when bundled lib64 GD.so is missing or ABI-incompatible with this perl.
if ! perl -MGD -e 1 2>/dev/null; then
  echo "note: bundled/system GD not loadable; will cpanm GD GD::Graph"
  pkgs+=(GD GD::Graph)
fi

echo "Installing into ${PERL5_ROOT}: ${pkgs[*]}"
cpanm --local-lib="$PERL5_ROOT" "${pkgs[@]}"

echo "Verifying..."
perl -MDB_File -MXML::Parser -MGD -e 'print "DB_File + XML::Parser + GD OK\n"'
"$REPO_ROOT/scripts/check-xquest-perl.pl" --xquest-root "$XQUEST_ROOT"
echo "Done. Slurm jobs pick up \$HOME/perl5 via scripts/run.sh (PERL5OPT)."
