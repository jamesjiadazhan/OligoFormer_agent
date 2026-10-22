#!/usr/bin/env bash
# Recreate the OligoFormer off-target tool prefix (perl + BioPerl + Statistics::Lite +
# ViennaRNA RNAplfold + GNU sed/grep) from the exact conda lock.
#
# Usage: deps/create_offtarget_env.sh [prefix]     (default: <PROJECT_ROOT>/deps/offtarget-env)
# Env:   CONDA_EXE  conda binary (default: conda on PATH; required)
#
# The lock (deps/offtarget-env.lock, `conda list --explicit --md5`) is platform osx-64
# (x86_64; runs under Rosetta 2 on Apple Silicon). Top-level packages it resolves from:
#   conda create -p <prefix> -c conda-forge -c bioconda \
#     perl=5.32.1 perl-bioperl=1.7.8 perl-statistics-lite=3.62 viennarna=2.7.2 sed=4.10 grep=3.12
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PREFIX="${1:-$SCRIPT_DIR/offtarget-env}"
LOCK="$SCRIPT_DIR/offtarget-env.lock"
CONDA_EXE="${CONDA_EXE:-$(command -v conda || true)}"
[[ -n "$CONDA_EXE" && -x "$CONDA_EXE" ]] || { echo "[offtarget-env] conda not found; install Miniconda/Miniforge or set CONDA_EXE" >&2; exit 1; }

# Exercise every tool the off-target pipeline actually runs; -x alone is not enough (a freshly
# created prefix can have an executable RNAplfold whose first invocation still aborts in dyld).
smoke_ok() {
  [[ -x "$PREFIX/bin/perl" ]] || return 1
  "$PREFIX/bin/perl" -e 'use Statistics::Lite; use Bio::TreeIO;' 2>/dev/null || return 1
  "$PREFIX/bin/RNAplfold" --version >/dev/null 2>&1 || return 1
  "$PREFIX/bin/sed" --version 2>/dev/null | grep -q 'GNU sed' || return 1
  "$PREFIX/bin/grep" --version 2>/dev/null | grep -q 'GNU grep' || return 1
  return 0
}

if smoke_ok; then
  echo "[offtarget-env] already complete: $PREFIX"
  exit 0
fi

echo "[offtarget-env] creating $PREFIX from $LOCK"
CONDA_SUBDIR=osx-64 "$CONDA_EXE" create -y -p "$PREFIX" --file "$LOCK"

# The first call into a just-created prefix can abort while dyld resolves @rpath libraries
# against the package cache; a retry after that resolves it. Fail only if it persists.
if ! smoke_ok; then
  echo "[offtarget-env] first smoke check failed; retrying once"
  sleep 2
  smoke_ok || { echo "[offtarget-env] prefix $PREFIX is not usable (perl/RNAplfold/sed/grep check failed)" >&2; exit 1; }
fi
"$PREFIX/bin/perl" -e 'print "perl modules OK\n"'
"$PREFIX/bin/RNAplfold" --version | head -1
"$PREFIX/bin/sed" --version | head -1
"$PREFIX/bin/grep" --version | head -1
echo "[offtarget-env] ready: $PREFIX"
