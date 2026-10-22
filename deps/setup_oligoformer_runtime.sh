#!/usr/bin/env bash
# Build a self-contained, runnable OligoFormer runtime copy (pinned source + RNA-FM + off-target
# tooling) WITHOUT modifying the pristine checkout under repo/OligoFormer.
#
# Usage: deps/setup_oligoformer_runtime.sh <target_dir>
#
# Result (<target_dir>):
#   - git checkout of lulab/OligoFormer @ e2f53ad63387bbe166bf123949151e2bc9bf6ec3
#   - RNA-FM/            ml4bio/RNA-FM @ 348951516e0963d22bbb33b3c9fc18c89081d38e (+ CPU-default patch),
#                        weights at RNA-FM/redevelop/pretrained/RNA-FM_pretrained.pth (sha256 verified),
#                        RNA-FM/torch_home/hub/checkpoints/RNA-FM_pretrained.pth -> that weights file
#   - off-target/pita    arm64 RNAduplex/RNAddG4 staged + `make install` (EXE_BASE_DIR -> absolute path)
#   - <target_dir>.env.sh  variables to source before running scripts/main.py
#
# Idempotent: re-running resets off-target/pita to the pinned state and re-applies the install,
# re-applies the RNA-FM patch from a clean file, and skips verified weights.
# Moving <target_dir> requires re-running (make install bakes in the absolute path).
#
# Optional environment variables:
#   OLIGOFORMER_SRC  local git clone to clone from (must contain the pinned commit; HEAD verified)
#   RNAFM_SRC        local git clone of ml4bio/RNA-FM (must contain the pinned commit)
#   RNAFM_WEIGHTS    existing RNA-FM_pretrained.pth (sha256 verified; hard-linked, else copied)
#   PROJECT_ENV      Python env whose bin/ provides `python` (default <PROJECT_ROOT>/OligoFormer-env)
#   OFFTARGET_ENV    conda prefix from deps/create_offtarget_env.sh (default <PROJECT_ROOT>/deps/offtarget-env)
#   PITA_BIN_DIR     dir with prebuilt arm64 RNAduplex + RNAddG4 (default <PROJECT_ROOT>/deps/pita-bin-arm64);
#                    if missing, or PITA_FORCE_BUILD=1, they are rebuilt from the bundled ViennaRNA-1.6 source.
#   SKIP_OFFTARGET=1 skip the PITA staging / make install step.
set -euo pipefail

OLIGOFORMER_URL="https://github.com/lulab/OligoFormer"
OLIGOFORMER_COMMIT="e2f53ad63387bbe166bf123949151e2bc9bf6ec3"
RNAFM_URL="https://github.com/ml4bio/RNA-FM.git"
RNAFM_COMMIT="348951516e0963d22bbb33b3c9fc18c89081d38e"
RNAFM_WEIGHTS_URL="https://huggingface.co/cuhkaih/rnafm/resolve/main/RNA-FM_pretrained.pth"
RNAFM_WEIGHTS_SHA256="5b5d7d87b37c291ef42c140ef9edf7aea29f255fa2a4fd435f776c52e93d5e99"
RNADUPLEX_ARM64_SHA256="b01a67068a5168f87df18a967c72398d0ea4d6ff6847f055ee6b743aea2c0df8"
RNADDG4_ARM64_SHA256="ada0a5af757ae3a086adda5dc5968d1b7703bcb352c706f3f61e94e64f28c607"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
PATCH_FILE="$SCRIPT_DIR/rnafm_cpu_default.patch"
PROJECT_ENV="${PROJECT_ENV:-$PROJECT_ROOT/OligoFormer-env}"
OFFTARGET_ENV="${OFFTARGET_ENV:-$PROJECT_ROOT/deps/offtarget-env}"
PITA_BIN_DIR="${PITA_BIN_DIR:-$PROJECT_ROOT/deps/pita-bin-arm64}"

log() { printf '[setup-runtime] %s\n' "$*" >&2; }
die() { printf '[setup-runtime] ERROR: %s\n' "$*" >&2; exit 1; }
sha256() { shasum -a 256 "$1" | awk '{print $1}'; }

[[ $# -eq 1 ]] || die "usage: $0 <target_dir>"
mkdir -p "$(dirname "$1")"
TARGET="$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
[[ -x "$PROJECT_ENV/bin/python" ]] || die "PROJECT_ENV python not found: $PROJECT_ENV/bin/python"

# ---------------------------------------------------------------- 1. OligoFormer source
if [[ -d "$TARGET/.git" ]]; then
  log "reusing existing checkout $TARGET"
else
  [[ ! -e "$TARGET" ]] || die "$TARGET exists but is not a git checkout"
  SRC="${OLIGOFORMER_SRC:-$OLIGOFORMER_URL}"
  if [[ -n "${OLIGOFORMER_SRC:-}" ]]; then
    head="$(git -C "$OLIGOFORMER_SRC" rev-parse HEAD)"
    [[ "$head" == "$OLIGOFORMER_COMMIT" ]] || die "OLIGOFORMER_SRC HEAD $head != pinned $OLIGOFORMER_COMMIT"
  fi
  log "cloning OligoFormer from $SRC"
  git clone --quiet --no-hardlinks "$SRC" "$TARGET"
  git -C "$TARGET" remote set-url origin "$OLIGOFORMER_URL"
fi
git -C "$TARGET" checkout --quiet --detach "$OLIGOFORMER_COMMIT"
[[ "$(git -C "$TARGET" rev-parse HEAD)" == "$OLIGOFORMER_COMMIT" ]] || die "target HEAD is not pinned commit"

# ---------------------------------------------------------------- 2. RNA-FM (+ CPU patch)
RNAFM_DIR="$TARGET/RNA-FM"
if [[ ! -d "$RNAFM_DIR/.git" ]]; then
  [[ ! -e "$RNAFM_DIR" ]] || die "$RNAFM_DIR exists but is not a git checkout"
  log "cloning RNA-FM from ${RNAFM_SRC:-$RNAFM_URL}"
  git clone --quiet --no-hardlinks "${RNAFM_SRC:-$RNAFM_URL}" "$RNAFM_DIR"
  git -C "$RNAFM_DIR" remote set-url origin "$RNAFM_URL"
fi
git -C "$RNAFM_DIR" checkout --quiet --detach "$RNAFM_COMMIT"
git -C "$RNAFM_DIR" checkout --quiet -- redevelop/launch/predict.py
git -C "$RNAFM_DIR" apply "$PATCH_FILE"
log "RNA-FM @ $(git -C "$RNAFM_DIR" rev-parse HEAD) with $(basename "$PATCH_FILE") applied"

# ---------------------------------------------------------------- 3. RNA-FM weights
WEIGHTS_DST="$RNAFM_DIR/redevelop/pretrained/RNA-FM_pretrained.pth"
if [[ -f "$WEIGHTS_DST" && "$(sha256 "$WEIGHTS_DST")" == "$RNAFM_WEIGHTS_SHA256" ]]; then
  log "weights already present and verified"
else
  rm -f "$WEIGHTS_DST"
  if [[ -n "${RNAFM_WEIGHTS:-}" ]]; then
    [[ "$(sha256 "$RNAFM_WEIGHTS")" == "$RNAFM_WEIGHTS_SHA256" ]] || die "RNAFM_WEIGHTS sha256 mismatch: $RNAFM_WEIGHTS"
    ln "$RNAFM_WEIGHTS" "$WEIGHTS_DST" 2>/dev/null || cp "$RNAFM_WEIGHTS" "$WEIGHTS_DST"
  else
    log "downloading weights from $RNAFM_WEIGHTS_URL"
    curl -L --fail -sS -o "$WEIGHTS_DST.part" "$RNAFM_WEIGHTS_URL"
    mv "$WEIGHTS_DST.part" "$WEIGHTS_DST"
  fi
  [[ "$(sha256 "$WEIGHTS_DST")" == "$RNAFM_WEIGHTS_SHA256" ]] || { rm -f "$WEIGHTS_DST"; die "weights sha256 mismatch"; }
  log "weights verified ($RNAFM_WEIGHTS_SHA256)"
fi
# fm.pretrained.rna_fm_t12() is called without a path, so RNA-FM loads its checkpoint through
# torch.hub's cache: $TORCH_HOME/hub/checkpoints/RNA-FM_pretrained.pth (else it downloads it).
mkdir -p "$RNAFM_DIR/torch_home/hub/checkpoints"
ln -sfn "../../../redevelop/pretrained/RNA-FM_pretrained.pth" "$RNAFM_DIR/torch_home/hub/checkpoints/RNA-FM_pretrained.pth"

# ---------------------------------------------------------------- 4. off-target (PITA)
PITA_DIR="$TARGET/off-target/pita"
PROGS_REL="Bin/ViennaRNA/ViennaRNA-1.6/Progs"
build_pita_bins() {  # $1 = output dir; builds from the pinned bundled ViennaRNA-1.6 source
  local out="$1" build
  build="$(mktemp -d "${TMPDIR:-/tmp}/pita-vienna.XXXXXX")"
  git -C "$TARGET" archive "$OLIGOFORMER_COMMIT" off-target/pita/Bin/ViennaRNA/ViennaRNA-1.6 | tar -x -C "$build"
  build="$build/off-target/pita/Bin/ViennaRNA/ViennaRNA-1.6"
  ( cd "$build"
    find . \( -name '*.o' -o -name '*.a' \) -delete
    for f in Progs/* lib/*; do if file "$f" | grep -q ELF; then rm -f "$f"; fi; done
    chmod +x configure config.guess config.sub install-sh 2>/dev/null || true
    ./configure CFLAGS="-O2 -std=gnu89 -Wno-implicit-function-declaration -Wno-implicit-int -Wno-return-type" \
      --without-perl >/dev/null
    make -C lib >/dev/null
    make -C Progs RNAduplex >/dev/null
    cc -O2 -std=gnu89 -Wno-implicit-function-declaration -Wno-implicit-int -Wno-return-type \
      -I. -I.. -IH -o Progs/RNAddG4 Progs/RNAddG4.c lib/libRNA.a -lm 2>/dev/null )
  mkdir -p "$out"
  cp "$build/Progs/RNAduplex" "$build/Progs/RNAddG4" "$out/"
  log "built PITA binaries into $out"
}

if [[ "${SKIP_OFFTARGET:-0}" != 1 ]]; then
  for t in perl RNAplfold sed grep; do [[ -x "$OFFTARGET_ENV/bin/$t" ]] || die "missing $OFFTARGET_ENV/bin/$t (run deps/create_offtarget_env.sh)"; done
  "$OFFTARGET_ENV/bin/perl" -e 'use Statistics::Lite; use Bio::TreeIO;' || die "perl modules missing in $OFFTARGET_ENV"
  "$OFFTARGET_ENV/bin/sed" --version | grep -q 'GNU sed' || die "GNU sed required"
  "$OFFTARGET_ENV/bin/grep" --version | grep -q 'GNU grep' || die "GNU grep required"

  # Prepare the host binaries BEFORE touching the checkout, so a build/verification failure
  # cannot leave the PITA tree reset (unusable Linux binaries, unsubstituted EXE_BASE_DIR).
  if [[ "$(uname -s)" == Darwin ]]; then
    if [[ "${PITA_FORCE_BUILD:-0}" == 1 || ! -f "$PITA_BIN_DIR/RNAduplex" || ! -f "$PITA_BIN_DIR/RNAddG4" ]]; then
      build_pita_bins "$PITA_BIN_DIR"
    elif [[ "$(uname -m)" == arm64 ]]; then
      # Reference hashes identify the binaries published with this package. A locally built
      # RNAduplex/RNAddG4 legitimately differs (ad-hoc signature, LC_UUID), so a mismatch means
      # "not ours", not "corrupt": rebuild from the pinned source instead of failing.
      if [[ "$(sha256 "$PITA_BIN_DIR/RNAduplex")" != "$RNADUPLEX_ARM64_SHA256" \
         || "$(sha256 "$PITA_BIN_DIR/RNAddG4")" != "$RNADDG4_ARM64_SHA256" ]]; then
        log "PITA binaries in $PITA_BIN_DIR are not the published ones; rebuilding from source"
        build_pita_bins "$PITA_BIN_DIR"
      fi
    fi
    for b in RNAduplex RNAddG4; do
      [[ -s "$PITA_BIN_DIR/$b" ]] || die "missing $PITA_BIN_DIR/$b after staging/build"
    done
  fi

  # reset PITA tree to the pinned state (undo previous install/staging/run debris)
  git -C "$TARGET" checkout --quiet -- off-target/pita
  git -C "$TARGET" clean -fdq -- off-target/pita

  if [[ "$(uname -s)" == Darwin ]]; then
    cp "$PITA_BIN_DIR/RNAduplex" "$PITA_BIN_DIR/RNAddG4" "$PITA_DIR/$PROGS_REL/"
    chmod 755 "$PITA_DIR/$PROGS_REL/RNAduplex" "$PITA_DIR/$PROGS_REL/RNAddG4"
    log "staged macOS PITA binaries"
  else
    log "non-macOS host: keeping shipped Linux x86-64 PITA binaries"
  fi
  ( cd "$PITA_DIR" && PATH="$OFFTARGET_ENV/bin:$PATH" make install >/dev/null )
  ! grep -rq EXE_BASE_DIR "$PITA_DIR/pita_prediction.pl" "$PITA_DIR/lib/RNAddG_compute.pl" || die "make install did not substitute EXE_BASE_DIR"
  log "PITA make install done (EXE_BASE_DIR -> $PITA_DIR)"
fi

# ---------------------------------------------------------------- 5. env file
ENV_FILE="$TARGET.env.sh"
cat > "$ENV_FILE" <<EOF
# Source before running OligoFormer from the runtime copy:
#   source "$ENV_FILE" && cd "\$OLIGOFORMER_RUNTIME" && python scripts/main.py ...
export OLIGOFORMER_RUNTIME="$TARGET"
# PROJECT_ENV first: scripts/RNA-FM.sh calls bare 'python'. OFFTARGET_ENV next: perl (BioPerl,
# Statistics::Lite), RNAplfold, GNU grep (-P, required by TargetScan) and GNU sed.
export PATH="$PROJECT_ENV/bin:$OFFTARGET_ENV/bin:\$PATH"
# Required: RNA-FM loads its checkpoint via torch.hub's cache under TORCH_HOME (no download).
export TORCH_HOME="$RNAFM_DIR/torch_home"
# Not needed: the conda perl finds BioPerl/Statistics::Lite in its own @INC; an inherited
# PERL5LIB could shadow them with incompatible modules, so it is cleared.
unset PERL5LIB
EOF
log "wrote $ENV_FILE"
log "done: $TARGET"
