#!/usr/bin/env bash
# Stage 2 driver for assignment infer_offtarget_toxicity (README "3.2 Inference with off-target").
# Runs the real upstream CLI (scripts/main.py) from this assignment's private runtime copy.
# Usage: run.sh [all|runs|changed|all_human|replay]   (default: all)
set -euo pipefail

N="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$N/runtime/OligoFormer.env.sh"
RT="$OLIGOFORMER_RUNTIME"
LOGS="$N/logs"; OUT="$N/outputs"; DATA="$N/data"
mkdir -p "$LOGS" "$OUT" "$DATA"
CMDLOG="$LOGS/commands.tsv"
[ -f "$CMDLOG" ] || printf 'label\tstart_utc\texit_status\twall_seconds\tcwd\tcommand\n' > "$CMDLOG"
RECORDS=(RNA0 RNA1)   # record names in data/example.fa

# run_main <label> <expect_fail:0|1> <main.py args...>
# Runs main.py with --output_dir <OUT>/<label>/, records exit status + wall time, copies the
# per-record intermediate files (pita.tab, targetscan.tab, siRNA.fa) out of the fixed ./data/infer/.
run_main() {
  local label="$1" expect_fail="$2"; shift 2
  local od="$OUT/$label"
  rm -rf "$od"; mkdir -p "$od/result"      # result/ for the mismatch module's fixed sub-path
  local cmd=(python scripts/main.py "$@" --output_dir "$od/")
  local t0 t1 st start
  start="$(date -u +%FT%TZ)"; t0=$(python -c 'import time;print(time.time())')
  set +e
  ( cd "$RT" && "${cmd[@]}" ) > "$LOGS/$label.stdout.log" 2> "$LOGS/$label.stderr.log"
  st=$?
  set -e
  t1=$(python -c 'import time;print(time.time())')
  local wall; wall=$(python -c "print(round($t1-$t0,2))")
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$label" "$start" "$st" "$wall" "$RT" "${cmd[*]}" >> "$CMDLOG"
  echo "[$label] exit=$st wall=${wall}s"
  for r in "${RECORDS[@]}"; do
    if [ -d "$RT/data/infer/$r" ]; then
      mkdir -p "$od/infer/$r"
      for f in pita.tab targetscan.tab siRNA.fa top_n_siRNA.fa; do
        [ -f "$RT/data/infer/$r/$f" ] && cp -p "$RT/data/infer/$r/$f" "$od/infer/$r/" || true
      done
    fi
  done
  if [ "$expect_fail" = 0 ] && [ "$st" -ne 0 ]; then
    echo "[$label] FAILED (exit $st); see $LOGS/$label.stderr.log" >&2; exit "$st"
  fi
  return 0
}

save_inputs() {
  for f in data/example.fa off-target/ref/100_human_UTR.txt off-target/ref/100_human_ORF.txt toxicity/cell_viability.txt; do
    cp -p "$RT/$f" "$DATA/"
  done
  ( cd "$DATA" && shasum -a 256 example.fa 100_human_UTR.txt 100_human_ORF.txt cell_viability.txt > SHA256SUMS )
}

mode="${1:-all}"
case "$mode" in all|runs|changed|all_human|replay) ;; *) echo "unknown mode $mode" >&2; exit 2 ;; esac
want() { [ "$mode" = all ] || [ "$mode" = "$1" ]; }   # (bash 3.2 compatible; no ';;&')
if want runs; then
  save_inputs
  for i in 1 2 3; do run_main "run$i" 0 --infer 1 -i1 ./data/example.fa -off -tox; done
fi
if want changed; then
  run_main changed_params 0 --infer 1 -i1 ./data/example.fa -off -tox --pita_threshold -5 --targetscan_threshold 0.2
fi
if want all_human; then
  # Literal README command (with -a). Expected to fail in the pinned code; recorded as evidence.
  run_main readme_all_human 1 --infer 1 -i1 ./data/example.fa -off -tox -a
fi
if want replay; then
  # Replay from the saved inputs in notebooks/<ID>/data (cell_viability path is hard-coded upstream;
  # its sha256 is compared instead).
  ( cd "$DATA" && shasum -a 256 -c SHA256SUMS )
  run_main replay 0 --infer 1 -i1 "$DATA/example.fa" --utr "$DATA/100_human_UTR.txt" --orf "$DATA/100_human_ORF.txt" -off -tox
fi
