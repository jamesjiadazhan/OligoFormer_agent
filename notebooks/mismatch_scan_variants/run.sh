#!/usr/bin/env bash
# Driver for assignment mismatch_scan_variants (README Docker usage `oligoformer -m 2`).
# Usage: run.sh [cases...]   cases: run1 run2 run3 changed_input depth3 err_T err_18nt replay
# Each case: printf '<base>\n' | python scripts/main.py -m <M> --output_dir <case_dir>
# mismatch.py writes Args.output_dir + '/result/mismatch.txt', so <case_dir>/result is pre-created.
set -euo pipefail

NB="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$NB/runtime/OligoFormer.env.sh"
cd "$OLIGOFORMER_RUNTIME"
mkdir -p "$NB/logs" "$NB/outputs"
STATUS_FILE="$NB/logs/status.tsv"
[ -f "$STATUS_FILE" ] || printf 'case\tmode\tinput\texit\twall_s\tcommand\n' > "$STATUS_FILE"

DEPTH3_TIMEOUT="${DEPTH3_TIMEOUT:-900}"   # seconds

run_case() {  # name mode base_input
  local name="$1" mode="$2" base="$3"
  local out="$NB/outputs/$name"
  rm -rf "$out"; mkdir -p "$out/result"
  local cmd="printf '$base\\n' | python scripts/main.py -m $mode --output_dir $out"
  echo "[run.sh] $name: $cmd (cwd=$PWD)"
  local t0 t1 rc
  t0=$(python -c 'import time;print(time.time())')
  set +e
  printf '%s\n' "$base" | python scripts/main.py -m "$mode" --output_dir "$out" \
    > "$NB/logs/$name.stdout" 2> "$NB/logs/$name.stderr"
  rc=$?
  set -e
  t1=$(python -c 'import time;print(time.time())')
  local wall; wall=$(python -c "print(round($t1-$t0,2))")
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$name" "$mode" "$base" "$rc" "$wall" "$cmd" >> "$STATUS_FILE"
  echo "[run.sh] $name exit=$rc wall=${wall}s"
  return 0
}

run_timeout_case() {  # name mode base timeout_s  (kills after timeout; exit 124 recorded)
  local name="$1" mode="$2" base="$3" tmo="$4"
  local out="$NB/outputs/$name"
  rm -rf "$out"; mkdir -p "$out/result"
  local cmd="printf '$base\\n' | python scripts/main.py -m $mode --output_dir $out  [timeout ${tmo}s]"
  echo "[run.sh] $name: $cmd"
  local t0 t1 rc
  t0=$(python -c 'import time;print(time.time())')
  set +e
  printf '%s\n' "$base" | python scripts/main.py -m "$mode" --output_dir "$out" \
    > "$NB/logs/$name.stdout" 2> "$NB/logs/$name.stderr" &
  local pid=$!
  local waited=0
  while kill -0 "$pid" 2>/dev/null && [ "$waited" -lt "$tmo" ]; do sleep 5; waited=$((waited+5)); done
  if kill -0 "$pid" 2>/dev/null; then
    kill "$pid"; wait "$pid"; rc=124
  else
    wait "$pid"; rc=$?
  fi
  set -e
  t1=$(python -c 'import time;print(time.time())')
  local wall; wall=$(python -c "print(round($t1-$t0,2))")
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$name" "$mode" "$base" "$rc" "$wall" "$cmd" >> "$STATUS_FILE"
  echo "[run.sh] $name exit=$rc wall=${wall}s"
}

BASE="$(head -1 "$NB/data/base_input.txt")"
CHANGED="$(head -1 "$NB/data/changed_input.txt")"

for c in "${@:-run1 run2 run3 changed_input err_T err_18nt}"; do
  for case in $c; do
    case "$case" in
      run1|run2|run3) run_case "$case" 2 "$BASE" ;;
      replay)         run_case replay 2 "$(head -1 "$NB/data/base_input.txt")" ;;
      changed_input)  run_case changed_input 2 "$CHANGED" ;;
      err_T)          run_case err_T 2 "TCGCAAATTAAATTGAACC" ;;
      err_18nt)       run_case err_18nt 2 "UCGCAAAUUAAAUUGAAC" ;;
      depth3)         run_timeout_case depth3 3 "$BASE" "$DEPTH3_TIMEOUT" ;;
      *) echo "unknown case $case" >&2; exit 2 ;;
    esac
  done
done
