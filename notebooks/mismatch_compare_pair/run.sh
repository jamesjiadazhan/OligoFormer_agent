#!/usr/bin/env bash
# Driver for execution assignment mismatch_compare_pair:
#   printf '<base>\n<mismatch>\n' | python scripts/main.py -m 1     (cwd = runtime copy)
# Usage: run.sh [out_root] [data_dir]
#   out_root defaults to <this dir>/outputs, data_dir to <this dir>/data.
# Writes, per case: <out>/stdout.txt, <out>/stderr.txt, <out>/meta.tsv (exit, wall seconds,
# files created/modified in the runtime copy during the run).
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT_ROOT="${1:-$HERE/outputs}"
DATA="${2:-$HERE/data}"
# shellcheck disable=SC1091
source "$HERE/runtime/OligoFormer.env.sh"
cd "$OLIGOFORMER_RUNTIME"

PY="$(command -v python)"
echo "[run] python=$PY ($("$PY" -V 2>&1)) cwd=$PWD head=$(git rev-parse HEAD)"

# run_case <out_dir> <stdin_file> <expect: ok|fail>
run_case() {
  local out="$1" stdin="$2" expect="$3"
  mkdir -p "$out"
  local marker; marker="$(mktemp "${TMPDIR:-/tmp}/mm_marker.XXXXXX")"
  sleep 1
  local t0 t1 status
  t0=$(python -c 'import time;print(time.time())')
  set +e
  python scripts/main.py -m 1 < "$stdin" > "$out/stdout.txt" 2> "$out/stderr.txt"
  status=$?
  set -e
  t1=$(python -c 'import time;print(time.time())')
  # anything new/modified in the runtime copy (excluding .git and RNA-FM)
  local changed
  changed=$(find . -path ./.git -prune -o -path ./RNA-FM -prune -o -newer "$marker" -print | grep -v '^\.$' | tr '\n' ' ' || true)
  rm -f "$marker"
  printf 'stdin\t%s\nexit\t%s\nwall_s\t%s\nchanged_paths\t%s\ngit_status_lines\t%s\n' \
    "$stdin" "$status" "$(python -c "print(round($t1-$t0,3))")" "$changed" \
    "$(git status --short | wc -l | tr -d ' ')" > "$out/meta.tsv"
  echo "[run] $(basename "$out"): exit=$status expect=$expect stdout_last=$(tail -n1 "$out/stdout.txt")"
  if [[ "$expect" == ok && "$status" -ne 0 ]]; then echo "[run] UNEXPECTED FAILURE in $out" >&2; return 1; fi
  if [[ "$expect" == fail && "$status" -eq 0 ]]; then echo "[run] UNEXPECTED SUCCESS in $out" >&2; return 1; fi
}

for i in 1 2 3; do
  run_case "$OUT_ROOT/run$i" "$DATA/reference.stdin.txt" ok
done
for c in a_identical b_2mm c_5mm d_other_base_1mm; do
  run_case "$OUT_ROOT/changed_input/$c" "$DATA/$c.stdin.txt" ok
done
for c in e_len_19_vs_18 f1_lowercase f2_T_containing; do
  run_case "$OUT_ROOT/error_cases/$c" "$DATA/$c.stdin.txt" fail
done
echo "[run] all cases finished as expected"
