#!/usr/bin/env bash
# Driver for execution assignment infer_specific_sirnas
# Source: repo/OligoFormer/README.md, Usage 3.1 Option 2 (pinned e2f53ad63387bbe166bf123949151e2bc9bf6ec3)
#   python scripts/main.py --infer 1 -i1 data/example.fa -i2 data/example_siRNA.fa --output_dir <dir>/
# Runs: run1..run3 (documented command, runtime-copy data paths), replay (same command on saved
# notebooks/<ID>/data copies), changed_input, error_cases/{len20, zero_match_RNA1}.
set -euo pipefail

N="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$N/runtime/OligoFormer.env.sh"
R="$OLIGOFORMER_RUNTIME"
OUT="$N/outputs"; LOG="$N/logs"; DATA="$N/data"
mkdir -p "$OUT" "$LOG"
CMDLOG="$LOG/commands.tsv"
printf 'label\texit_status\twall_seconds\tcwd\tcommand\n' > "$CMDLOG"

now() { python -c 'import time; print(f"{time.time():.3f}")'; }

# run_case <label> <outdir> <expect: ok|fail> <args...>
run_case() {
  local label="$1" od="$2" expect="$3"; shift 3
  rm -rf "$od"; mkdir -p "$od"
  # main.py writes intermediates to a fixed ./data/infer/<record>/ path; start clean
  rm -rf "$R/data/infer"
  local cmd=(python scripts/main.py --infer 1 "$@" --output_dir "$od/")
  local t0 t1 rc
  t0=$(now)
  set +e
  ( cd "$R" && "${cmd[@]}" ) > "$LOG/$label.stdout" 2> "$LOG/$label.stderr"
  rc=$?
  set -e
  t1=$(now)
  local wall; wall=$(python -c "print(round($t1-$t0,3))")
  printf '%s\t%s\t%s\t%s\t%s\n' "$label" "$rc" "$wall" "$R" "${cmd[*]}" >> "$CMDLOG"
  echo "[$label] exit=$rc wall=${wall}s"
  # keep small intermediates (fasta inputs to RNA-FM + listing of embeddings), not the .npy files
  if [ -d "$R/data/infer" ]; then
    for d in "$R"/data/infer/*/; do
      [ -d "$d" ] || continue
      local rec; rec=$(basename "$d")
      mkdir -p "$od/infer_intermediate/$rec"
      for f in siRNA.fa mRNA.fa; do [ -f "$d/$f" ] && cp "$d/$f" "$od/infer_intermediate/$rec/"; done
      ( cd "$d" && find . -type f | sort ) > "$od/infer_intermediate/$rec/file_listing.txt"
    done
  fi
  if [ "$expect" = ok ] && [ "$rc" -ne 0 ]; then
    echo "FAILED: $label exited $rc (see $LOG/$label.stderr)" >&2; exit "$rc"
  fi
  if [ "$expect" = fail ] && [ "$rc" -eq 0 ]; then
    echo "UNEXPECTED SUCCESS: $label" >&2; exit 1
  fi
  return 0
}

for i in 1 2 3; do
  run_case "run$i" "$OUT/run$i" ok -i1 data/example.fa -i2 data/example_siRNA.fa
done
run_case replay "$OUT/replay" ok -i1 "$DATA/example.fa" -i2 "$DATA/example_siRNA.fa"
run_case changed_input "$OUT/changed_input" ok -i1 "$DATA/example.fa" -i2 "$DATA/changed_siRNA.fa"
run_case error_len20 "$OUT/error_cases/len20" fail -i1 "$DATA/example.fa" -i2 "$DATA/error_20nt_siRNA.fa"
# zero-match record: exit status recorded whatever it is (behaviour under investigation)
run_case error_zero_match_RNA1 "$OUT/error_cases/zero_match_RNA1" any -i1 "$DATA/example.fa" -i2 "$DATA/zero_match_RNA1_siRNA.fa"
rm -rf "$R/data/infer"
echo "all done"
