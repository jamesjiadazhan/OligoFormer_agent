#!/usr/bin/env bash
# Stage-2 driver for infer_offtarget_top_n (README 3.2 second command, -top N variant).
# Usage: run.sh <label> <top_n> <input_fasta_abs> [expect_fail]
set -euo pipefail
N="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LABEL="$1"; TOPN="$2"; INPUT="$3"; EXPECT_FAIL="${4:-0}"
source "$N/runtime/OligoFormer.env.sh"
cd "$OLIGOFORMER_RUNTIME"
OUT="$N/outputs/$LABEL"
rm -rf "$OUT"; mkdir -p "$OUT/result"
LOG="$N/logs/$LABEL"
CMD=(python scripts/main.py --infer 1 -i1 "$INPUT" -off -tox -top "$TOPN" --output_dir "$OUT/")
echo "cwd=$PWD" > "$LOG.cmd"; printf '%q ' "${CMD[@]}" >> "$LOG.cmd"; echo >> "$LOG.cmd"
# record pre-existing top_n_siRNA.fa (accumulation check)
for r in RNA0 RNA1; do f="data/infer/$r/top_n_siRNA.fa"; if [ -f "$f" ]; then echo "pre_$r=$(wc -l < "$f")"; else echo "pre_$r=absent"; fi; done >> "$LOG.cmd"
START=$(python -c 'import time;print(time.time())')
set +e
"${CMD[@]}" > "$LOG.stdout" 2> "$LOG.stderr"
RC=$?
set -e
END=$(python -c 'import time;print(time.time())')
echo "exit_status=$RC" >> "$LOG.cmd"
python -c "print('wall_seconds=%.1f' % ($END-$START))" >> "$LOG.cmd"
# keep relevant native intermediates (no .npy)
for r in RNA0 RNA1; do
  if [ -d "data/infer/$r" ]; then
    mkdir -p "$OUT/data_infer/$r"
    for f in top_n_siRNA.fa siRNA.fa mRNA.fa pita.tab targetscan.tab; do
      [ -f "data/infer/$r/$f" ] && cp -p "data/infer/$r/$f" "$OUT/data_infer/$r/"
    done
  fi
done
cat "$LOG.cmd"
if [ "$EXPECT_FAIL" = "1" ]; then
  [ "$RC" -ne 0 ] && { echo "expected failure observed"; exit 0; }
  echo "expected failure NOT observed"; exit 1
fi
exit "$RC"
