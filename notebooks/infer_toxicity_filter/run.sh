#!/usr/bin/env bash
# Driver for infer_toxicity_filter: README `oligoformer -tox` example, run natively in a private runtime copy.
# Usage: run.sh <label> <input_fasta (abs)> [extra args...]   -> outputs/<label>/, logs/<label>.{out,err,meta}
set -euo pipefail
N="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$N/runtime/OligoFormer.env.sh"
LABEL="$1"; INFA="$2"; shift 2
OUT="$N/outputs/$LABEL"
rm -rf "$OUT"; mkdir -p "$OUT/result" "$N/logs/$(dirname "$LABEL")"
cd "$OLIGOFORMER_RUNTIME"
CMD=(python scripts/main.py --infer 1 -i1 "$INFA" -tox --output_dir "$OUT/" "$@")
META="$N/logs/$LABEL.meta"
{ echo "cwd=$PWD"; echo "python=$(command -v python)"; python --version; printf 'cmd='; printf '%q ' "${CMD[@]}"; echo; echo "start=$(date -u +%FT%TZ)"; } > "$META"
t0=$(python -c 'import time;print(time.time())')
set +e
"${CMD[@]}" > "$N/logs/$LABEL.out" 2> "$N/logs/$LABEL.err"
rc=$?
set -e
t1=$(python -c 'import time;print(time.time())')
{ echo "exit=$rc"; python -c "print('wall_s=%.1f'%($t1-$t0))"; } >> "$META"
# keep relevant intermediate files (skip large .npy embeddings)
for rec in $(grep '^>' "$INFA" | sed 's/^>//; s/ /_@_/g'); do
  mkdir -p "$OUT/data_infer/$rec"; cp "data/infer/$rec/siRNA.fa" "data/infer/$rec/mRNA.fa" "$OUT/data_infer/$rec/" 2>/dev/null || true
done
cat "$META"
exit $rc
