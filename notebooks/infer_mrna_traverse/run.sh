#!/usr/bin/env bash
# Driver for execution assignment infer_mrna_traverse
# Source: repo/OligoFormer/README.md, Usage 3.1 Option 1 (traverse mRNA, 19-nt windows)
# README command: python scripts/main.py --infer 1 -i1 data/example.fa
# cwd for every command = the private runtime copy (runtime/OligoFormer).
set -euo pipefail

NB="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$NB/runtime/OligoFormer.env.sh"
RT="$OLIGOFORMER_RUNTIME"
OUT="$NB/outputs"; LOG="$NB/logs"
mkdir -p "$OUT" "$LOG" "$NB/data"
TSV="$LOG/commands.tsv"
[ -f "$TSV" ] || printf 'label\tcommand\texit\twall_s\n' > "$TSV"

# ---- saved inputs ----
cp -f "$RT/data/example.fa" "$NB/data/example.fa"
# changed input: data/fasta/Hu_mRNA.fa record RNA3 (full 57-nt sequence, <=120 nt), header renamed
if [ ! -f "$NB/data/changed_input_Hu_mRNA_RNA3.fa" ]; then
  python - "$RT/data/fasta/Hu_mRNA.fa" "$NB/data/changed_input_Hu_mRNA_RNA3.fa" <<'PY'
import sys
src,dst=sys.argv[1],sys.argv[2]
seqs={};name=None
for line in open(src):
    line=line.strip()
    if line.startswith('>'): name=line[1:]; seqs[name]=''
    elif name: seqs[name]+=line
s=seqs['RNA3'][:120]
open(dst,'w').write('>HuMRNA3\n%s\n'%s)
print('HuMRNA3 len', len(s))
PY
fi

run_cmd () { # label, then command words
  local label="$1"; shift
  local t0 t1 rc
  t0=$(date +%s)
  set +e
  ( cd "$RT" && "$@" ) >"$LOG/${label}.stdout.log" 2>"$LOG/${label}.stderr.log"
  rc=$?
  set -e
  t1=$(date +%s)
  printf '%s\t%s\t%s\t%s\n' "$label" "$*" "$rc" "$((t1-t0))" >> "$TSV"
  echo "[$label] exit=$rc wall=$((t1-t0))s"
  return $rc
}

collect () { # label dest_dir  -- copy result tables + window FASTAs
  local dest="$1"; shift
  mkdir -p "$dest/data_infer"
  for rec in "$@"; do
    for d in "$RT/data/infer/$rec"; do
      [ -d "$d" ] || continue
      mkdir -p "$dest/data_infer/$rec"
      for f in siRNA.fa mRNA.fa; do
        [ -f "$d/$f" ] && cp -f "$d/$f" "$dest/data_infer/$rec/$f"
      done
    done
  done
}

# ---- 1. literal README command (writes ./result/ inside the runtime copy) ----
rm -rf "$RT/result"
run_cmd readme python scripts/main.py --infer 1 -i1 data/example.fa
mkdir -p "$OUT/readme_literal"
cp -f "$RT"/result/*.txt "$OUT/readme_literal/"
collect "$OUT/readme_literal" RNA0 RNA1

# ---- 2. three recorded runs with --output_dir ----
for i in 1 2 3; do
  d="$OUT/run$i"; rm -rf "$d"; mkdir -p "$d"
  run_cmd "run$i" python scripts/main.py --infer 1 -i1 data/example.fa --output_dir "$d/"
  collect "$d" RNA0 RNA1
done

# ---- 3. replay from the saved input copy ----
d="$OUT/replay"; rm -rf "$d"; mkdir -p "$d"
run_cmd replay python scripts/main.py --infer 1 -i1 "$NB/data/example.fa" --output_dir "$d/"
collect "$d" RNA0 RNA1

# ---- 4. changed input ----
d="$OUT/changed_input"; rm -rf "$d"; mkdir -p "$d"
run_cmd changed_input python scripts/main.py --infer 1 -i1 "$NB/data/changed_input_Hu_mRNA_RNA3.fa" --output_dir "$d/"
collect "$d" HuMRNA3

# ---- hashes ----
( cd "$NB" && find data outputs -type f | sort | xargs shasum -a 256 ) > "$LOG/sha256.txt"
echo "driver done"
