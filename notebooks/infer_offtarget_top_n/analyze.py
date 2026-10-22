"""Compare top_n_siRNA.fa against RESULT tables, filter codes, run-to-run spread. Reads saved outputs only."""
import json, sys, hashlib, pathlib
import pandas as pd, numpy as np
from scipy.stats import spearmanr
N = pathlib.Path(__file__).resolve().parent
runs = ['run1', 'run2', 'run3', 'replay']
out = {}
def fasta(p):
    L = p.read_text().split('\n'); L = [x for x in L if x]
    return [(L[i][1:], L[i+1]) for i in range(0, len(L), 2)]
tables = {}
for r in runs:
    d = N/'outputs'/r
    if not (d/'RNA0.txt').exists(): continue
    out[r] = {}
    for rec in ['RNA0', 'RNA1']:
        t = pd.read_csv(d/f'{rec}.txt', sep='\t'); tables[(r, rec)] = t
        rk = pd.read_csv(d/f'{rec}_ranked.txt', sep='\t')
        fa = fasta(d/'data_infer'/rec/'top_n_siRNA.fa')
        true_top = [int(p)-1 for p in rk['pos'][:10]]
        hdr_idx = [int(h[3:]) for h, _ in fa]
        seqs = [s for _, s in fa]
        seq_rows = [int(t.index[t['siRNA'] == s][0]) for s in seqs]
        out[r][rec] = dict(
            n_rows=len(t), top_n_fa_entries=len(fa),
            headers_equal_true_top10=hdr_idx == true_top,
            sequence_rows_written=seq_rows,
            sequences_equal_first10_positional=seq_rows == list(range(10)),
            n_header_seq_mismatch=sum(h != s for h, s in zip(hdr_idx, seq_rows)),
            true_top10_idx=true_top,
            off_target_filter_counts={str(k): int(v) for k, v in t['off_target_filter'].value_counts().items()},
            rows_with_pita_score=int(t['pita_score'].notna().sum()),
            rows_with_targetscan_score=int(t['targetscan_score'].notna().sum()),
            scored_rows_idx=sorted(int(i) for i in t.index[t['off_target_filter'] != -5]),
            filter_eq0_rows=int((t['filter'] == 0).sum()),
            ranked_filtered_rows=len(pd.read_csv(d/f'{rec}_ranked_filtered.txt', sep='\t')),
            columns=list(t.columns),
            sha256={f.name: hashlib.sha256(f.read_bytes()).hexdigest() for f in sorted(d.rglob('*')) if f.is_file()},
        )
comp = {}
base = [r for r in ['run1', 'run2', 'run3'] if r in out]
for rec in ['RNA0', 'RNA1']:
    c = {}
    det = ['pos', 'sense', 'siRNA', 'func_filter', 'Seed', 'cell_viability', 'toxicity_filter']
    for a in base + (['replay'] if 'replay' in out else []):
        for b in base:
            if a >= b and a != 'replay': continue
            if a == b: continue
            ta, tb = tables[(a, rec)], tables[(b, rec)]
            key = f'{a}_vs_{b}'
            sa = set(out[a][rec]['true_top10_idx']); sb = set(out[b][rec]['true_top10_idx'])
            c[key] = dict(
                deterministic_cols_identical=bool(ta[det].equals(tb[det])),
                efficacy_max_abs_diff=float(np.max(np.abs(ta['efficacy'] - tb['efficacy']))),
                efficacy_max_rel_diff=float(np.max(np.abs(ta['efficacy'] - tb['efficacy']) / np.abs(tb['efficacy']))),
                efficacy_spearman=float(spearmanr(ta['efficacy'], tb['efficacy']).correlation),
                top10_overlap=len(sa & sb), top10_order_identical=out[a][rec]['true_top10_idx'] == out[b][rec]['true_top10_idx'],
                top_n_fa_identical=(N/'outputs'/a/'data_infer'/rec/'top_n_siRNA.fa').read_bytes() == (N/'outputs'/b/'data_infer'/rec/'top_n_siRNA.fa').read_bytes(),
                pita_tab_identical=(N/'outputs'/a/'data_infer'/rec/'pita.tab').read_bytes() == (N/'outputs'/b/'data_infer'/rec/'pita.tab').read_bytes(),
                targetscan_tab_identical=(N/'outputs'/a/'data_infer'/rec/'targetscan.tab').read_bytes() == (N/'outputs'/b/'data_infer'/rec/'targetscan.tab').read_bytes(),
                pita_score_identical=bool(ta['pita_score'].equals(tb['pita_score'])),
                targetscan_score_identical=bool(ta['targetscan_score'].equals(tb['targetscan_score'])),
                off_target_filter_identical=bool(ta['off_target_filter'].equals(tb['off_target_filter'])),
                filter_identical=bool(ta['filter'].equals(tb['filter'])),
            )
    if len(base) == 3:
        sets = [set(out[r][rec]['true_top10_idx']) for r in base]
        c['top10_intersection_all3'] = len(sets[0] & sets[1] & sets[2])
        c['top10_union_all3'] = len(sets[0] | sets[1] | sets[2])
    comp[rec] = c
json.dump(dict(per_run=out, comparisons=comp), sys.stdout, indent=1)
