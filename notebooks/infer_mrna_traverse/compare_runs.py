"""Compare OligoFormer infer_mrna_traverse outputs across runs (no recomputation of the model)."""
import json, sys, os
import pandas as pd
from scipy.stats import spearmanr
NB = os.path.dirname(os.path.abspath(__file__))
O = os.path.join(NB, 'outputs')
runs = ['run1', 'run2', 'run3', 'replay', 'readme_literal']
det = ['pos', 'sense', 'siRNA', 'func_filter', 'filter']
res = {'per_record': {}}
for rec in ['RNA0', 'RNA1']:
    t = {r: pd.read_csv(f'{O}/{r}/{rec}.txt', sep='\t') for r in runs}
    rk = {r: pd.read_csv(f'{O}/{r}/{rec}_ranked.txt', sep='\t') for r in runs}
    rf = {r: pd.read_csv(f'{O}/{r}/{rec}_ranked_filtered.txt', sep='\t') for r in runs}
    base = t['run1']
    d = {'rows': int(base.shape[0]), 'columns': list(base.columns),
         'ranked_filtered_rows': {r: int(rf[r].shape[0]) for r in runs},
         'func_filter_counts': {str(k): int(v) for k, v in base['func_filter'].value_counts().sort_index().items()},
         'deterministic_columns_identical': {r: bool(t[r][det].equals(base[det])) for r in runs},
         'filtered_pos_set_identical': {r: sorted(rf[r]['pos']) == sorted(rf['run1']['pos']) for r in runs},
         'efficacy_range_run1': [float(base.efficacy.min()), float(base.efficacy.max())],
         'pairwise': {}}
    for i, a in enumerate(runs):
        for b in runs[i+1:]:
            ea, eb = t[a].efficacy, t[b].efficacy
            d['pairwise'][f'{a}~{b}'] = {
                'max_abs_diff_efficacy': float((ea - eb).abs().max()),
                'max_rel_diff_efficacy': float(((ea - eb).abs() / eb.abs()).max()),
                'spearman_efficacy': float(spearmanr(ea, eb).correlation),
                'ranked_order_identical': rk[a]['pos'].tolist() == rk[b]['pos'].tolist(),
                'ranked_filtered_order_identical': rf[a]['pos'].tolist() == rf[b]['pos'].tolist(),
                'top5_ranked_pos': [rk[a]['pos'].tolist()[:5], rk[b]['pos'].tolist()[:5]]}
    # window FASTAs identical across runs
    d['window_fastas_identical'] = all(
        open(f'{O}/{r}/data_infer/{rec}/{f}').read() == open(f'{O}/run1/data_infer/{rec}/{f}').read()
        for r in runs for f in ['siRNA.fa', 'mRNA.fa'])
    d['n_windows_in_siRNA_fa'] = open(f'{O}/run1/data_infer/{rec}/siRNA.fa').read().count('>')
    res['per_record'][rec] = d
c = pd.read_csv(f'{O}/changed_input/HuMRNA3.txt', sep='\t')
cf = pd.read_csv(f'{O}/changed_input/HuMRNA3_ranked_filtered.txt', sep='\t')
res['changed_input'] = {'rows': int(c.shape[0]), 'columns': list(c.columns),
    'ranked_filtered_rows': int(cf.shape[0]),
    'func_filter_counts': {str(k): int(v) for k, v in c['func_filter'].value_counts().sort_index().items()},
    'efficacy_range': [float(c.efficacy.min()), float(c.efficacy.max())],
    'top3_ranked': pd.read_csv(f'{O}/changed_input/HuMRNA3_ranked.txt', sep='\t').head(3)[['pos','siRNA','efficacy']].to_dict('records')}
json.dump(res, open(os.path.join(NB, 'logs', 'comparison.json'), 'w'), indent=1)
print(json.dumps(res, indent=1))
