import pandas as pd, numpy as np, json, sys, hashlib, glob, os
from scipy.stats import spearmanr
N=sys.argv[1]; O=N+'/outputs'
tox=pd.read_csv(N+'/data/cell_viability.txt',sep='\t')
res={}
det=['pos','sense','siRNA','func_filter','Seed','cell_viability','toxicity_filter','filter']
for rec in ['RNA0','RNA1']:
    R={l:pd.read_csv(f'{O}/{l}/{rec}.txt',sep='\t') for l in ['run1','run2','run3','replay']}
    b=R['run1']; r={'rows':len(b),'columns':list(b.columns)}
    r['filter_eq_func_plus_tox']=bool((b['filter']==b['func_filter']+b['toxicity_filter']).all())
    r['seed_is_siRNA_2_7']=bool((b['Seed']==b['siRNA'].str.slice(1,7)).all())
    m=b.merge(tox,on='Seed',suffixes=('','_tab'))
    r['cell_viability_matches_table']=bool(len(m)==len(b) and (m['cell_viability']==m['cell_viability_tab']).all())
    r['toxicity_flagged']=int(b['toxicity_filter'].sum()); r['cell_viability_nan']=int(b['cell_viability'].isna().sum())
    r['func_flagged']=int((b['func_filter']>0).sum()); r['filter_pass']=int((b['filter']==0).sum())
    r['tox_flagged_rows']=b.loc[b['toxicity_filter']==1,['pos','siRNA','Seed','cell_viability','func_filter','filter']].to_dict('records')
    rf=pd.read_csv(f'{O}/run1/{rec}_ranked_filtered.txt',sep='\t'); r['ranked_filtered_rows']=len(rf)
    cmp={}
    for l in ['run2','run3','replay']:
        o=R[l]; c={'deterministic_cols_identical':bool(o[det].equals(b[det])),
          'efficacy_max_abs_diff':float(np.abs(o.efficacy-b.efficacy).max()),
          'efficacy_max_rel_diff':float((np.abs(o.efficacy-b.efficacy)/np.abs(b.efficacy)).max()),
          'spearman':float(spearmanr(o.efficacy,b.efficacy).correlation),
          'ranked_order_identical':bool(pd.read_csv(f'{O}/{l}/{rec}_ranked.txt',sep='\t').pos.tolist()==pd.read_csv(f'{O}/run1/{rec}_ranked.txt',sep='\t').pos.tolist()),
          'ranked_filtered_order_identical':bool(pd.read_csv(f'{O}/{l}/{rec}_ranked_filtered.txt',sep='\t').pos.tolist()==rf.pos.tolist())}
        cmp[l+'_vs_run1']=c
    r['comparisons']=cmp
    t=pd.read_csv(f'{O}/changed_params/thr80/{rec}.txt',sep='\t')
    r['thr80']={'columns':list(t.columns),'toxicity_flagged':int(t.toxicity_filter.sum()),
      'toxicity_filter_eq_cv_lt_80':bool((t.toxicity_filter==(t.cell_viability<80).astype(int)).all()),
      'filter_eq_func_plus_tox':bool((t['filter']==t.func_filter+t.toxicity_filter).all()),
      'filter_pass':int((t['filter']==0).sum()),'cols_other_than_tox_filter_efficacy_identical':bool(t[['pos','siRNA','func_filter','Seed','cell_viability']].equals(b[['pos','siRNA','func_filter','Seed','cell_viability']]))}
    f=pd.read_csv(f'{O}/changed_params/nf/{rec}.txt',sep='\t')
    r['nf']={'columns':list(f.columns),'has_func_filter':'func_filter' in f.columns,'toxicity_flagged':int(f.toxicity_filter.sum()),
      'filter_eq_tox':bool((f['filter']==f.toxicity_filter).all()),'filter_pass':int((f['filter']==0).sum()),
      'tox_cols_identical_to_run1':bool(f[['pos','siRNA','Seed','cell_viability','toxicity_filter']].equals(b[['pos','siRNA','Seed','cell_viability','toxicity_filter']])),
      'ranked_filtered_rows':len(pd.read_csv(f'{O}/changed_params/nf/{rec}_ranked_filtered.txt',sep='\t'))}
    res[rec]=r
h={os.path.relpath(p,N):hashlib.sha256(open(p,'rb').read()).hexdigest() for p in sorted(glob.glob(O+'/**/*.*',recursive=True)) if os.path.isfile(p)}
json.dump({'per_record':res,'output_sha256':h},open(N+'/outputs/comparison.json','w'),indent=1)
print(json.dumps(res,indent=1,default=str)[:6000])
