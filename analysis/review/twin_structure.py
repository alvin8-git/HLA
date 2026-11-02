"""Are within-BMDP 10/10 twins independent draws or duplicate/related records?
Test: sample_id index gap between twins vs a permutation null (random pairs).
Usage: python twin_structure.py <hla_clean.csv>
"""
import sys
import numpy as np, pandas as pd
d = pd.read_csv(sys.argv[1], dtype=str)
d = d[d.source == 'BMDP_OUT']
d['pair'] = [a + '+' + b if a <= b else b + '+' + a for a, b in zip(d.allele1, d.allele2)]
w = d.pivot_table(index=['sample_id', 'ethnicity'], columns='locus', values='pair', aggfunc='first').dropna().reset_index()
w['idx'] = w.sample_id.str.rsplit('_', n=1).str[-1].astype(int)
w['g'] = w[['HLA-A', 'HLA-B', 'HLA-C', 'DRB1', 'DQB1']].agg('|'.join, axis=1)
rng = np.random.default_rng(0)
for eth in ['Chinese', 'Malay', 'Indian', 'Others']:
    e = w[w.ethnicity == eth].sort_values('idx')
    gaps = []
    for g, grp in e.groupby('g'):
        if len(grp) > 1:
            ix = grp.idx.to_numpy(); gaps += list(np.diff(ix))   # consecutive members of a twin cluster
    gaps = np.array(gaps)
    null = np.abs(np.diff(rng.choice(e.idx.to_numpy(), size=(len(gaps), 2)), axis=1)).ravel()
    sizes = e.g.value_counts(); sizes = sizes[sizes > 1]
    print(f"{eth:8s} twin-adjacent pairs={len(gaps):5d}  gap<=10: {np.mean(gaps<=10):.3f} (null {np.mean(null<=10):.4f})  "
          f"gap<=1000: {np.mean(gaps<=1000):.3f} (null {np.mean(null<=1000):.4f})  median gap {np.median(gaps):.0f} (null {np.median(null):.0f})  "
          f"cluster sizes: {sizes.value_counts().sort_index().head(6).to_dict()}")
