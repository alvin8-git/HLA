"""For BMDP 2-field 10/10 twin pairs, are they also identical at the raw
high-resolution typing strings (G-groups, 3-field, suffix pattern)?
Usage: python dup_resolution.py <HLA Data combined.xlsx> <HLA Data.cleaned.xlsx>
"""
import sys
import numpy as np, pandas as pd
raw = pd.read_excel(sys.argv[1], sheet_name='BMDP', dtype=str).fillna('')
cln = pd.read_excel(sys.argv[2], sheet_name='BMDP.out', dtype=str).fillna('')
assert len(raw) == len(cln)
cols = ['A1', 'A2', 'B1', 'B2', 'C1', 'C2', 'DRB11', 'DRB12', 'DQB11', 'DQB12']
def key2(r):
    return '|'.join('+'.join(sorted((r[c1], r[c2]))) for c1, c2 in zip(cols[::2], cols[1::2]))
cln['k2'] = cln.apply(key2, axis=1)
raw['kraw'] = raw[cols].agg('|'.join, axis=1)
cln['kraw'] = raw['kraw'].values
cln['eth'] = cln.Ethnicity
cln['idx'] = np.arange(len(cln))
rows = []
for k, grp in cln.groupby('k2'):
    if len(grp) < 2: continue
    ix = grp.idx.to_numpy(); kr = grp.kraw.to_numpy(); e = grp.eth.iloc[0]
    for a in range(len(ix) - 1):
        rows.append(dict(eth=e, gap=ix[a + 1] - ix[a], raw_identical=kr[a] == kr[a + 1]))
p = pd.DataFrame(rows)
p['near'] = np.where(p.gap <= 10, 'gap<=10', 'gap>10')
print(p.groupby(['eth', 'near']).raw_identical.agg(['size', 'mean']).round(3).to_string())
print('\nall:', p.groupby('near').raw_identical.agg(['size', 'mean']).round(3).to_string())
print('\nfully identical raw rows (keep=False):', raw.duplicated(subset=cols, keep=False).sum(),
      '| removable duplicates (keep=first):', raw.duplicated(subset=cols + ['Ethnicity']).sum(), 'of', len(raw))
