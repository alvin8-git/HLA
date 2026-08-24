#!/usr/bin/env python
"""Reconcile v2.15 / v2.16 / observed reality at ONE point: N ~ registry size.

The Chinese registry holds 44,400 five-locus-typed donors, and the v2.15 model
puts N*(95%) at 41,727 -- close enough to compare a model prediction against a
directly observable fact. That fact is the leave-one-out twin rate: the share
of donors whose exact 10/10 genotype occurs at least twice. It assumes nothing
-- no EM, no HWE, no floor -- and is unbiased at N = n-1.

v2.15 renormalised after applying the floor. Multiplying its coverage back by
the retained mass undoes that and asks whether the underlying arithmetic was
sound or whether the whole model was.
"""
import pandas as pd

fc = pd.read_csv('paper_BMT_workdir/floor_curve_full.csv')
em = pd.read_csv('paper_BMT_workdir/empirical_match.csv')
obs = em[em.match_level == '10of10'].set_index('ethnicity')

print('Typing resolution (drives whether an exact-genotype match is credible):')
d = pd.read_csv('analysis/data/hla_clean.csv', usecols=['locus', 'allele1'])
d['fields'] = d.allele1.astype(str).str.count(':') + 1
print(d.groupby('locus', observed=True).fields.value_counts().unstack().fillna(0)
        .astype(int).to_string(), '\n')

rows = []
for _, r in fc.iterrows():
    e = r.ethnicity
    if e not in obs.index:
        continue
    o = obs.loc[e]
    rows.append(dict(
        ethnicity=e, floor=r.floor, mass=round(r.mass_retained, 4),
        # as published: coverage of the RETAINED subpopulation only
        cov50k_published=r.cov_50k,
        # denominator restored: coverage of ALL patients
        cov50k_unconditional=round(r.cov_50k * r.mass_retained, 1),
        observed_twin_rate=round(o.pct_donors_with_a_twin, 1),
        n_donors=int(o.n_donors),
        N95_published=int(r.N95),
    ))
t = pd.DataFrame(rows)
t.to_csv('paper_BMT_workdir/reconcile.csv', index=False)
pd.set_option('display.width', 200)
print(t.to_string(index=False))

print('\nChinese, the one group where n (44,400) ~ v2.15 N*95 (41,727):')
c = t[t.ethnicity == 'Chinese']
for _, r in c.iterrows():
    print(f'  floor {r.floor:<8} published {r.cov50k_published:5.1f}%  '
          f'x mass {r.mass:6.4f} -> {r.cov50k_unconditional:5.1f}%  '
          f'| observed {r.observed_twin_rate}%')
