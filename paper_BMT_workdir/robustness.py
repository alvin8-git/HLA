#!/usr/bin/env python
"""Are the paper's COMPARATIVE conclusions floor-robust?

Absolute N* is not (it moves 2,000-fold between 1e-3 and 1e-6). But the paper's
clinical claims are ratios, not absolutes:
  (a) cross-ethnic donors are harder to find than same-ethnicity ones
  (b) relaxing the match helps more than expanding the registry
A ratio can survive a bias that cancels in numerator and denominator. This
checks whether these two do, by recomputing each ratio under both floors.
"""
import io
import subprocess

import numpy as np
import pandas as pd

V215 = '16bd246'   # last commit before the 1e-6 re-run


def load(ref=None):
    if ref is None:
        return pd.read_csv('analysis/data/registry_size_targets.csv')
    txt = subprocess.run(['git', 'show', f'{ref}:analysis/data/registry_size_targets.csv'],
                         capture_output=True, text=True).stdout
    return pd.read_csv(io.StringIO(txt))


def penalty(d, level='10of10', tgt=0.95):
    """Cross-ethnic N* divided by same-ethnicity N*, per group."""
    s = d[(d.match_level == level) & (d.target_coverage == tgt)]
    p = s.pivot_table(index='ethnicity', columns='model_variant', values='registry_size')
    p = p[p.index != 'Combined']
    p['ratio'] = p.cross_ethnic / p.same_ethnicity
    return p


print('=' * 74)
print('(a) CROSS-ETHNIC PENALTY: N*(cross) / N*(same), 10/10 @ 95%')
print('=' * 74)
a, b = penalty(load(V215)), penalty(load())
cmp = pd.DataFrame({
    'N*same 1e-3': a.same_ethnicity, 'N*cross 1e-3': a.cross_ethnic, 'ratio 1e-3': a.ratio,
    'N*same 1e-6': b.same_ethnicity, 'N*cross 1e-6': b.cross_ethnic, 'ratio 1e-6': b.ratio,
})
cmp['ratio drift'] = (cmp['ratio 1e-6'] / cmp['ratio 1e-3']).round(2)
print(cmp.to_string(float_format=lambda x: f'{x:,.2f}'))
print('\n  direction identical in every group? ',
      bool(((cmp['ratio 1e-3'] > 1) == (cmp['ratio 1e-6'] > 1)).all()))

print('\n' + '=' * 74)
print('(b) MMUD: is one mismatch worth more than registry expansion?')
print('=' * 74)
mc = pd.read_csv('paper_BMT_workdir/partial_match_mc_all.csv')
mc = mc[mc.framework == '10-locus']
cov = mc[mc.coverage == 'cov@50k'].set_index(['ethnicity', 'min_match']).N_star

rows = []
for eth in cov.index.get_level_values(0).unique():
    c10, c9, c8 = cov[eth, '10/10'], cov[eth, '≥9/10'], cov[eth, '≥8/10']
    # How large must a 10/10-only registry be to reach what >=9/10 gives at 50k?
    n = mc[(mc.ethnicity == eth) & (mc.min_match == '10/10') &
           (mc.coverage != 'cov@50k')]
    x = n.coverage.astype(float).values * 100
    y = np.log10(n.N_star.values.astype(float))
    need9 = 10 ** np.interp(c9, x, y) if c9 <= x.max() else np.nan
    rows.append(dict(ethnicity=eth, cov10=c10, cov9=c9, cov8=c8,
                     gain_1mm=round(c9 - c10, 1), gain_2mm=round(c8 - c10, 1),
                     donors_for_same_gain=need9, expansion_x=need9 / 50_000))
r = pd.DataFrame(rows)
print(r.to_string(index=False, float_format=lambda x: f'{x:,.1f}'))
print('\n  cov* = % of patients matched by a 50,000-donor registry.')
print('  donors_for_same_gain = 10/10-only registry size reaching the >=9/10 figure')
print('  (log-interpolated on that group\'s own 10/10 curve; extrapolation beyond')
print('   the 75-95% anchors is not attempted).')
