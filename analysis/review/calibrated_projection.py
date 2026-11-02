"""Bias-corrected 2x/3x coverage projections from sample_size_projection.csv.

The optimism of a k-fold extrapolation is measured by fitting on n/k donors and predicting
coverage at n, against the SCBB probe truth at n. Assuming the error depends on the
extrapolation factor rather than the absolute scale, corrected C(kn) = model C(kn) - err(k).
err(3) is interpolated in log k between err(2) and err(4). The interval combines the
Clopper-Pearson CI of the probe truth (which every err shares) with the rep spread of err(k).
Usage: python calibrated_projection.py <outdir>
"""
import sys
import numpy as np, pandas as pd
from scipy.stats import beta
OUT = sys.argv[1]
r = pd.read_csv(f'{OUT}/sample_size_projection.csv')
pt = pd.read_csv(f'{OUT}/dedup_em_results.csv').query("scenario=='near'").set_index('ethnicity')
r['k'] = 1 / r.frac
r['err'] = r.pred_at_n - r.truth_at_n
rows = []
for eth, g in r.groupby('ethnicity'):
    full = g[g.frac == 1].iloc[0]
    # probe count at full n, recovered from the rate and the SCBB group size
    m_probe = {'Chinese': 1028, 'Malay': 94, 'Indian': 93, 'Others': 235}[eth]
    hit = round(full.truth_at_n * m_probe)
    lo, hi = beta.ppf(0.025, hit, m_probe - hit + 1), beta.ppf(0.975, hit + 1, m_probe - hit)
    e = g.groupby('k').err.median(); sd = g.groupby('k').err.std().fillna(0)
    e3 = np.interp(np.log(3), np.log([2, 4]), [e[2], e[4]]); sd3 = np.interp(np.log(3), np.log([2, 4]), [sd[2], sd[4]])
    for k, ek, sk, model in [(1, e[1], 0, full.pred_at_n), (2, e[2], sd[2], pt.loc[eth, 'cov_2n']),
                             (3, e3, sd3, pt.loc[eth, 'cov_3n'])]:
        c = model - ek
        # probe CI shifts the corrected value one-for-one: err = pred - truth
        clo = c - (full.truth_at_n - lo) - 1.96 * sk; chi = c + (hi - full.truth_at_n) + 1.96 * sk
        rows.append(dict(ethnicity=eth, size=f'{k}n', N=int(full.n * k), model=model, err_pts=100 * ek,
                         corrected=c, lo=max(clo, 0), hi=chi))
t = pd.DataFrame(rows); t.to_csv(f'{OUT}/calibrated_projection.csv', index=False)
pd.set_option('display.width', 200)
print(t.to_string(index=False, float_format=lambda v: f'{v:.3f}'))
# observed (model-free) coverage gain per doubling of BMDP, from probes against nested subsamples
o = r.groupby(['ethnicity', 'frac']).truth_at_m.median().unstack()
print('\nmodel-free probe coverage at n/8, n/4, n/2, n:\n', (100 * o).round(1).to_string())
