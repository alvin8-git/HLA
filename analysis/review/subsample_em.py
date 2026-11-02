"""Half-sampling (without replacement) uncertainty for the dedup 1e-6 EM model.

Each replicate takes half of each group's deduplicated BMDP donors without replacement and
reruns EM. That avoids the duplicate-individual artefact of the with-replacement bootstrap,
which inflates coverage (dedup_em_bootstrap.csv). For m = n/2 without replacement, the
half-sample variance approximately equals the full-sample variance (finite-population factor
(n/m)(1 - m/n) = 1). So the interval is the full-sample point estimate plus the centred
half-sample quantiles. The median shift (half minus full) is also reported: it shows how far
the estimator drifts with sample size.
Statistics are evaluated at the full-sample N (n, 2n, 3n), the same as the point estimate.
Usage: python subsample_em.py <repo> <outdir> <B> <workers>
"""
import sys
import numpy as np, pandas as pd
from multiprocessing import Pool
REPO, OUT, B, W = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
sys.path.insert(0, f'{REPO}/analysis')
from hwe_test import run_em_haplotypes
COLS = ['A1', 'A2', 'B1', 'B2', 'C1', 'C2', 'DRB11', 'DRB12', 'DQB11', 'DQB12']


def coverage(f, N):
    f = np.sort(f)[::-1]; tot = 0.0
    for i in range(len(f)):
        d = 2.0 * f[i] * f[i:]; d[0] = f[i] ** 2
        tot += np.sum(d * -np.expm1(N * np.log1p(-d)))
    return tot


def nstar(f, C, lo=1e2, hi=1e13):
    if coverage(f, hi) < C: return np.inf
    for _ in range(30):
        mid = np.sqrt(lo * hi)
        lo, hi = (mid, hi) if coverage(f, mid) < C else (lo, mid)
    return hi


raw = pd.read_excel(f'{REPO}/HLA Data combined.xlsx', sheet_name='BMDP', dtype=str).fillna('')
cln = pd.read_excel(f'{REPO}/HLA Data.cleaned.xlsx', sheet_name='BMDP.out', dtype=str).fillna('')
cln.columns = cln.columns.str.strip()
key = raw[COLS].agg('|'.join, axis=1).values + '#' + cln.Ethnicity.values
last, near = {}, np.zeros(len(cln), bool)
for i, k in enumerate(key):
    j = last.get(k); near[i] = j is not None and i - j <= 10; last[k] = i
df = pd.read_csv(f'{REPO}/analysis/data/hla_clean.csv')
D = df[(df.source == 'BMDP_OUT') & ~df.sample_id.isin({f'BMDP_OUT_{i}' for i in np.flatnonzero(near)})]
pt = pd.read_csv(f'{OUT}/dedup_em_results.csv').query("scenario=='near'").set_index('ethnicity')


def rep(b):
    rng = np.random.default_rng(1000 + b)
    keep = [rng.choice(g.sample_id.unique(), len(g.sample_id.unique()) // 2, replace=False)
            for _, g in D.groupby('ethnicity')]
    h = run_em_haplotypes(D[D.sample_id.isin(np.concatenate(keep))], freq_threshold=1e-6)
    out = []
    for eth, g in h.groupby('ethnicity'):
        f = g.frequency.to_numpy(); n = pt.loc[eth, 'n']
        out.append(dict(b=b, ethnicity=eth, cov_n=coverage(f, n), cov_2n=coverage(f, 2 * n),
                        cov_3n=coverage(f, 3 * n), N50=nstar(f, 0.5), N75=nstar(f, 0.75)))
    print('rep', b, flush=True)
    return out


if __name__ == '__main__':
    with Pool(W) as p:
        res = pd.DataFrame([r for x in p.map(rep, range(B)) for r in x])
    res.to_csv(f'{OUT}/subsample_em.csv', index=False)
    rows = []
    for eth, g in res.groupby('ethnicity'):
        for s in ['cov_n', 'cov_2n', 'cov_3n', 'N50', 'N75']:
            x = g[s].to_numpy(); med = np.median(x); p0 = pt.loc[eth, s]
            lo, hi = p0 + np.quantile(x, [0.025, 0.975]) - med
            rows.append(dict(ethnicity=eth, stat=s, point=p0, lo=lo, hi=hi, half_median=med,
                             drift_pct=100 * (med / p0 - 1)))
    t = pd.DataFrame(rows); t.to_csv(f'{OUT}/subsample_em_ci.csv', index=False)
    print(t.to_string(index=False, float_format=lambda v: f'{v:,.4g}'))
