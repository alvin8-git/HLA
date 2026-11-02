"""Floor-free registry coverage from unphased genotypes.

A 10/10 (8/8) match is identity of the 5-locus (4-locus) genotype, so each
distinct genotype is a 'species' and match probability at registry size N is
sample coverage C(N) = sum_g p_g [1-(1-p_g)^N]. No EM, no frequency floor.

  interpolation  N <= n : Chao et al. 2014 rarefied coverage (exact, unbiased)
  extrapolation  N >  n : (a) Chao & Jost 2012 coverage extrapolation
                          (b) Pitman-Yor posterior predictive (MLE theta, sigma)
  validation           : fit on a subsample of size n/k, predict C at full n,
                          compare with leave-one-out observed 1 - f1/n.
Usage: python genotype_coverage.py <hla_clean.csv> <outdir>
"""
import sys, json
import numpy as np, pandas as pd
from scipy.special import gammaln
from scipy.optimize import minimize, brentq

SRC, OUT = sys.argv[1], sys.argv[2]
LOCI10 = ['HLA-A', 'HLA-B', 'HLA-C', 'DRB1', 'DQB1']
LOCI8 = LOCI10[:4]
TARGETS = [0.75, 0.85, 0.90, 0.95]
ETH = ['Chinese', 'Malay', 'Indian', 'Others']
rng = np.random.default_rng(42)

d = pd.read_csv(SRC, dtype=str)
d = d[d.source.isin(['BMDP_OUT', 'SCBB_OUT'])]
ok = d.allele1.str.fullmatch(r'\d+:\d+[NLQS]?') & d.allele2.str.fullmatch(r'\d+:\d+[NLQS]?')
bad_ids = set(d.loc[~ok, 'sample_id'])           # drop donors with ambiguity codes
d = d[~d.sample_id.isin(bad_ids)]
d['pair'] = [a + '+' + b if a <= b else b + '+' + a for a, b in zip(d.allele1, d.allele2)]
w = d.pivot_table(index=['sample_id', 'ethnicity'], columns='locus', values='pair', aggfunc='first').dropna().reset_index()
w['g10'] = w[LOCI10].agg('|'.join, axis=1)
w['g8'] = w[LOCI8].agg('|'.join, axis=1)


def lchoose(a, b):
    return gammaln(a + 1) - gammaln(b + 1) - gammaln(a - b + 1)


def rarefied_cov(counts, m):
    """Chao et al. 2014 eq: expected coverage of a size-m subsample (m < n)."""
    n = counts.sum()
    ks, fk = np.unique(counts, return_counts=True)
    ks_ok = ks[n - ks >= m]; fk_ok = fk[n - ks >= m]
    t = np.exp(lchoose(n - ks_ok, m) - lchoose(n - 1, m))
    return 1 - np.sum(fk_ok * ks_ok / n * t)


def cj_params(counts):
    n = counts.sum(); f1 = (counts == 1).sum(); f2 = (counts == 2).sum()
    A = (n - 1) * f1 / ((n - 1) * f1 + 2 * f2) if f2 > 0 else (n - 1) * (f1 - 1) / ((n - 1) * (f1 - 1) + 2)
    return n, f1, A


def cj_cov(counts, N):
    n, f1, A = cj_params(counts)
    return 1 - f1 / n * A ** (N - n + 1)


def cj_nstar(counts, C):
    n, f1, A = cj_params(counts)
    m1 = np.log((1 - C) * n / f1) / np.log(A)
    return n + m1 - 1


def py_fit(counts):
    n = counts.sum(); K = len(counts)
    def nll(x):
        s = 1 / (1 + np.exp(-x[0])); th = np.exp(x[1]) - s
        if th <= -s: return 1e18
        ll = np.sum(np.log(th + np.arange(1, K) * s)) - (gammaln(th + n) - gammaln(th + 1)) \
             + np.sum(gammaln(counts - s) - gammaln(1 - s))
        return -ll
    best = min((minimize(nll, [a, b], method='Nelder-Mead', options={'xatol': 1e-8, 'fatol': 1e-8, 'maxiter': 4000})
                for a in (-1, 0, 1) for b in (2, 6, 9)), key=lambda r: r.fun)
    s = 1 / (1 + np.exp(-best.x[0])); th = np.exp(best.x[1]) - s
    return th, s


def py_cov(counts, th, s, N):
    """Posterior predictive P(next draw is an already-seen type) after N total draws."""
    n = counts.sum(); K = len(counts); m = N - n
    if m <= 0:
        return np.nan
    # E[new types in m further draws | data] (Favaro et al. 2009)
    ratio = np.exp(gammaln(th + n + s + m) - gammaln(th + n + s) - gammaln(th + n + m) + gammaln(th + n))
    Knew = (K + th / s) * (ratio - 1)
    return 1 - (th + s * (K + Knew)) / (th + N)


def py_nstar(counts, th, s, C, hi=1e13):
    n = counts.sum()
    f = lambda logN: py_cov(counts, th, s, np.exp(logN)) - C
    lo = np.log(n + 1)
    if f(lo) >= 0: return float(n)
    if f(np.log(hi)) < 0: return np.inf
    return float(np.exp(brentq(f, lo, np.log(hi))))


rows, curves, valid = [], [], []
for level, key in (('10of10', 'g10'), ('8of8', 'g8')):
    for eth in ETH:
        g = w.loc[w.ethnicity == eth, key]
        counts = g.value_counts().to_numpy()
        n = counts.sum(); K = len(counts); f1 = int((counts == 1).sum()); f2 = int((counts == 2).sum())
        loo = 1 - f1 / n
        th, s = py_fit(counts)
        r = dict(match_level=level, ethnicity=eth, n=int(n), K=K, f1=f1, f2=f2,
                 observed_loo_cov=loo, py_theta=th, py_sigma=s)
        for C in TARGETS:
            cjn = cj_nstar(counts, C) if C > loo else np.nan
            pyn = py_nstar(counts, th, s, C) if C > loo else np.nan
            r[f'cj_N{int(C*100)}'] = cjn; r[f'py_N{int(C*100)}'] = pyn
        rows.append(r)
        # curves: rarefaction up to n, both extrapolations to 1e9
        for N in np.unique(np.r_[np.geomspace(100, n - 1, 40).astype(int), n - 1]):
            curves.append(dict(match_level=level, ethnicity=eth, N=int(N), method='rarefaction', coverage=rarefied_cov(counts, N)))
        for N in np.geomspace(n + 1, 1e9, 60):
            curves.append(dict(match_level=level, ethnicity=eth, N=float(N), method='chao_jost', coverage=cj_cov(counts, N)))
            curves.append(dict(match_level=level, ethnicity=eth, N=float(N), method='pitman_yor', coverage=py_cov(counts, th, s, N)))
        # validation: subsample at 1/2 and 1/4, extrapolate to full n
        if level == '10of10':
            vals = g.to_numpy()
            for frac in (0.5, 0.25):
                for rep in range(20):
                    sub = rng.choice(vals, size=int(n * frac), replace=False)
                    sc = pd.Series(sub).value_counts().to_numpy()
                    th_s, s_s = py_fit(sc)
                    valid.append(dict(ethnicity=eth, frac=frac, rep=rep, n_sub=int(sc.sum()), target_n=int(n),
                                      truth=loo, cj_pred=cj_cov(sc, n - 1), py_pred=py_cov(sc, th_s, s_s, n - 1)))
        print(f"{level} {eth}: n={n} K={K} f1={f1} f2={f2} LOO={loo:.4f} PY(theta={th:.1f}, sigma={s:.3f})", flush=True)

pd.DataFrame(rows).to_csv(f'{OUT}/genotype_coverage_targets.csv', index=False)
pd.DataFrame(curves).to_csv(f'{OUT}/genotype_coverage_curves.csv', index=False)
v = pd.DataFrame(valid); v.to_csv(f'{OUT}/genotype_coverage_validation.csv', index=False)
summ = v.groupby(['ethnicity', 'frac']).agg(truth=('truth', 'first'), cj_mean=('cj_pred', 'mean'), cj_sd=('cj_pred', 'std'),
                                            py_mean=('py_pred', 'mean'), py_sd=('py_pred', 'std')).reset_index()
summ['cj_err_pts'] = 100 * (summ.cj_mean - summ.truth); summ['py_err_pts'] = 100 * (summ.py_mean - summ.truth)
summ.to_csv(f'{OUT}/genotype_coverage_validation_summary.csv', index=False)
print('dropped donors with ambiguity codes:', len(bad_ids))
print(summ.round(4).to_string())
print(pd.DataFrame(rows).round(3).to_string())
