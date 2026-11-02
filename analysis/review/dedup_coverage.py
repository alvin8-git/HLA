"""Genotype-coverage estimators under three duplicate-handling scenarios,
validated against SCBB cord units as independent probes of the BMDP registry.

  none  : records as ingested
  near  : drop a record if its raw high-resolution typing string equals that of a
          record <=10 rows earlier (near-certain duplicate entry). Lower bound.
  all   : drop every raw-string repeat (keep first). Upper bound: also removes
          some genuine independent twins.
Truth for validation: fraction of SCBB units (same ethnicity) with a 2-field 10/10
twin in BMDP. Dedup cannot change it (it depends only on the set of BMDP genotypes),
so it scores each scenario's coverage estimate at N = that scenario's BMDP size.
Usage: python dedup_coverage.py <combined.xlsx> <cleaned.xlsx> <outdir>
"""
import sys
import numpy as np, pandas as pd
from scipy.special import gammaln
from scipy.optimize import minimize, brentq
from scipy.stats import beta

RAW, CLN, OUT = sys.argv[1:4]
COLS = ['A1', 'A2', 'B1', 'B2', 'C1', 'C2', 'DRB11', 'DRB12', 'DQB11', 'DQB12']
ETH = {'C': 'Chinese', 'M': 'Malay', 'I': 'Indian', 'O': 'Others'}
TARGETS = [0.75, 0.85, 0.90, 0.95]


def keys(df, cols):
    pairs = ['+'.join(sorted((a, b))) for c1, c2 in zip(cols[::2], cols[1::2]) for a, b in [(None, None)]] if False else None
    out = []
    for r in df[cols].itertuples(index=False):
        out.append('|'.join('+'.join(sorted((r[i], r[i + 1]))) for i in range(0, len(cols), 2)))
    return out


raw_b = pd.read_excel(RAW, sheet_name='BMDP', dtype=str).fillna('')
cln_b = pd.read_excel(CLN, sheet_name='BMDP.out', dtype=str).fillna('')
cln_s = pd.read_excel(CLN, sheet_name='SCBB.out', dtype=str).fillna('')
for df in (cln_b, cln_s):
    df.columns = df.columns.str.strip()
    ok = np.all([df[c].str.fullmatch(r'\d+:\d+[NLQS]?') for c in COLS], axis=0)
    df['ok'] = ok
    df['k10'] = keys(df, COLS)
    df['eth'] = df.Ethnicity.str.strip().map(ETH)
cln_b['raw'] = raw_b[COLS].agg('|'.join, axis=1).values + '#' + cln_b.Ethnicity
cln_b['idx'] = np.arange(len(cln_b))

# near-duplicate flag: same raw string seen within the previous 10 rows
last_seen = {}
near = np.zeros(len(cln_b), bool)
for i, k in enumerate(cln_b.raw):
    j = last_seen.get(k)
    near[i] = j is not None and i - j <= 10
    last_seen[k] = i
cln_b['dup_near'] = near
cln_b['dup_all'] = cln_b.raw.duplicated()
B = cln_b[cln_b.ok & cln_b.eth.notna()]
S = cln_s[cln_s.ok & cln_s.eth.notna()]


def lchoose(a, b): return gammaln(a + 1) - gammaln(b + 1) - gammaln(a - b + 1)


def cj(counts):
    n = counts.sum(); f1 = (counts == 1).sum(); f2 = (counts == 2).sum()
    A = (n - 1) * f1 / ((n - 1) * f1 + 2 * f2)
    return n, f1, A


def cj_cov(counts, N):
    n, f1, A = cj(counts); return 1 - f1 / n * A ** (N - n + 1)


def cj_nstar(counts, C):
    n, f1, A = cj(counts); return n + np.log((1 - C) * n / f1) / np.log(A) - 1


def py_fit(counts):
    n = counts.sum(); K = len(counts)
    def nll(x):
        s = 1 / (1 + np.exp(-x[0])); th = np.exp(x[1]) - s
        return -(np.sum(np.log(th + np.arange(1, K) * s)) - (gammaln(th + n) - gammaln(th + 1))
                 + np.sum(gammaln(counts - s) - gammaln(1 - s)))
    r = min((minimize(nll, [a, b], method='Nelder-Mead', options={'xatol': 1e-8, 'fatol': 1e-8, 'maxiter': 4000})
             for a in (-1, 0, 1) for b in (2, 6, 9)), key=lambda r: r.fun)
    s = 1 / (1 + np.exp(-r.x[0])); return np.exp(r.x[1]) - s, s


def py_cov(counts, th, s, N):
    n = counts.sum(); K = len(counts); m = N - n
    ratio = np.exp(gammaln(th + n + s + m) - gammaln(th + n + s) - gammaln(th + n + m) + gammaln(th + n))
    return 1 - (th + s * (K + (K + th / s) * (ratio - 1))) / (th + N)


def py_nstar(counts, th, s, C, hi=1e13):
    n = counts.sum(); f = lambda x: py_cov(counts, th, s, np.exp(x)) - C
    if f(np.log(n + 1)) >= 0: return float(n)
    if f(np.log(hi)) < 0: return np.inf
    return float(np.exp(brentq(f, np.log(n + 1), np.log(hi))))


def rarefied(counts, m):
    n = counts.sum(); ks, fk = np.unique(counts, return_counts=True)
    keep = n - ks >= m
    return 1 - np.sum(fk[keep] * ks[keep] / n * np.exp(lchoose(n - ks[keep], m) - lchoose(n - 1, m)))


rows, curves = [], []
rng = np.random.default_rng(7)
for scen, mask in (('none', np.ones(len(B), bool)), ('near', ~B.dup_near.to_numpy()), ('all', ~B.dup_all.to_numpy())):
    Bs = B[mask]
    for eth in ETH.values():
        g_b = Bs.loc[Bs.eth == eth, 'k10']
        cb = g_b.value_counts(); counts = cb.to_numpy(); n = int(counts.sum())
        probes = S.loc[S.eth == eth, 'k10']
        hit = int(probes.isin(cb.index).sum()); m = len(probes)
        lo, hi = beta.ppf(0.025, hit, m - hit + 1), beta.ppf(0.975, hit + 1, m - hit)
        th, s = py_fit(counts)
        loo = 1 - (counts == 1).sum() / n
        # half-sample -> full extrapolation test (20 reps), scored against the probe truth
        halves_cj, halves_py = [], []
        for _ in range(20):
            sub = pd.Series(rng.choice(g_b.to_numpy(), n // 2, replace=False)).value_counts().to_numpy()
            t2, s2 = py_fit(sub)
            halves_cj.append(cj_cov(sub, n)); halves_py.append(py_cov(sub, t2, s2, n))
        r = dict(scenario=scen, ethnicity=eth, n_bmdp=n, K=len(counts), f1=int((counts == 1).sum()), f2=int((counts == 2).sum()),
                 loo_cov=loo, probe_n=m, probe_cov=hit / m, probe_lo=lo, probe_hi=hi,
                 half_cj_pred=np.mean(halves_cj), half_py_pred=np.mean(halves_py), py_theta=th, py_sigma=s)
        for C in TARGETS:
            r[f'cj_N{int(C*100)}'] = cj_nstar(counts, C) if C > loo else np.nan
            r[f'py_N{int(C*100)}'] = py_nstar(counts, th, s, C) if C > loo else np.nan
        rows.append(r)
        for N in np.unique(np.r_[np.geomspace(200, n - 1, 30).astype(int), n - 1]):
            curves.append(dict(scenario=scen, ethnicity=eth, N=float(N), method='rarefaction', coverage=rarefied(counts, N)))
        for N in np.geomspace(n + 1, 1e9, 50):
            curves.append(dict(scenario=scen, ethnicity=eth, N=float(N), method='chao_jost', coverage=cj_cov(counts, N)))
            curves.append(dict(scenario=scen, ethnicity=eth, N=float(N), method='pitman_yor', coverage=py_cov(counts, th, s, N)))
        print(f"{scen:4s} {eth:8s} n={n:6d} LOO={loo:.3f} probe={hit}/{m}={hit/m:.3f} [{lo:.3f},{hi:.3f}] "
              f"halfCJ={np.mean(halves_cj):.3f} halfPY={np.mean(halves_py):.3f} PY(th={th:.0f},s={s:.2f})", flush=True)

out = pd.DataFrame(rows)
out.to_csv(f'{OUT}/dedup_coverage_targets.csv', index=False)
pd.DataFrame(curves).to_csv(f'{OUT}/dedup_coverage_curves.csv', index=False)
print('\nflagged near-duplicates:', int(B.dup_near.sum()), '| all raw repeats:', int(B.dup_all.sum()), '| BMDP usable:', len(B))
pd.set_option('display.width', 250)
print(out[['scenario', 'ethnicity', 'n_bmdp', 'loo_cov', 'probe_cov', 'cj_N75', 'cj_N90', 'cj_N95', 'py_N75', 'py_N90', 'py_N95']].round(3).to_string(index=False))
