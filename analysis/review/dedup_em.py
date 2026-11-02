"""EM re-run on deduplicated BMDP donors (1e-6 floor, no pooling, unconditional coverage).

Fits haplotypes on BMDP only, so the SCBB cord units are a genuinely out-of-sample probe
(the live haplo_freqs_em.csv was fitted on BMDP+SCBB, so the earlier model-vs-probe check
leaked). Then reports coverage at n, 2n, 3n and N* for 75/90/95%.
Usage: python dedup_em.py <repo> <outdir> [B_bootstrap]
"""
import sys, time
import numpy as np, pandas as pd
REPO, OUT = sys.argv[1], sys.argv[2]
B = int(sys.argv[3]) if len(sys.argv) > 3 else 0
sys.path.insert(0, f'{REPO}/analysis')
from hwe_test import run_em_haplotypes

COLS = ['A1', 'A2', 'B1', 'B2', 'C1', 'C2', 'DRB11', 'DRB12', 'DQB11', 'DQB12']
ETH = {'C': 'Chinese', 'M': 'Malay', 'I': 'Indian', 'O': 'Others'}
probe = {'Chinese': (261, 1000), 'Malay': (7, 90), 'Indian': (4, 91), 'Others': (7, 220)}  # dedup_coverage.log


def coverage(f, N):
    f = np.sort(f)[::-1]; tot = 0.0
    for i in range(len(f)):
        d = 2.0 * f[i] * f[i:]; d[0] = f[i] ** 2
        tot += np.sum(d * -np.expm1(N * np.log1p(-d)))
    return tot


def nstar(f, C, lo=1e2, hi=1e13):
    if coverage(f, hi) < C: return np.inf
    for _ in range(40):  # bisection on log N
        mid = np.sqrt(lo * hi)
        lo, hi = (mid, hi) if coverage(f, mid) < C else (lo, mid)
    return hi


# near-duplicate flags, same rule as dedup_coverage.py
raw = pd.read_excel(f'{REPO}/HLA Data combined.xlsx', sheet_name='BMDP', dtype=str).fillna('')
cln = pd.read_excel(f'{REPO}/HLA Data.cleaned.xlsx', sheet_name='BMDP.out', dtype=str).fillna('')
cln.columns = cln.columns.str.strip()
key = raw[COLS].agg('|'.join, axis=1).values + '#' + cln.Ethnicity.values
last, near = {}, np.zeros(len(cln), bool)
for i, k in enumerate(key):
    j = last.get(k); near[i] = j is not None and i - j <= 10; last[k] = i
drop = {'near': {f'BMDP_OUT_{i}' for i in np.flatnonzero(near)},
        'all': {f'BMDP_OUT_{i}' for i in np.flatnonzero(pd.Series(key).duplicated())}}

df = pd.read_csv(f'{REPO}/analysis/data/hla_clean.csv')
bm = df[df.source == 'BMDP_OUT']
rows, rng = [], np.random.default_rng(11)
for scen, ids in drop.items():
    d = bm[~bm.sample_id.isin(ids)]
    t0 = time.time(); h = run_em_haplotypes(d, freq_threshold=1e-6); dt = time.time() - t0
    h.to_csv(f'{OUT}/dedup_em_haplo_{scen}.csv', index=False)
    for eth, g in h.groupby('ethnicity'):
        f = g.frequency.to_numpy()
        n = d[(d.ethnicity == eth) & d.allele1.notna()].groupby('sample_id').locus.nunique().eq(5).sum()
        r = dict(scenario=scen, ethnicity=eth, n=n, K=len(f), mass=f.sum(), em_s=round(dt),
                 probe=probe[eth][0] / probe[eth][1],
                 cov_n=coverage(f, n), cov_2n=coverage(f, 2 * n), cov_3n=coverage(f, 3 * n),
                 **{f'N{int(c*100)}': nstar(f, c) for c in (0.5, 0.75, 0.9, 0.95)})
        rows.append(r); print(r, flush=True)
res = pd.DataFrame(rows); res.to_csv(f'{OUT}/dedup_em_results.csv', index=False)
print(res.to_string())

# donor-level non-parametric bootstrap, EM re-run per replicate ('near' scenario)
if B:
    d = bm[~bm.sample_id.isin(drop['near'])]
    bs = []
    for b in range(B):
        parts = []
        for eth, g in d.groupby('ethnicity'):
            ids = g.sample_id.unique(); pick = rng.choice(ids, len(ids))
            m = pd.DataFrame({'sample_id': pick, 'new': [f'{eth}_{b}_{k}' for k in range(len(pick))]})
            parts.append(m.merge(g, on='sample_id').drop(columns='sample_id').rename(columns={'new': 'sample_id'}))
        h = run_em_haplotypes(pd.concat(parts), freq_threshold=1e-6)
        for eth, g in h.groupby('ethnicity'):
            f = g.frequency.to_numpy(); n = res.query("scenario=='near' and ethnicity==@eth").n.iat[0]
            bs.append(dict(b=b, ethnicity=eth, cov_n=coverage(f, n), cov_2n=coverage(f, 2 * n),
                           cov_3n=coverage(f, 3 * n), N50=nstar(f, 0.5), N75=nstar(f, 0.75)))
        pd.DataFrame(bs).to_csv(f'{OUT}/dedup_em_bootstrap.csv', index=False)
        print('boot', b, flush=True)
    q = pd.DataFrame(bs).groupby('ethnicity').quantile([0.025, 0.5, 0.975]).drop(columns='b')
    print(q.to_string())
