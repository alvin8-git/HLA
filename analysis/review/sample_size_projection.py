"""How do the 1e-6 model's outputs depend on how many donors were typed?

For fractions 1/8, 1/4, 1/2 and 1 of the near-deduplicated BMDP donors (without replacement,
R reps per fraction), EM is refitted. For each fit it records:
  - N* at 50/75/90/95% (does the target settle down as the typed sample grows?)
  - predicted coverage at the subsample's own size m, and at the full BMDP size n
    (an extrapolation of n/m-fold)
  - the truth for both: the fraction of SCBB cord units with a 10/10 twin in the subsample,
    and in full BMDP
The second pair is an out-of-sample test of extrapolation: fitting on n/8 and predicting at n
is the same 8-fold step as fitting on n and predicting at 8n.
Usage: python sample_size_projection.py <repo> <outdir> <reps> <workers>
"""
import sys
import numpy as np, pandas as pd
from multiprocessing import Pool
REPO, OUT, R, W = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
sys.path.insert(0, f'{REPO}/analysis')
from hwe_test import run_em_haplotypes
COLS = ['A1', 'A2', 'B1', 'B2', 'C1', 'C2', 'DRB11', 'DRB12', 'DQB11', 'DQB12']
FRACS = [0.125, 0.25, 0.5, 1.0]
TARGETS = [0.5, 0.75, 0.9, 0.95]


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
df = df[df.source.isin(['BMDP_OUT', 'SCBB_OUT'])]
df = df[~df.sample_id.isin({f'BMDP_OUT_{i}' for i in np.flatnonzero(near)})]
# 10/10 genotype key per sample: unordered allele pair per locus, loci in fixed order
df['pair'] = [a + '+' + b if a <= b else b + '+' + a for a, b in zip(df.allele1, df.allele2)]
G = df.pivot_table(index=['sample_id', 'source', 'ethnicity'], columns='locus', values='pair', aggfunc='first').dropna()
G['g'] = G[['HLA-A', 'HLA-B', 'HLA-C', 'DRB1', 'DQB1']].agg('|'.join, axis=1)
G = G.reset_index()[['sample_id', 'source', 'ethnicity', 'g']]
BM, SC = G[G.source == 'BMDP_OUT'], G[G.source == 'SCBB_OUT']
D = df[df.source == 'BMDP_OUT']
NFULL = BM.groupby('ethnicity').size()


def probe_rate(eth, bm_ids):
    s = set(BM.loc[BM.sample_id.isin(bm_ids) & (BM.ethnicity == eth), 'g'])
    p = SC.loc[SC.ethnicity == eth, 'g']
    return p.isin(s).mean(), len(p)


def job(args):
    frac, r = args
    rng = np.random.default_rng(int(frac * 1000) * 100 + r)
    keep = np.concatenate([g.sample_id.to_numpy() if frac == 1 else
                           rng.choice(g.sample_id.to_numpy(), int(len(g) * frac), replace=False)
                           for _, g in BM.groupby('ethnicity')])
    h = run_em_haplotypes(D[D.sample_id.isin(keep)], freq_threshold=1e-6)
    out = []
    for eth, g in h.groupby('ethnicity'):
        f = g.frequency.to_numpy(); n = NFULL[eth]; m = int(BM[(BM.ethnicity == eth) & BM.sample_id.isin(keep)].shape[0])
        t_m, np_ = probe_rate(eth, keep)
        out.append(dict(frac=frac, rep=r, ethnicity=eth, m=m, n=n, K=len(f),
                        pred_at_m=coverage(f, m), truth_at_m=t_m,
                        pred_at_n=coverage(f, n), pred_at_2n=coverage(f, 2 * n), pred_at_3n=coverage(f, 3 * n),
                        **{f'N{int(c * 100)}': nstar(f, c) for c in TARGETS}))
    print('done', frac, r, flush=True)
    return out


if __name__ == '__main__':
    truth_n = {e: probe_rate(e, BM.sample_id) for e in NFULL.index}
    print('probe truth at full n:', {e: (round(t, 4), k) for e, (t, k) in truth_n.items()}, flush=True)
    jobs = [(f, r) for f in FRACS for r in (range(1) if f == 1 else range(R))]
    with Pool(W) as p:
        res = pd.DataFrame([x for o in p.map(job, jobs) for x in o])
    res['truth_at_n'] = res.ethnicity.map({e: t for e, (t, _) in truth_n.items()})
    res.to_csv(f'{OUT}/sample_size_projection.csv', index=False)
    s = res.groupby(['ethnicity', 'frac']).median(numeric_only=True).drop(columns='rep')
    s['extrap_err_pts'] = 100 * (s.pred_at_n - s.truth_at_n)
    s['own_err_pts'] = 100 * (s.pred_at_m - s.truth_at_m)
    pd.set_option('display.width', 250)
    print(s[['m', 'K', 'truth_at_m', 'pred_at_m', 'own_err_pts', 'truth_at_n', 'pred_at_n', 'extrap_err_pts',
             'N50', 'N75', 'N90', 'N95']].to_string(float_format=lambda v: f'{v:,.4g}'))
    # elasticity of N* in typed sample size between successive fractions
    for c in ['N50', 'N75', 'N95']:
        e = np.log(s[c]).groupby(level=0).diff() / np.log(2)
        print(c, 'elasticity d log N* / d log m:', e.dropna().round(2).to_dict())
