#!/usr/bin/env python
"""Model-free empirical coverage: no EM, no HWE, no frequency floor.

Adjudicates the v2.15-vs-v2.16 dispute. Both report C(N) from an EM haplotype
table expanded under HWE; they disagree because the floor renormalises away
half the mass. This script bypasses the model entirely: it counts how often
two real donors in the registry are an exact 10/10 (or 8/8) genotype match.

  C(N) = sum_g p_g * [1 - (1 - p_g)^N]      p_g = observed genotype frequency

Same estimator, but p_g is a raw count, not an HWE-expanded EM product. This
is the ground truth the model must reproduce at N ~ n before its extrapolation
to N >> n means anything.

Known bias, stated rather than hidden: observed p_g is quantised at 1/n, so a
genotype seen once is scored 1/n when its true frequency may be far lower.
That makes empirical C(N) OPTIMISTIC -- an upper bound on real coverage. Any
model predicting more than this at N <= n is definitively wrong.
"""
import numpy as np
import pandas as pd

LOCI_10 = ['HLA-A', 'HLA-B', 'HLA-C', 'DRB1', 'DQB1']
LOCI_8 = ['HLA-A', 'HLA-B', 'HLA-C', 'DRB1']
MAIN_SOURCES = ['BMDP_OUT', 'SCBB_OUT']


def genotype_keys(df, loci):
    """One row per sample, value = 'A*01:01+A*02:01|B*...' across `loci`.
    Alleles sorted within locus so phase/order never splits a true match."""
    d = df[df.locus.isin(loci)].copy()
    lo = np.minimum(d.allele1.values.astype(str), d.allele2.values.astype(str))
    hi = np.maximum(d.allele1.values.astype(str), d.allele2.values.astype(str))
    d['g'] = d.locus.astype(str) + ':' + lo + '+' + hi
    w = d.pivot_table(index='sample_id', columns='locus', values='g',
                      aggfunc='first')
    w = w.dropna()                       # 5-locus (or 4-locus) complete only
    return w[loci].agg('|'.join, axis=1)


def coverage(p, N):
    """C(N) = sum p*[1-(1-p)^N], computed in log space for tiny p * huge N."""
    return float(np.sum(p * -np.expm1(N * np.log1p(-p))))


def nstar(p, target, hi=10**10):
    lo = 1
    if coverage(p, hi) < target:
        return np.inf
    for _ in range(60):
        mid = np.sqrt(lo * hi)
        if coverage(p, mid) < target:
            lo = mid
        else:
            hi = mid
    return hi


def main():
    df = pd.read_csv('analysis/data/hla_clean.csv')
    df = df[df.source.isin(MAIN_SOURCES)]

    NS = [10_000, 50_000, 100_000, 500_000, 1_000_000, 10_000_000]
    rows = []
    for level, loci in [('10of10', LOCI_10), ('8of8', LOCI_8)]:
        for eth, sub in df.groupby('ethnicity', observed=True):
            keys = genotype_keys(sub, loci)
            n = len(keys)
            counts = keys.value_counts().values.astype(float)
            p = counts / n
            singletons = int((counts == 1).sum())
            rows.append(dict(
                match_level=level, ethnicity=eth, n_donors=n,
                n_distinct_genotypes=len(counts),
                pct_singleton_genotypes=100 * singletons / len(counts),
                pct_donors_with_a_twin=100 * float(counts[counts > 1].sum()) / n,
                **{f'cov_{N}': round(100 * coverage(p, N), 1) for N in NS},
                N90=nstar(p, 0.90), N95=nstar(p, 0.95),
            ))

    out = pd.DataFrame(rows)
    out.to_csv('paper_BMT_workdir/empirical_match.csv', index=False)
    pd.set_option('display.width', 200, 'display.max_columns', 40)
    print(out.to_string(index=False))
    return out


if __name__ == '__main__':
    o = main()
    # ponytail: one runnable check -- coverage must be monotone in N and the
    # estimator must reproduce a hand-computed two-genotype case exactly.
    p = np.array([0.5, 0.5])
    assert abs(coverage(p, 1) - 0.5) < 1e-12
    assert abs(coverage(p, 2) - 0.75) < 1e-12
    c = o[[c for c in o.columns if c.startswith('cov_')]].values
    assert np.all(np.diff(c, axis=1) >= -1e-9), 'coverage not monotone in N'
    print('\nself-check OK')
