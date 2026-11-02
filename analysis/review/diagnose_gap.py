"""Why does the HWE-haplotype model under-predict genotype repetition?
  (1) per-locus: observed genotype twin rate vs HWE prediction from allele freqs
      (no phasing involved -> isolates HWE departure)
  (2) per-source: twin rate within BMDP, within SCBB, and BMDP<->SCBB cross
      (cross-source twins cannot be within-registry duplicate records)
Usage: python diagnose_gap.py <hla_clean.csv>
"""
import sys
import numpy as np, pandas as pd

d = pd.read_csv(sys.argv[1], dtype=str)
d = d[d.source.isin(['BMDP_OUT', 'SCBB_OUT'])]
ok = d.allele1.str.fullmatch(r'\d+:\d+[NLQS]?') & d.allele2.str.fullmatch(r'\d+:\d+[NLQS]?')
d = d[~d.sample_id.isin(set(d.loc[~ok, 'sample_id']))]
d['pair'] = [a + '+' + b if a <= b else b + '+' + a for a, b in zip(d.allele1, d.allele2)]
LOCI = ['HLA-A', 'HLA-B', 'HLA-C', 'DRB1', 'DQB1']


def hwe_cov(alleles, N):
    p = pd.Series(alleles).value_counts(normalize=True).to_numpy()
    i, j = np.triu_indices(len(p))
    g = np.where(i == j, p[i] ** 2, 2 * p[i] * p[j])
    return np.sum(g * -np.expm1(N * np.log1p(-g)))


print("(1) per-locus genotype twin rate (LOO) vs HWE")
rows = []
for eth in ['Chinese', 'Malay', 'Indian', 'Others']:
    for loc in LOCI:
        s = d[(d.ethnicity == eth) & (d.locus == loc)]
        c = s.pair.value_counts(); n = c.sum()
        obs = 1 - (c == 1).sum() / n
        homo_obs = (s.allele1 == s.allele2).mean()
        p = pd.concat([s.allele1, s.allele2]).value_counts(normalize=True)
        homo_exp = (p ** 2).sum()
        rows.append(dict(eth=eth, locus=loc, obs_twin=round(obs, 4), hwe_twin=round(hwe_cov(pd.concat([s.allele1, s.allele2]), n - 1), 4),
                         homo_obs=round(homo_obs, 4), homo_exp=round(homo_exp, 4), F=round((homo_obs - homo_exp) / (1 - homo_exp), 4)))
print(pd.DataFrame(rows).to_string(index=False))

print("\n(2) 10/10 twins by source")
w = d.pivot_table(index=['sample_id', 'ethnicity', 'source'], columns='locus', values='pair', aggfunc='first').dropna().reset_index()
w['g'] = w[LOCI].agg('|'.join, axis=1)
for eth in ['Chinese', 'Malay', 'Indian', 'Others']:
    e = w[w.ethnicity == eth]
    b, s = e[e.source == 'BMDP_OUT'], e[e.source == 'SCBB_OUT']
    cb, cs = b.g.value_counts(), s.g.value_counts()
    within_b = (b.g.map(cb) > 1).mean()
    within_s = (s.g.map(cs) > 1).mean() if len(s) else np.nan
    # cross: fraction of SCBB units with a 10/10 twin among BMDP donors, vs expected if same population
    cross_s = s.g.isin(cb.index).mean() if len(s) else np.nan
    # expectation: rarefied within-BMDP coverage at N=len(b) is ~ within_b ; compare
    print(f"{eth:8s} BMDP n={len(b):6d} within={within_b:.3f} | SCBB n={len(s):5d} within={within_s:.3f} | SCBB unit has BMDP twin={cross_s:.3f}")
