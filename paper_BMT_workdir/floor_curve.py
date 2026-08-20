"""Bias-vs-floor curve: is the rare tail real signal or EM phase-enumeration noise?

Answers Reviewer 1's challenge that 57,934 distinct haplotypes from <=10,000
haplotype copies must contain phase-ambiguity artefact. If the artefact dominates,
N* should be unstable across floors; if there is real structure between 1e-4 and
1e-3, N* should move sharply there and stabilise below it.
"""
import sys, os, pandas as pd
HERE = '/data/alvin/HLA/analysis'
sys.path.insert(0, HERE)
os.chdir(HERE)
from hwe_test import run_em_haplotypes
from registry_model import get_diplotype_frequencies, find_registry_size, compute_coverage

h = pd.read_csv('data/hla_clean.csv')
h = h[h.ethnicity.isin(['Chinese', 'Others'])]

em_full = run_em_haplotypes(h, freq_threshold=0.0, cap=5000)
rows = []
for floor in (0.0, 1e-5, 1e-4, 5e-4, 1e-3):
    for eth in ['Chinese', 'Others']:
        s = em_full[em_full.ethnicity == eth][['haplotype', 'frequency']].copy()
        s = s[s.frequency >= floor]
        mass = float(s.frequency.sum())
        s['frequency'] = s.frequency / mass
        dip = get_diplotype_frequencies(s)
        rec = dict(ethnicity=eth, floor=floor, n_haps=len(s),
                   mass_retained=round(mass, 4),
                   cov_50k=round(100 * compute_coverage(dip, 50_000), 1),
                   cov_500k=round(100 * compute_coverage(dip, 500_000), 1),
                   N90=find_registry_size(dip, 0.90))
        rows.append(rec)
        print(rec, flush=True)

df = pd.DataFrame(rows)
df.to_csv('/data/alvin/HLA/paper_BMT_workdir/floor_curve.csv', index=False)
print("DONE")
