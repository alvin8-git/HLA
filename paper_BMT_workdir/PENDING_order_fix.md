# Pending fix: order-sensitive cross-ethnic merge (found 2026-08-20)

**Not yet applied** — held until the freq_threshold=1e-4 run completes, so the two
changes stay separately attributable.

## Defect

`registry_model.get_diplotype_frequencies` labels each diplotype
`(haplotype1, haplotype2)` in the order the two haplotypes appear in *that
population's* frequency ranking. `04_registry_model.compute_coverage_cross` then
merges patient and donor frames on `['haplotype1','haplotype2']`.

Because the patient population and the combined donor pool rank haplotypes
differently, the same unordered pair {X, Y} is stored as (X,Y) in one frame and
(Y,X) in the other, and the ordered merge misses it — silently assigning
`donor_freq = 0` via the `fillna(0.0)`.

## Measured impact (Malay patients vs combined pool, 1e-3 data)

| Quantity | Value |
|---|---|
| Patient diplotype pairs | 8,911 |
| Unmatched by ordered merge | 5,528 (62.0%) |
| Unmatched after canonicalising pair order | 4,060 (45.6%) — these are genuine absences |
| **Patient frequency mass wrongly given donor_freq = 0** | **30.4%** |
| Cross-ethnic coverage at N=1e6 | 0.6353 (as-implemented) → 0.7757 (order-corrected) |

Direction: the defect **understates** cross-ethnic coverage, making cross-ethnic
matching look less feasible than the model actually implies. This affects the
manuscript's Table 3 and Recommendation 1 ("cross-ethnic matching is infeasible for
Malay, Indian and Others at any registry size").

Note the qualitative conclusion may well survive — 45.6% of pairs are genuinely
absent from the donor pool even after correction — but the reported N* values and
the ">10,000,000" entries must be recomputed before being restated.

## Fix

In `get_diplotype_frequencies`, emit the pair in canonical (lexicographic) order
instead of frequency-rank order. `diplotype_freq` is unchanged; only the two label
columns change, so every merge aligns:

```python
a = hap_arr[ii[order]]
b = hap_arr[jj[order]]
swap = a > b
h1 = np.where(swap, b, a)
h2 = np.where(swap, a, b)
```

Then re-run `04_registry_model.py` and compare the cross-ethnic rows against the
uncorrected run.
