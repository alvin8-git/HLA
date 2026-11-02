# Registry-size review scripts

These scripts back `docs/review-registry-size-estimation.md`. They are standalone and are not
part of `run_all.sh`. They read the raw Excel sources and `analysis/data/hla_clean.csv`, and write
only to `analysis/review/out/`. Run them from the repo root with the pinned interpreter.

```bash
PY=/home/alvin/miniconda3/bin/python; O=analysis/review/out; D=analysis/data
# 1. duplicates and the SCBB probe benchmark
$PY analysis/review/dup_resolution.py "HLA Data combined.xlsx" "HLA Data.cleaned.xlsx" > $O/dup_resolution.log
$PY analysis/review/twin_structure.py $D/hla_clean.csv > $O/twin_structure.log
$PY analysis/review/dedup_coverage.py "HLA Data combined.xlsx" "HLA Data.cleaned.xlsx" $O > $O/dedup_coverage.log
# 2. genotype-level estimators and the floor sweep (live 1e-6 EM)
$PY analysis/review/genotype_coverage.py $D/hla_clean.csv $O > $O/run1.log
$PY analysis/review/floor_sweep.py $D/haplo_freqs_em.csv $O/genotype_coverage_targets.csv $O/floor_sweep.csv > $O/floor_sweep.log
$PY analysis/review/diagnose_gap.py $D/hla_clean.csv > $O/diagnose_gap.log
$PY analysis/review/model_vs_probe.py $D/haplo_freqs_em.csv $O/dedup_coverage_targets.csv > $O/model_vs_probe.log
$PY analysis/review/publishable_range.py $D/haplo_freqs_em.csv $O/dedup_coverage_targets.csv $D/registry_size_targets.csv analysis/snapshot_1e-3/data/registry_size_targets.csv > $O/publishable_range.log
# 3. EM refit on deduplicated BMDP (about 10 min), half-sample intervals, sample-size projection
$PY analysis/review/dedup_em.py . $O > $O/dedup_em.log
$PY analysis/review/subsample_em.py . $O 100 12 > $O/subsample_em.log
$PY analysis/review/sample_size_projection.py . $O 12 20 > $O/sample_size_projection.log
$PY analysis/review/calibrated_projection.py $O > $O/calibrated_projection.log
```

`dedup_em.py . $O 100` also runs the with-replacement bootstrap, which takes about 3.5 h and
writes `dedup_boot.log`. It is kept only as evidence that this bootstrap is biased; do not use its
intervals. The `out/` logs are the exact outputs that the review quotes.
