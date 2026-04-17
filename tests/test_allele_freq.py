"""tests/test_allele_freq.py"""
import pathlib, importlib.util
import numpy as np
import pandas as pd

_spec = importlib.util.spec_from_file_location(
    "af_mod",
    pathlib.Path(__file__).parent.parent / "analysis" / "02_allele_freq.py"
)
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)

compute_allele_frequencies = _mod.compute_allele_frequencies
compare_frequencies        = _mod.compare_frequencies


def _make_clean_df():
    return pd.DataFrame([
        {"sample_id":"S1","source":"BMDP_OUT","ethnicity":"Chinese","locus":"HLA-A",
         "allele1":"11:01","allele2":"24:02"},
        {"sample_id":"S2","source":"BMDP_OUT","ethnicity":"Chinese","locus":"HLA-A",
         "allele1":"11:01","allele2":"33:03"},
        {"sample_id":"S3","source":"BMDP_OUT","ethnicity":"Chinese","locus":"HLA-A",
         "allele1":"11:01","allele2":np.nan},
    ])


def test_compute_allele_frequencies_counts_correctly():
    df = _make_clean_df()
    freq = compute_allele_frequencies(df)
    row = freq[(freq["ethnicity"]=="Chinese") & (freq["locus"]=="HLA-A") &
               (freq["allele"]=="11:01")]
    assert len(row) == 1
    assert abs(row["frequency"].values[0] - 3/5) < 1e-9


def test_compute_allele_frequencies_excludes_nan():
    df = pd.DataFrame([
        {"sample_id":"S1","source":"BMDP_OUT","ethnicity":"Chinese","locus":"HLA-A",
         "allele1":"11:01","allele2": np.nan},
        {"sample_id":"S2","source":"BMDP_OUT","ethnicity":"Chinese","locus":"HLA-A",
         "allele1": np.nan,"allele2":"24:02"},
    ])
    freq = compute_allele_frequencies(df)
    total = freq[(freq["ethnicity"]=="Chinese") & (freq["locus"]=="HLA-A")]["frequency"].sum()
    assert abs(total - 1.0) < 1e-9


def test_compute_allele_frequencies_excludes_hsa_sources():
    df = pd.DataFrame([
        {"sample_id":"S1","source":"BMDP_OUT","ethnicity":"Chinese","locus":"HLA-A",
         "allele1":"11:01","allele2":"24:02"},
        {"sample_id":"S2","source":"HSA-Donor","ethnicity":"Chinese","locus":"HLA-A",
         "allele1":"33:03","allele2": np.nan},
    ])
    freq = compute_allele_frequencies(df)
    alleles = set(freq[freq["ethnicity"]=="Chinese"]["allele"])
    assert "33:03" not in alleles


def test_compare_frequencies_flags_large_discrepancy():
    observed = pd.DataFrame([
        {"ethnicity":"Chinese","locus":"HLA-A","allele":"11:01","frequency":0.30}
    ])
    published = pd.DataFrame([
        {"ethnicity":"Chinese","locus":"HLA-A","allele":"11:01","pub_frequency":0.25}
    ])
    result = compare_frequencies(observed, published, threshold=0.005)
    assert result["flagged"].values[0] == True
    assert abs(result["difference"].values[0] - 0.05) < 1e-9


def test_compare_frequencies_does_not_flag_small_discrepancy():
    observed = pd.DataFrame([
        {"ethnicity":"Chinese","locus":"HLA-A","allele":"11:01","frequency":0.251}
    ])
    published = pd.DataFrame([
        {"ethnicity":"Chinese","locus":"HLA-A","allele":"11:01","pub_frequency":0.250}
    ])
    result = compare_frequencies(observed, published, threshold=0.005)
    assert result["flagged"].values[0] == False
