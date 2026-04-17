"""tests/test_ingest.py"""
import sys
import pathlib
import importlib.util
import numpy as np
import pandas as pd
import tempfile, os

_spec = importlib.util.spec_from_file_location(
    "ingest_mod",
    pathlib.Path(__file__).parent.parent / "analysis" / "01_ingest.py"
)
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)

normalize_allele      = _mod.normalize_allele
map_ethnicity         = _mod.map_ethnicity
detect_allele_columns = _mod.detect_allele_columns
load_txt_haplotypes   = _mod.load_txt_haplotypes


def test_normalize_allele_2field():
    assert normalize_allele("11:01") == "11:01"

def test_normalize_allele_strips_G_suffix():
    assert normalize_allele("03:01:01G") == "03:01"

def test_normalize_allele_missing_dash():
    assert pd.isna(normalize_allele("-"))

def test_normalize_allele_empty_string():
    assert pd.isna(normalize_allele(""))

def test_normalize_allele_4field():
    assert normalize_allele("11:01:02:03") == "11:01"

def test_map_ethnicity_single_codes():
    assert map_ethnicity("C") == "Chinese"
    assert map_ethnicity("M") == "Malay"
    assert map_ethnicity("I") == "Indian"
    assert map_ethnicity("O") == "Others"

def test_map_ethnicity_full_names_case_insensitive():
    assert map_ethnicity("Chinese") == "Chinese"
    assert map_ethnicity("MALAY") == "Malay"
    assert map_ethnicity("indian") == "Indian"

def test_detect_allele_columns_standard_names():
    df = pd.DataFrame(columns=[
        "A1","A2","B1","B2","C1","C2","DRB1_1","DRB1_2","DQB1_1","DQB1_2","Ethnicity"
    ])
    m = detect_allele_columns(df)
    assert m["A1"]     == ("HLA-A", 1)
    assert m["A2"]     == ("HLA-A", 2)
    assert m["DRB1_1"] == ("DRB1",  1)
    assert m["DQB1_2"] == ("DQB1",  2)
    assert "Ethnicity" not in m

def test_detect_allele_columns_no_underscore():
    df = pd.DataFrame(columns=["A1","A2","B1","B2","C1","C2","DRB11","DRB12","DQB11","DQB12","Ethnicity"])
    m = detect_allele_columns(df)
    assert m["DRB11"] == ("DRB1", 1)
    assert m["DRB12"] == ("DRB1", 2)
    assert m["DQB11"] == ("DQB1", 1)
    assert m["DQB12"] == ("DQB1", 2)
    assert "Ethnicity" not in m

def test_detect_allele_columns_with_hla_prefix():
    df = pd.DataFrame(columns=["HLA-A1","HLA-A2","HLA-B1","HLA-B2"])
    m = detect_allele_columns(df)
    assert m["HLA-A1"] == ("HLA-A", 1)
    assert m["HLA-B2"] == ("HLA-B", 2)

def test_load_txt_haplotypes_one_row_five_loci():
    content = "M\t11:01\t15:02\t08:01\t12:02\t03:01\n"
    with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
        f.write(content)
        fname = f.name
    try:
        df = load_txt_haplotypes(pathlib.Path(fname), "TEST")
        assert len(df) == 5
        assert set(df["locus"]) == {"HLA-A","HLA-B","HLA-C","DRB1","DQB1"}
        assert df.loc[df["locus"]=="HLA-A","allele1"].values[0] == "11:01"
        assert pd.isna(df["allele2"].values[0])
        assert df["ethnicity"].values[0] == "Malay"
    finally:
        os.unlink(fname)

def test_load_txt_haplotypes_missing_allele():
    content = "C\t11:01\t15:02\t08:01\t12:02\t-\n"
    with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
        f.write(content)
        fname = f.name
    try:
        df = load_txt_haplotypes(pathlib.Path(fname), "TEST")
        dqb1_row = df[df["locus"] == "DQB1"]
        assert pd.isna(dqb1_row["allele1"].values[0])
    finally:
        os.unlink(fname)
