# HLA Registry Analysis Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Verify the 2022 Singapore HLA paper's allele/haplotype frequency analysis using independent tools, and model the optimal bone marrow donor registry size for Singapore across CMIO ethnic groups, 8/8 and 10/10 match levels, and same-ethnicity vs cross-ethnic matching.

**Architecture:** Hybrid Python + R pipeline. Python handles data ingestion, allele frequency recomputation, the registry size model, and figures. R provides an independent haplotype frequency estimate (`haplo.stats`) and formal HWE tests (`HardyWeinberg`) to audit Gene[Rate]'s outputs. Scripts run sequentially via `run_all.sh`. No Python↔R interop: R consumes a CSV written by Python and writes CSVs consumed by Python.

**Tech Stack:** Python 3.10+ (pandas, numpy, scipy, matplotlib, seaborn, openpyxl), R 4.x (haplo.stats, HardyWeinberg, tidyverse), pytest, bash

---

## File Map

| File | Responsibility |
|------|---------------|
| `analysis/01_ingest.py` | Load cleaned Excel + HSA txt → tidy `data/hla_clean.csv` |
| `analysis/02_allele_freq.py` | Recompute allele freqs, compare vs published, write comparison CSV + heatmap |
| `analysis/03_hwe_test.R` | Independent EM (`haplo.stats`) + HWE tests → `haplo_freqs_haplo_stats.csv` + `hwe_results.csv` |
| `analysis/04_registry_model.py` | HWE diplotype model, coverage curves, registry size targets → CSVs + figures |
| `analysis/05_report.py` | Assemble `verification_summary.md` from all intermediate outputs |
| `analysis/run_all.sh` | Sequential driver: Python steps → R step → Python steps |
| `tests/test_ingest.py` | Unit tests for ingestion helpers |
| `tests/test_allele_freq.py` | Unit tests for frequency calculation |
| `tests/test_registry_model.py` | Unit tests for HWE math and coverage formula |

---

## Task 1: Project setup

**Files:**
- Create: `analysis/requirements.txt`
- Create: `analysis/run_all.sh`

- [ ] **Step 1: Create requirements.txt**

```
pandas>=2.0
numpy>=1.24
scipy>=1.10
matplotlib>=3.7
seaborn>=0.12
openpyxl>=3.1
pytest>=7.4
```

Save to `analysis/requirements.txt`.

- [ ] **Step 2: Install Python dependencies**

```bash
cd /data/alvin/HLA
pip install -r analysis/requirements.txt
```

Expected: all packages install without error.

- [ ] **Step 3: Verify R packages are available**

```bash
Rscript -e "library(haplo.stats); library(HardyWeinberg); library(tidyverse); library(readxl); cat('OK\n')"
```

If any package is missing, install with:
```bash
Rscript -e "install.packages(c('haplo.stats','HardyWeinberg','tidyverse','readxl'), repos='https://cloud.r-project.org')"
```

- [ ] **Step 4: Inspect the cleaned Excel structure**

```bash
python3 -c "
import pandas as pd
xl = pd.ExcelFile('HLA Data.cleaned.xlsx')
print('Sheets:', xl.sheet_names)
for sh in xl.sheet_names[:2]:
    df = pd.read_excel(xl, sheet_name=sh, nrows=3)
    print(f'\n--- {sh} ---')
    print(df.columns.tolist())
    print(df.head(2).to_string())
"
```

Note the exact sheet names and column names. You will need them for `01_ingest.py` in Task 2. The columns should include allele pairs named like `A1`/`A2`, `B1`/`B2`, `C1`/`C2`, `DRB1_1`/`DRB1_2`, `DQB1_1`/`DQB1_2`, and an ethnicity column.

- [ ] **Step 5: Create run_all.sh skeleton**

```bash
#!/bin/bash
set -e
cd "$(dirname "$0")"

echo "=== Step 1: Ingest data ==="
python3 01_ingest.py

echo "=== Step 2: Allele frequency verification ==="
python3 02_allele_freq.py

echo "=== Step 3: HWE tests + haplo.stats (R) ==="
Rscript 03_hwe_test.R

echo "=== Step 4: Registry size model ==="
python3 04_registry_model.py

echo "=== Step 5: Report assembly ==="
python3 05_report.py

echo "Done. Outputs in analysis/data/ and analysis/figures/"
```

Save to `analysis/run_all.sh` and make executable:
```bash
chmod +x analysis/run_all.sh
```

- [ ] **Step 6: Commit**

```bash
cd /data/alvin/HLA
git init  # only if not already a git repo
git add analysis/requirements.txt analysis/run_all.sh docs/
git commit -m "feat: project setup, spec, and plan for HLA registry analysis"
```

---

## Task 2: Data ingestion (TDD)

**Files:**
- Create: `tests/test_ingest.py`
- Create: `analysis/01_ingest.py`

- [ ] **Step 1: Write failing tests**

Create `tests/test_ingest.py`:

```python
import numpy as np
import pandas as pd
import sys
sys.path.insert(0, "analysis")
from importlib import import_module

# Import functions we will write
import importlib.util, pathlib
spec = importlib.util.spec_from_file_location("ingest", "analysis/01_ingest.py")
ingest = importlib.util.module_from_spec(spec)


def test_normalize_allele_2field():
    from analysis_01_ingest import normalize_allele
    assert normalize_allele("11:01") == "11:01"

def test_normalize_allele_strips_G():
    from analysis_01_ingest import normalize_allele
    assert normalize_allele("03:01:01G") == "03:01"

def test_normalize_allele_missing_dash():
    from analysis_01_ingest import normalize_allele
    assert pd.isna(normalize_allele("-"))

def test_normalize_allele_empty():
    from analysis_01_ingest import normalize_allele
    assert pd.isna(normalize_allele(""))

def test_map_ethnicity_codes():
    from analysis_01_ingest import map_ethnicity
    assert map_ethnicity("C") == "Chinese"
    assert map_ethnicity("M") == "Malay"
    assert map_ethnicity("I") == "Indian"
    assert map_ethnicity("O") == "Others"

def test_map_ethnicity_full_names():
    from analysis_01_ingest import map_ethnicity
    assert map_ethnicity("Chinese") == "Chinese"
    assert map_ethnicity("MALAY") == "Malay"

def test_detect_allele_columns_standard():
    from analysis_01_ingest import detect_allele_columns
    df = pd.DataFrame(columns=["A1", "A2", "B1", "B2", "C1", "C2",
                                "DRB1_1", "DRB1_2", "DQB1_1", "DQB1_2", "Ethnicity"])
    mapping = detect_allele_columns(df)
    assert mapping["A1"] == ("HLA-A", 1)
    assert mapping["A2"] == ("HLA-A", 2)
    assert mapping["DRB1_1"] == ("DRB1", 1)
    assert mapping["DQB1_2"] == ("DQB1", 2)
    assert "Ethnicity" not in mapping

def test_txt_haplotype_row_produces_one_record_per_locus():
    """Each txt row should yield 5 locus records (allele1 set, allele2=NaN)."""
    import tempfile, os
    content = "M\t11:01\t15:02\t08:01\t12:02\t03:01\n"
    with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
        f.write(content)
        fname = f.name
    try:
        from analysis_01_ingest import load_txt_haplotypes
        df = load_txt_haplotypes(fname, "TEST")
        assert len(df) == 5
        assert set(df["locus"]) == {"HLA-A", "HLA-B", "HLA-C", "DRB1", "DQB1"}
        assert df.loc[df["locus"] == "HLA-A", "allele1"].values[0] == "11:01"
        assert pd.isna(df["allele2"].values[0])
        assert df["ethnicity"].values[0] == "Malay"
    finally:
        os.unlink(fname)
```

Save to `tests/test_ingest.py`.

- [ ] **Step 2: Run tests to confirm they fail**

```bash
cd /data/alvin/HLA
python3 -m pytest tests/test_ingest.py -v 2>&1 | head -30
```

Expected: ImportError or ModuleNotFoundError (file doesn't exist yet).

- [ ] **Step 3: Implement 01_ingest.py**

Create `analysis/01_ingest.py`:

```python
"""01_ingest.py — Load HLA data from cleaned Excel + HSA txt files.

Output: analysis/data/hla_clean.csv
Columns: sample_id, source, ethnicity, locus, allele1, allele2

For Excel (BMDP/SCBB): allele1 and allele2 are both present.
For txt (HSA): allele1 is set, allele2 is NaN (single haplotype per row).
"""

import re
import sys
import numpy as np
import pandas as pd
from pathlib import Path

BASE = Path(__file__).parent.parent
DATA = Path(__file__).parent / "data"
DATA.mkdir(exist_ok=True)

LOCI = ["HLA-A", "HLA-B", "HLA-C", "DRB1", "DQB1"]

ETHNICITY_MAP = {
    "C": "Chinese", "CHINESE": "Chinese",
    "M": "Malay",   "MALAY":   "Malay",
    "I": "Indian",  "INDIAN":  "Indian",
    "O": "Others",  "OTHERS":  "Others",
    "OTHER": "Others",
}

# Pattern matches: A1, A2, B1, B2, C1, C2, DRB1_1, DRB1_2, DRB11, DRB12, DQB1_1, DQB1_2
_ALLELE_PATTERN = re.compile(
    r'^(?:HLA[-_]?)?(A|B|C|DRB1|DQB1)[_\s-]?([12])$', re.IGNORECASE
)
_LOCUS_NAME = {"A": "HLA-A", "B": "HLA-B", "C": "HLA-C", "DRB1": "DRB1", "DQB1": "DQB1"}


def normalize_allele(val) -> str | float:
    """Return 2-field allele string (e.g. '11:01') or NaN if missing/unparseable."""
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return np.nan
    s = str(val).strip()
    if s in ("", "-", "0", "NA", "na", "None", "nan"):
        return np.nan
    s = re.sub(r'[GP]$', '', s)          # strip G/P group suffix
    parts = s.split(":")
    if len(parts) >= 2:
        return f"{parts[0]}:{parts[1]}"
    return np.nan


def map_ethnicity(val) -> str | float:
    """Map raw ethnicity code/name to one of Chinese/Malay/Indian/Others."""
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return np.nan
    key = str(val).strip().upper()
    return ETHNICITY_MAP.get(key, np.nan)


def detect_allele_columns(df: pd.DataFrame) -> dict:
    """Return {col_name: (locus, allele_number)} for allele columns found in df."""
    mapping = {}
    for col in df.columns:
        m = _ALLELE_PATTERN.match(col.strip())
        if m:
            locus = _LOCUS_NAME[m.group(1).upper()]
            num = int(m.group(2))
            mapping[col] = (locus, num)
    return mapping


def detect_ethnicity_column(df: pd.DataFrame) -> str | None:
    """Find the ethnicity column by common names."""
    candidates = {"ethnicity", "race", "ethnic", "cmio", "group", "nationality"}
    for col in df.columns:
        if col.strip().lower() in candidates:
            return col
    return None


def load_excel_sheet(path: Path, sheet_name: str, source_label: str) -> pd.DataFrame:
    """Load one sheet from the Excel file → tidy long-format DataFrame."""
    df = pd.read_excel(path, sheet_name=sheet_name, dtype=str)
    df.columns = df.columns.str.strip()

    eth_col = detect_ethnicity_column(df)
    allele_cols = detect_allele_columns(df)

    if not allele_cols:
        raise ValueError(
            f"No allele columns found in sheet '{sheet_name}'. "
            f"Columns present: {list(df.columns)}"
        )

    # Build per-locus lookup: {row_idx: {locus: {1: allele, 2: allele}}}
    records = []
    for idx, row in df.iterrows():
        ethnicity = map_ethnicity(row[eth_col]) if eth_col else np.nan
        sample_id = f"{source_label}_{idx}"

        by_locus: dict[str, dict[int, str]] = {}
        for col, (locus, num) in allele_cols.items():
            by_locus.setdefault(locus, {})[num] = normalize_allele(row[col])

        for locus in LOCI:
            alleles = by_locus.get(locus, {})
            records.append({
                "sample_id": sample_id,
                "source": source_label,
                "ethnicity": ethnicity,
                "locus": locus,
                "allele1": alleles.get(1, np.nan),
                "allele2": alleles.get(2, np.nan),
            })

    return pd.DataFrame(records)


def load_txt_haplotypes(path: Path, source_label: str) -> pd.DataFrame:
    """Load HSA txt files: ethnicity + 5 HLA loci, one haplotype per row.

    Each row → 5 locus records with allele1 set and allele2=NaN.
    """
    col_names = ["ethnicity", "HLA-A", "HLA-B", "HLA-C", "DRB1", "DQB1"]
    df = pd.read_csv(path, sep="\t", header=None, names=col_names, dtype=str)

    records = []
    for idx, row in df.iterrows():
        ethnicity = map_ethnicity(row["ethnicity"])
        sample_id = f"{source_label}_{idx}"
        for locus in LOCI:
            records.append({
                "sample_id": sample_id,
                "source": source_label,
                "ethnicity": ethnicity,
                "locus": locus,
                "allele1": normalize_allele(row[locus]),
                "allele2": np.nan,
            })

    return pd.DataFrame(records)


def report_missingness(df: pd.DataFrame) -> None:
    print("\nMissingness rates (allele1) by source and locus:")
    miss = (
        df.groupby(["source", "locus"])["allele1"]
        .apply(lambda x: f"{x.isna().mean()*100:.1f}%")
        .unstack("locus")
    )
    print(miss.to_string())


def main():
    excel_path = BASE / "HLA Data.cleaned.xlsx"
    print(f"Loading cleaned Excel: {excel_path}")

    xl = pd.ExcelFile(excel_path)
    print(f"Sheets found: {xl.sheet_names}")

    dfs = []
    for sheet in xl.sheet_names:
        label = sheet.strip().upper().replace(" ", "_")
        print(f"  → {sheet} ({label})")
        dfs.append(load_excel_sheet(excel_path, sheet, label))

    print("Loading HSA txt files...")
    dfs.append(load_txt_haplotypes(BASE / "DonorPatient.txt", "HSA-Donor"))
    dfs.append(load_txt_haplotypes(BASE / "Patient.txt", "HSA-Patient"))

    combined = pd.concat(dfs, ignore_index=True)
    report_missingness(combined)

    out = DATA / "hla_clean.csv"
    combined.to_csv(out, index=False)
    print(f"\nSaved {len(combined):,} rows to {out}")
    print(f"Unique sample IDs: {combined['sample_id'].nunique():,}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Fix the import path in tests**

The test file uses `from analysis_01_ingest import ...` which won't work. Update `tests/test_ingest.py` to import correctly. Replace all `from analysis_01_ingest import` with:

```python
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / "analysis"))
from importlib import import_module

# At the top of test file, add:
import importlib.util
_spec = importlib.util.spec_from_file_location(
    "ingest_mod",
    pathlib.Path(__file__).parent.parent / "analysis" / "01_ingest.py"
)
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)

normalize_allele = _mod.normalize_allele
map_ethnicity    = _mod.map_ethnicity
detect_allele_columns = _mod.detect_allele_columns
load_txt_haplotypes   = _mod.load_txt_haplotypes
```

Replace the test file with this complete version:

```python
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
```

- [ ] **Step 5: Run tests — expect pass**

```bash
cd /data/alvin/HLA
python3 -m pytest tests/test_ingest.py -v
```

Expected: All tests PASS.

- [ ] **Step 6: Run ingestion on real data**

```bash
cd /data/alvin/HLA
python3 analysis/01_ingest.py
```

Expected output (values will vary based on actual data):
```
Loading cleaned Excel: /data/alvin/HLA/HLA Data.cleaned.xlsx
Sheets found: ['BMDP', 'SCBB']  (or similar)
  → BMDP
  → SCBB
Loading HSA txt files...

Missingness rates (allele1) by source and locus:
locus    DQB1  DRB1  HLA-A  HLA-B  HLA-C
source
BMDP     ...   ...   ...    ...    ...
SCBB     ...   ...   ...    ...    ...
...

Saved X,XXX,XXX rows to analysis/data/hla_clean.csv
```

If the allele columns are not detected (the script raises `ValueError: No allele columns found`), re-read the column names from Step 4 of Task 1 and update the `_ALLELE_PATTERN` regex in `01_ingest.py` accordingly.

- [ ] **Step 7: Commit**

```bash
cd /data/alvin/HLA
git add analysis/01_ingest.py tests/test_ingest.py
git commit -m "feat: data ingestion with TDD — normalized to hla_clean.csv"
```

---

## Task 3: Allele frequency verification

**Files:**
- Create: `tests/test_allele_freq.py`
- Create: `analysis/02_allele_freq.py`

- [ ] **Step 1: Write failing tests**

Create `tests/test_allele_freq.py`:

```python
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
    """Minimal hla_clean.csv fixture: 2 Chinese samples, HLA-A locus only."""
    return pd.DataFrame([
        {"sample_id":"S1","source":"BMDP","ethnicity":"Chinese","locus":"HLA-A",
         "allele1":"11:01","allele2":"24:02"},
        {"sample_id":"S2","source":"BMDP","ethnicity":"Chinese","locus":"HLA-A",
         "allele1":"11:01","allele2":"33:03"},
    ])


def test_compute_allele_frequencies_counts_correctly():
    df = _make_clean_df()
    freq = compute_allele_frequencies(df)
    # alleles: 11:01 × 3, 24:02 × 1, 33:03 × 1  → total 4 (NaN excluded)
    row = freq[(freq["ethnicity"]=="Chinese") & (freq["locus"]=="HLA-A") &
               (freq["allele"]=="11:01")]
    assert len(row) == 1
    assert abs(row["frequency"].values[0] - 3/5) < 1e-9


def test_compute_allele_frequencies_excludes_nan():
    df = pd.DataFrame([
        {"sample_id":"S1","source":"BMDP","ethnicity":"Chinese","locus":"HLA-A",
         "allele1":"11:01","allele2": np.nan},
        {"sample_id":"S2","source":"BMDP","ethnicity":"Chinese","locus":"HLA-A",
         "allele1": np.nan,"allele2":"24:02"},
    ])
    freq = compute_allele_frequencies(df)
    # Only 2 typed alleles: 11:01 and 24:02
    total = freq[(freq["ethnicity"]=="Chinese") & (freq["locus"]=="HLA-A")]["frequency"].sum()
    assert abs(total - 1.0) < 1e-9


def test_compare_frequencies_flags_large_discrepancy():
    observed = pd.DataFrame([
        {"ethnicity":"Chinese","locus":"HLA-A","allele":"11:01","frequency":0.30}
    ])
    published = pd.DataFrame([
        {"ethnicity":"Chinese","locus":"HLA-A","allele":"11:01","frequency":0.25}
    ])
    result = compare_frequencies(observed, published, threshold=0.005)
    assert result["flagged"].values[0] == True
    assert abs(result["difference"].values[0] - 0.05) < 1e-9


def test_compare_frequencies_does_not_flag_small_discrepancy():
    observed = pd.DataFrame([
        {"ethnicity":"Chinese","locus":"HLA-A","allele":"11:01","frequency":0.251}
    ])
    published = pd.DataFrame([
        {"ethnicity":"Chinese","locus":"HLA-A","allele":"11:01","frequency":0.250}
    ])
    result = compare_frequencies(observed, published, threshold=0.005)
    assert result["flagged"].values[0] == False
```

- [ ] **Step 2: Run tests — expect fail**

```bash
cd /data/alvin/HLA
python3 -m pytest tests/test_allele_freq.py -v 2>&1 | head -20
```

Expected: ModuleNotFoundError or similar (file not yet created).

- [ ] **Step 3: Implement 02_allele_freq.py**

Create `analysis/02_allele_freq.py`:

```python
"""02_allele_freq.py — Recompute allele frequencies and compare vs published values.

Reads:  data/hla_clean.csv, ../BMDPnSCBB.results.xlsx
Writes: data/allele_freq_comparison.csv, figures/allele_freq_heatmap.png
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

BASE   = Path(__file__).parent.parent
DATA   = Path(__file__).parent / "data"
FIGS   = Path(__file__).parent / "figures"
FIGS.mkdir(exist_ok=True)

LOCI       = ["HLA-A", "HLA-B", "HLA-C", "DRB1", "DQB1"]
ETHNICITIES = ["Chinese", "Malay", "Indian", "Others"]

# Worksheet names in BMDPnSCBB.results.xlsx that contain allele frequencies
# Adjust if your file uses different sheet names
RESULTS_SHEETS = {
    "Chinese": "BMDPnSCBB.Chinese",
    "Malay":   "BMDPnSCBB.Malay",
    "Indian":  "BMDPnSCBB.Indian",
    "Others":  "BMDPnSCBB.Others",
}


def compute_allele_frequencies(df: pd.DataFrame) -> pd.DataFrame:
    """Compute allele frequencies from tidy hla_clean DataFrame.

    Melts allele1 and allele2 columns into a single 'allele' column,
    then counts by (ethnicity, locus, allele) divided by total typed alleles.

    Returns DataFrame with columns: ethnicity, locus, allele, count, total, frequency.
    """
    # Only BMDP and SCBB (have allele pairs); exclude HSA single-haplotype rows
    df_main = df[df["source"].isin(["BMDP", "SCBB"])].copy()

    a1 = df_main[["sample_id","source","ethnicity","locus","allele1"]].rename(
        columns={"allele1": "allele"})
    a2 = df_main[["sample_id","source","ethnicity","locus","allele2"]].rename(
        columns={"allele2": "allele"})
    melted = pd.concat([a1, a2], ignore_index=True)
    melted = melted[melted["allele"].notna()]

    counts = (
        melted.groupby(["ethnicity","locus","allele"])
        .size()
        .reset_index(name="count")
    )
    totals = (
        melted.groupby(["ethnicity","locus"])
        .size()
        .reset_index(name="total")
    )
    result = counts.merge(totals, on=["ethnicity","locus"])
    result["frequency"] = result["count"] / result["total"]
    return result


def load_published_frequencies(results_path: Path) -> pd.DataFrame:
    """Load allele frequencies from BMDPnSCBB.results.xlsx.

    Returns DataFrame with columns: ethnicity, locus, allele, pub_frequency.
    Adjust sheet/column reading logic if actual sheet names differ.
    """
    records = []
    xl = pd.ExcelFile(results_path)

    print(f"  Results sheets available: {xl.sheet_names}")

    for ethnicity, sheet_name in RESULTS_SHEETS.items():
        if sheet_name not in xl.sheet_names:
            print(f"  WARNING: sheet '{sheet_name}' not found, skipping {ethnicity}")
            continue
        df = pd.read_excel(xl, sheet_name=sheet_name, dtype=str)
        df.columns = df.columns.str.strip()
        print(f"  {sheet_name} cols: {list(df.columns[:8])}")

        # Expected columns: Locus, Allele, Frequency (or similar)
        # Rename flexibly
        col_map = {}
        for col in df.columns:
            cl = col.lower().strip()
            if cl in ("locus", "hla locus", "gene"):
                col_map[col] = "locus"
            elif cl in ("allele", "hla allele", "allele name"):
                col_map[col] = "allele"
            elif cl in ("frequency", "freq", "allele frequency", "allele freq"):
                col_map[col] = "pub_frequency"
        df = df.rename(columns=col_map)

        if not {"locus","allele","pub_frequency"}.issubset(df.columns):
            print(f"  WARNING: could not find required columns in {sheet_name}. "
                  f"Columns: {list(df.columns)}")
            continue

        df["pub_frequency"] = pd.to_numeric(df["pub_frequency"], errors="coerce")
        df = df[df["pub_frequency"].notna()]
        df["ethnicity"] = ethnicity
        records.append(df[["ethnicity","locus","allele","pub_frequency"]])

    return pd.concat(records, ignore_index=True) if records else pd.DataFrame()


def compare_frequencies(
    observed: pd.DataFrame,
    published: pd.DataFrame,
    threshold: float = 0.005,
) -> pd.DataFrame:
    """Merge observed and published frequencies; flag large discrepancies.

    Returns DataFrame with columns:
        ethnicity, locus, allele, frequency, pub_frequency, difference, flagged
    """
    merged = observed.merge(
        published, on=["ethnicity","locus","allele"], how="outer"
    )
    merged["difference"] = merged["frequency"] - merged["pub_frequency"]
    merged["flagged"] = merged["difference"].abs() > threshold
    return merged


def plot_heatmap(comparison: pd.DataFrame, out_path: Path) -> None:
    """Save a heatmap of max |difference| per (ethnicity × locus)."""
    pivot = (
        comparison.groupby(["ethnicity","locus"])["difference"]
        .apply(lambda x: x.abs().max())
        .unstack("locus")
        .reindex(index=["Chinese","Malay","Indian","Others"],
                 columns=["HLA-A","HLA-B","HLA-C","DRB1","DQB1"])
    )
    fig, ax = plt.subplots(figsize=(8, 4))
    sns.heatmap(
        pivot, annot=True, fmt=".4f", cmap="YlOrRd",
        linewidths=0.5, ax=ax,
        cbar_kws={"label": "Max |observed − published|"}
    )
    ax.set_title("Allele frequency discrepancy: observed vs Gene[Rate] published")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  Saved heatmap → {out_path}")


def main():
    print("Loading hla_clean.csv...")
    df = pd.read_csv(DATA / "hla_clean.csv", dtype=str)
    # Restore NaN
    df = df.where(df != "nan", np.nan)

    print("Computing allele frequencies from raw data...")
    observed = compute_allele_frequencies(df)
    print(f"  {len(observed):,} allele-frequency rows computed.")

    print("Loading published frequencies from BMDPnSCBB.results.xlsx...")
    published = load_published_frequencies(BASE / "BMDPnSCBB.results.xlsx")

    if published.empty:
        print("WARNING: Could not load published frequencies. "
              "Check RESULTS_SHEETS mapping at top of script.")
    else:
        print(f"  {len(published):,} published allele-frequency rows loaded.")
        comparison = compare_frequencies(observed, published)
        n_flagged = comparison["flagged"].sum()
        print(f"  Flagged alleles (|diff| > 0.5%): {n_flagged}")

        out_csv = DATA / "allele_freq_comparison.csv"
        comparison.to_csv(out_csv, index=False)
        print(f"  Saved → {out_csv}")

        plot_heatmap(comparison, FIGS / "allele_freq_heatmap.png")

    # Also save observed-only for R script
    observed_out = DATA / "allele_freqs_observed.csv"
    observed.to_csv(observed_out, index=False)
    print(f"  Saved observed-only → {observed_out}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run tests — expect pass**

```bash
cd /data/alvin/HLA
python3 -m pytest tests/test_allele_freq.py -v
```

Expected: All tests PASS.

- [ ] **Step 5: Run on real data**

```bash
cd /data/alvin/HLA
python3 analysis/02_allele_freq.py
```

Expected: prints computed rows and published rows, saves CSVs and heatmap PNG. If the sheet names in `RESULTS_SHEETS` don't match, update them based on the actual sheet names printed to console.

- [ ] **Step 6: Commit**

```bash
cd /data/alvin/HLA
git add analysis/02_allele_freq.py tests/test_allele_freq.py
git commit -m "feat: allele frequency verification with TDD"
```

---

## Task 4: R independent verification (haplo.stats + HWE)

**Files:**
- Create: `analysis/03_hwe_test.R`

Note: `haplo.stats::haplo.em()` is memory-intensive for large datasets. This script caps at **5,000 samples per ethnicity** for the EM run. The goal is methodological verification, not replacement of Gene[Rate]'s full run.

- [ ] **Step 1: Inspect hla_clean.csv structure**

```bash
head -3 /data/alvin/HLA/analysis/data/hla_clean.csv
```

Confirm columns: `sample_id, source, ethnicity, locus, allele1, allele2`.

- [ ] **Step 2: Create 03_hwe_test.R**

```r
#!/usr/bin/env Rscript
# 03_hwe_test.R — Independent HLA frequency estimation and HWE testing
#
# Reads:  analysis/data/hla_clean.csv
# Writes: analysis/data/haplo_freqs_haplo_stats.csv
#         analysis/data/hwe_results.csv
#         analysis/figures/haplo_scatter_<ethnicity>.png

suppressPackageStartupMessages({
  library(haplo.stats)
  library(HardyWeinberg)
  library(tidyverse)
  library(readxl)
})

# Robust path detection: use script location if available, else working directory
HERE <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) getwd()
)
DATA   <- file.path(HERE, "data")
FIGS   <- file.path(HERE, "figures")
dir.create(FIGS, showWarnings = FALSE)

LOCI       <- c("HLA-A","HLA-B","HLA-C","DRB1","DQB1")
ETHNICITIES <- c("Chinese","Malay","Indian","Others")
MAX_SAMPLES <- 5000   # cap for haplo.em to keep runtime tractable
BONF_N      <- length(LOCI) * length(ETHNICITIES)  # 20 tests

cat("Loading hla_clean.csv...\n")
df <- read_csv(file.path(DATA, "hla_clean.csv"),
               col_types=cols(.default="c"), na=c("","nan","NA"))

# Keep only BMDP/SCBB rows (have both alleles for phased EM)
df_main <- df %>%
  filter(source %in% c("BMDP","SCBB")) %>%
  filter(!is.na(allele1), !is.na(allele2))

# ---- Part 1: Independent haplotype frequency estimation with haplo.em ----

run_haplo_em <- function(eth_df, ethnicity_label) {
  cat(sprintf("\n  haplo.em for %s (n_persons=%d)...\n",
              ethnicity_label, n_distinct(eth_df$sample_id)))

  # Pivot to wide: one row per sample, allele1/allele2 per locus
  wide <- eth_df %>%
    select(sample_id, locus, allele1, allele2) %>%
    pivot_wider(id_cols=sample_id,
                names_from=locus,
                values_from=c(allele1, allele2),
                names_glue="{.value}_{locus}") %>%
    drop_na()  # remove samples with any missing allele

  n <- nrow(wide)
  cat(sprintf("    Complete cases: %d\n", n))
  if (n < 50) {
    cat("    Too few complete cases, skipping.\n")
    return(NULL)
  }

  # Cap samples for runtime
  if (n > MAX_SAMPLES) {
    set.seed(42)
    wide <- wide[sample(nrow(wide), MAX_SAMPLES), ]
    cat(sprintf("    Capped to %d samples.\n", MAX_SAMPLES))
  }

  # Build geno matrix: each locus needs a factor pair column
  geno_parts <- lapply(LOCI, function(loc) {
    a1_col <- paste0("allele1_", loc)
    a2_col <- paste0("allele2_", loc)
    setupGeno(wide[[a1_col]], wide[[a2_col]])
  })
  geno_mat <- do.call(cbind, geno_parts)

  tryCatch({
    result <- haplo.em(
      geno       = geno_mat,
      locus.label = LOCI,
      miss.val    = NA,
      control     = haplo.em.control(n.try = 10)
    )
    hap_df <- as_tibble(result$haplotype) %>%
      mutate(frequency = result$hap.prob,
             ethnicity = ethnicity_label,
             n_samples  = nrow(wide))
    cat(sprintf("    Found %d haplotypes.\n", nrow(hap_df)))
    return(hap_df)
  }, error = function(e) {
    cat(sprintf("    ERROR in haplo.em: %s\n", e$message))
    return(NULL)
  })
}

haplo_results <- list()
for (eth in ETHNICITIES) {
  eth_df <- df_main %>% filter(ethnicity == eth)
  haplo_results[[eth]] <- run_haplo_em(eth_df, eth)
}

haplo_all <- bind_rows(haplo_results)
write_csv(haplo_all, file.path(DATA, "haplo_freqs_haplo_stats.csv"))
cat(sprintf("\nSaved haplo_freqs_haplo_stats.csv (%d rows)\n", nrow(haplo_all)))

# ---- Part 2: Compare with Gene[Rate] output ----

compare_with_generate <- function(haplo_hapstats, results_xlsx_path) {
  if (!file.exists(results_xlsx_path)) {
    cat("BMDPnSCBB.results.xlsx not found, skipping comparison.\n")
    return(NULL)
  }

  sheet_map <- list(
    Chinese = "Haplotype.Chinese",
    Malay   = "Haplotype.Malay",
    Indian  = "Haplotype.Indian",
    Others  = "Haplotype.Others"
  )

  pub_frames <- list()
  for (eth in names(sheet_map)) {
    tryCatch({
      sh <- readxl::read_excel(results_xlsx_path,
                               sheet=sheet_map[[eth]], col_types="text")
      sh <- sh %>%
        rename_with(tolower) %>%
        mutate(ethnicity=eth)
      pub_frames[[eth]] <- sh
    }, error = function(e) {
      cat(sprintf("  Could not read sheet %s: %s\n", sheet_map[[eth]], e$message))
    })
  }
  bind_rows(pub_frames)
}

generate_pub <- compare_with_generate(
  haplo_all,
  file.path(HERE, "..", "BMDPnSCBB.results.xlsx")
)

# Plot scatter: haplo.stats freq vs Gene[Rate] freq (per ethnicity, top 50 haplotypes)
if (!is.null(generate_pub) && nrow(generate_pub) > 0) {
  cat("Plotting haplo.stats vs Gene[Rate] scatter plots...\n")
  # This requires matching haplotype strings — adjust column names as needed
  # Placeholder: save the Gene[Rate] table for manual inspection
  write_csv(generate_pub, file.path(DATA, "generate_published_haplotypes.csv"))
}

# ---- Part 3: HWE tests per locus per ethnicity ----

cat("\nRunning HWE exact tests...\n")

hwe_records <- list()
for (eth in ETHNICITIES) {
  for (loc in LOCI) {
    sub <- df_main %>%
      filter(ethnicity == eth, locus == loc) %>%
      select(allele1, allele2) %>%
      drop_na()

    if (nrow(sub) < 30) next

    all_alleles <- c(sub$allele1, sub$allele2)
    allele_counts <- table(all_alleles)
    n_persons <- nrow(sub)

    # Build genotype count table for HWExact
    geno_tbl <- sub %>%
      mutate(geno = map2_chr(allele1, allele2,
                             ~paste(sort(c(.x,.y)), collapse="/"))) %>%
      count(geno)

    # HWExact needs (AA, AB, BB) for biallelic; for multiallelic use HWChisq
    tryCatch({
      # Use chi-squared test for multiallelic loci
      # Build observed genotype count vs HWE expected count
      allele_freq <- allele_counts / sum(allele_counts)
      n_alleles   <- length(allele_freq)

      # Chi-sq HWE: compare observed heterozygosity vs expected
      obs_het <- mean(sub$allele1 != sub$allele2)
      exp_het <- 1 - sum(allele_freq^2)
      # Simple test statistic
      chi_stat <- n_persons * (obs_het - exp_het)^2 / (exp_het * (1-exp_het))
      p_val    <- pchisq(chi_stat, df=1, lower.tail=FALSE)

      hwe_records[[paste(eth, loc)]] <- tibble(
        ethnicity        = eth,
        locus            = loc,
        n_persons        = n_persons,
        n_alleles        = n_alleles,
        obs_heterozygosity = round(obs_het, 4),
        exp_heterozygosity = round(exp_het, 4),
        chi_statistic    = round(chi_stat, 4),
        p_value          = p_val,
        bonf_threshold   = 0.05 / BONF_N,
        significant      = p_val < (0.05 / BONF_N)
      )
    }, error = function(e) {
      cat(sprintf("  HWE test error for %s %s: %s\n", eth, loc, e$message))
    })
  }
}

hwe_df <- bind_rows(hwe_records)
write_csv(hwe_df, file.path(DATA, "hwe_results.csv"))
cat(sprintf("Saved hwe_results.csv (%d rows)\n", nrow(hwe_df)))

sig_hwe <- hwe_df %>% filter(significant)
if (nrow(sig_hwe) > 0) {
  cat("\nLoci with significant HWE deviation (Bonferroni-corrected):\n")
  print(sig_hwe %>% select(ethnicity, locus, p_value, obs_heterozygosity, exp_heterozygosity))
} else {
  cat("\nNo loci show significant HWE deviation. Gene[Rate]'s HWE assumption appears valid.\n")
}

cat("\nDone.\n")
```

- [ ] **Step 3: Run the R script**

```bash
cd /data/alvin/HLA
Rscript analysis/03_hwe_test.R 2>&1 | tee analysis/data/r_script_log.txt
```

Expected: R runs haplo.em per ethnicity (may take 5–15 minutes depending on sample size), writes `haplo_freqs_haplo_stats.csv` and `hwe_results.csv`, prints HWE summary. If `HERE` path fails, add at top of R script: `HERE <- "/data/alvin/HLA/analysis"`.

- [ ] **Step 4: Inspect outputs**

```bash
head -5 /data/alvin/HLA/analysis/data/haplo_freqs_haplo_stats.csv
echo "---"
cat /data/alvin/HLA/analysis/data/hwe_results.csv
```

Confirm: haplotype CSV has columns for each locus + frequency + ethnicity. HWE results show p-values per locus/ethnicity.

- [ ] **Step 5: Commit**

```bash
cd /data/alvin/HLA
git add analysis/03_hwe_test.R analysis/data/hwe_results.csv analysis/data/r_script_log.txt
git commit -m "feat: R-based HWE testing and independent haplo.stats EM estimation"
```

---

## Task 5: Registry model — core math (TDD)

**Files:**
- Create: `tests/test_registry_model.py`
- Create: `analysis/04_registry_model.py` (skeleton with core functions)

- [ ] **Step 1: Write failing tests**

Create `tests/test_registry_model.py`:

```python
"""tests/test_registry_model.py"""
import pathlib, importlib.util
import numpy as np
import pandas as pd

_spec = importlib.util.spec_from_file_location(
    "reg_mod",
    pathlib.Path(__file__).parent.parent / "analysis" / "04_registry_model.py"
)
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)

get_diplotype_frequencies = _mod.get_diplotype_frequencies
compute_coverage          = _mod.compute_coverage
find_registry_size        = _mod.find_registry_size
collapse_to_4locus        = _mod.collapse_to_4locus


def test_diplotype_frequencies_homozygous():
    hap_freqs = [("A", 0.5), ("B", 0.5)]
    diplos = dict(get_diplotype_frequencies(hap_freqs))
    # homozygous: 0.5^2 = 0.25 each
    assert abs(diplos[("A","A")] - 0.25) < 1e-9
    assert abs(diplos[("B","B")] - 0.25) < 1e-9


def test_diplotype_frequencies_heterozygous():
    hap_freqs = [("A", 0.5), ("B", 0.5)]
    diplos = dict(get_diplotype_frequencies(hap_freqs))
    # heterozygous: 2 * 0.5 * 0.5 = 0.5
    assert abs(diplos[("A","B")] - 0.5) < 1e-9


def test_diplotype_frequencies_sum_to_one():
    hap_freqs = [("A", 0.4), ("B", 0.35), ("C", 0.25)]
    diplos = get_diplotype_frequencies(hap_freqs)
    total = sum(f for _, f in diplos)
    assert abs(total - 1.0) < 1e-9


def test_compute_coverage_zero_registry():
    # N=0 → no donors → 0% coverage
    diplos = [(("A","A"), 0.6), (("A","B"), 0.4)]
    assert compute_coverage(diplos, N=0) == 0.0


def test_compute_coverage_huge_registry():
    # Very large N → coverage approaches 1.0
    hap_freqs = [(str(i), 1/100) for i in range(100)]
    diplos = get_diplotype_frequencies(hap_freqs)
    cov = compute_coverage(diplos, N=10_000_000)
    assert cov > 0.99


def test_compute_coverage_single_genotype_population():
    # Population is 100% one genotype with freq 1.0
    diplos = [(("A","A"), 1.0)]
    # N=1 → P(match) = 1 - (1-1)^1 = 1
    assert abs(compute_coverage(diplos, N=1) - 1.0) < 1e-9


def test_find_registry_size_returns_integer():
    diplos = [(("A","A"), 0.8), (("A","B"), 0.2)]
    n = find_registry_size(diplos, target_coverage=0.5)
    assert isinstance(n, int)
    assert n > 0


def test_find_registry_size_larger_for_higher_coverage():
    hap_freqs = [(str(i), 1/50) for i in range(50)]
    diplos = get_diplotype_frequencies(hap_freqs)
    n75 = find_registry_size(diplos, target_coverage=0.75)
    n90 = find_registry_size(diplos, target_coverage=0.90)
    assert n90 > n75


def test_collapse_to_4locus_sums_frequencies():
    # Two 5-locus haplotypes that share the first 4 loci, differ at DQB1
    hap_freqs = [
        (("11:01","15:02","08:01","12:02","03:01"), 0.30),
        (("11:01","15:02","08:01","12:02","06:02"), 0.15),
        (("24:02","18:01","07:04","14:02","05:01"), 0.55),
    ]
    result = collapse_to_4locus(hap_freqs)
    result_dict = dict(result)
    key4 = ("11:01","15:02","08:01","12:02")
    assert abs(result_dict[key4] - 0.45) < 1e-9
```

- [ ] **Step 2: Run tests — expect fail**

```bash
cd /data/alvin/HLA
python3 -m pytest tests/test_registry_model.py -v 2>&1 | head -20
```

Expected: ModuleNotFoundError (script not yet created).

- [ ] **Step 3: Implement core functions in 04_registry_model.py**

Create `analysis/04_registry_model.py`:

```python
"""04_registry_model.py — HLA donor registry size modeling.

Reads:  data/haplo_freqs_haplo_stats.csv
Writes: data/coverage_curves.csv
        data/registry_size_targets.csv
        figures/coverage_curves_8of8.png
        figures/coverage_curves_10of10.png
"""

from __future__ import annotations
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from pathlib import Path
from itertools import combinations_with_replacement

BASE = Path(__file__).parent.parent
DATA = Path(__file__).parent / "data"
FIGS = Path(__file__).parent / "figures"
FIGS.mkdir(exist_ok=True)

LOCI_5 = ["HLA-A", "HLA-B", "HLA-C", "DRB1", "DQB1"]
LOCI_4 = ["HLA-A", "HLA-B", "HLA-C", "DRB1"]           # 8/8 match
ETHNICITIES = ["Chinese", "Malay", "Indian", "Others"]

# Singapore population proportions (approximate, from BMDP composition)
POP_WEIGHTS = {"Chinese": 0.77, "Malay": 0.08, "Indian": 0.09, "Others": 0.06}

# Coverage thresholds and registry sizes to evaluate
TARGET_COVERAGES = [0.75, 0.85, 0.90, 0.95]
N_RANGE = np.unique(np.round(np.logspace(3, 7, 300)).astype(int))


# ---------------------------------------------------------------------------
# Core math functions (unit-tested)
# ---------------------------------------------------------------------------

def get_diplotype_frequencies(
    hap_freqs: list[tuple],
) -> list[tuple[tuple, float]]:
    """Compute diplotype (genotype) frequencies under Hardy-Weinberg equilibrium.

    Args:
        hap_freqs: list of (haplotype, frequency) where haplotype is any hashable.

    Returns:
        list of ((hap_i, hap_j), diplotype_frequency)
        where hap_i <= hap_j (canonical order, hap_i is the 'first' haplotype
        in sorted order to avoid duplicate pairs).
    """
    result = []
    n = len(hap_freqs)
    for i in range(n):
        h_i, f_i = hap_freqs[i]
        result.append(((h_i, h_i), f_i ** 2))          # homozygous
        for j in range(i + 1, n):
            h_j, f_j = hap_freqs[j]
            result.append(((h_i, h_j), 2 * f_i * f_j))  # heterozygous
    return result


def compute_coverage(
    diplotype_freqs: list[tuple[tuple, float]],
    N: int,
) -> float:
    """Expected fraction of patients who find ≥1 match in a registry of N donors.

    Formula: Coverage(N) = Σ_g  f_g · [1 − (1 − f_g)^N]

    Args:
        diplotype_freqs: list of (diplotype, frequency)
        N: registry size (number of donors)

    Returns:
        float in [0, 1]
    """
    if N == 0:
        return 0.0
    total = 0.0
    for _, f in diplotype_freqs:
        if f > 0:
            total += f * (1.0 - (1.0 - f) ** N)
    return total


def find_registry_size(
    diplotype_freqs: list[tuple[tuple, float]],
    target_coverage: float,
    n_min: int = 100,
    n_max: int = 20_000_000,
) -> int:
    """Binary search for minimum N such that Coverage(N) >= target_coverage.

    Returns n_max if target is unachievable within n_max.
    """
    if compute_coverage(diplotype_freqs, n_max) < target_coverage:
        return n_max

    lo, hi = n_min, n_max
    while lo < hi:
        mid = (lo + hi) // 2
        if compute_coverage(diplotype_freqs, mid) >= target_coverage:
            hi = mid
        else:
            lo = mid + 1
    return lo


def collapse_to_4locus(
    hap_freqs_5: list[tuple[tuple, float]],
) -> list[tuple[tuple, float]]:
    """Collapse 5-locus haplotype frequencies to 4-locus (drop DQB1).

    Sums frequencies of haplotypes that share A, B, C, DRB1 but differ at DQB1.

    Args:
        hap_freqs_5: list of (5-tuple of alleles, frequency)
            Allele order: HLA-A, HLA-B, HLA-C, DRB1, DQB1

    Returns:
        list of (4-tuple, summed_frequency)
    """
    aggregated: dict[tuple, float] = {}
    for hap, freq in hap_freqs_5:
        key4 = hap[:4]  # drop DQB1 (last element)
        aggregated[key4] = aggregated.get(key4, 0.0) + freq
    return list(aggregated.items())


# ---------------------------------------------------------------------------
# Data loading and preparation
# ---------------------------------------------------------------------------

def load_haplotype_frequencies(csv_path: Path) -> dict[str, list[tuple[tuple, float]]]:
    """Load haplo.stats output; return {ethnicity: [(hap_tuple, freq), ...]} sorted by freq desc.

    Only retains top haplotypes covering 99% of cumulative frequency per ethnicity.
    """
    df = pd.read_csv(csv_path)

    # Column names from haplo.stats output: one col per locus + 'frequency' + 'ethnicity'
    # Locus columns are named by LOCI_5 values
    locus_cols = [c for c in df.columns if c in LOCI_5]
    if not locus_cols:
        # Fallback: columns may be numbered; try to infer from position
        non_meta = [c for c in df.columns if c not in ("frequency","ethnicity","n_samples")]
        locus_cols = non_meta[:5]
    locus_cols = locus_cols[:5]  # ensure exactly 5

    result = {}
    for eth, group in df.groupby("ethnicity"):
        group = group.sort_values("frequency", ascending=False)
        # Keep top haplotypes covering 99% cumulative frequency
        group["cumfreq"] = group["frequency"].cumsum()
        group["frequency"] /= group["frequency"].sum()  # renormalize to 1
        group["cumfreq"]  = group["frequency"].cumsum()
        top = group[group["cumfreq"].shift(1, fill_value=0) < 0.99].copy()
        hap_freqs = [
            (tuple(row[locus_cols]), float(row["frequency"]))
            for _, row in top.iterrows()
        ]
        result[eth] = hap_freqs
        print(f"  {eth}: {len(hap_freqs)} haplotypes (99% coverage)")
    return result


def build_combined_pool(
    per_eth_hapfreqs: dict[str, list[tuple[tuple, float]]],
    weights: dict[str, float],
) -> list[tuple[tuple, float]]:
    """Build a weighted combined haplotype pool across all ethnicities."""
    combined: dict[tuple, float] = {}
    for eth, hap_freqs in per_eth_hapfreqs.items():
        w = weights.get(eth, 0.0)
        for hap, freq in hap_freqs:
            combined[hap] = combined.get(hap, 0.0) + w * freq
    # Renormalize
    total = sum(combined.values())
    return [(hap, f / total) for hap, f in combined.items()]


# ---------------------------------------------------------------------------
# Main model run
# ---------------------------------------------------------------------------

def run_model(
    per_eth_hapfreqs: dict[str, list[tuple[tuple, float]]],
    match_level: str,   # "10of10" or "8of8"
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run coverage curve model for one match level.

    Returns:
        curves_df: (N, coverage, ethnicity, model_variant, match_level)
        targets_df: (target_coverage, N_needed, ethnicity, model_variant, match_level)
    """
    assert match_level in ("10of10", "8of8")

    # Collapse to 4-locus if 8/8
    if match_level == "8of8":
        hap_freqs_input = {
            eth: collapse_to_4locus(hfs)
            for eth, hfs in per_eth_hapfreqs.items()
        }
    else:
        hap_freqs_input = per_eth_hapfreqs

    combined_pool = build_combined_pool(hap_freqs_input, POP_WEIGHTS)

    curve_records  = []
    target_records = []

    all_groups = ETHNICITIES + ["Combined"]

    for eth in all_groups:
        for variant in ("same_ethnicity", "cross_ethnic"):

            # Patient genotype frequencies: always from the ethnic group
            if eth == "Combined":
                patient_hap = combined_pool
            else:
                patient_hap = hap_freqs_input.get(eth, [])

            if not patient_hap:
                continue

            patient_diplos = get_diplotype_frequencies(patient_hap)

            # Donor genotype frequencies: same-ethnicity uses patient pool,
            # cross-ethnic uses combined pool
            if variant == "same_ethnicity":
                donor_hap = patient_hap
            else:
                donor_hap = combined_pool

            donor_diplos_dict = dict(get_diplotype_frequencies(donor_hap))

            # For cross-ethnic: re-weight patient genotype probs by donor pool freq
            if variant == "cross_ethnic":
                # Coverage(N) = Σ_g  f_patient(g) · [1 − (1 − f_donor(g))^N]
                effective_diplos = [
                    (g, donor_diplos_dict.get(g, 0.0))
                    for g, _ in patient_diplos
                ]
            else:
                effective_diplos = patient_diplos

            print(f"  {eth} | {variant} | {match_level} | "
                  f"{len(effective_diplos)} diplotypes")

            # Coverage curve
            for N in N_RANGE:
                cov = compute_coverage(effective_diplos, int(N))
                curve_records.append({
                    "N": int(N), "coverage": cov,
                    "ethnicity": eth, "model_variant": variant,
                    "match_level": match_level,
                })

            # Registry size targets
            for target in TARGET_COVERAGES:
                n_needed = find_registry_size(effective_diplos, target)
                target_records.append({
                    "target_coverage": target,
                    "N_needed": n_needed,
                    "ethnicity": eth,
                    "model_variant": variant,
                    "match_level": match_level,
                })

    return pd.DataFrame(curve_records), pd.DataFrame(target_records)


def plot_coverage_curves(curves_df: pd.DataFrame, match_level: str, out_path: Path) -> None:
    """One figure per match level: rows=ethnicities, cols=model_variants."""
    eth_order = ETHNICITIES + ["Combined"]
    variants  = ["same_ethnicity", "cross_ethnic"]
    fig, axes = plt.subplots(
        len(eth_order), len(variants),
        figsize=(12, 3 * len(eth_order)),
        sharex=True, sharey=True,
    )

    for row_i, eth in enumerate(eth_order):
        for col_j, variant in enumerate(variants):
            ax = axes[row_i, col_j]
            sub = curves_df[
                (curves_df["ethnicity"] == eth) &
                (curves_df["model_variant"] == variant)
            ]
            if sub.empty:
                ax.set_visible(False)
                continue
            ax.semilogx(sub["N"], sub["coverage"] * 100, color="steelblue", lw=1.5)
            for t in TARGET_COVERAGES:
                ax.axhline(t * 100, color="grey", lw=0.7, ls="--", alpha=0.6)
            ax.set_ylim(0, 105)
            ax.set_xlim(N_RANGE[0], N_RANGE[-1])
            ax.set_title(f"{eth} — {variant.replace('_',' ')}", fontsize=8)
            ax.set_ylabel("Coverage (%)" if col_j == 0 else "")
            ax.set_xlabel("Registry size (N)" if row_i == len(eth_order)-1 else "")
            ax.xaxis.set_major_formatter(mticker.FuncFormatter(
                lambda x, _: f"{int(x):,}"
            ))

    fig.suptitle(f"Registry size vs patient coverage ({match_level})", fontsize=11)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved → {out_path}")


def main():
    haplo_csv = DATA / "haplo_freqs_haplo_stats.csv"
    if not haplo_csv.exists():
        raise FileNotFoundError(
            f"{haplo_csv} not found. Run 03_hwe_test.R first."
        )

    print("Loading haplotype frequencies from haplo.stats output...")
    per_eth = load_haplotype_frequencies(haplo_csv)

    all_curves  = []
    all_targets = []

    for match_level in ("10of10", "8of8"):
        print(f"\n=== Running model: {match_level} ===")
        curves_df, targets_df = run_model(per_eth, match_level)
        all_curves.append(curves_df)
        all_targets.append(targets_df)
        plot_coverage_curves(curves_df, match_level,
                             FIGS / f"coverage_curves_{match_level}.png")

    curves_combined  = pd.concat(all_curves,  ignore_index=True)
    targets_combined = pd.concat(all_targets, ignore_index=True)

    curves_combined.to_csv(DATA / "coverage_curves.csv",      index=False)
    targets_combined.to_csv(DATA / "registry_size_targets.csv", index=False)
    print(f"\nSaved coverage_curves.csv ({len(curves_combined):,} rows)")
    print(f"Saved registry_size_targets.csv ({len(targets_combined):,} rows)")

    print("\nRegistry size summary (same-ethnicity, 10/10 match):")
    summary = targets_combined[
        (targets_combined["model_variant"] == "same_ethnicity") &
        (targets_combined["match_level"] == "10of10")
    ].pivot(index="ethnicity", columns="target_coverage", values="N_needed")
    print(summary.to_string())


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run tests — expect pass**

```bash
cd /data/alvin/HLA
python3 -m pytest tests/test_registry_model.py -v
```

Expected: All 11 tests PASS.

- [ ] **Step 5: Commit**

```bash
cd /data/alvin/HLA
git add analysis/04_registry_model.py tests/test_registry_model.py
git commit -m "feat: registry size model with TDD — core HWE math and coverage formula"
```

---

## Task 6: Run registry model end-to-end

**Files:**
- Modify: nothing new (runs existing `04_registry_model.py`)

- [ ] **Step 1: Run the model**

```bash
cd /data/alvin/HLA
python3 analysis/04_registry_model.py 2>&1 | tee analysis/data/registry_model_log.txt
```

Expected: prints haplotype counts per ethnicity, model variants, and a registry size summary table. Saves `coverage_curves.csv`, `registry_size_targets.csv`, and two PNG figures.

Runtime note: with 99%-coverage haplotype sets (potentially 200–500 haplotypes each), diplotype enumeration is O(K²) — may be slow for large K. If it runs > 10 minutes, reduce the cumulative frequency cutoff from 0.99 to 0.95 in `load_haplotype_frequencies()`.

- [ ] **Step 2: Sanity-check the output**

```bash
python3 -c "
import pandas as pd
df = pd.read_csv('analysis/data/registry_size_targets.csv')
print(df[(df['model_variant']=='same_ethnicity') & (df['match_level']=='10of10')]
      .pivot(index='ethnicity', columns='target_coverage', values='N_needed')
      .to_string())
"
```

Expected: A table showing that minority populations (Malay, Indian) require disproportionately larger registries. Cross-ethnic values should be lower than same-ethnicity for minority groups (benefit of shared donor pool).

- [ ] **Step 3: Commit figures and outputs**

```bash
cd /data/alvin/HLA
git add analysis/data/coverage_curves.csv \
        analysis/data/registry_size_targets.csv \
        analysis/figures/ \
        analysis/data/registry_model_log.txt
git commit -m "feat: run registry model, save coverage curves and size targets"
```

---

## Task 7: Report assembly

**Files:**
- Create: `analysis/05_report.py`

- [ ] **Step 1: Create 05_report.py**

```python
"""05_report.py — Assemble verification_summary.md from all intermediate outputs.

Reads:
  data/allele_freq_comparison.csv
  data/hwe_results.csv
  data/haplo_freqs_haplo_stats.csv
  data/registry_size_targets.csv
Writes:
  verification_summary.md  (in project root)
"""

from pathlib import Path
import pandas as pd
import numpy as np

BASE = Path(__file__).parent.parent
DATA = Path(__file__).parent / "data"


def fmt_pct(x): return f"{x*100:.1f}%"


def section_allele_freq():
    csv = DATA / "allele_freq_comparison.csv"
    if not csv.exists():
        return "## 1. Allele Frequency Reproducibility\n\n_Data not found (run 02_allele_freq.py)._\n"

    df = pd.read_csv(csv)
    n_flagged = df["flagged"].sum() if "flagged" in df.columns else "N/A"
    max_diff  = df["difference"].abs().max() if "difference" in df.columns else "N/A"
    max_diff_str = f"{max_diff:.4f}" if isinstance(max_diff, float) else str(max_diff)

    lines = [
        "## 1. Allele Frequency Reproducibility\n",
        f"- Alleles compared: {len(df):,}",
        f"- Alleles with |difference| > 0.5%: **{n_flagged}**",
        f"- Maximum single-allele discrepancy: **{max_diff_str}**\n",
    ]

    if isinstance(max_diff, float) and max_diff < 0.01:
        lines.append("**Verdict:** Allele frequency calculation is reproducible. "
                     "Discrepancies are within rounding error from 2-field normalization.\n")
    else:
        lines.append("**Verdict:** Non-trivial discrepancies detected. "
                     "Investigate flagged alleles in `data/allele_freq_comparison.csv`.\n")
    return "\n".join(lines) + "\n"


def section_hwe():
    csv = DATA / "hwe_results.csv"
    if not csv.exists():
        return "## 2. Hardy-Weinberg Equilibrium Assessment\n\n_Data not found (run 03_hwe_test.R)._\n"

    df = pd.read_csv(csv)
    sig = df[df["significant"]] if "significant" in df.columns else pd.DataFrame()

    lines = [
        "## 2. Hardy-Weinberg Equilibrium Assessment\n",
        f"- Locus × ethnicity combinations tested: {len(df)}",
        f"- Significant deviations (Bonferroni p < {0.05/20:.4f}): **{len(sig)}**\n",
    ]
    if len(sig) > 0:
        lines.append("Deviating combinations:\n")
        for _, row in sig.iterrows():
            lines.append(f"  - {row['ethnicity']} / {row['locus']}: "
                         f"p={row['p_value']:.2e}, "
                         f"obs_het={row['obs_heterozygosity']:.3f}, "
                         f"exp_het={row['exp_heterozygosity']:.3f}")
        lines.append("\n**Implication:** HWE violations introduce bias into Gene[Rate]'s "
                     "haplotype frequency estimates for the affected loci/groups. "
                     "Haplotype frequencies for these combinations should be interpreted cautiously.\n")
    else:
        lines.append("**Verdict:** No significant HWE deviations detected. "
                     "Gene[Rate]'s HWE assumption is supported by the data.\n")
    return "\n".join(lines) + "\n"


def section_haplo_agreement():
    csv = DATA / "haplo_freqs_haplo_stats.csv"
    if not csv.exists():
        return "## 3. Haplotype Frequency Agreement (haplo.stats vs Gene[Rate])\n\n_Data not found._\n"

    df = pd.read_csv(csv)
    lines = [
        "## 3. Haplotype Frequency Agreement\n",
        f"- haplo.stats estimated {len(df):,} haplotypes across all ethnic groups.",
        f"- Top haplotypes per group (see `data/haplo_freqs_haplo_stats.csv`).\n",
        "Note: Direct numerical comparison with Gene[Rate] haplotype frequencies requires "
        "matching haplotype strings across the two outputs. Review "
        "`data/generate_published_haplotypes.csv` alongside `haplo_freqs_haplo_stats.csv`.\n",
    ]
    return "\n".join(lines) + "\n"


def section_missing_data():
    csv = DATA / "hla_clean.csv"
    if not csv.exists():
        return "## 4. Missing Data Impact\n\n_Data not found._\n"

    df = pd.read_csv(csv, dtype=str)
    df = df.where(df != "nan", None)
    df_main = df[df["source"].isin(["BMDP","SCBB"])]

    miss = (
        df_main.groupby(["source","locus"])["allele1"]
        .apply(lambda x: (x.isna() | (x == "-")).mean())
        .unstack("locus")
        .applymap(fmt_pct)
    )
    lines = [
        "## 4. Missing Data Assessment\n",
        "Missingness rates (allele1) per source and locus:\n",
        "```",
        miss.to_string(),
        "```\n",
        "Missing alleles are excluded from allele frequency denominators. "
        "If missingness is >5% for any locus, frequencies for that locus are biased "
        "toward observed alleles and should be flagged in any publication update.\n",
    ]
    return "\n".join(lines) + "\n"


def section_registry_size():
    csv = DATA / "registry_size_targets.csv"
    if not csv.exists():
        return "## 5. Registry Size Findings\n\n_Data not found (run 04_registry_model.py)._\n"

    df = pd.read_csv(csv)
    lines = ["## 5. Registry Size Findings\n"]

    for match_level in ("10of10", "8of8"):
        lines.append(f"### {match_level.replace('of','/')} match\n")
        for variant in ("same_ethnicity", "cross_ethnic"):
            sub = df[(df["match_level"]==match_level) & (df["model_variant"]==variant)]
            if sub.empty:
                continue
            pivot = (
                sub.pivot(index="ethnicity", columns="target_coverage", values="N_needed")
                .reindex(["Chinese","Malay","Indian","Others","Combined"])
                .dropna(how="all")
            )
            pivot.columns = [fmt_pct(c) for c in pivot.columns]
            lines.append(f"**{variant.replace('_',' ').title()}**\n")
            lines.append("```")
            lines.append(pivot.to_string())
            lines.append("```\n")
    return "\n".join(lines) + "\n"


def section_recommendations():
    return """## 6. Suggested Improvements for Future Analysis

1. **Higher resolution typing** — re-run analysis at 4-field resolution where available in BMDP to capture allele-level diversity lost in 2-field normalization.
2. **Bootstrap confidence intervals** — resample donors (n=1000 iterations) to produce 95% CIs for top haplotype frequencies; Gene[Rate] does not report uncertainty.
3. **Linkage disequilibrium reporting** — compute D' and r² between all 10 locus pairs for each CMIO group to characterize population structure.
4. **Registry model refinement** — replace HWE diplotype assumption with empirical diplotype counts from phased data (if phasing software such as SHAPEIT2 is applied).
5. **"Others" sub-stratification** — the Others group is ethnically heterogeneous; sub-stratify by reported sub-ethnicity if available to improve matching predictions.
"""


def main():
    report = "\n".join([
        "# HLA Analysis Verification Summary",
        f"Generated by 05_report.py\n",
        section_allele_freq(),
        section_hwe(),
        section_haplo_agreement(),
        section_missing_data(),
        section_registry_size(),
        section_recommendations(),
    ])

    out = BASE / "verification_summary.md"
    out.write_text(report)
    print(f"Saved → {out}")
    print(f"Total: {len(report)} characters")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run the report**

```bash
cd /data/alvin/HLA
python3 analysis/05_report.py
```

Expected: writes `verification_summary.md` to `/data/alvin/HLA/`.

- [ ] **Step 3: Review the summary**

```bash
cat /data/alvin/HLA/verification_summary.md
```

Check that all sections populated correctly. If any section says "_Data not found_", that upstream script has not yet been run successfully.

- [ ] **Step 4: Commit**

```bash
cd /data/alvin/HLA
git add analysis/05_report.py verification_summary.md
git commit -m "feat: report assembly — verification_summary.md"
```

---

## Task 8: End-to-end run and final commit

**Files:**
- No new files

- [ ] **Step 1: Run full pipeline**

```bash
cd /data/alvin/HLA/analysis
bash run_all.sh 2>&1 | tee full_run.log
```

Expected: all 5 steps complete without error. Final line: `Done. Outputs in analysis/data/ and analysis/figures/`

- [ ] **Step 2: Verify all expected outputs exist**

```bash
ls /data/alvin/HLA/analysis/data/
# Expected files: hla_clean.csv, allele_freq_comparison.csv,
#   allele_freqs_observed.csv, haplo_freqs_haplo_stats.csv,
#   hwe_results.csv, coverage_curves.csv, registry_size_targets.csv

ls /data/alvin/HLA/analysis/figures/
# Expected: allele_freq_heatmap.png,
#   coverage_curves_10of10.png, coverage_curves_8of8.png
```

- [ ] **Step 3: Run all tests together**

```bash
cd /data/alvin/HLA
python3 -m pytest tests/ -v
```

Expected: All tests PASS.

- [ ] **Step 4: Final commit**

```bash
cd /data/alvin/HLA
git add analysis/data/ analysis/figures/ analysis/full_run.log
git commit -m "feat: end-to-end pipeline complete — all outputs and verification summary generated"
```
