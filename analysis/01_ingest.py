"""01_ingest.py — Load HLA data from cleaned Excel + HSA txt files.

Output: analysis/data/hla_clean.csv
Columns: sample_id, source, ethnicity, locus, allele1, allele2

For Excel (BMDP/SCBB): allele1 and allele2 are both present.
For txt (HSA): allele1 is set, allele2 is NaN (single haplotype per row).
"""

import re
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

# Matches: A1, A2, B1, B2, C1, C2, DRB1_1, DRB1_2, DRB11, DRB12, DQB11, DQB12
# Also handles HLA-A1, HLA-B2 prefixed variants
_ALLELE_PATTERN = re.compile(
    r'^(?:HLA[-_]?)?(A|B|C|DRB1|DQB1)[_\s-]?([12])$', re.IGNORECASE
)
_LOCUS_NAME = {"A": "HLA-A", "B": "HLA-B", "C": "HLA-C", "DRB1": "DRB1", "DQB1": "DQB1"}


def normalize_allele(val) -> str:
    """Return 2-field allele string (e.g. '11:01') or NaN if missing/unparseable."""
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return np.nan
    s = str(val).strip()
    if s in ("", "-", "0", "NA", "na", "None", "nan"):
        return np.nan
    s = re.sub(r'[GP]$', '', s)
    parts = s.split(":")
    if len(parts) >= 2:
        return f"{parts[0]}:{parts[1]}"
    return np.nan


def map_ethnicity(val) -> str:
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
    """Load HSA txt files: ethnicity + 5 HLA loci, one haplotype per row."""
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
        label = sheet.strip().upper().replace(" ", "_").replace(".", "_")
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
