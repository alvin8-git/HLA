"""02_allele_freq.py — Recompute allele frequencies and compare vs published values.

Reads:  analysis/data/hla_clean.csv, ../BMDPnSCBB.results.xlsx
Writes: analysis/data/allele_freq_comparison.csv
        analysis/data/allele_freqs_observed.csv
        analysis/figures/allele_freq_heatmap.png
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

BASE  = Path(__file__).parent.parent
DATA  = Path(__file__).parent / "data"
FIGS  = Path(__file__).parent / "figures"
FIGS.mkdir(exist_ok=True)

LOCI        = ["HLA-A", "HLA-B", "HLA-C", "DRB1", "DQB1"]
ETHNICITIES = ["Chinese", "Malay", "Indian", "Others"]

LOCUS_COLS = {
    "HLA-A": ("HLA-A", "%"),
    "HLA-B": ("HLA-B", "%.4"),
    "HLA-C": ("HLA-C", "%.8"),
    "DRB1":  ("DRB1",  "%.12"),
    "DQB1":  ("DQB1",  "%.16"),
}

SHEET_MAP = {
    "Chinese": "BMDPnSCBB.Chinese",
    "Malay":   "BMDPnSCBB.Malay",
    "Indian":  "BMDPnSCBB.Indian",
    "Others":  "BMDPnSCBB.Others",
}

MAIN_SOURCES = {"BMDP_OUT", "SCBB_OUT"}


def compute_allele_frequencies(df: pd.DataFrame) -> pd.DataFrame:
    """Compute allele frequencies from tidy hla_clean DataFrame.

    Uses only BMDP_OUT and SCBB_OUT sources (have allele pairs).
    Returns DataFrame: ethnicity, locus, allele, count, total, frequency
    """
    df_main = df[df["source"].isin(MAIN_SOURCES)].copy()

    a1 = df_main[["sample_id","source","ethnicity","locus","allele1"]].rename(
        columns={"allele1": "allele"})
    a2 = df_main[["sample_id","source","ethnicity","locus","allele2"]].rename(
        columns={"allele2": "allele"})
    melted = pd.concat([a1, a2], ignore_index=True)
    melted = melted[melted["allele"].notna() & (melted["allele"].astype(str) != "nan")]

    counts = (
        melted.groupby(["ethnicity","locus","allele"])
        .size().reset_index(name="count")
    )
    totals = (
        melted.groupby(["ethnicity","locus"])
        .size().reset_index(name="total")
    )
    result = counts.merge(totals, on=["ethnicity","locus"])
    result["frequency"] = result["count"] / result["total"]
    return result


def load_published_frequencies(results_path: Path) -> pd.DataFrame:
    """Extract allele frequencies from BMDPnSCBB.results.xlsx wide format.

    Returns DataFrame: ethnicity, locus, allele, pub_frequency
    """
    records = []
    xl = pd.ExcelFile(results_path)

    for ethnicity, sheet_name in SHEET_MAP.items():
        if sheet_name not in xl.sheet_names:
            print(f"  WARNING: sheet '{sheet_name}' not found")
            continue

        df = pd.read_excel(xl, sheet_name=sheet_name, dtype=str)
        df.columns = df.columns.str.strip()

        for locus, (allele_col, freq_col) in LOCUS_COLS.items():
            if allele_col not in df.columns or freq_col not in df.columns:
                print(f"  WARNING: cols '{allele_col}'/'{freq_col}' not found in {sheet_name}")
                print(f"  Available: {list(df.columns)[:10]}")
                continue

            sub = df[[allele_col, freq_col]].copy()
            sub.columns = ["allele", "pub_frequency"]
            sub["allele"] = sub["allele"].astype(str).str.strip()
            sub = sub[sub["allele"].notna() & (sub["allele"] != "") & (sub["allele"] != "nan")]
            sub["pub_frequency"] = pd.to_numeric(sub["pub_frequency"], errors="coerce")
            sub = sub[sub["pub_frequency"].notna()]
            sub["ethnicity"] = ethnicity
            sub["locus"] = locus
            records.append(sub[["ethnicity","locus","allele","pub_frequency"]])

    return pd.concat(records, ignore_index=True) if records else pd.DataFrame()


def compare_frequencies(
    observed: pd.DataFrame,
    published: pd.DataFrame,
    threshold: float = 0.005,
) -> pd.DataFrame:
    """Merge observed and published; flag |difference| > threshold."""
    merged = observed.merge(published, on=["ethnicity","locus","allele"], how="outer")
    merged["difference"] = merged["frequency"] - merged["pub_frequency"]
    merged["flagged"] = merged["difference"].abs() > threshold
    return merged


def plot_heatmap(comparison: pd.DataFrame, out_path: Path) -> None:
    pivot = (
        comparison.groupby(["ethnicity","locus"])["difference"]
        .apply(lambda x: x.abs().max())
        .unstack("locus")
        .reindex(index=ETHNICITIES, columns=LOCI)
    )
    fig, ax = plt.subplots(figsize=(9, 4))
    sns.heatmap(
        pivot.astype(float), annot=True, fmt=".4f", cmap="YlOrRd",
        linewidths=0.5, ax=ax,
        cbar_kws={"label": "Max |observed − published|"},
    )
    ax.set_title("Allele frequency discrepancy: observed vs published (Gene[Rate])")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  Heatmap → {out_path}")


def main():
    print("Loading hla_clean.csv...")
    df = pd.read_csv(DATA / "hla_clean.csv", dtype=str).where(lambda x: x != "nan", np.nan)

    print("Computing allele frequencies (BMDP_OUT + SCBB_OUT)...")
    observed = compute_allele_frequencies(df)
    print(f"  {len(observed):,} rows computed")
    observed.to_csv(DATA / "allele_freqs_observed.csv", index=False)

    print("Loading published frequencies...")
    published = load_published_frequencies(BASE / "BMDPnSCBB.results.xlsx")
    if published.empty:
        print("WARNING: No published frequencies loaded.")
        return
    print(f"  {len(published):,} published rows loaded")

    comparison = compare_frequencies(observed, published)
    n_flagged = int(comparison["flagged"].sum())
    max_diff  = comparison["difference"].abs().max()
    print(f"  Flagged (|diff|>0.5%): {n_flagged}  |  Max diff: {max_diff:.4f}")

    comparison.to_csv(DATA / "allele_freq_comparison.csv", index=False)
    plot_heatmap(comparison, FIGS / "allele_freq_heatmap.png")


if __name__ == "__main__":
    main()
