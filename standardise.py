"""
Standardise the Run2 CC (CESC) and HNC datasets into a uniform format.

- Unifies column names: sample id -> 'Sample_ID', label -> 'Label'
- Unifies label values: positive / HPV+ -> 'HPV_positive'
                        negative / HPV- -> 'HPV_negative'
- Saves CC_data.csv and HNC_data.csv
- Reports the per-gene mean and SD to check whether the gene values are
  standardised (z-scored: mean ~0, SD ~1).
"""

import pandas as pd

CC_IN = "CC CESC_HPV_VST_52genes (1).csv"
HNC_IN = "HNC VST_genes_of_interest_HPVpos_vs_HPVneg.xlsx"

# Map every raw label spelling onto one canonical value
LABEL_MAP = {
    "positive": "HPV_positive",
    "HPV+": "HPV_positive",
    "negative": "HPV_negative",
    "HPV-": "HPV_negative",
}


def standardise(df):
    """Rename the first two columns and normalise the label values."""
    df = df.rename(columns={df.columns[0]: "Sample_ID", df.columns[1]: "Label"})
    df["Label"] = df["Label"].map(LABEL_MAP)
    if df["Label"].isna().any():
        raise ValueError("Unmapped label value found")
    return df


def summarise(datasets):
    """Print a summary table: samples, genes and per-class sample counts."""
    classes = sorted({c for df in datasets.values() for c in df["Label"].unique()})
    header = f"{'Dataset':<10}{'Samples':>9}{'Genes':>7}" + "".join(f"{c:>16}" for c in classes)
    print("\n" + "=" * len(header))
    print("DATASET SUMMARY")
    print("=" * len(header))
    print(header)
    print("-" * len(header))
    for name, df in datasets.items():
        n_samples = len(df)
        n_genes = df.shape[1] - 2  # exclude Sample_ID and Label
        counts = df["Label"].value_counts()
        row = f"{name:<10}{n_samples:>9}{n_genes:>7}"
        for c in classes:
            n = int(counts.get(c, 0))
            pct = 100 * n / n_samples if n_samples else 0
            row += f"{f'{n} ({pct:.1f}%)':>16}"
        print(row)
    print("=" * len(header))


def check_standardisation(df, name):
    """Compute the per-gene mean and SD across the gene columns."""
    genes = df.columns[2:]
    stats = pd.DataFrame({
        "mean": df[genes].mean(),
        "sd": df[genes].std(),
    })
    print(f"\n=== {name}: per-gene mean / SD ({len(genes)} genes) ===")
    print(stats.round(3).to_string())
    print(f"Overall mean of gene means: {stats['mean'].mean():.4f}")
    print(f"Overall mean of gene SDs:   {stats['sd'].mean():.4f}")
    standardised = stats["mean"].abs().max() < 0.01 and abs(stats["sd"].mean() - 1) < 0.05
    print(f"Values appear z-standardised (mean~0, SD~1): {standardised}")
    return stats


def main():
    cc = standardise(pd.read_csv(CC_IN))
    hnc = standardise(pd.read_excel(HNC_IN))

    cc.to_csv("CC_data.csv", index=False)
    hnc.to_csv("HNC_data.csv", index=False)
    print(f"Saved CC_data.csv  {cc.shape}")
    print(f"Saved HNC_data.csv {hnc.shape}")

    summarise({"CC": cc, "HNC": hnc})

    check_standardisation(cc, "CC")
    check_standardisation(hnc, "HNC")


if __name__ == "__main__":
    main()
