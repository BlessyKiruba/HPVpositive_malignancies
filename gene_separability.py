######################################################################
## Model-free per-gene separability of HPV-positive vs HPV-negative. ##
##                                                                  ##
## For each cohort (CC, HNC) and each gene, computes:               ##
##   - ROC AUC (oriented to HPV_positive as the positive class)     ##
##   - Mann-Whitney / Wilcoxon rank-sum p-value, then BH-FDR        ##
##   - log2 fold-change (VST values are already on a log2 scale,    ##
##     so mean(HPV+) - mean(HPV-) is the log2 fold-change)          ##
##                                                                  ##
## Then requires CONSISTENCY across both cohorts (same direction +  ##
## significant in both) and ranks by the weaker cohort's            ##
## separability. This is the primary, model-free result.            ##
######################################################################

import os
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from sklearn.metrics import roc_auc_score

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
POS_LABEL = "HPV_positive"
NEG_LABEL = "HPV_negative"
FDR_ALPHA = 0.05

DATASETS = {
    "CC": os.path.join(BASE_DIR, "CC_data.csv"),
    "HNC": os.path.join(BASE_DIR, "HNC_data.csv"),
}


def bh_fdr(pvals):
    """Benjamini-Hochberg FDR correction (no statsmodels dependency)."""
    p = np.asarray(pvals, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    # enforce monotonicity from the largest p-value downward
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty(n)
    q[order] = np.clip(ranked, 0, 1)
    return q


def analyse_cohort(df, cohort):
    """Per-gene AUC, Wilcoxon p (+FDR) and log2FC for one cohort."""
    genes = [c for c in df.columns if c not in ("Sample_ID", "Label")]
    y = (df["Label"] == POS_LABEL).astype(int).values
    n_pos, n_neg = int(y.sum()), int((1 - y).sum())
    print(f"{cohort}: {len(df)} samples ({n_pos} {POS_LABEL}, {n_neg} {NEG_LABEL}), {len(genes)} genes")

    rows = []
    for g in genes:
        x = df[g].values.astype(float)
        pos, neg = x[y == 1], x[y == 0]

        # AUC oriented so HPV_positive is the positive class
        auc = roc_auc_score(y, x)
        # Mann-Whitney U (Wilcoxon rank-sum)
        try:
            _, p = mannwhitneyu(pos, neg, alternative="two-sided")
        except ValueError:
            p = 1.0  # e.g. all values identical
        log2fc = pos.mean() - neg.mean()  # VST is log2 scale

        rows.append({
            "Gene": g,
            "AUC": auc,
            "Separability": max(auc, 1 - auc),  # direction-agnostic 0.5..1
            "log2FC": log2fc,
            "Direction": "Up in HPV+" if log2fc > 0 else "Down in HPV+",
            "p_value": p,
        })

    res = pd.DataFrame(rows)
    res["FDR"] = bh_fdr(res["p_value"].values)
    res["Significant"] = res["FDR"] < FDR_ALPHA
    return res.sort_values("Separability", ascending=False).reset_index(drop=True)


def main():
    out_dir = os.path.join(BASE_DIR, "Separability_Results")
    os.makedirs(out_dir, exist_ok=True)

    per_cohort = {}
    for cohort, path in DATASETS.items():
        res = analyse_cohort(pd.read_csv(path), cohort)
        res.to_csv(os.path.join(out_dir, f"separability_{cohort}.csv"), index=False)
        per_cohort[cohort] = res

    cc, hnc = per_cohort["CC"], per_cohort["HNC"]

    # Merge the two cohorts on gene and assess consistency
    merged = cc.merge(hnc, on="Gene", suffixes=("_CC", "_HNC"))
    merged["Direction_Consistent"] = np.sign(merged["log2FC_CC"]) == np.sign(merged["log2FC_HNC"])
    merged["Significant_Both"] = merged["Significant_CC"] & merged["Significant_HNC"]
    merged["Consistent"] = merged["Direction_Consistent"] & merged["Significant_Both"]
    # rank by the weaker cohort: a consistent gene is only as good as its worst cohort
    merged["Min_Separability"] = merged[["Separability_CC", "Separability_HNC"]].min(axis=1)
    merged["Mean_Separability"] = merged[["Separability_CC", "Separability_HNC"]].mean(axis=1)

    merged = merged.sort_values(
        ["Consistent", "Min_Separability"], ascending=[False, False]
    ).reset_index(drop=True)
    merged["Rank"] = np.arange(1, len(merged) + 1)

    cols = [
        "Rank", "Gene", "Consistent", "Direction_Consistent", "Significant_Both",
        "Direction_CC", "Direction_HNC",
        "AUC_CC", "AUC_HNC", "Separability_CC", "Separability_HNC",
        "Min_Separability", "Mean_Separability",
        "log2FC_CC", "log2FC_HNC", "FDR_CC", "FDR_HNC",
    ]
    merged = merged[cols]
    merged.to_csv(os.path.join(out_dir, "separability_combined.csv"), index=False)

    n_consistent = int(merged["Consistent"].sum())
    print("\n" + "=" * 100)
    print(f"CONSISTENT GENES (significant FDR<{FDR_ALPHA} AND same direction in BOTH cohorts): {n_consistent}/{len(merged)}")
    print("=" * 100)
    show = merged.head(max(n_consistent, 15)).copy()
    with pd.option_context("display.width", 200, "display.max_columns", None):
        print(show[[
            "Rank", "Gene", "Consistent", "Direction_CC",
            "AUC_CC", "AUC_HNC", "Min_Separability",
            "log2FC_CC", "log2FC_HNC", "FDR_CC", "FDR_HNC",
        ]].round(4).to_string(index=False))

    print(f"\nResults saved to: {os.path.relpath(out_dir, BASE_DIR)}/")
    print("  - separability_CC.csv, separability_HNC.csv (per-cohort, all genes)")
    print("  - separability_combined.csv (merged + consistency, ranked)")


if __name__ == "__main__":
    main()
