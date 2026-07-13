######################################################################
## Reduce the 52-gene panel to a compact (5-6 gene) HPV+/- signature ##
## using MULTIVARIATE feature selection (not univariate screening).  ##
##                                                                    ##
## Pipeline:                                                          ##
##   1. Elastic-Net logistic regression + stability selection        ##
##      (bootstrapped) per cohort  -> selection frequency per gene    ##
##   2. mRMR (min-redundancy max-relevance) per cohort  -> ranking    ##
##   3. Composite cross-cohort rank  -> final 5-6 gene panel          ##
##   4. Cross-cohort validation (train CC -> test HNC, and reverse)   ##
##   5. Sanity-check vs the model-free separability ranking           ##
##                                                                    ##
## Rationale: the top univariate genes are highly correlated          ##
## (proliferation module); a panel needs genes that are predictive,   ##
## NON-redundant, STABLE under the class imbalance, and that          ##
## TRANSFER across both cancers.                                      ##
######################################################################

import os
import warnings
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn.feature_selection import mutual_info_classif
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.metrics import balanced_accuracy_score, roc_auc_score, f1_score
from sklearn.utils import resample
from sklearn.exceptions import ConvergenceWarning

warnings.filterwarnings("ignore", category=ConvergenceWarning)
warnings.filterwarnings("ignore", category=FutureWarning)  # sklearn 1.8 penalty/l1_ratio deprecation

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
POS_LABEL = "HPV_positive"
SEED = 19
N_BOOTSTRAP = 200
L1_RATIO = 0.5          # elastic-net mix (0=ridge, 1=lasso)
PANEL_SIZE = 6
COEF_EPS = 1e-6         # |coef| above this counts as "selected"

# If True, restrict the candidate pool to genes that are cross-cohort
# consistent (FDR-significant AND same direction in BOTH cohorts) per
# gene_separability.py. Drops individually-weak / cohort-specific genes
# (e.g. SBSPON) before multivariate selection runs.
REQUIRE_CONSISTENT = True

DATASETS = {
    "CC": os.path.join(BASE_DIR, "CC_data.csv"),
    "HNC": os.path.join(BASE_DIR, "HNC_data.csv"),
}


def load(path):
    df = pd.read_csv(path)
    genes = [c for c in df.columns if c not in ("Sample_ID", "Label")]
    X = df[genes]
    y = (df["Label"] == POS_LABEL).astype(int).values
    return X, y, genes


def consistent_candidates(all_genes):
    """Genes flagged Consistent by gene_separability.py; falls back to all."""
    path = os.path.join(BASE_DIR, "Separability_Results", "separability_combined.csv")
    if not os.path.exists(path):
        print("  WARNING: separability_combined.csv not found "
              "(run gene_separability.py first). Using ALL genes.")
        return list(all_genes)
    sep = pd.read_csv(path)
    is_consistent = sep["Consistent"].astype(str).str.lower() == "true"
    consistent = [g for g in sep.loc[is_consistent, "Gene"] if g in all_genes]
    return consistent


# ------------------------------------------------------------------ #
# 1. Elastic-Net stability selection                                  #
# ------------------------------------------------------------------ #
def stability_selection(X, y, cohort):
    """Bootstrap an elastic-net logit; return per-gene selection frequency."""
    # Choose C once via CV optimising balanced accuracy (class-weighted)
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
    selector = make_pipeline(
        StandardScaler(),
        LogisticRegressionCV(
            penalty="elasticnet", solver="saga", l1_ratios=[L1_RATIO],
            Cs=10, cv=cv, class_weight="balanced",
            scoring="balanced_accuracy", max_iter=5000, random_state=SEED, n_jobs=-1,
        ),
    )
    selector.fit(X, y)
    best_C = float(selector[-1].C_[0])

    counts = np.zeros(X.shape[1])
    for b in range(N_BOOTSTRAP):
        Xb, yb = resample(X, y, replace=True, n_samples=len(y),
                          stratify=y, random_state=SEED + b)
        clf = make_pipeline(
            StandardScaler(),
            LogisticRegression(
                penalty="elasticnet", solver="saga", l1_ratio=L1_RATIO,
                C=best_C, class_weight="balanced", max_iter=5000,
                random_state=SEED + b,
            ),
        )
        clf.fit(Xb, yb)
        counts += np.abs(clf[-1].coef_[0]) > COEF_EPS

    freq = pd.Series(counts / N_BOOTSTRAP, index=X.columns)
    print(f"  {cohort}: best C={best_C:.4g}, "
          f"{int((freq > 0.5).sum())} genes selected in >50% of bootstraps")
    return freq


# ------------------------------------------------------------------ #
# 2. mRMR (minimum redundancy, maximum relevance)                     #
# ------------------------------------------------------------------ #
def mrmr_rank(X, y):
    """Return genes ordered by mRMR (MID: relevance - mean redundancy)."""
    relevance = pd.Series(
        mutual_info_classif(X, y, random_state=SEED), index=X.columns
    )
    corr = X.corr().abs()

    ordered = [relevance.idxmax()]
    remaining = [g for g in X.columns if g not in ordered]
    while remaining:
        best_gene, best_score = None, -np.inf
        for g in remaining:
            redundancy = corr.loc[g, ordered].mean()
            score = relevance[g] - redundancy
            if score > best_score:
                best_gene, best_score = g, score
        ordered.append(best_gene)
        remaining.remove(best_gene)
    return {g: i + 1 for i, g in enumerate(ordered)}  # gene -> rank (1=best)


# ------------------------------------------------------------------ #
# 4. Cross-cohort validation                                          #
# ------------------------------------------------------------------ #
def cross_cohort(name, genes, data):
    rows = []
    pairs = [("CC", "HNC"), ("HNC", "CC")]
    for tr, te in pairs:
        Xtr, ytr, _ = data[tr]
        Xte, yte, _ = data[te]
        clf = make_pipeline(
            StandardScaler(),
            LogisticRegression(class_weight="balanced", max_iter=5000, random_state=SEED),
        )
        clf.fit(Xtr[genes], ytr)
        yp = clf.predict(Xte[genes])
        proba = clf.predict_proba(Xte[genes])[:, 1]
        rows.append({
            "Panel": name, "N_genes": len(genes), "Train": tr, "Test": te,
            "Balanced_Acc": balanced_accuracy_score(yte, yp),
            "AUC": roc_auc_score(yte, proba),
            "F1_Macro": f1_score(yte, yp, average="macro", zero_division=0),
        })
    return rows


def within_cohort_cv(genes, data):
    rows = []
    for cohort, (X, y, _) in data.items():
        clf = make_pipeline(
            StandardScaler(),
            LogisticRegression(class_weight="balanced", max_iter=5000, random_state=SEED),
        )
        cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
        scores = cross_val_score(clf, X[genes], y, cv=cv, scoring="balanced_accuracy")
        rows.append({"Cohort": cohort, "CV_Balanced_Acc": scores.mean(), "CV_Std": scores.std()})
    return rows


def get_candidates(genes):
    """Candidate pool for selection (optionally consistency-filtered)."""
    if REQUIRE_CONSISTENT:
        candidates = consistent_candidates(genes)
        print(f"REQUIRE_CONSISTENT=True -> candidate pool restricted to "
              f"{len(candidates)}/{len(genes)} cross-cohort-consistent genes")
    else:
        candidates = list(genes)
        print(f"REQUIRE_CONSISTENT=False -> using all {len(genes)} genes as candidates")
    return candidates


def composite_ranking(data, candidates):
    """Stability selection + mRMR per cohort, fused into one Borda ranking."""
    # ---- 1. stability selection per cohort (on candidate pool) ----
    print("Elastic-Net stability selection:")
    freq = {c: stability_selection(X[candidates], y, c) for c, (X, y, _) in data.items()}

    # ---- 2. mRMR per cohort (on candidate pool) ----
    print("mRMR ranking per cohort...")
    mrmr = {c: mrmr_rank(X[candidates], y) for c, (X, y, _) in data.items()}

    # ---- 3. composite cross-cohort ranking ----
    tbl = pd.DataFrame(index=candidates)
    tbl["StabFreq_CC"] = freq["CC"]
    tbl["StabFreq_HNC"] = freq["HNC"]
    tbl["Min_StabFreq"] = tbl[["StabFreq_CC", "StabFreq_HNC"]].min(axis=1)
    tbl["mRMR_CC"] = pd.Series(mrmr["CC"])
    tbl["mRMR_HNC"] = pd.Series(mrmr["HNC"])
    # Borda-style: average of 4 within-signal ranks (lower = better)
    r_stab_cc = tbl["StabFreq_CC"].rank(ascending=False, method="average")
    r_stab_hnc = tbl["StabFreq_HNC"].rank(ascending=False, method="average")
    r_mrmr_cc = tbl["mRMR_CC"].rank(ascending=True, method="average")
    r_mrmr_hnc = tbl["mRMR_HNC"].rank(ascending=True, method="average")
    tbl["Mean_Rank"] = (r_stab_cc + r_stab_hnc + r_mrmr_cc + r_mrmr_hnc) / 4
    tbl = tbl.sort_values(["Mean_Rank"]).reset_index().rename(columns={"index": "Gene"})
    tbl["Final_Rank"] = np.arange(1, len(tbl) + 1)
    return tbl


def main():
    out_dir = os.path.join(BASE_DIR, "Panel_Selection_Results")
    os.makedirs(out_dir, exist_ok=True)

    data = {c: load(p) for c, p in DATASETS.items()}
    genes = data["CC"][2]  # full 52-gene list (used for the baseline)

    candidates = get_candidates(genes)
    tbl = composite_ranking(data, candidates)
    tbl.to_csv(os.path.join(out_dir, "panel_ranking.csv"), index=False)

    panel = tbl["Gene"].head(PANEL_SIZE).tolist()

    print("\n" + "=" * 92)
    print(f"COMPOSITE GENE RANKING (top {PANEL_SIZE} = selected panel)")
    print("=" * 92)
    with pd.option_context("display.width", 200, "display.max_columns", None):
        print(tbl.head(15)[[
            "Final_Rank", "Gene", "StabFreq_CC", "StabFreq_HNC",
            "Min_StabFreq", "mRMR_CC", "mRMR_HNC", "Mean_Rank",
        ]].round(3).to_string(index=False))

    print(f"\n>>> SELECTED PANEL ({PANEL_SIZE} genes): {', '.join(panel)}\n")

    # ---- 4. validation ----
    val = cross_cohort("Selected_Panel", panel, data)
    val += cross_cohort("All_52_genes", genes, data)
    val_df = pd.DataFrame(val)
    val_df.to_csv(os.path.join(out_dir, "cross_cohort_validation.csv"), index=False)

    print("=" * 92)
    print("CROSS-COHORT VALIDATION (train on one cancer, test on the other)")
    print("=" * 92)
    print(val_df.round(4).to_string(index=False))

    wc = pd.DataFrame(within_cohort_cv(panel, data))
    print("\nWithin-cohort 5-fold CV of the selected panel:")
    print(wc.round(4).to_string(index=False))

    # ---- 5. sanity-check vs separability ----
    sep_path = os.path.join(BASE_DIR, "Separability_Results", "separability_combined.csv")
    if os.path.exists(sep_path):
        sep = pd.read_csv(sep_path)[["Gene", "Rank", "AUC_CC", "AUC_HNC", "Direction_CC"]]
        sep = sep.rename(columns={"Rank": "Separability_Rank"})
        check = pd.DataFrame({"Gene": panel}).merge(sep, on="Gene", how="left")
        print("\nPanel genes vs model-free separability ranking:")
        print(check.round(4).to_string(index=False))
    else:
        print("\n(Run gene_separability.py first to cross-check against separability ranking.)")

    print(f"\nResults saved to: {os.path.relpath(out_dir, BASE_DIR)}/")
    print("  - panel_ranking.csv (all genes, composite ranked)")
    print("  - cross_cohort_validation.csv (panel vs all-52 baseline)")


if __name__ == "__main__":
    main()
