######################################################################
## Sweep the panel size from 1 to ALL genes and see which performs   ##
## best. Genes are ordered ONCE by the composite stability+mRMR       ##
## ranking (from gene_panel_selection.py); then for each k we take    ##
## the top-k genes and measure:                                       ##
##   - cross-cohort balanced accuracy & AUC (CC->HNC and HNC->CC)     ##
##   - within-cohort 5-fold CV balanced accuracy                      ##
##                                                                    ##
## "Best" is judged by mean cross-cohort balanced accuracy (the       ##
## average of both train/test directions) - the metric that reflects ##
## a panel transferring across both cancers.                          ##
######################################################################

import os
import numpy as np
import pandas as pd

import gene_panel_selection as gps  # reuse load / ranking / validation

# Sweep spans 1..all genes, so it always ranks the FULL gene set
# (regardless of gps.REQUIRE_CONSISTENT, which only governs the fixed panel).


def main():
    out_dir = os.path.join(gps.BASE_DIR, "Panel_Selection_Results")
    os.makedirs(out_dir, exist_ok=True)

    data = {c: gps.load(p) for c, p in gps.DATASETS.items()}
    genes = data["CC"][2]
    candidates = list(genes)

    print(f"Sweeping panel size 1..{len(candidates)} over the full gene set\n")
    tbl = gps.composite_ranking(data, candidates)
    ordered = tbl["Gene"].tolist()

    rows = []
    for k in range(1, len(ordered) + 1):
        panel = ordered[:k]
        val = gps.cross_cohort(f"top{k}", panel, data)
        cc2hnc = next(v for v in val if v["Train"] == "CC")
        hnc2cc = next(v for v in val if v["Train"] == "HNC")
        wc = gps.within_cohort_cv(panel, data)
        cv_mean = float(np.mean([w["CV_Balanced_Acc"] for w in wc]))

        rows.append({
            "k": k,
            "Added_Gene": panel[-1],
            "CC2HNC_BalAcc": cc2hnc["Balanced_Acc"],
            "HNC2CC_BalAcc": hnc2cc["Balanced_Acc"],
            "Mean_CrossCohort_BalAcc": (cc2hnc["Balanced_Acc"] + hnc2cc["Balanced_Acc"]) / 2,
            "Mean_CrossCohort_AUC": (cc2hnc["AUC"] + hnc2cc["AUC"]) / 2,
            "WithinCohort_CV_BalAcc": cv_mean,
            "Panel_Genes": ", ".join(panel),
        })

    sweep = pd.DataFrame(rows)
    sweep.to_csv(os.path.join(out_dir, "panel_size_sweep.csv"), index=False)

    # ---- best & parsimonious choices ----
    primary = "Mean_CrossCohort_BalAcc"
    best_k = int(sweep.loc[sweep[primary].idxmax(), "k"])
    best_val = sweep[primary].max()
    # smallest k within 1% (0.01) of the best - the parsimonious panel
    within = sweep[sweep[primary] >= best_val - 0.01]
    parsimonious_k = int(within["k"].min())

    print("\n" + "=" * 84)
    print("PANEL-SIZE SWEEP (ranked genes added one at a time)")
    print("=" * 84)
    show = sweep[[
        "k", "Added_Gene", "CC2HNC_BalAcc", "HNC2CC_BalAcc",
        "Mean_CrossCohort_BalAcc", "Mean_CrossCohort_AUC", "WithinCohort_CV_BalAcc",
    ]].copy()
    with pd.option_context("display.width", 200, "display.max_rows", None):
        print(show.round(4).to_string(index=False))

    print("\n" + "=" * 84)
    print(f"BEST panel size by mean cross-cohort balanced accuracy: k={best_k} "
          f"({best_val:.4f})")
    print(f"PARSIMONIOUS choice (smallest k within 0.01 of best): k={parsimonious_k}")
    print("=" * 84)
    for k in sorted({best_k, parsimonious_k}):
        gset = sweep.loc[sweep["k"] == k, "Panel_Genes"].iloc[0]
        print(f"  k={k}: {gset}")
    print(f"\nFull sweep saved to: Panel_Selection_Results/panel_size_sweep.csv")


if __name__ == "__main__":
    main()
