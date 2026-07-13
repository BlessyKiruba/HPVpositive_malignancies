######################################################################
## Plot the panel-size sweep: k vs cross-cohort balanced accuracy,   ##
## cross-cohort AUC, and within-cohort CV balanced accuracy.         ##
## Marks the best balanced-accuracy size and the AUC-plateau elbow.  ##
## Reads Panel_Selection_Results/panel_size_sweep.csv.               ##
######################################################################

import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # headless / file output
import matplotlib.pyplot as plt

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
RES_DIR = os.path.join(BASE_DIR, "Panel_Selection_Results")
AUC_TOL = 0.005  # elbow = smallest k within this of the peak AUC


def main():
    df = pd.read_csv(os.path.join(RES_DIR, "panel_size_sweep.csv"))
    k = df["k"]

    # Best balanced-accuracy size, and the AUC-plateau elbow
    best_bal_k = int(df.loc[df["Mean_CrossCohort_BalAcc"].idxmax(), "k"])
    best_bal = df["Mean_CrossCohort_BalAcc"].max()
    auc_peak = df["Mean_CrossCohort_AUC"].max()
    elbow_k = int(df.loc[df["Mean_CrossCohort_AUC"] >= auc_peak - AUC_TOL, "k"].min())

    fig, ax = plt.subplots(figsize=(9.5, 5.8))

    ax.plot(k, df["Mean_CrossCohort_AUC"], "-o", ms=3.5, lw=1.8,
            color="#1f77b4", label="Cross-cohort AUC")
    ax.plot(k, df["Mean_CrossCohort_BalAcc"], "-s", ms=3.5, lw=1.8,
            color="#d62728", label="Cross-cohort balanced accuracy")
    ax.plot(k, df["WithinCohort_CV_BalAcc"], "-^", ms=3.5, lw=1.8,
            color="#2ca02c", label="Within-cohort 5-fold CV balanced accuracy")

    ax.set_xlabel("Panel size (number of genes, k)", fontsize=11)
    ax.set_ylabel("Performance", fontsize=11)
    ax.set_title("Performance vs Panel Size", fontsize=12, fontweight="bold")
    ax.set_xlim(0.5, len(df) + 0.5)
    ax.set_ylim(0.70, 1.0)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="lower right", fontsize=9, framealpha=0.95)
    fig.tight_layout()

    png = os.path.join(RES_DIR, "panel_size_sweep.png")
    pdf = os.path.join(RES_DIR, "panel_size_sweep.pdf")
    fig.savefig(png, dpi=300)
    fig.savefig(pdf)  # vector for manuscripts
    print(f"Saved:\n  {os.path.relpath(png, BASE_DIR)}\n  {os.path.relpath(pdf, BASE_DIR)}")
    print(f"Best balanced-accuracy size: k={best_bal_k} ({best_bal:.4f})")
    print(f"AUC plateau elbow: k={elbow_k} (peak AUC {auc_peak:.4f})")


if __name__ == "__main__":
    main()
