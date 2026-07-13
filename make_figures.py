######################################################################
## Manuscript figures for the 6-gene HPV+/- panel.                   ##
##   1. Cross-cohort ROC + Precision-Recall curves                   ##
##   2. Per-gene expression boxplots (HPV+ vs - in each cohort)      ##
##   3. Bootstrap selection-frequency bars                           ##
##   4. Confusion matrices (consolidated): within-cohort 5-fold CV,  ##
##      cross-cohort transfer, and external GEO cohort               ##
## Outputs (png 300dpi + pdf) -> Manuscript_Results/figures/         ##
######################################################################

import os
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import (roc_curve, auc, precision_recall_curve,
                             average_precision_score, confusion_matrix,
                             balanced_accuracy_score)

import gene_panel_selection as gps

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(BASE_DIR, "Manuscript_Results", "figures")
GEO_DIR = os.path.join(BASE_DIR, "GEO_test")
POS, NEG = "HPV_positive", "HPV_negative"
POS_C, NEG_C = "#d62728", "#1f77b4"          # HPV+ / HPV-
COH_C = {"CC": "#9467bd", "HNC": "#8c564b"}  # cohort colours

os.makedirs(OUT, exist_ok=True)


def save(fig, name):
    fig.savefig(os.path.join(OUT, name + ".png"), dpi=300, bbox_inches="tight")
    fig.savefig(os.path.join(OUT, name + ".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(f"  {name}.png / .pdf")


def fit_predict(train, test, genes):
    ytr = (train["Label"] == POS).astype(int)
    yte = (test["Label"] == POS).astype(int)
    clf = make_pipeline(StandardScaler(),
                        LogisticRegression(class_weight="balanced",
                                           max_iter=5000, random_state=gps.SEED))
    clf.fit(train[genes], ytr)
    return yte.values, clf.predict_proba(test[genes])[:, 1]


# ------------------------------------------------------------------ #
def fig_roc_pr(data, panel):
    fig, (axr, axp) = plt.subplots(1, 2, figsize=(11, 5))
    for (tr, te), color in [(("CC", "HNC"), "#1f77b4"), (("HNC", "CC"), "#ff7f0e")]:
        y, p = fit_predict(data[tr], data[te], panel)
        fpr, tpr, _ = roc_curve(y, p)
        axr.plot(fpr, tpr, color=color, lw=2,
                 label=f"{tr}→{te}  (AUC={auc(fpr, tpr):.3f})")
        prec, rec, _ = precision_recall_curve(y, p)
        ap = average_precision_score(y, p)
        axp.plot(rec, prec, color=color, lw=2, label=f"{tr}→{te}  (AP={ap:.3f})")
        axp.axhline(y.mean(), color=color, ls=":", lw=1)  # prevalence baseline

    axr.plot([0, 1], [0, 1], "k--", lw=0.8)
    axr.set(xlabel="False positive rate", ylabel="True positive rate",
            title="Cross-cohort ROC", xlim=(0, 1), ylim=(0, 1.02))
    axr.legend(loc="lower right", fontsize=9)
    axp.set(xlabel="Recall", ylabel="Precision",
            title="Cross-cohort Precision-Recall", xlim=(0, 1), ylim=(0, 1.02))
    axp.legend(loc="lower left", fontsize=9)
    fig.tight_layout()
    save(fig, "fig1_roc_pr")


def fig_boxplots(data, panel):
    fig, axes = plt.subplots(2, 3, figsize=(12, 7))
    rng = np.random.RandomState(gps.SEED)
    for ax, gene in zip(axes.flat, panel):
        groups, colors = [], []
        for c in ("CC", "HNC"):
            for lab in (NEG, POS):
                groups.append(data[c].loc[data[c]["Label"] == lab, gene].values)
                colors.append(NEG_C if lab == NEG else POS_C)
        pos = [1, 2, 3.6, 4.6]
        bp = ax.boxplot(groups, positions=pos, widths=0.7, patch_artist=True,
                        showfliers=False, medianprops=dict(color="black"))
        for patch, c in zip(bp["boxes"], colors):
            patch.set_facecolor(c); patch.set_alpha(0.55)
        for x, g, c in zip(pos, groups, colors):  # jittered points
            ax.scatter(x + rng.uniform(-0.18, 0.18, len(g)), g, s=6,
                       color=c, alpha=0.5, edgecolors="none", zorder=3)
        ax.set_xticks([1.5, 4.1]); ax.set_xticklabels(["CC", "HNC"])
        ax.set_title(gene, fontsize=11, fontweight="bold")
        ax.set_ylabel("VST expression")
    handles = [Patch(facecolor=NEG_C, alpha=0.55, label="HPV-"),
               Patch(facecolor=POS_C, alpha=0.55, label="HPV+")]
    fig.legend(handles=handles, loc="upper right", fontsize=10)
    fig.suptitle("Panel gene expression by HPV status", fontsize=13, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    save(fig, "fig2_gene_boxplots")


def fig_stability(panel):
    r = pd.read_csv(os.path.join(BASE_DIR, "Panel_Selection_Results", "panel_ranking.csv"))
    r = r.head(20).iloc[::-1]  # top 20, best on top
    y = np.arange(len(r))
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.barh(y - 0.2, r["StabFreq_CC"], height=0.4, color=COH_C["CC"], label="CC")
    ax.barh(y + 0.2, r["StabFreq_HNC"], height=0.4, color=COH_C["HNC"], label="HNC")
    ax.set_yticks(y); ax.set_yticklabels(r["Gene"])
    for tick, gene in zip(ax.get_yticklabels(), r["Gene"]):
        if gene in panel:
            tick.set_color("#d62728"); tick.set_fontweight("bold")
    ax.axvline(0.5, color="gray", ls="--", lw=1)
    ax.set(xlabel="Bootstrap selection frequency", xlim=(0, 1.02),
           title="Elastic-net stability selection (top 20)\npanel genes in red")
    ax.legend(loc="lower right")
    fig.tight_layout()
    save(fig, "fig3_stability")


def _lr():
    return make_pipeline(StandardScaler(),
                         LogisticRegression(class_weight="balanced",
                                            max_iter=5000, random_state=gps.SEED))


def _draw_cm(ax, y_true, y_pred, title):
    cm = confusion_matrix(y_true, y_pred, labels=[0, 1])
    cm_norm = cm / cm.sum(axis=1, keepdims=True)  # row-normalised (recall)
    ax.imshow(cm_norm, cmap="Blues", vmin=0, vmax=1)
    for i in range(2):
        for j in range(2):
            ax.text(j, i, f"{cm[i, j]}\n({cm_norm[i, j]:.2f})", ha="center",
                    va="center", fontsize=10,
                    color="white" if cm_norm[i, j] > 0.5 else "black")
    ax.set_xticks([0, 1]); ax.set_xticklabels(["HPV-", "HPV+"])
    ax.set_yticks([0, 1]); ax.set_yticklabels(["HPV-", "HPV+"])
    ax.set_xlabel("Predicted"); ax.set_ylabel("True")
    ax.set_title(f"{title}\nbal. acc = {balanced_accuracy_score(y_true, y_pred):.3f}",
                 fontsize=10)


def fig_confusion(data, panel):
    """One consolidated confusion-matrix figure (2 x 3).

    Rows    = target/test cohort (CC, HNC).
    Columns = validation regime:
      within-cohort 5-fold CV | cross-cohort transfer | external GEO cohort.
    External-GEO panels are read from the per-sample predictions written by
    geo_external_test.py so they match the fitted model exactly.
    """
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=gps.SEED)

    # within-cohort out-of-fold 5-fold CV predictions (honest)
    within = {}
    for c in ("CC", "HNC"):
        y = (data[c]["Label"] == POS).astype(int)
        within[c] = (y.values, cross_val_predict(_lr(), data[c][panel], y, cv=skf))

    # cross-cohort: train one cohort, test the other; keyed by TEST cohort
    cross = {}
    for tr, te in [("HNC", "CC"), ("CC", "HNC")]:
        ytr = (data[tr]["Label"] == POS).astype(int)
        yte = (data[te]["Label"] == POS).astype(int)
        clf = _lr().fit(data[tr][panel], ytr)
        cross[te] = (yte.values, clf.predict(data[te][panel]), f"{tr}→{te}")

    # external GEO: per-sample panel predictions from geo_external_test.py
    geo_pred = pd.read_csv(os.path.join(GEO_DIR, "geo_panel_predictions.csv"))
    geo = {}
    for c in ("CC", "HNC"):
        sub = geo_pred[geo_pred["Cohort"] == c]
        yt = (sub["True_Label"] == POS).astype(int).values
        geo[c] = (yt, sub["Pred_HPV_positive"].values)

    fig, axes = plt.subplots(2, 3, figsize=(13.5, 9))
    for r, c in enumerate(("CC", "HNC")):
        _draw_cm(axes[r, 0], *within[c], f"{c} · within-cohort (5-fold CV)")
        yte, yp, name = cross[c]
        _draw_cm(axes[r, 1], yte, yp, f"{name} · cross-cohort")
        _draw_cm(axes[r, 2], *geo[c], f"{c} · external GEO")
    fig.suptitle("6-gene panel confusion matrices: within-cohort, cross-cohort "
                 "and external GEO\ncounts; row-normalised (recall) in parentheses",
                 fontsize=13, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    save(fig, "fig4_confusion")


def main():
    # drop any stale manuscript figures from earlier numbering schemes
    for f in glob.glob(os.path.join(OUT, "fig*.png")) + glob.glob(os.path.join(OUT, "fig*.pdf")):
        os.remove(f)

    data = {c: pd.read_csv(os.path.join(BASE_DIR, f"{c}_data.csv")) for c in ("CC", "HNC")}
    ranking = pd.read_csv(os.path.join(BASE_DIR, "Panel_Selection_Results", "panel_ranking.csv"))
    panel = ranking.head(gps.PANEL_SIZE)["Gene"].tolist()
    print(f"Panel: {', '.join(panel)}\nWriting figures to Manuscript_Results/figures/")

    fig_roc_pr(data, panel)
    fig_boxplots(data, panel)
    fig_stability(panel)
    fig_confusion(data, panel)


if __name__ == "__main__":
    main()
