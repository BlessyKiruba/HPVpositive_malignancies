######################################################################
## External validation on independent GEO cohorts.                   ##
##                                                                    ##
## For each cancer the TCGA data (CC_data / HNC_data) is split into   ##
## train/test exactly as before (stratified, test_size=0.3,          ##
## random_state=SEED). The model is trained ONLY on the TCGA train    ##
## split, then evaluated on:                                          ##
##   - the held-out TCGA internal test split, and                     ##
##   - the matched external GEO cohort (held out entirely).           ##
##     CC  <- GSE151666 (cervical)                                    ##
##     HNC <- GSE74927  (head & neck)                                 ##
##                                                                    ##
## Model = the 6-gene panel logistic regression (all-52 for           ##
## reference). Scaler is fit on the TCGA train split only and applied ##
## to every test set (no leakage). Results -> GEO_test/               ##
##                                                                    ##
## Writes metrics (geo_external_results.csv) and per-sample panel     ##
## predictions (geo_panel_predictions.csv). The manuscript confusion  ##
## matrices are drawn from those predictions by make_figures.py.      ##
######################################################################

import os
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.metrics import (balanced_accuracy_score, roc_auc_score, f1_score,
                             confusion_matrix)

import gene_panel_selection as gps
import standardise as std  # reuse column/label standardisation

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(BASE_DIR, "GEO_test")
POS = "HPV_positive"

# cohort -> (TCGA file, GEO file, GEO reader)
COHORTS = {
    "CC":  ("CC_data.csv",  "CC GSE151666_HPV_VST_52genes.csv", pd.read_csv),
    "HNC": ("HNC_data.csv", "HNC VST_genes_of_interest_GSE74927_HPVpos_vs_HPVneg.xlsx", pd.read_excel),
}


def lr():
    # SimpleImputer (fit on TCGA train) neutralises any gene missing in an
    # external cohort -> train mean -> ~0 after scaling. No-op when no NaNs,
    # so the panel/internal results are identical to a plain pipeline.
    return make_pipeline(SimpleImputer(strategy="mean"), StandardScaler(),
                         LogisticRegression(class_weight="balanced",
                                            max_iter=5000, random_state=gps.SEED))


def evaluate(clf, X, y):
    yp = clf.predict(X)
    proba = clf.predict_proba(X)[:, 1]
    tn, fp, fn, tp = confusion_matrix(y, yp, labels=[0, 1]).ravel()
    metrics = {
        "N": len(y), "N_pos": int(y.sum()), "N_neg": int((1 - y).sum()),
        "Balanced_Acc": balanced_accuracy_score(y, yp),
        "AUC": roc_auc_score(y, proba) if len(np.unique(y)) > 1 else np.nan,
        "F1_Macro": f1_score(y, yp, average="macro", zero_division=0),
        "Sensitivity_HPVpos": tp / (tp + fn) if (tp + fn) else np.nan,
        "Specificity_HPVneg": tn / (tn + fp) if (tn + fp) else np.nan,
    }
    return metrics, (y.values, yp)


def main():
    os.makedirs(OUT, exist_ok=True)
    panel = pd.read_csv(os.path.join(BASE_DIR, "Panel_Selection_Results",
                                     "panel_ranking.csv")).head(gps.PANEL_SIZE)["Gene"].tolist()
    print(f"Panel ({gps.PANEL_SIZE}): {', '.join(panel)}\n")

    rows, preds = [], []
    for cohort, (tcga_f, geo_f, reader) in COHORTS.items():
        tcga = pd.read_csv(os.path.join(BASE_DIR, tcga_f))
        geo = std.standardise(reader(os.path.join(BASE_DIR, geo_f)))
        geo.to_csv(os.path.join(OUT, f"{cohort}_GEO_standardised.csv"), index=False)

        missing_panel = [g for g in panel if g not in geo.columns]
        if missing_panel:
            raise ValueError(f"{cohort} GEO missing PANEL genes: {missing_panel}")
        all52 = [c for c in tcga.columns if c not in ("Sample_ID", "Label")]
        absent = [g for g in all52 if g not in geo.columns]
        nan_genes = [g for g in all52 if g in geo.columns and geo[g].isna().all()]
        unusable = sorted(set(absent + nan_genes))
        if unusable:
            print(f"  {cohort} GEO: {len(unusable)} of 52 genes unavailable "
                  f"(imputed for All_52; none are panel genes): {', '.join(unusable)}")

        y_all = (tcga["Label"] == POS).astype(int)
        # SAME split as before: stratified 70/30, fixed seed
        idx_tr, idx_te = train_test_split(
            tcga.index, test_size=0.3, random_state=gps.SEED, stratify=y_all)
        y_geo = (geo["Label"] == POS).astype(int)

        for fs_name, feats in [("Panel_6", panel), ("All_52", all52)]:
            clf = lr().fit(tcga.loc[idx_tr, feats], y_all.loc[idx_tr])

            m_int, _ = evaluate(clf, tcga.loc[idx_te, feats], y_all.loc[idx_te])
            m_geo, _ = evaluate(clf, geo[feats], y_geo)
            for tset, m in [("Internal_test", m_int), ("GEO_external", m_geo)]:
                rows.append({"Cohort": cohort, "Feature_Set": fs_name,
                             "Test_Set": tset, **m})

            if fs_name == "Panel_6":
                # per-sample GEO predictions (also feed the manuscript
                # confusion matrices drawn by make_figures.py)
                proba = clf.predict_proba(geo[panel])[:, 1]
                preds.append(pd.DataFrame({
                    "Cohort": cohort, "Sample_ID": geo["Sample_ID"],
                    "True_Label": geo["Label"],
                    "Pred_HPV_positive": clf.predict(geo[panel]),
                    "Proba_HPV_positive": proba.round(4),
                }))

    res = pd.DataFrame(rows)
    res_round = res.copy()
    for c in ["Balanced_Acc", "AUC", "F1_Macro", "Sensitivity_HPVpos", "Specificity_HPVneg"]:
        res_round[c] = res_round[c].round(4)
    res_round.to_csv(os.path.join(OUT, "geo_external_results.csv"), index=False)
    pd.concat(preds, ignore_index=True).to_csv(
        os.path.join(OUT, "geo_panel_predictions.csv"), index=False)

    print(res_round.to_string(index=False))
    print(f"\nResults saved to: GEO_test/")
    print("  geo_external_results.csv, geo_panel_predictions.csv,")
    print("  {CC,HNC}_GEO_standardised.csv")


if __name__ == "__main__":
    main()
