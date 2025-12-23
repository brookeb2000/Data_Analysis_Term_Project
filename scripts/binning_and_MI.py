
"""
binning_and_MI.py

Purpose:
    Perform mutual information analysis between binned microbial
    abundances and Crohn's disease status, including permutation-based
    p-values and Benjamini–Hochberg FDR correction.

Inputs (from data/cleaned/):
    - joined_MI_ready.csv : one row per tube_id with
      diagnosis label, Final_Metabolomics_Weight and relative abundances for 9 selected taxa.

Outputs (written to results/):
    - tables/MI_results.csv          : MI (bits), raw p-values,
                                       FDR-adjusted p-values,
                                       significance flags per taxon
    - figures/MI_dotplot.png         : dot plot of MI vs bacteria,
                                       highlighting significant taxa

Statistical overview:
    - For each taxon, abundance is discretized into 3 bins:
        0 = absent, 1 = low, 2 = high (based on data distribution).
    - Mutual information is computed between binned abundance and
      Crohn's vs healthy status.
    - Significance is assessed via permutation testing by shuffling
      disease labels many times to approximate the null distribution
      of MI.
    - P-values are adjusted using the Benjamini–Hochberg procedure
      to control false discovery rate across the 9 tests.
"""

import pandas as pd
import numpy as np
from sklearn.metrics import mutual_info_score
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt

# ============================================================
# 1. Load cleaned dataset
#    - Read joined microbiome data with relative abundances
# ============================================================
df = pd.read_csv("data/cleaned/joined_MI_ready.csv")

# ============================================================
# 2. Bin microbial abundances into categorical variables
#    - 0 = absent
#    - 1 = low (0 < abundance ≤ median of non-zero)
#    - 2 = high (abundance > median of non-zero)
#    Note: Binning strategy was justified using distribution
#          exploration in the cleaning script.
# ============================================================

binned_df = df.copy()

# Identify taxon columns: everything except ID / label / metadata
taxa_cols = [
    c for c in binned_df.columns
    if c not in ["tube_id", "Diagnosis", "Final_Metabolomics_Weight"]
]

for taxon in taxa_cols:
    vals = binned_df[taxon].astype(float)

    # median among non-zero abundances
    nonzero = vals[vals > 0]
    if nonzero.empty:
        print(f"Warning: {taxon} has no non-zero values; skipping binning.")
        continue

    cutoff = nonzero.median()

    bin_col = f"{taxon}_bin"

    # 0 = absent, 1 = low, 2 = high
    binned_df[bin_col] = np.where(
        vals == 0,
        0,
        np.where(vals <= cutoff, 1, 2)
    )

# Make sure bin columns are integer-typed
bin_cols = [c for c in binned_df.columns if c.endswith("_bin")]
binned_df[bin_cols] = binned_df[bin_cols].astype("Int64")

# Keep only metadata + binned columns (drop raw abundance columns)
keep_cols = ["tube_id", "Diagnosis", "Final_Metabolomics_Weight"] + bin_cols
binned_df = binned_df[keep_cols]

print("\nColumns retained for MI analysis:")
print(list(binned_df.columns))

print("\nData types after binning and filtering:")
print(binned_df.info())

# ============================================================
# Bin count summary per taxon
#    - Count how many samples fall into 0,1,2 for each microbe
# ============================================================

bin_summary_rows = []

for col in bin_cols:
    taxon = col.replace("_bin", "")

    # Count 0, 1, 2 — and fill missing categories with 0
    counts = binned_df[col].value_counts().reindex([0, 1, 2], fill_value=0)

    bin_summary_rows.append({
        "taxon": taxon,
        "n_zero": counts[0],
        "n_low": counts[1],
        "n_high": counts[2],
    })

bin_summary = pd.DataFrame(bin_summary_rows).sort_values("taxon")

print("\n================ Bin Counts Per Taxon ================")
print(bin_summary.to_string(index=False))
print("======================================================\n")


# ============================================================
# 3. Compute Mutual Information (MI)
#    - MI between binned abundance and disease status
#    - One MI value per microbe
# ============================================================

# Encode diagnosis as integers (e.g. 0 = Crohn, 1 = Healthy; order doesn't matter for MI)
diagnosis_cat = binned_df["Diagnosis"].astype("category")
y_all = diagnosis_cat.cat.codes.to_numpy()   # global label array

mi_results = []

for col in bin_cols:
    taxon = col.replace("_bin", "")

    # Drop any rows with missing bin or diagnosis (just in case)
    mask = binned_df[col].notna() & binned_df["Diagnosis"].notna()
    x = binned_df.loc[mask, col].to_numpy()
    y = y_all[mask.to_numpy()]

    # mutual_info_score returns MI in nats; convert to bits by dividing by log(2)
    mi_nats = mutual_info_score(x, y)
    mi_bits = mi_nats / np.log(2.0)

    mi_results.append({
        "taxon": taxon,
        "MI_bits": mi_bits,
        "n_samples": mask.sum()
    })

mi_df = pd.DataFrame(mi_results).sort_values("taxon")

print("\n================ MI (bits) Per Taxon ================")
print(mi_df.to_string(index=False))
print("=====================================================\n")

# ============================================================
# 4. Permutation test for statistical significance
#    - Shuffle disease labels many times (e.g., 10,000)
#    - Build null distribution of MI for each microbe
#    - Compute raw p-values as tail probability
# ============================================================

#  we create a distribution that shuffles the disease and healthy labels many times and computes MI and we see where our MI value lies
n_permutations = 10000
rng = np.random.default_rng(42)

perm_pvals = []

for idx, row in mi_df.iterrows():
    taxon = row["taxon"]
    col = f"{taxon}_bin"

    # same mask as before, just in case
    mask = binned_df[col].notna() & binned_df["Diagnosis"].notna()
    x = binned_df.loc[mask, col].to_numpy()
    y = y_all[mask.to_numpy()]

    mi_obs = row["MI_bits"]   # observed MI in bits

    null_mi_bits = np.empty(n_permutations)

    for i in range(n_permutations):
        y_perm = rng.permutation(y)
        mi_nats_perm = mutual_info_score(x, y_perm)
        null_mi_bits[i] = mi_nats_perm / np.log(2.0)

    # permutation p-value: how often null MI ≥ observed MI
    p_val = (np.sum(null_mi_bits >= mi_obs) + 1) / (n_permutations + 1)
    perm_pvals.append(p_val)

mi_df["p_perm"] = perm_pvals

print("\n================ Permutation p-values ================")
print(mi_df[["taxon", "MI_bits", "p_perm"]].to_string(index=False))
print("======================================================\n")

# ============================================================
# 5. Benjamini–Hochberg FDR correction (we will use FDR = .05, of the microbes we call significant, we expect about 5% to be false positives)
#    - Adjust p-values across the 10 microbes
#    - Flag taxa with significant adjusted p-values
# ============================================================
reject, p_fdr, _, _ = multipletests(mi_df["p_perm"], alpha=0.05, method="fdr_bh")
mi_df["p_fdr"] = p_fdr
mi_df["significant_FDR_0.05"] = reject

print("\n============ FDR-adjusted p-values ============")
print(mi_df[["taxon", "MI_bits", "p_perm", "p_fdr", "significant_FDR_0.05"]]
      .to_string(index=False))
print("===============================================\n")

# ============================================================
# 6. Save results to CSV
#    - Taxon name, MI (bits), raw p-value, FDR p-value,
#      and significance indicator
# ============================================================
import os

# Ensure results/tables directory exists
os.makedirs("results/tables", exist_ok=True)

# Add a column that records the multiple-testing correction used
mi_df["correction"] = "Benjamini–Hochberg FDR (alpha=0.05)"

# Order rows nicely (e.g., by decreasing MI)
mi_results_out = (
    mi_df
    .sort_values("MI_bits", ascending=False)
    [["taxon",
      "MI_bits",
      "p_perm",
      "p_fdr",
      "correction",
      "significant_FDR_0.05",
      "n_samples"]]
)

# Save to CSV
out_path = "results/tables/MI_results.csv"
mi_results_out.to_csv(out_path, index=False)

print("\n================ Final MI Results ================")
print(mi_results_out.to_string(index=False))
print(f"\nSaved MI results to: {out_path}")
print("==================================================\n")


# ============================================================
# 7. Visualize MI results
#    - Dot plot of MI vs. taxa
#    - Side-By-Side Boxplots for top 3 MI bacteria)
#    - Highlight significant microbes
#    - Export to results/figures/
# ============================================================
import os
from matplotlib.lines import Line2D

os.makedirs("results/figures", exist_ok=True)

# ---------- 7a. Dot plot of MI vs taxa (highlight FDR-significant) ----------

# Sort taxa by MI for nicer plotting
mi_plot_df = mi_df.sort_values("MI_bits", ascending=False).reset_index(drop=True)

x_positions = np.arange(len(mi_plot_df))
mi_values = mi_plot_df["MI_bits"].to_numpy()
sig_flags = mi_plot_df["significant_FDR_0.05"].to_numpy()

# Color significant vs non-significant
colors = ["red" if sig else "gray" for sig in sig_flags]

plt.figure(figsize=(8, 5))
plt.scatter(x_positions, mi_values, s=60, c=colors)

# Horizontal reference line at 0
plt.axhline(0, color="black", linewidth=0.8, linestyle="--")

plt.xticks(x_positions, mi_plot_df["taxon"], rotation=45, ha="right")
plt.ylabel("Mutual Information (bits)")
plt.title("Mutual Information between Binned Abundance and Crohn Status")

# Legend for significance
legend_elements = [
    Line2D([0], [0], marker='o', linestyle='None',
           label='FDR-significant', markerfacecolor='red'),
    Line2D([0], [0], marker='o', linestyle='None',
           label='Not significant', markerfacecolor='gray'),
]
plt.legend(handles=legend_elements, loc="upper right")

plt.tight_layout()
dotplot_path = "results/figures/MI_dotplot.png"
plt.savefig(dotplot_path, dpi=300)
plt.close()
print(f"Saved MI dot plot to: {dotplot_path}")

# ---------- 7b. Side-by-side boxplots for top 3 MI taxa ----------

# Use original relative abundances from df (not binned_df)
top3 = mi_plot_df["taxon"].head(3).tolist()
print("\nTop 3 taxa by MI (for boxplots):", top3)

for taxon in top3:
    crohn_vals = df.loc[df["Diagnosis"] == "Crohn", taxon].dropna()
    healthy_vals = df.loc[df["Diagnosis"] == "Healthy", taxon].dropna()

    plt.figure(figsize=(5, 5))
    plt.boxplot(
        [crohn_vals, healthy_vals],
        labels=["Crohn", "Healthy"],
        showfliers=False
    )

    plt.ylabel("Relative abundance")
    plt.title(f"{taxon}: abundance by diagnosis")

    plt.tight_layout()
    boxplot_path = f"results/figures/boxplot_{taxon}.png"
    plt.savefig(boxplot_path, dpi=300)
    plt.close()
    print(f"Saved boxplot for {taxon} to: {boxplot_path}")

# ---------- 7c. Side-by-side boxplots for ALL taxa ----------

# Use the same taxa order as in the MI table
all_taxa = mi_plot_df["taxon"].tolist()
print("\nAll taxa for boxplots:", all_taxa)

for taxon in all_taxa:
    crohn_vals = df.loc[df["Diagnosis"] == "Crohn", taxon].dropna()
    healthy_vals = df.loc[df["Diagnosis"] == "Healthy", taxon].dropna()

    # Skip taxa with too few observations in either group
    if len(crohn_vals) < 3 or len(healthy_vals) < 3:
        print(f"Skipping {taxon}: too few observations for boxplot.")
        continue

    plt.figure(figsize=(5, 5))
    plt.boxplot(
        [crohn_vals, healthy_vals],
        labels=["Crohn", "Healthy"],
        showfliers=False
    )

    plt.ylabel("Relative abundance")
    plt.title(f"{taxon}: abundance by diagnosis")

    plt.tight_layout()
    boxplot_path = f"results/figures/boxplot_all_{taxon}.png"
    plt.savefig(boxplot_path, dpi=300)
    plt.close()
    print(f"Saved boxplot for {taxon} to: {boxplot_path}")
