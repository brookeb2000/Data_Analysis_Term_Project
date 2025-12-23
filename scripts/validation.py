"""
validation.py

Purpose:
    Perform a simple biological/clinical validation of the top microbes
    identified by mutual information, by relating their abundance to an
    independent disease-related measure (Final_Metabolomics_Weight).

Inputs:
    - data/cleaned/joined_MI_ready.csv
        * continuous relative abundances for the 9 taxa
        * diagnosis label
        * Final_Metabolomics_Weight (validation measure)
    - results/tables/MI_results.csv
        * MI values and FDR-adjusted p-values for the 9 taxa

Outputs:
    - figures/validation_scatter_<taxon>.png
        * scatter plots of microbe abundance vs Final_Metabolomics_Weight for the
          top 3 taxa in Crohn's disease patients

Validation idea:
    - Use the MI results to identify the most informative taxa.
    - For the top taxa, examine whether their abundance varies with
      Final_Metabolomics_Weight (a disease severity proxy).
    - For a pathogenic taxon, we expect higher abundance in patients
      with worse severity; for a protective taxon, the opposite trend.
"""

import os
import pandas as pd
import matplotlib.pyplot as plt

# ============================================================
# 1. Load data
# ============================================================
df = pd.read_csv("data/cleaned/joined_MI_ready.csv")

# ============================================================
# 2. Define top taxa to validate (from MI results)
# ============================================================
top_taxa = [
    "Escherichia_coli",
    "Faecalibacterium_prausnitzii",
    "Faecalibacterium",
]

# Sanity check: make sure all columns exist
missing_taxa = [t for t in top_taxa if t not in df.columns]
if missing_taxa:
    raise ValueError(f"Missing taxa in joined_MI_ready.csv: {missing_taxa}")

if "Final_Metabolomics_Weight" not in df.columns:
    raise ValueError("Missing Final_Metabolomics_Weight in joined_MI_ready.csv")

# Drop rows missing the validation variable
df = df.dropna(subset=["Final_Metabolomics_Weight"])

print(f"Number of patients with metabolomics data: {len(df)}")

# ============================================================
# 3. Make scatter plots: abundance (x) vs metabolomics weight (y)
#    - One plot per taxon, all diagnoses
#    - Color points by Diagnosis
# ============================================================
os.makedirs("results/figures", exist_ok=True)

# Simple color mapping; anything not listed falls back to gray
color_map = {
    "Healthy": "tab:blue",
    "Crohn": "tab:red",
    "colitis": "tab:green",
}

for taxon in top_taxa:
    subset = df[[taxon, "Final_Metabolomics_Weight", "Diagnosis"]].dropna()

    plt.figure(figsize=(6, 4))

    # Plot each diagnosis separately so we can color + label
    for diagnosis in subset["Diagnosis"].unique():
        diag_data = subset[subset["Diagnosis"] == diagnosis]
        plt.scatter(
            diag_data[taxon],
            diag_data["Final_Metabolomics_Weight"],
            alpha=0.6,
            label=diagnosis,
            color=color_map.get(diagnosis, "gray"),
        )

    plt.xlabel(f"{taxon} relative abundance")
    plt.ylabel("Final_Metabolomics_Weight")
    plt.title(f"{taxon} vs metabolomic weight (all diagnoses)")
    plt.legend(title="Diagnosis")

    plt.tight_layout()
    out_path = f"results/figures/validation_scatter_all_{taxon}.png"
    plt.savefig(out_path, dpi=300)
    plt.close()

    print(f"Saved scatter plot for {taxon} (all diagnoses) to: {out_path}")



