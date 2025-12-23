"""
clean_and_join.py

Purpose:
    Prepare a single analysis-ready dataset linking Crohn's disease status
    to microbial relative abundances and selected clinical metadata.

Inputs (expected in data/raw/):
    - patient_metadata.csv     : patient-level phenotype and clinical variables
    - bacteria_table.csv     : microbial relative abundances per patient

Outputs (written to data/cleaned/):
    - joined_microbiome_MI_ready.csv : one row per tube_id with:
        * tube_id
        * Diagnosis (Crohn's vs healthy label)
        * Relative abundance for 9 selected taxa
        * Final_Metabolomics_Weight (validation variable)

Notes:
    - This script is intended to be fully reproducible: running it should
      overwrite any existing cleaned dataset using the current raw inputs.
"""

import pandas as pd
import numpy as np

# ============================================================
# 1. Load raw data
# ============================================================

raw_patient_metadata = pd.read_csv(
    "data/raw/patient_metadata.tsv",
    sep="\t"
)

raw_bacteria_table = pd.read_csv(
    "data/raw/bacteria_table.tsv",
    sep="\t"
)

# ============================================================
# 2. Basic inspection and sanity checks
#   - Shape (rows/columns)
#   - Validate patient counts by diagnosis group
#   - Detect duplicate sample_ids
# ============================================================

print("raw patient metadata shape:", raw_patient_metadata.shape)
print("raw bacteria table shape:", raw_bacteria_table.shape)

# Restructure bacteria table → one row per tube_id

# Rename the first column to something meaningful (taxon names)
raw_bacteria_table.rename(
    columns={raw_bacteria_table.columns[0]: "taxon"},
    inplace=True
)

# Transpose so that samples become rows and taxa become columns
bact_by_sample = (
    raw_bacteria_table
    .set_index("taxon")   # index = taxon names
    .T                    # now index = tube_id, columns = taxa
    .reset_index()        # bring tube_id out of index
    .rename(columns={"index": "tube_id"})
)

# Confirm structure
print("Restructured bacteria table:", bact_by_sample.shape)
print(bact_by_sample.head())

#check count of patients by diagnosis group for raw_patient_metadata
print("Patient counts by diagnosis group:")
print(raw_patient_metadata['Diagnosis'].value_counts())

#look for duplicate tupe_ids in raw_patient_metadata
duplicate_tube_ids = raw_patient_metadata[raw_patient_metadata.duplicated(subset=['tube_id'], keep=False)]
if not duplicate_tube_ids.empty:
    print("Duplicate tube_ids found in patient metadata:")
    print(duplicate_tube_ids)

# we found a sample that is the same record copied twice with a typo, so we will drop one
raw_patient_metadata = raw_patient_metadata.drop_duplicates(
    subset="tube_id",
    keep="first"
)

#validate we have no more duplicates
duplicate_tube_ids = raw_patient_metadata[raw_patient_metadata.duplicated(subset=['tube_id'], keep=False)]
if duplicate_tube_ids.empty:
    print("No duplicate tube_ids found in patient metadata after cleanup.")

#rename dataframes for clarity
patient_metadata_cleaned = raw_patient_metadata
bacteria_table_cleaned = bact_by_sample

print("patient_metadata_cleaned shape:", patient_metadata_cleaned.shape)
print("bacteria_table_cleaned shape:", bacteria_table_cleaned.shape)

# ============================================================
# 3. Standardize column names, types, and labels
#   - Convert types (e.g., timestamps, numeric abundance)
#   - Ensure consistent disease categories (e.g. "Crohn's" vs "CD")
#   - For each sampkle, convert raw bacterial counts to relative abundances
# ============================================================

# fix types for patient_metadata_cleaned for tube_id, Diagnosis and Final_Metabolomics_Weight
patient_metadata_cleaned['tube_id'] = patient_metadata_cleaned['tube_id'].astype(str)   
patient_metadata_cleaned['Final_Metabolomics_Weight'] = pd.to_numeric(
    patient_metadata_cleaned['Final_Metabolomics_Weight'],
    errors='coerce'
)

# fix data types for bacteria_table_cleaned for tube_id and taxa columns
bacteria_table_cleaned['tube_id'] = bacteria_table_cleaned['tube_id'].astype(str)
for col in bacteria_table_cleaned.columns[1:]:
    bacteria_table_cleaned[col] = pd.to_numeric(
        bacteria_table_cleaned[col],
        errors='coerce'
    )

# standardize Diagnosis labels
patient_metadata_cleaned["Diagnosis"] = (
    patient_metadata_cleaned["Diagnosis"]
    .str.strip()
    .replace({"Healthy_control": "Healthy", "CD": "Crohn"})
    .astype("category")
)
#verify data types are correct
print("Patient metadata data types:")
print(patient_metadata_cleaned.dtypes)
print("Bacteria table data types:")
print(bacteria_table_cleaned.dtypes)


# ------------------------------------------------------------
# Convert raw bacterial counts → relative abundance per sample
# ------------------------------------------------------------

# Identify microbe columns (all except metadata)
microbe_cols = [
    c for c in bacteria_table_cleaned.columns
    if c not in ["tube_id"]  # only tube_id is metadata in bacteria table
]

# Normalize each row so abundances sum to 1
row_sums = bacteria_table_cleaned[microbe_cols].sum(axis=1)
bacteria_table_cleaned[microbe_cols] = (
    bacteria_table_cleaned[microbe_cols].div(row_sums, axis=0)
)

# Confirm normalization worked: each sample should sum to 1
check_sums = bacteria_table_cleaned[microbe_cols].sum(axis=1).round(5)
print("\nRow sums after normalization (should be 1.00000):")
print(check_sums.head())

# Optional: ensure no row sums deviate too far from 1
if not np.allclose(check_sums, 1.0, rtol=1e-5):
    raise ValueError("Normalization failed: not all rows sum to 1")
else:
    print("✓ Relative abundance normalization successful!")

# ============================================================
# 4. Join metadata and microbial abundance tables
#   - Inner join on tube_id
#   - Confirm joined patient counts match expectations
# ============================================================
metadata_cols = ["tube_id", "Diagnosis", "Final_Metabolomics_Weight"]

patient_subset = patient_metadata_cleaned[metadata_cols].copy()

# Perform the join (only matching tube_ids kept)
joined_df = bacteria_table_cleaned.merge(
    patient_subset,
    on="tube_id",
    how="inner"
)
# Confirm joined shape, (some tube_ids didnt have matching patient records, some patient records didnt have tube_ids)
print("Joined dataset shape:", joined_df.shape)

print("Diagnosis counts in joined_df:")
print(joined_df["Diagnosis"].value_counts())



# ============================================================
# 5. Explore distributions of microbial abundances
#   - Histograms to choose meaningful binning strategy
# ============================================================

#inspect distributions of microbial abundances
import numpy as np
import matplotlib.pyplot as plt

bacteria_cols = [
    c for c in joined_df.columns
    if c not in ["tube_id", "Diagnosis", "Final_Metabolomics_Weight"]
]

for taxon in bacteria_cols:
    vals = joined_df[taxon].dropna().values
    
    # Compute quartiles
    q1 = np.percentile(vals, 25)
    q2 = np.percentile(vals, 50)  # median
    q3 = np.percentile(vals, 75)

    plt.figure(figsize=(8,4))

    bins = np.linspace(0, 1, 30)
    plt.hist(vals, bins=bins, edgecolor="black", alpha=0.7)

    # Overlay quartile cut markers
    for q, label in zip([q1, q2, q3], ["Q1", "Median", "Q3"]):
        plt.axvline(q, color="red", linestyle="--", linewidth=1)
        plt.text(q, plt.ylim()[1]*0.9, label,
                 rotation=90, verticalalignment="top",
                 color="red", fontsize=9, fontweight="bold")
    
    plt.title(f"Distribution of {taxon} with quartiles")
    plt.xlabel("Relative abundance")
    plt.ylabel("Number of samples")
    plt.xlim(0, 1)
    plt.xticks(np.arange(0, 1.01, 0.1))
    plt.tight_layout()
    plt.show()

# Binning strategy that makes sense based on these histograms:
# Microbial abundances are zero-inflated and right-skewed.
# Treat true zeros as a separate category (absent),
# then split non-zero values at their median into low and high abundance
summary_rows = []

for taxon in bacteria_cols:
    vals = joined_df[taxon].dropna()

    n_total = vals.shape[0]
    n_zero = (vals == 0).sum()
    prop_zero = n_zero / n_total if n_total > 0 else np.nan

    nonzero_vals = vals[vals > 0]
    median_nonzero = nonzero_vals.median() if not nonzero_vals.empty else np.nan

    summary_rows.append({
        "taxon": taxon,
        "median_nonzero": median_nonzero,
        "n_zero": n_zero,
        "prop_zero": prop_zero,
    })

bacteria_summary = pd.DataFrame(summary_rows).sort_values("taxon")

print("\nBinning summary per taxon:")
print(bacteria_summary)

# Note:
# Because zero-inflation was modest across taxa (mostly 2–15%), 
# presence vs. absence is not the dominant driver of variation in microbial abundance. 
# Therefore, splitting the non-zero values into low and high (using the median among non-zeros) 
# preserves the biologically relevant abundance differences that are likely to contribute most to the MI signal.

# ============================================================
# 6. Final cleanup and export
#   - Write cleaned dataset to data/cleaned/
# ============================================================

joined_df.to_csv("data/cleaned/joined_MI_ready.csv", index=False)
print("\n✓ Cleaned and joined dataset written to data/cleaned/joined_MI_ready.csv")








