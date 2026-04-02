#!/usr/bin/env python3
"""
Compare significant pathways between inflamed and uninflamed UC meta-analyses.
Usage:
    python ibd/scripts/S_compare_pathways.py

Outputs:
  1. pathway_comparison.csv - full outer join with per-pathway annotations
  2. Summary stats printed to stdout
"""

import pandas as pd
import sys

# --- Config ---
FILE_INFLAMED = "results/paper/uc_vs_hc/pathway_meta_significant.csv"
FILE_UNINFLAMED = "results/paper/uninf_vs_hc/pathway_meta_significant.csv"
OUTPUT = "results/paper/pathway_comparison.csv"

KEEP_COLS = [
    "term", "combined_nes", "direction", "q_value",
    "I2", "direction_ratio", "n_up", "n_down", "datasets",
    "lead_genes_core", "lead_genes_core_size",
]

# --- Load ---
inf = pd.read_csv(FILE_INFLAMED)
uni = pd.read_csv(FILE_UNINFLAMED)

# Subset to useful columns
inf_sub = inf[KEEP_COLS].copy()
uni_sub = uni[KEEP_COLS].copy()

# --- Merge (outer join on term) ---
merged = inf_sub.merge(
    uni_sub,
    on="term",
    how="outer",
    suffixes=("_inflamed", "_uninflamed"),
    indicator=True,
)

# --- Classify presence ---
presence_map = {
    "both": "common",
    "left_only": "inflamed_only",
    "right_only": "uninflamed_only",
}
merged["presence"] = merged["_merge"].map(presence_map).astype(str)
merged.drop(columns=["_merge"], inplace=True)

# --- Direction agreement (only for common pathways) ---
mask_common = merged["presence"] == "common"

merged["direction_agreement"] = pd.Series(dtype="object")
agree_series = (
    merged.loc[mask_common, "direction_inflamed"]
    == merged.loc[mask_common, "direction_uninflamed"]
).map({True: "agree", False: "disagree"})
merged.loc[agree_series.index, "direction_agreement"] = agree_series.values

# Classify the concordance pattern for common pathways
def classify_pattern(row):
    if row["presence"] != "common":
        return None
    d_inf = row["direction_inflamed"]
    d_uni = row["direction_uninflamed"]
    if d_inf == "up" and d_uni == "up":
        return "both_up"
    elif d_inf == "down" and d_uni == "down":
        return "both_down"
    elif d_inf == "up" and d_uni == "down":
        return "inflamed_up_uninflamed_down"
    elif d_inf == "down" and d_uni == "up":
        return "inflamed_down_uninflamed_up"
    return None

merged["direction_pattern"] = merged.apply(classify_pattern, axis=1)

# NES difference (inflamed – uninflamed) for common pathways
merged["nes_diff"] = merged["combined_nes_inflamed"] - merged["combined_nes_uninflamed"]

# For unique pathways, carry over direction into a single column for convenience
merged["direction_any"] = merged["direction_inflamed"].fillna(merged["direction_uninflamed"])

# --- Sort: common first (by nes_diff), then inflamed_only, then uninflamed_only ---
merged["_sort"] = merged["presence"].map({"common": 0, "inflamed_only": 1, "uninflamed_only": 2}).astype(int)
merged = merged.sort_values(
    ["_sort", "nes_diff", "combined_nes_inflamed"],
    ascending=[True, True, True],
    na_position="last",
).drop(columns=["_sort"]).reset_index(drop=True)

# --- Save ---
merged.to_csv(OUTPUT, index=False)

# --- Summary (recompute mask after sort) ---
n_inf = len(inf)
n_uni = len(uni)
mask_common = merged["presence"] == "common"
n_common = mask_common.sum()
n_inf_only = (merged["presence"] == "inflamed_only").sum()
n_uni_only = (merged["presence"] == "uninflamed_only").sum()

common = merged[mask_common]
n_agree = (common["direction_agreement"] == "agree").sum()
n_disagree = (common["direction_agreement"] == "disagree").sum()

patterns = common["direction_pattern"].value_counts()

print("=" * 60)
print("PATHWAY COMPARISON SUMMARY")
print("=" * 60)
print(f"Significant pathways – inflamed:    {n_inf}")
print(f"Significant pathways – uninflamed:  {n_uni}")
print(f"Common (in both):                   {n_common}")
print(f"Inflamed only:                      {n_inf_only}")
print(f"Uninflamed only:                    {n_uni_only}")
print()
print(f"Direction agreement (common):       {n_agree} / {n_common}  "
      f"({100*n_agree/n_common:.1f}%)" if n_common else "")
print(f"Direction disagreement (common):    {n_disagree} / {n_common}  "
      f"({100*n_disagree/n_common:.1f}%)" if n_common else "")
print()
print("Direction patterns among common pathways:")
for pat, cnt in patterns.items():
    print(f"  {pat:40s} {cnt:4d}")
print()

# Direction breakdown for each group
for grp_name, grp_label in [("inflamed_only", "Inflamed-only"),
                             ("uninflamed_only", "Uninflamed-only")]:
    sub = merged[merged["presence"] == grp_name]
    if len(sub):
        up = (sub["direction_any"] == "up").sum()
        down = (sub["direction_any"] == "down").sum()
        print(f"{grp_label} pathways: {len(sub)} total  ({up} up, {down} down)")

print()
print("Disagreeing pathways (direction flips):")
disagree = common[common["direction_agreement"] == "disagree"].sort_values("nes_diff")
if len(disagree):
    for _, row in disagree.iterrows():
        print(f"  {row['term'][:70]:70s}  "
              f"inf={row['direction_inflamed']:4s}  "
              f"uni={row['direction_uninflamed']:4s}  "
              f"NES_inf={row['combined_nes_inflamed']:+.2f}  "
              f"NES_uni={row['combined_nes_uninflamed']:+.2f}")
else:
    print("  (none)")

print(f"\nOutput saved to: {OUTPUT}")
