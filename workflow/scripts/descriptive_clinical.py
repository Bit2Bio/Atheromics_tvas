#!/usr/bin/env python3
import pandas as pd
import numpy as np
from pathlib import Path

CLINICAL = "data/processed/clinical_tvas.csv"
GUIDE    = "data/processed/clinical_guide.csv"
OUT_DIR  = Path("results/clinics/descriptive")

df    = pd.read_csv(CLINICAL, index_col=0)
guide = pd.read_csv(GUIDE)

continuous   = guide.loc[guide["type"] == "continuous",   "harmonized_name"].tolist()
categorical  = guide.loc[guide["type"] == "categorical",  "harmonized_name"].tolist()

# --------------------------------------------------------------------------
# Continuous
# --------------------------------------------------------------------------
rows = []
for var in continuous:
    col = df[var]
    rows.append({
        "variable":  var,
        "n_total":   len(col),
        "n_missing": col.isna().sum(),
        "n_valid":   col.notna().sum(),
        "mean":      col.mean(),
        "sd":        col.std(),
        "median":    col.median(),
        "q25":       col.quantile(0.25),
        "q75":       col.quantile(0.75),
    })
cont_df = pd.DataFrame(rows)

# --------------------------------------------------------------------------
# Categorical
# --------------------------------------------------------------------------
rows = []
for var in categorical:
    col = df[var]
    n_total   = len(col)
    n_missing = col.isna().sum()
    counts    = col.value_counts(dropna=True).sort_index()
    for cat, n in counts.items():
        rows.append({
            "variable":  var,
            "n_total":   n_total,
            "n_missing": n_missing,
            "category":  cat,
            "n":         n,
            "pct":       round(n / (n_total - n_missing) * 100, 1),
        })
cat_df = pd.DataFrame(rows)

# --------------------------------------------------------------------------
# Save
# --------------------------------------------------------------------------
OUT_DIR.mkdir(parents=True, exist_ok=True)
cont_df.to_csv(OUT_DIR / "continuous.csv",   index=False)
cat_df.to_csv( OUT_DIR / "categorical.csv",  index=False)

print(f"continuous.csv:  {len(cont_df)} variabili")
print(f"categorical.csv: {len(cat_df)} righe ({len(categorical)} variabili)")
