#!/usr/bin/env python3
"""
Preprocessing clinical data from TVAS sheet.

Output:
  data/processed/clinical_tvas.csv
  data/processed/clinical_guide.csv
"""

import pandas as pd
import numpy as np
from pathlib import Path

EXCEL    = "data/raw/metabolon e glicemia 090926_rev_LA.xlsx"
SHEET    = "TVAS"
MET_FILE = "data/processed/metabolomics_tvas.csv"
OUT_DIR  = Path("data/processed")

VARIABLES = [
    (v["excel_name"], v["name"], v["type"])
    for v in snakemake.config["variables"]
]

# --------------------------------------------------------------------------
# Load
# --------------------------------------------------------------------------
df = pd.read_excel(EXCEL, sheet_name=SHEET)

# CLIENT_IDENTIFIER
df["CLIENT_IDENTIFIER"] = (
    "TVAS_" + df["n. placche TVAS "].astype(float).astype(int).astype(str)
)

# --------------------------------------------------------------------------
# Build output dataframe
# --------------------------------------------------------------------------
out = pd.DataFrame(index=df.index)
out["CLIENT_IDENTIFIER"] = df["CLIENT_IDENTIFIER"]

for orig, name, _ in VARIABLES:
    out[name] = df[orig]

out = out.set_index("CLIENT_IDENTIFIER")

# --------------------------------------------------------------------------
# Generic numeric conversion for continuous variables
# --------------------------------------------------------------------------
continuous = [name for _, name, typ in VARIABLES if typ == "continuous"]
for col in continuous:
    out[col] = pd.to_numeric(out[col], errors="coerce")

# --------------------------------------------------------------------------
# Replace invalid values with NA for categorical variables
# --------------------------------------------------------------------------
categorical = [name for _, name, typ in VARIABLES if typ == "categorical"]
invalid = {"?", "1?", "/", "", " "}
for col in categorical:
    out[col] = out[col].apply(lambda x: np.nan if str(x).strip() in invalid else x)

# --------------------------------------------------------------------------
# Specific fixes
# --------------------------------------------------------------------------
# sesso: 2 → 0
out["gender"] = out["gender"].replace(2, 0)

# hsPCR: strip '<' and convert
out["hscrp"] = (
    out["hscrp"].astype(str)
    .str.replace("<", "", regex=False)
    .str.replace(",", ".", regex=False)
    .str.strip()
)
out["hscrp"] = pd.to_numeric(out["hscrp"], errors="coerce")

# plaque_type: 2 → 1 (instabile)
out["plaque_type"] = out["plaque_type"].replace(2, 1)

# --------------------------------------------------------------------------
# Guide
# --------------------------------------------------------------------------
guide = pd.DataFrame(
    [(name, typ) for _, name, typ in VARIABLES],
    columns=["harmonized_name", "type"]
)

# --------------------------------------------------------------------------
# Align to metabolomics sample order
# --------------------------------------------------------------------------
met_ids = pd.read_csv(MET_FILE, index_col=0).index
out = out.reindex(met_ids)

# --------------------------------------------------------------------------
# Save
# --------------------------------------------------------------------------
OUT_DIR.mkdir(parents=True, exist_ok=True)
out.to_csv(OUT_DIR / "clinical_tvas.csv")
guide.to_csv(OUT_DIR / "clinical_guide.csv", index=False)

print(f"clinical_tvas.csv:  {out.shape[0]} campioni x {out.shape[1]} variabili")
print(f"clinical_guide.csv: {len(guide)} variabili")
print("\nNA per variabile:")
for col in out.columns:
    na = out[col].isna().sum()
    if na > 0:
        print(f"  {col:<20} NA={na}")
