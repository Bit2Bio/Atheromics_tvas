#!/usr/bin/env python3
"""
Preprocessing metabolomics data: filter TVAS samples only.

Output:
  data/processed/metabolomics_tvas.csv
"""

import pandas as pd
from pathlib import Path

LOG_DATA   = "data/raw/metabolomics/log_transformed_data.csv"
CHEM_META  = "data/raw/metabolomics/chemical_metadata.csv"
OUT_DIR    = Path("data/processed")

# --------------------------------------------------------------------------
# Load
# --------------------------------------------------------------------------
met  = pd.read_csv(LOG_DATA, index_col=0)
chem = pd.read_csv(CHEM_META, index_col=0)

# --------------------------------------------------------------------------
# Filter TVAS
# --------------------------------------------------------------------------
tvas = met[met.index.str.startswith("TVAS_")]

# --------------------------------------------------------------------------
# Save
# --------------------------------------------------------------------------
OUT_DIR.mkdir(parents=True, exist_ok=True)
tvas.to_csv(OUT_DIR / "metabolomics_tvas.csv")

print(f"metabolomics_tvas.csv: {tvas.shape[0]} campioni x {tvas.shape[1]} features")
