#!/usr/bin/env python3
"""
Preprocessing proteomics for TVAS samples.

Steps:
  1. Extract NPX from raw Olink file (SAMPLE + assay rows only)
  2. Filter to TVAS_ samples only
  3. LOD filter: keep proteins with >25% samples above 5th percentile
  4. Variance filter: remove bottom 20th percentile
"""

import pandas as pd
from pathlib import Path

OLINK_RAW = "data/raw/proteomics/FedericiAllPlates_Extended_NPX_2026-03-11.csv"
OUT_DIR   = Path("data/processed")

LOD_QUANTILE      = 0.05   # NPX percentile used as limit-of-detection proxy
LOD_MIN_FRACTION  = 0.25   # min fraction of samples above LOD to keep a protein
VAR_QUANTILE      = 0.20   # bottom variance percentile removed

# --------------------------------------------------------------------------
# Extract NPX from raw Olink file
# --------------------------------------------------------------------------
df = pd.read_csv(OLINK_RAW, sep=";", usecols=["SampleID", "SampleType", "Assay", "AssayType", "NPX"])

mat = (
    df.loc[(df["SampleType"] == "SAMPLE") & (df["AssayType"] == "assay")]
    .pivot(index="Assay", columns="SampleID", values="NPX")
)

print(f"After extraction: {mat.shape[0]} proteins x {mat.shape[1]} samples")

# --------------------------------------------------------------------------
# Save full matrix (all samples)
# --------------------------------------------------------------------------
OUT_DIR.mkdir(parents=True, exist_ok=True)
mat.to_csv(OUT_DIR / "proteomics_all.csv")
print(f"proteomics_all.csv:  {mat.shape[0]} proteine x {mat.shape[1]} campioni")

# --------------------------------------------------------------------------
# Filter to TVAS samples only
# --------------------------------------------------------------------------
tvas_cols = [c for c in mat.columns if c.startswith("TVAS_")]
mat_tvas = mat[tvas_cols].astype(float)
print(f"After TVAS filter: {mat_tvas.shape[0]} proteine x {mat_tvas.shape[1]} campioni")

# --------------------------------------------------------------------------
# LOD filter: keep proteins with >25% of samples above the 5th percentile
# --------------------------------------------------------------------------
lod_threshold = mat_tvas.stack().quantile(LOD_QUANTILE)
above_lod = (mat_tvas > lod_threshold).where(mat_tvas.notna())
keep_lod  = above_lod.mean(axis=1) >= LOD_MIN_FRACTION
mat_lod   = mat_tvas.loc[keep_lod]
print(f"After LOD filter (>{LOD_QUANTILE:.0%} pct={lod_threshold:.3f} in >={LOD_MIN_FRACTION:.0%} samples): "
      f"{mat_lod.shape[0]} proteine (rimosse {mat_tvas.shape[0] - mat_lod.shape[0]})")

# --------------------------------------------------------------------------
# Variance filter: remove bottom 20th percentile
# --------------------------------------------------------------------------
variance     = mat_lod.var(axis=1)
var_threshold = variance.quantile(VAR_QUANTILE)
keep_var     = variance > var_threshold
mat_filt     = mat_lod.loc[keep_var]
print(f"After variance filter (>{VAR_QUANTILE:.0%} pct={var_threshold:.3f}): "
      f"{mat_filt.shape[0]} proteine (rimosse {mat_lod.shape[0] - mat_filt.shape[0]})")

# --------------------------------------------------------------------------
# Save
# --------------------------------------------------------------------------
mat_filt.to_csv(OUT_DIR / "proteomics_tvas.csv")
print(f"proteomics_tvas.csv: {mat_filt.shape[0]} proteine x {mat_filt.shape[1]} campioni")
