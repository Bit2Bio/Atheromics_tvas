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
# Filter to TVAS samples only and save
# --------------------------------------------------------------------------
tvas_cols = [c for c in mat.columns if c.startswith("TVAS_")]
mat_tvas = mat[tvas_cols].astype(float)

mat_tvas.to_csv(OUT_DIR / "proteomics_tvas.csv")
print(f"proteomics_tvas.csv: {mat_tvas.shape[0]} proteine x {mat_tvas.shape[1]} campioni")
