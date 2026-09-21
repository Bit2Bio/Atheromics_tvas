#!/usr/bin/env python3
"""
Preprocessing transcriptomics (RNA-seq) for TVAS samples.

Steps:
  1. Load harmonized, QC-filtered featureCounts matrix (genes x samples, raw counts)
  2. Filter to TVAS_ samples only
  3. Align sample order to metabolomics (subset — not all TVAS samples have RNA-seq)

Output:
  data/processed/transcriptomics_tvas.csv
"""

import pandas as pd
from pathlib import Path

COUNTS_RAW = "data/raw/transcriptomics/featureCounts_napoli_harmonized_qcfilt.tsv"
MET_FILE   = "data/processed/metabolomics_tvas.csv"
OUT_DIR    = Path("data/processed")

# --------------------------------------------------------------------------
# Load
# --------------------------------------------------------------------------
counts = pd.read_csv(COUNTS_RAW, sep="\t", index_col=0)

# --------------------------------------------------------------------------
# Filter TVAS
# --------------------------------------------------------------------------
tvas_cols = [c for c in counts.columns if c.startswith("TVAS_")]
counts_tvas = counts[tvas_cols]
print(f"After TVAS filter: {counts_tvas.shape[0]} geni x {counts_tvas.shape[1]} campioni")

# --------------------------------------------------------------------------
# Align sample order to metabolomics (subset: not all TVAS have RNA-seq)
# --------------------------------------------------------------------------
met_ids = pd.read_csv(MET_FILE, index_col=0).index
common_ids = [s for s in met_ids if s in counts_tvas.columns]
missing = set(met_ids) - set(counts_tvas.columns)
print(f"Campioni TVAS in metabolomica senza RNA-seq: {len(missing)}")

counts_tvas = counts_tvas[common_ids]

# --------------------------------------------------------------------------
# Save
# --------------------------------------------------------------------------
OUT_DIR.mkdir(parents=True, exist_ok=True)
counts_tvas.to_csv(OUT_DIR / "transcriptomics_tvas.csv")

print(f"transcriptomics_tvas.csv: {counts_tvas.shape[0]} geni x {counts_tvas.shape[1]} campioni")
