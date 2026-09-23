#!/usr/bin/env python3
"""
Consolidate all per-variable DESeq2 CSVs (transcriptomics) into a single
Excel workbook, one sheet per variable (or variable+contrast for categorical
variables with multiple pairwise contrasts).

Same logic as consolidate_excel.py, adapted for DESeq2's file naming
("{variable}_results.csv" for continuous variables instead of
"contrast_*.csv") and column names (FDR/PValue instead of adj.P.Val/P.Value).

Each sheet is capped: rows with FDR <= FDR_CAP, or the top N_CAP rows by
PValue, whichever is more (union) — a null result on ~79k genes still
shows its top hits, and a strong result isn't truncated below its
significant set.

Usage:
    python rnaseq_consolidate_excel.py --results results/transcriptomics/deseq2/ --out results/transcriptomics/transcriptomics_deseq2_results.xlsx
"""
import argparse
from pathlib import Path
import pandas as pd

FDR_CAP = 0.25
N_CAP   = 3000

parser = argparse.ArgumentParser()
parser.add_argument("--results", required=True, help="results/transcriptomics/deseq2/ directory")
parser.add_argument("--out",     required=True, help="Output .xlsx path")
args = parser.parse_args()

RESULTS = Path(args.results)
OUT     = Path(args.out)

OUT.parent.mkdir(parents=True, exist_ok=True)

seen_names = set()
n_sheets = 0
with pd.ExcelWriter(OUT, engine="openpyxl") as writer:
    for var_dir in sorted(RESULTS.iterdir()):
        if not var_dir.is_dir():
            continue
        variable = var_dir.name
        for csv in sorted(var_dir.glob("*.csv")):
            contrast = csv.stem.replace("contrast_", "").replace(f"{variable}_results", variable)
            sheet_name = variable if contrast == variable else f"{variable}_{contrast}"
            sheet_name = sheet_name[:31]
            if sheet_name in seen_names:
                sheet_name = f"{sheet_name[:28]}_{len(seen_names)}"
            seen_names.add(sheet_name)

            df = pd.read_csv(csv).sort_values("PValue")
            n_sig = (df["FDR"] <= FDR_CAP).sum()
            df = df.head(max(n_sig, N_CAP))
            df.to_excel(writer, sheet_name=sheet_name, index=False)
            n_sheets += 1

print(f"Saved: {OUT} ({n_sheets} sheets)")
