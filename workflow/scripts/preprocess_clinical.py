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

EXCEL  = "data/raw/metabolon e glicemia 090926_rev_LA.xlsx"
SHEET  = "TVAS"
OUT_DIR = Path("data/processed")

VARIABLES = [
    # (nome_originale, nome_armonizzato, tipo)
    ("Età",                        "age",              "continuous"),
    ("sesso",                      "gender",           "categorical"),
    ("ADA con IGM (IFG+IGT)",      "ada_igm",          "categorical"),
    ("Familiarità diabete",        "familiarity_dm",   "categorical"),
    ("Familiarità CVD",            "familiarity_cvd",  "categorical"),
    ("Fumo (0=no 1=si 2=ex)",      "smoking",          "categorical"),
    ("ATS carotidea",              "ats_carotidea",    "categorical"),
    ("sintomi Ats carot",          "ats_symptoms",     "categorical"),
    ("Peso",                       "weight",           "continuous"),
    ("Altezza",                    "height",           "continuous"),
    ("BMI",                        "bmi",              "continuous"),
    ("vita",                       "waist",            "continuous"),
    ("PA syst",                    "pa_syst",          "continuous"),
    ("PA diast",                   "pa_diast",         "continuous"),
    ("colesterolo tot",            "cholesterol",      "continuous"),
    ("HDL",                        "hdl",              "continuous"),
    ("LDL",                        "ldl",              "continuous"),
    ("trigliceridi",               "triglycerides",    "continuous"),
    ("AST",                        "ast",              "continuous"),
    ("ALT",                        "alt",              "continuous"),
    ("azotemia",                   "azotemia",         "continuous"),
    ("creatinina",                 "creatinine",       "continuous"),
    ("GFR",                        "gfr",              "continuous"),
    ("FPG 0' mg/dl",               "fpg",              "continuous"),
    ("120' mg/dl",                 "gluc_120",         "continuous"),
    ("AUC OGTT GLIC",              "auc_ogtt_gluc",    "continuous"),
    ("Insulinemia 0'",             "insulinemia",      "continuous"),
    ("120'",                       "ins_120",          "continuous"),
    ("AUC OGTT INS ",              "auc_ogtt_ins",     "continuous"),
    ("HbA1c",                      "hba1c",            "continuous"),
    ("IGF1",                       "igf1",             "continuous"),
    ("Peptide C",                  "peptide_c",        "continuous"),
    ("hsPCR",                      "hscrp",            "continuous"),
    ("PCR",                        "crp",              "continuous"),
    ("fibrinogeno",                "fibrinogen",       "continuous"),
    ("WBC",                        "wbc",              "continuous"),
    ("placca stabile/instabile.1", "plaque_type",      "categorical"),
    ("HOMA IR",                    "homa_ir",          "continuous"),
    ("Tyg Index",                  "tyg_index",        "continuous"),
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
