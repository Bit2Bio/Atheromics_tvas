#!/usr/bin/env python3
"""
Cross-variable heatmap — how the dm2-associated metabolites (restricted to
the metabolites in the 7 significant enrichment sub_pathways) behave across
all the OTHER clinical variables already analyzed.

Rows: 68 metabolites (Lipid/Amino Acid sub_pathways significant for dm2),
      grouped by super_pathway then sub_pathway.
Columns: every other variable/contrast (dm2 itself excluded — tautological).
         Continuous variables: logFC x SD(X) ("per-SD" effect, comparable
         across variables with different units). SD taken from
         clinical_tvas_transformed.csv, i.e. on whatever scale (raw or log)
         that variable was actually modeled on.
         Categorical contrasts: raw logFC (group-mean difference) — a
         different kind of quantity, kept in a visually separated block.
Cell color: diverging scale centered at 0 (robust-clipped at 1st/99th pct).
Cell marker: a black dot overlay where adj.P.Val < 0.05 for that specific
         metabolite x variable/contrast combination.
Row annotation: two colored strips (super_pathway, sub_pathway).

Excluded entirely: hba1c, fibrinogen, homa_ir, crp (pending data review,
same exclusion as config.yaml) — their results/ folders are stale leftovers
from before that exclusion.

Standalone, one-off script (not wired into Snakemake) — run manually:
    python workflow/scripts/dm2/metabolomics_cross_variable_heatmap.py

Output:
    results/dm2/metabolomics_cross_variable_heatmap.png
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from pathlib import Path

DM2_CSV      = "results/metabolomics/limma/dm2/contrast_G1_vs_G0.csv"
SUB_ENRICH   = "results/dm2/metabolomics_enrichment_sub_pathway.csv"
LIMMA_ROOT   = Path("results/metabolomics/limma")
CLINICAL_TRANSFORMED = "data/processed/clinical_tvas_transformed.csv"
OUT          = Path("results/dm2/metabolomics_cross_variable_heatmap.png")
EXCLUDE_VARS = {"dm2", "hba1c", "fibrinogen", "homa_ir", "crp"}
FDR_SIG      = 0.05

# --------------------------------------------------------------------------
# Row set: metabolites in the significant sub_pathways
# --------------------------------------------------------------------------
sub_enrich = pd.read_csv(SUB_ENRICH)
sig_subs = set(sub_enrich.loc[sub_enrich["FDR"] < FDR_SIG, "sub_pathway"])

dm2 = pd.read_csv(DM2_CSV)
rows = dm2[(dm2["sub_pathway"].isin(sig_subs)) & (dm2["adj.P.Val"] < FDR_SIG)].copy()
rows["feature"] = rows["feature"].astype(str)
rows = rows.sort_values(["super_pathway", "sub_pathway", "plot_name"]).reset_index(drop=True)
print(f"Rows: {len(rows)} metabolites, {rows['super_pathway'].nunique()} super_pathway, "
      f"{rows['sub_pathway'].nunique()} sub_pathway")

feature_ids = rows["feature"].tolist()

# --------------------------------------------------------------------------
# SD lookup (scale each variable was actually modeled on)
# --------------------------------------------------------------------------
clinical = pd.read_csv(CLINICAL_TRANSFORMED, index_col=0)
sd_lookup = clinical.std(numeric_only=True)

# --------------------------------------------------------------------------
# Columns: every other variable/contrast file
# --------------------------------------------------------------------------
col_values = {}   # column label -> Series indexed by feature (logFC, scaled if continuous)
col_sig    = {}   # column label -> Series indexed by feature (bool, adj.P.Val<FDR_SIG)
col_is_categorical = {}

for var_dir in sorted(LIMMA_ROOT.iterdir()):
    if not var_dir.is_dir() or var_dir.name in EXCLUDE_VARS:
        continue
    variable = var_dir.name
    for csv in sorted(var_dir.glob("*.csv")):
        df = pd.read_csv(csv)
        df["feature"] = df["feature"].astype(str)
        df = df.set_index("feature").reindex(feature_ids)

        is_categorical = csv.stem.startswith("contrast_")
        if is_categorical:
            contrast = csv.stem.replace("contrast_", "").replace("_vs_", " vs ")
            label = f"{variable} ({contrast})"
            logfc = df["logFC"]
        else:
            label = variable
            sd = sd_lookup.get(variable, np.nan)
            logfc = df["logFC"] * sd

        col_values[label] = logfc
        col_sig[label]    = df["adj.P.Val"] < FDR_SIG
        col_is_categorical[label] = is_categorical

cont_cols = sorted([c for c, is_cat in col_is_categorical.items() if not is_cat])
cat_cols  = sorted([c for c, is_cat in col_is_categorical.items() if is_cat])
col_order = cont_cols + cat_cols
n_cont = len(cont_cols)

mat = pd.DataFrame(col_values)[col_order].loc[feature_ids]
sig = pd.DataFrame(col_sig)[col_order].loc[feature_ids]
print(f"Columns: {len(col_order)} ({n_cont} continuous, {len(cat_cols)} categorical contrasts)")

# --------------------------------------------------------------------------
# Robust symmetric color scale (clip at 1st/99th percentile)
# --------------------------------------------------------------------------
vals = mat.values.flatten()
vals = vals[~np.isnan(vals)]
vmax = np.percentile(np.abs(vals), 99)
norm = mcolors.Normalize(vmin=-vmax, vmax=vmax)
cmap = plt.colormaps["RdBu_r"]

# --------------------------------------------------------------------------
# Row annotation colors
# --------------------------------------------------------------------------
super_cats = rows["super_pathway"].unique().tolist()
sub_cats   = rows["sub_pathway"].unique().tolist()
super_palette = dict(zip(super_cats, plt.colormaps["Set2"].colors))
sub_palette   = dict(zip(sub_cats, plt.colormaps["tab20"].colors))

# --------------------------------------------------------------------------
# Plot
# --------------------------------------------------------------------------
n_rows, n_cols = mat.shape
fig_w = max(10, n_cols * 0.28 + 4)
fig_h = max(8, n_rows * 0.22 + 1.5)
fig = plt.figure(figsize=(fig_w, fig_h))

ax_super = fig.add_axes([0.02, 0.08, 0.015, 0.82])
ax_sub   = fig.add_axes([0.04, 0.08, 0.015, 0.82])
ax_main  = fig.add_axes([0.16, 0.08, 0.78, 0.82])
ax_cbar  = fig.add_axes([0.96, 0.35, 0.015, 0.3])

for i, (_, row) in enumerate(rows.iterrows()):
    ax_super.add_patch(plt.Rectangle((0, n_rows - i - 1), 1, 1, color=super_palette[row["super_pathway"]]))
    ax_sub.add_patch(plt.Rectangle((0, n_rows - i - 1), 1, 1, color=sub_palette[row["sub_pathway"]]))
for ax in (ax_super, ax_sub):
    ax.set_xlim(0, 1); ax.set_ylim(0, n_rows); ax.axis("off")

im = ax_main.imshow(mat.values, aspect="auto", cmap=cmap, norm=norm, origin="upper")

# significance dots
for r in range(n_rows):
    for c in range(n_cols):
        if sig.values[r, c]:
            ax_main.plot(c, r, marker="o", markersize=2.2, color="black")

# block separator between continuous and categorical
ax_main.axvline(n_cont - 0.5, color="black", linewidth=1.2)

ax_main.set_xticks(range(n_cols))
ax_main.set_xticklabels(col_order, rotation=90, fontsize=6)
ax_main.set_yticks(range(n_rows))
ax_main.set_yticklabels(rows["plot_name"], fontsize=6)
ax_main.set_ylim(n_rows - 0.5, -0.5)
ax_main.set_title("Metabolites (dm2 significant sub_pathways) vs all other variables\n"
                   "continuous: logFC x SD  |  categorical: raw logFC (group diff)  |  "
                   "● = FDR<0.05 for that cell", fontsize=9)

cbar = fig.colorbar(im, cax=ax_cbar)
cbar.set_label("effect size", fontsize=8)

# legends for row annotations
handles_super = [plt.Rectangle((0, 0), 1, 1, color=super_palette[c]) for c in super_cats]
fig.legend(handles_super, super_cats, loc="lower left", bbox_to_anchor=(0.0, 0.91),
           fontsize=6, title="super_pathway", title_fontsize=7, ncol=len(super_cats))

handles_sub = [plt.Rectangle((0, 0), 1, 1, color=sub_palette[c]) for c in sub_cats]
fig.legend(handles_sub, sub_cats, loc="lower left", bbox_to_anchor=(0.0, 0.955),
           fontsize=6, title="sub_pathway", title_fontsize=7, ncol=4)

OUT.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT, dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"Saved: {OUT}")
