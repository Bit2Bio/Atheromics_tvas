#!/usr/bin/env python3
"""
Dotplot of significant sub_pathway enrichment for dm2 metabolomics.

Single panel, all significant sub_pathways together, sorted by fold
enrichment. Sub_pathway name (y-tick label) colored by its parent
super_pathway.

Standalone, one-off script (not wired into Snakemake) — run manually:
    python workflow/scripts/dm2/metabolomics_enrichment_plot.py

Output:
    results/dm2/metabolomics_enrichment_plot.png
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

SUB_CSV = "results/dm2/metabolomics_enrichment_sub_pathway.csv"
DM2_CSV = "results/metabolomics/limma/dm2/contrast_G1_vs_G0.csv"
OUT     = Path("results/dm2/metabolomics_enrichment_plot.png")
FDR_THR = 0.05

enrich_sub = pd.read_csv(SUB_CSV)

sub_sig = enrich_sub[(enrich_sub["FDR"] < FDR_THR) & (enrich_sub["k"] > 0)].copy()
if sub_sig.empty:
    raise SystemExit("No significant sub_pathway — nothing to plot.")

# sub_pathway -> super_pathway lookup, taken from the annotated dm2 results
sub_to_super = (
    pd.read_csv(DM2_CSV)[["super_pathway", "sub_pathway"]]
    .dropna()
    .drop_duplicates(subset="sub_pathway")
)

plot_df = sub_sig.merge(sub_to_super, on="sub_pathway", how="left")
plot_df["neglog10_fdr"] = -np.log10(plot_df["FDR"])
plot_df = plot_df.sort_values("enrichment").reset_index(drop=True)

# --------------------------------------------------------------------------
# One color per super_pathway, used for the sub_pathway name text
# --------------------------------------------------------------------------
super_cats = plot_df["super_pathway"].unique().tolist()
palette = dict(zip(super_cats, plt.colormaps["Dark2"].colors))

fig, ax = plt.subplots(figsize=(7, max(2.5, len(plot_df) * 0.5 + 1)))

vmax = plot_df["neglog10_fdr"].max()
norm = plt.Normalize(vmin=0, vmax=vmax)
cmap = plt.colormaps["Blues"]
size_scale = lambda k: 80 + 500 * (k / plot_df["k"].max())

y_pos = np.arange(len(plot_df))
sc = ax.scatter(
    plot_df["enrichment"], y_pos,
    s=size_scale(plot_df["k"]), c=plot_df["neglog10_fdr"], cmap=cmap, norm=norm,
    edgecolors="#333333", linewidths=0.5, zorder=3,
)
ax.set_yticks(y_pos)
ax.set_yticklabels(plot_df["sub_pathway"], fontsize=9)
for tick, super_name in zip(ax.get_yticklabels(), plot_df["super_pathway"]):
    tick.set_color(palette[super_name])
    tick.set_fontweight("bold")

ax.grid(axis="x", color="#eeeeee", linewidth=0.8, zorder=0)
for spine in ("top", "right"):
    ax.spines[spine].set_visible(False)
ax.set_ylim(-0.7, len(plot_df) - 0.3)
ax.set_xlabel("Fold enrichment")
ax.set_title("Metabolomics pathway enrichment — dm2 (G1 vs G0)", fontsize=12, fontweight="bold")

cbar = fig.colorbar(sc, ax=ax, pad=0.02)
cbar.set_label("-log10(FDR)")

# size legend
for k_val in sorted(plot_df["k"].unique()):
    ax.scatter([], [], s=size_scale(k_val), color="gray", edgecolors="#333333",
               linewidths=0.5, label=f"k={k_val}")
size_legend = ax.legend(title="Overlap", frameon=False, loc="lower right", fontsize=7, title_fontsize=8)
ax.add_artist(size_legend)

# super_pathway color legend
handles = [plt.Line2D([0], [0], marker="s", linestyle="", color=palette[s], markersize=8)
           for s in super_cats]
ax.legend(handles, super_cats, title="super_pathway", frameon=False,
          loc="upper left", bbox_to_anchor=(0, 1), fontsize=7, title_fontsize=8)

OUT.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT, dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"Saved: {OUT}")
