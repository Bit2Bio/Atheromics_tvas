#!/usr/bin/env python3
"""
Dotplot of significant sub_pathway enrichment for dm2 metabolomics,
faceted by super_pathway — matplotlib port of Atheromics'
metabolomics_enrichment_ora.R::make_enrichment_plot().

Only sub_pathways whose parent super_pathway is ALSO significant are
shown (same filter as the original), so this reads the two CSVs already
produced by metabolomics_enrichment.py.

Standalone, one-off script (not wired into Snakemake) — run manually:
    python workflow/scripts/dm2/metabolomics_enrichment_plot.py

Output:
    results/dm2/metabolomics_enrichment_plot.png
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

SUPER_CSV = "results/dm2/metabolomics_enrichment_super_pathway.csv"
SUB_CSV   = "results/dm2/metabolomics_enrichment_sub_pathway.csv"
DM2_CSV   = "results/metabolomics/limma/dm2/contrast_G1_vs_G0.csv"
OUT       = Path("results/dm2/metabolomics_enrichment_plot.png")
FDR_THR   = 0.05

enrich_super = pd.read_csv(SUPER_CSV)
enrich_sub   = pd.read_csv(SUB_CSV)

sub_sig = enrich_sub[(enrich_sub["FDR"] < FDR_THR) & (enrich_sub["k"] > 0)].copy()
if sub_sig.empty:
    raise SystemExit("No significant sub_pathway — nothing to plot.")

super_sig_names = set(enrich_super.loc[enrich_super["FDR"] < FDR_THR, "super_pathway"])

# sub_pathway -> super_pathway lookup, taken from the annotated dm2 results
sub_to_super = (
    pd.read_csv(DM2_CSV)[["super_pathway", "sub_pathway"]]
    .dropna()
    .drop_duplicates(subset="sub_pathway")
)

plot_df = sub_sig.merge(sub_to_super, on="sub_pathway", how="left")
plot_df = plot_df[plot_df["super_pathway"].isin(super_sig_names)].copy()
if plot_df.empty:
    raise SystemExit("No sub_pathway left after filtering to significant super_pathway — nothing to plot.")

plot_df["neglog10_fdr"] = -np.log10(plot_df["FDR"])

# --------------------------------------------------------------------------
# One row of panels per super_pathway (facet), sub_pathway ordered by
# enrichment within each panel — same layout as facet_grid(super ~ .)
# --------------------------------------------------------------------------
super_groups = (
    plot_df.groupby("super_pathway")["enrichment"].max()
    .sort_values(ascending=False).index.tolist()
)
panel_heights = [max(1, plot_df[plot_df["super_pathway"] == s].shape[0]) for s in super_groups]

fig, axes = plt.subplots(
    nrows=len(super_groups), ncols=1, sharex=True,
    figsize=(7, max(2.5, sum(panel_heights) * 0.5 + 1)),
    gridspec_kw=dict(height_ratios=panel_heights),
)
axes = np.atleast_1d(axes)

vmax = plot_df["neglog10_fdr"].max()
norm = plt.Normalize(vmin=0, vmax=vmax)
cmap = plt.colormaps["Blues"]
size_scale = lambda k: 80 + 500 * (k / plot_df["k"].max())

sc = None
for ax, super_name in zip(axes, super_groups):
    sub = plot_df[plot_df["super_pathway"] == super_name].sort_values("enrichment")
    y_pos = np.arange(len(sub))
    sc = ax.scatter(
        sub["enrichment"], y_pos,
        s=size_scale(sub["k"]), c=sub["neglog10_fdr"], cmap=cmap, norm=norm,
        edgecolors="#333333", linewidths=0.5, zorder=3,
    )
    ax.set_yticks(y_pos)
    ax.set_yticklabels(sub["sub_pathway"], fontsize=9)
    ax.set_ylabel(super_name, fontsize=9, fontweight="bold", rotation=0, ha="right", va="center")
    ax.grid(axis="x", color="#eeeeee", linewidth=0.8, zorder=0)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.set_ylim(-0.7, len(sub) - 0.3)

axes[-1].set_xlabel("Fold enrichment")
fig.suptitle("Metabolomics pathway enrichment — dm2 (G1 vs G0)", fontsize=12, fontweight="bold")

cbar = fig.colorbar(sc, ax=axes, location="right", pad=0.02, shrink=0.6)
cbar.set_label("-log10(FDR)")

# size legend
for k_val in sorted(plot_df["k"].unique()):
    axes[0].scatter([], [], s=size_scale(k_val), color="gray", edgecolors="#333333",
                     linewidths=0.5, label=f"k={k_val}")
axes[0].legend(title="Overlap", frameon=False, loc="lower right", fontsize=7, title_fontsize=8)

OUT.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT, dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"Saved: {OUT}")
