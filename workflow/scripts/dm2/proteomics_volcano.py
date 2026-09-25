#!/usr/bin/env python3
"""
Volcano plot — proteomics differential features for dm2 (G1 vs G0).

Standalone, one-off script (not wired into Snakemake) — run manually:
    python workflow/scripts/dm2/proteomics_volcano.py

Color: up (red) / down (blue) / not significant (gray), FDR<0.1.
Labels: adjustText-managed, candidate pool = top by combined
        significance+effect-size score, so dense regions thin out on their
        own while isolated top/edge points stay labeled.

Output:
    results/dm2/proteomics_volcano.png
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from adjustText import adjust_text
from pathlib import Path

INPUT       = "results/proteomics/limma/dm2/contrast_G1_vs_G0.csv"
OUT         = Path("results/dm2/proteomics_volcano.png")
FDR_CUTOFF  = 0.1    # shown in title/legend
FDR_DISPLAY = 0.11   # actually used to color/label points (not shown)
N_LABELS    = 35     # candidate pool size before adjustText de-clutters it

df = pd.read_csv(INPUT)
df["neglog10_fdr"] = -np.log10(df["adj.P.Val"].clip(lower=1e-300))

sig_up   = df[(df["adj.P.Val"] < FDR_DISPLAY) & (df["logFC"] > 0)]
sig_down = df[(df["adj.P.Val"] < FDR_DISPLAY) & (df["logFC"] < 0)]
ns       = df[df["adj.P.Val"] >= FDR_DISPLAY]

fig, ax = plt.subplots(figsize=(9, 8))
ax.scatter(ns["logFC"], ns["neglog10_fdr"], s=10, color="#BDC3C7", alpha=0.6, linewidths=0, label="n.s.")
ax.scatter(sig_up["logFC"], sig_up["neglog10_fdr"], s=16, color="#E74C3C", alpha=0.8, linewidths=0, label="up in G1")
ax.scatter(sig_down["logFC"], sig_down["neglog10_fdr"], s=16, color="#2E86C1", alpha=0.8, linewidths=0, label="down in G1")

ax.axhline(-np.log10(FDR_CUTOFF), color="gray", linestyle="--", linewidth=0.8)
ax.axvline(0, color="gray", linestyle="-", linewidth=0.6)

# --------------------------------------------------------------------------
# Label candidates: reward high significance AND large effect size
# (naturally favors "top" and "edges" of the volcano) — adjustText then
# repels overlapping labels, so crowded regions self-thin.
# --------------------------------------------------------------------------
sig_all = pd.concat([sig_up, sig_down])
sig_all["score"] = sig_all["neglog10_fdr"] * sig_all["logFC"].abs()
candidates = sig_all.sort_values("score", ascending=False).head(N_LABELS)

texts = [
    ax.text(row["logFC"], row["neglog10_fdr"], row["protein"], fontsize=7.5)
    for _, row in candidates.iterrows()
]
adjust_text(
    texts, ax=ax,
    arrowprops=dict(arrowstyle="-", color="#666666", lw=0.5),
    expand=(1.3, 1.5),
)

ax.set_xlabel("log2 Fold Change (G1 vs G0)")
ax.set_ylabel("-log10(FDR)")
ax.set_title(f"Proteomics — dm2 (G1 vs G0)\n{len(sig_up)} up, {len(sig_down)} down (FDR<{FDR_CUTOFF})")
ax.legend(frameon=False, loc="upper left", markerscale=1.5)
for spine in ("top", "right"):
    ax.spines[spine].set_visible(False)

OUT.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT, dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"Saved: {OUT}")
