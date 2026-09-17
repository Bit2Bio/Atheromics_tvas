#!/usr/bin/env python3
import argparse
import re
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap

parser = argparse.ArgumentParser()
parser.add_argument("--omic",    required=True, help="metabolomics or proteomics")
parser.add_argument("--results", required=True, help="results/{omic}/limma/ directory")
parser.add_argument("--out",     required=True, help="Output PDF path")
parser.add_argument("--fdr",     type=float, default=0.05)
args = parser.parse_args()

RESULTS = Path(args.results)
OUT     = Path(args.out)
FDR     = args.fdr

# --------------------------------------------------------------------------
# Collect results
# --------------------------------------------------------------------------
rows = []
for var_dir in sorted(RESULTS.iterdir()):
    if not var_dir.is_dir():
        continue
    variable = var_dir.name
    for csv in sorted(var_dir.glob("*.csv")):
        df = pd.read_csv(csv)
        if "adj.P.Val" not in df.columns or "logFC" not in df.columns:
            continue
        contrast = csv.stem.replace("contrast_", "")
        sig   = df[df["adj.P.Val"] < FDR]
        n_up   = (sig["logFC"] > 0).sum()
        n_down = (sig["logFC"] < 0).sum()
        label  = variable if contrast == variable else f"{variable}  {contrast}"
        rows.append(dict(label=label, variable=variable,
                         n_up=n_up, n_down=n_down, n_total=n_up + n_down))

data = pd.DataFrame(rows).sort_values(["variable", "label"]).reset_index(drop=True)

# --------------------------------------------------------------------------
# Build matrix: 3 columns
# --------------------------------------------------------------------------
mat   = data[["n_up", "n_down", "n_total"]].values.astype(float)
ylabs = data["label"].tolist()
xlabs = ["Up", "Down", "Total"]

cmaps = [
    LinearSegmentedColormap.from_list("reds",   ["#fff5f0", "#cb181d"]),
    LinearSegmentedColormap.from_list("blues",  ["#f0f4ff", "#2171b5"]),
    LinearSegmentedColormap.from_list("purps",  ["#fcf0ff", "#7a0177"]),
]

n_rows, n_cols = mat.shape
cell_h = 0.42
cell_w = 1.6
fig_h  = max(5, n_rows * cell_h + 1.8)
fig_w  = n_cols * cell_w + 3.2

fig, ax = plt.subplots(figsize=(fig_w, fig_h))
fig.patch.set_facecolor("white")
ax.set_facecolor("white")

# draw cells
for ci in range(n_cols):
    col_vals = mat[:, ci]
    vmax = max(col_vals.max(), 1)
    norm = mcolors.Normalize(vmin=0, vmax=vmax)
    cmap = cmaps[ci]
    for ri, val in enumerate(col_vals):
        color = cmap(norm(val))
        ax.add_patch(plt.Rectangle((ci, ri), 1, 1,
                                   color=color, linewidth=0))
        txt = str(int(val)) if val > 0 else "·"
        brightness = 0.299*color[0] + 0.587*color[1] + 0.114*color[2]
        fc = "white" if brightness < 0.45 else "#222222"
        ax.text(ci + 0.5, ri + 0.5, txt,
                ha="center", va="center",
                fontsize=8.5, fontweight="bold" if val > 0 else "normal",
                color=fc, fontfamily="DejaVu Sans")

# column separators
for ci in range(1, n_cols):
    ax.axvline(ci, color="white", linewidth=2.5, zorder=3)

# row separators — light line every row, thicker between variables
prev_var = None
for ri, row in data.iterrows():
    if prev_var is not None and row["variable"] != prev_var:
        ax.axhline(ri, color="#cccccc", linewidth=1.2, zorder=3)
    else:
        ax.axhline(ri, color="#eeeeee", linewidth=0.5, zorder=2)
    prev_var = row["variable"]
ax.axhline(n_rows, color="#cccccc", linewidth=1.2)
ax.axhline(0,      color="#cccccc", linewidth=1.2)

# axes
ax.set_xlim(0, n_cols)
ax.set_ylim(0, n_rows)
ax.invert_yaxis()

ax.set_xticks([i + 0.5 for i in range(n_cols)])
ax.set_xticklabels(xlabs, fontsize=11, fontweight="bold",
                   fontfamily="DejaVu Sans")
ax.xaxis.set_ticks_position("top")
ax.xaxis.set_label_position("top")
ax.tick_params(axis="x", length=0, pad=6)

ax.set_yticks([i + 0.5 for i in range(n_rows)])
ax.set_yticklabels(ylabs, fontsize=8.5, fontfamily="DejaVu Sans")
ax.tick_params(axis="y", length=0, pad=6)

for spine in ax.spines.values():
    spine.set_visible(False)

omic_label = args.omic.capitalize()
ax.set_title(f"{omic_label}  —  FDR < {FDR}",
             fontsize=13, fontweight="bold", pad=18,
             fontfamily="DejaVu Sans", loc="left")

plt.tight_layout(rect=[0, 0, 1, 0.97])
OUT.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT, bbox_inches="tight", dpi=200)
plt.close(fig)
print(f"Saved: {OUT}")
