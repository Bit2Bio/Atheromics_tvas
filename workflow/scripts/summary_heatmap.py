#!/usr/bin/env python3
import argparse
import math
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap

parser = argparse.ArgumentParser()
parser.add_argument("--omic",    required=True, help="metabolomics or proteomics")
parser.add_argument("--results", required=True, help="results/{omic}/limma/ directory")
parser.add_argument("--out",     required=True, help="Output PNG path")
parser.add_argument("--fdr",     type=float, default=0.05)
parser.add_argument("--rows_per_panel", type=int, default=1000,
                     help="Wrap rows into side-by-side panels of this many rows max (default: effectively never — single tall panel)")
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

mat   = data[["n_up", "n_down", "n_total"]].values.astype(float)
ylabs = data["label"].tolist()
xlabs = ["Up", "Down", "Total"]

cmaps = [
    LinearSegmentedColormap.from_list("reds",   ["#fff5f0", "#cb181d"]),
    LinearSegmentedColormap.from_list("blues",  ["#f0f4ff", "#2171b5"]),
    LinearSegmentedColormap.from_list("purps",  ["#fcf0ff", "#7a0177"]),
]

n_rows, n_cols = mat.shape

# --------------------------------------------------------------------------
# Wrap rows into side-by-side panels for a slide-friendly (landscape) layout
# --------------------------------------------------------------------------
n_panels = max(1, math.ceil(n_rows / args.rows_per_panel))
row_chunks = np.array_split(np.arange(n_rows), n_panels)

cell_h = 0.25
cell_w = 1.6
panel_pad_w = 3.2   # space for row labels + margins per panel
panel_gap   = 0.6   # horizontal gap between panels

max_panel_rows = max(len(c) for c in row_chunks)
fig_h = max(5, max_panel_rows * cell_h + 1.8)
fig_w = n_panels * (n_cols * cell_w + panel_pad_w) + (n_panels - 1) * panel_gap

fig, axes = plt.subplots(ncols=n_panels, figsize=(fig_w, fig_h))
fig.patch.set_facecolor("white")
if n_panels == 1:
    axes = [axes]

# global color scale per column, shared across panels
col_vmax = [max(mat[:, ci].max(), 1) for ci in range(n_cols)]

for pi, idx in enumerate(row_chunks):
    ax = axes[pi]
    ax.set_facecolor("white")
    sub_mat   = mat[idx]
    sub_ylabs = [ylabs[i] for i in idx]
    sub_vars  = [data["variable"].iloc[i] for i in idx]
    n_sub     = len(idx)

    for ci in range(n_cols):
        norm = mcolors.Normalize(vmin=0, vmax=col_vmax[ci])
        cmap = cmaps[ci]
        for ri, val in enumerate(sub_mat[:, ci]):
            color = cmap(norm(val))
            ax.add_patch(plt.Rectangle((ci, ri), 1, 1, color=color, linewidth=0))
            txt = str(int(val)) if val > 0 else "·"
            brightness = 0.299 * color[0] + 0.587 * color[1] + 0.114 * color[2]
            fc = "white" if brightness < 0.45 else "#222222"
            ax.text(ci + 0.5, ri + 0.5, txt,
                    ha="center", va="center",
                    fontsize=8.5, fontweight="bold" if val > 0 else "normal",
                    color=fc, fontfamily="DejaVu Sans")

    for ci in range(1, n_cols):
        ax.axvline(ci, color="white", linewidth=2.5, zorder=3)

    prev_var = None
    for ri, var in enumerate(sub_vars):
        if prev_var is not None and var != prev_var:
            ax.axhline(ri, color="#cccccc", linewidth=1.2, zorder=3)
        else:
            ax.axhline(ri, color="#eeeeee", linewidth=0.5, zorder=2)
        prev_var = var
    ax.axhline(n_sub, color="#cccccc", linewidth=1.2)
    ax.axhline(0,     color="#cccccc", linewidth=1.2)

    ax.set_xlim(0, n_cols)
    ax.set_ylim(0, max_panel_rows)
    ax.invert_yaxis()

    ax.set_xticks([i + 0.5 for i in range(n_cols)])
    ax.set_xticklabels(xlabs, fontsize=11, fontweight="bold", fontfamily="DejaVu Sans")
    ax.xaxis.set_ticks_position("top")
    ax.xaxis.set_label_position("top")
    ax.tick_params(axis="x", length=0, pad=6)

    ax.set_yticks([i + 0.5 for i in range(n_sub)])
    ax.set_yticklabels(sub_ylabs, fontsize=8.5, fontfamily="DejaVu Sans")
    ax.tick_params(axis="y", length=0, pad=6)

    for spine in ax.spines.values():
        spine.set_visible(False)

omic_label = args.omic.capitalize()
fig.suptitle(f"{omic_label}  —  FDR < {FDR}",
             fontsize=13, fontweight="bold", x=0.01, ha="left",
             fontfamily="DejaVu Sans", y=0.99)

plt.tight_layout(rect=[0, 0, 1, 0.96], w_pad=3)
OUT.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT, bbox_inches="tight", dpi=200)
plt.close(fig)
print(f"Saved: {OUT}")
