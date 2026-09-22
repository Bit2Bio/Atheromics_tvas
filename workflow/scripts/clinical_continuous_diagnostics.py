#!/usr/bin/env python3
"""
Diagnostic grid for continuous clinical variables: raw vs log distribution,
side by side, with outliers labeled by sample ID.

Standalone QC script (not wired into the Snakefile) — run manually:
    python workflow/scripts/clinical_continuous_diagnostics.py

Output:
    results/clinics/qc/continuous_distributions.png
"""

import numpy as np
import pandas as pd
import yaml
import matplotlib.pyplot as plt
from pathlib import Path

CLINICAL = "data/processed/clinical_tvas.csv"
CONFIG   = "config/config.yaml"
OUT      = Path("results/clinics/qc/continuous_distributions.png")

OUTLIER_Z = 3.0

with open(CONFIG) as f:
    cfg = yaml.safe_load(f)

cont_vars = [v for v in cfg["variables"] if v["type"] == "continuous"]
meta = pd.read_csv(CLINICAL, index_col=0)

n_vars = len(cont_vars)
n_cols_vars = 3           # variables per row
n_rows = int(np.ceil(n_vars / n_cols_vars))
fig, axes = plt.subplots(n_rows, n_cols_vars * 2,
                          figsize=(n_cols_vars * 2 * 3.2, n_rows * 2.6))
axes = np.atleast_2d(axes)


def plot_panel(ax, values, title, active):
    ax.hist(values, bins=20, color="#4C72B0" if active else "#B0B0B0",
            edgecolor="white", linewidth=0.5)
    z = (values - values.mean()) / values.std()
    outliers = z[z.abs() > OUTLIER_Z]
    for idx in outliers.index:
        ax.axvline(values.loc[idx], color="#E74C3C", linewidth=1, linestyle="--")
        ax.text(values.loc[idx], ax.get_ylim()[1] * 0.95, str(idx),
                rotation=90, fontsize=6, color="#E74C3C",
                ha="right", va="top")
    border_color = "#4C72B0" if active else "#cccccc"
    for spine in ax.spines.values():
        spine.set_edgecolor(border_color)
        spine.set_linewidth(2.2 if active else 0.8)
    ax.set_title(title, fontsize=8, fontweight="bold" if active else "normal")
    ax.tick_params(labelsize=6)


for i, v in enumerate(cont_vars):
    row = i // n_cols_vars
    block = i % n_cols_vars
    col_raw, col_log = block * 2, block * 2 + 1

    name = v["name"]
    transform = v.get("transform", "raw")
    s = meta[name].dropna()

    plot_panel(axes[row, col_raw], s, f"{name}\nraw", active=(transform == "raw"))

    s_log = np.log(s)
    plot_panel(axes[row, col_log], s_log, f"{name}\nlog", active=(transform == "log"))

# hide unused panels if n_vars is not a multiple of n_cols_vars
for i in range(n_vars, n_rows * n_cols_vars):
    row = i // n_cols_vars
    block = i % n_cols_vars
    axes[row, block * 2].axis("off")
    axes[row, block * 2 + 1].axis("off")

fig.suptitle("Continuous clinical variables — raw vs log (blue border = current pipeline choice, "
             f"red dashed = |z|>{OUTLIER_Z:.0f} outlier)",
             fontsize=11, fontweight="bold")
plt.tight_layout(rect=[0, 0, 1, 0.97])
OUT.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT, dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"Saved: {OUT}")
