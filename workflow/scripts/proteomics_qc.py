#!/usr/bin/env python3
"""
Quality control for proteomics data (163 samples x 814 proteins, Olink NPX).

1. Per-sample boxplot — all 163 samples side by side, to visually confirm
   normalization worked (similar median/IQR across samples) and show the
   scale the data is compressed into.
2. PCA — features z-scored, PC1 vs PC2, 95% confidence ellipse (bivariate
   normal), samples outside it flagged and labeled by ID. Points colored by
   Olink PlateID (3 plates) — the main technical batch axis specific to
   Olink data, absent from the metabolomics (single-batch) run.

Output:
    results/proteomics/qc/sample_boxplot.png
    results/proteomics/qc/pca.png
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.stats import chi2
from adjustText import adjust_text
from pathlib import Path

INPUT       = "data/processed/proteomics_tvas.csv"
OLINK_RAW   = "data/raw/proteomics/FedericiAllPlates_Extended_NPX_2026-03-11.csv"
OUT_DIR     = Path("results/proteomics/qc")
CONF        = 0.95   # PCA confidence ellipse level

# proteins x samples -> samples x proteins
prot = pd.read_csv(INPUT, index_col=0).T
print(f"Matrix: {prot.shape[0]} samples x {prot.shape[1]} proteins")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# --------------------------------------------------------------------------
# Plate assignment (Olink batch)
# --------------------------------------------------------------------------
raw = pd.read_csv(OLINK_RAW, sep=";", usecols=["SampleID", "SampleType", "AssayType", "PlateID"])
plate = (
    raw[(raw["SampleType"] == "SAMPLE") & (raw["AssayType"] == "assay")]
    [["SampleID", "PlateID"]].drop_duplicates().set_index("SampleID")["PlateID"]
)
plate = plate.reindex(prot.index)
plate = plate.str.replace("_Layout_Federici_Reveal", "", regex=False)  # Plate1/2/3
print(plate.value_counts())

# --------------------------------------------------------------------------
# 1. Per-sample boxplot
# --------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(16, 5))
ax.boxplot(prot.values.T, showfliers=False, widths=0.7,
           patch_artist=True,
           boxprops=dict(facecolor="#4C72B0", alpha=0.6, linewidth=0.3),
           medianprops=dict(color="#E74C3C", linewidth=0.8),
           whiskerprops=dict(linewidth=0.3), capprops=dict(linewidth=0.3))
ax.set_xticks(range(1, prot.shape[0] + 1))
ax.set_xticklabels(prot.index, rotation=90, fontsize=4)
ax.set_xlabel(f"Samples ({prot.shape[0]}, original file order)")
ax.set_ylabel("NPX (Olink normalized protein expression, log2 scale)")
ax.set_title("Proteomics — per-sample distribution")
for spine in ("top", "right"):
    ax.spines[spine].set_visible(False)
fig.savefig(OUT_DIR / "sample_boxplot.png", dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"Saved: {OUT_DIR / 'sample_boxplot.png'}")

# --------------------------------------------------------------------------
# 2. PCA
# --------------------------------------------------------------------------
X = StandardScaler().fit_transform(prot.values)
pca = PCA(n_components=5)
scores = pca.fit_transform(X)
var_exp = pca.explained_variance_ratio_ * 100

pc_df = pd.DataFrame(scores[:, :2], columns=["PC1", "PC2"], index=prot.index)
pc_df["plate"] = plate.values

mean = pc_df[["PC1", "PC2"]].mean().values
cov = np.cov(pc_df[["PC1", "PC2"]].values, rowvar=False)
eigvals, eigvecs = np.linalg.eigh(cov)
order = eigvals.argsort()[::-1]
eigvals, eigvecs = eigvals[order], eigvecs[:, order]

chi2_val = chi2.ppf(CONF, df=2)
width, height = 2 * np.sqrt(eigvals * chi2_val)
angle = np.degrees(np.arctan2(*eigvecs[:, 0][::-1]))

inv_cov = np.linalg.inv(cov)
diffs = pc_df[["PC1", "PC2"]].values - mean
maha2 = np.einsum("ij,jk,ik->i", diffs, inv_cov, diffs)
pc_df["outlier"] = maha2 > chi2_val

n_out = pc_df["outlier"].sum()
print(f"PC1: {var_exp[0]:.1f}% var  |  PC2: {var_exp[1]:.1f}% var")
print(f"Samples outside {CONF:.0%} ellipse: {n_out}")
if n_out:
    print(pc_df[pc_df["outlier"]][["PC1", "PC2", "plate"]].to_string())

fig, ax = plt.subplots(figsize=(8, 7))
ax.add_patch(Ellipse(mean, width, height, angle=angle, facecolor="none",
                      edgecolor="gray", linestyle="--", linewidth=1))

plate_colors = dict(zip(sorted(pc_df["plate"].unique()), ["#4C72B0", "#DD8452", "#55A868"]))
for plate_name, color in plate_colors.items():
    sub = pc_df[pc_df["plate"] == plate_name]
    ax.scatter(sub["PC1"], sub["PC2"], s=25, color=color, alpha=0.75, label=plate_name)

outside = pc_df[pc_df["outlier"]]
ax.scatter(outside["PC1"], outside["PC2"], s=90, facecolors="none",
           edgecolors="black", linewidths=1.2, label="outlier", zorder=4)
texts = [ax.text(row["PC1"], row["PC2"], sid, fontsize=8, color="black")
         for sid, row in outside.iterrows()]
adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="#999999", lw=0.5))

ax.axhline(0, color="#dddddd", linewidth=0.6, zorder=0)
ax.axvline(0, color="#dddddd", linewidth=0.6, zorder=0)
ax.set_xlabel(f"PC1 ({var_exp[0]:.1f}%)")
ax.set_ylabel(f"PC2 ({var_exp[1]:.1f}%)")
ax.set_title(f"Proteomics PCA — sample QC\n{prot.shape[0]} samples, {CONF:.0%} confidence ellipse, {n_out} outlier(s), colored by Olink plate")
ax.legend(frameon=False)
for spine in ("top", "right"):
    ax.spines[spine].set_visible(False)
fig.savefig(OUT_DIR / "pca.png", dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"Saved: {OUT_DIR / 'pca.png'}")
