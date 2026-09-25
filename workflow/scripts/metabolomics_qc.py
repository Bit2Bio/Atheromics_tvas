#!/usr/bin/env python3
"""
Preliminary QC for metabolomics data (163 samples x 1738 features).

Per Metabolon's methodology: raw peak areas -> batch normalization ->
imputation of missing/below-detection-limit values -> natural log transform
("Log Transformed Data" sheet). No further transformation applied here.

1. Per-sample boxplot — all 163 samples side by side, to visually confirm
   normalization worked (similar median/IQR across samples) and show the
   scale the data is compressed into.
2. PCA — features z-scored, PC1 vs PC2, 95% confidence ellipse (bivariate
   normal), samples outside it flagged and labeled by ID.

Output:
    results/metabolomics/qc/sample_boxplot.png
    results/metabolomics/qc/pca.png
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

INPUT       = "data/processed/metabolomics_tvas.csv"
OUT_DIR     = Path("results/metabolomics/qc")
CONF        = 0.95   # PCA confidence ellipse level

met = pd.read_csv(INPUT, index_col=0)
feat = met.drop(columns=["PARENT_SAMPLE_NAME"])
print(f"Matrix: {feat.shape[0]} samples x {feat.shape[1]} features")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# --------------------------------------------------------------------------
# 1. Per-sample boxplot
# --------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(16, 5))
ax.boxplot(feat.values.T, showfliers=False, widths=0.7,
           patch_artist=True,
           boxprops=dict(facecolor="#4C72B0", alpha=0.6, linewidth=0.3),
           medianprops=dict(color="#E74C3C", linewidth=0.8),
           whiskerprops=dict(linewidth=0.3), capprops=dict(linewidth=0.3))
ax.set_xticks(range(1, feat.shape[0] + 1))
ax.set_xticklabels(feat.index, rotation=90, fontsize=4)
ax.set_xlabel(f"Samples ({feat.shape[0]}, original file order)")
ax.set_ylabel("Batch-normalized, imputed, log-transformed abundance (Metabolon)")
ax.set_title("Metabolomics — per-sample distribution")
for spine in ("top", "right"):
    ax.spines[spine].set_visible(False)
fig.savefig(OUT_DIR / "sample_boxplot.png", dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"Saved: {OUT_DIR / 'sample_boxplot.png'}")

# --------------------------------------------------------------------------
# 2. PCA
# --------------------------------------------------------------------------
X = StandardScaler().fit_transform(feat.values)
pca = PCA(n_components=5)
scores = pca.fit_transform(X)
var_exp = pca.explained_variance_ratio_ * 100

pc_df = pd.DataFrame(scores[:, :2], columns=["PC1", "PC2"], index=feat.index)

mean = pc_df.mean().values
cov = np.cov(pc_df.values, rowvar=False)
eigvals, eigvecs = np.linalg.eigh(cov)
order = eigvals.argsort()[::-1]
eigvals, eigvecs = eigvals[order], eigvecs[:, order]

chi2_val = chi2.ppf(CONF, df=2)
width, height = 2 * np.sqrt(eigvals * chi2_val)
angle = np.degrees(np.arctan2(*eigvecs[:, 0][::-1]))

inv_cov = np.linalg.inv(cov)
diffs = pc_df.values - mean
maha2 = np.einsum("ij,jk,ik->i", diffs, inv_cov, diffs)
pc_df["outlier"] = maha2 > chi2_val

n_out = pc_df["outlier"].sum()
print(f"PC1: {var_exp[0]:.1f}% var  |  PC2: {var_exp[1]:.1f}% var")
print(f"Samples outside {CONF:.0%} ellipse: {n_out}")
if n_out:
    print(pc_df[pc_df["outlier"]][["PC1", "PC2"]].to_string())

fig, ax = plt.subplots(figsize=(8, 7))
ax.add_patch(Ellipse(mean, width, height, angle=angle, facecolor="none",
                      edgecolor="gray", linestyle="--", linewidth=1))
inside = pc_df[~pc_df["outlier"]]
outside = pc_df[pc_df["outlier"]]
ax.scatter(inside["PC1"], inside["PC2"], s=25, color="#4C72B0", alpha=0.7, label="in range")
ax.scatter(outside["PC1"], outside["PC2"], s=35, color="#E74C3C", label="outlier")
texts = [ax.text(row["PC1"], row["PC2"], sid, fontsize=8, color="#E74C3C")
         for sid, row in outside.iterrows()]
adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="#999999", lw=0.5))
ax.axhline(0, color="#dddddd", linewidth=0.6, zorder=0)
ax.axvline(0, color="#dddddd", linewidth=0.6, zorder=0)
ax.set_xlabel(f"PC1 ({var_exp[0]:.1f}%)")
ax.set_ylabel(f"PC2 ({var_exp[1]:.1f}%)")
ax.set_title(f"Metabolomics PCA — sample QC\n{feat.shape[0]} samples, {CONF:.0%} confidence ellipse, {n_out} outlier(s)")
ax.legend(frameon=False)
for spine in ("top", "right"):
    ax.spines[spine].set_visible(False)
fig.savefig(OUT_DIR / "pca.png", dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"Saved: {OUT_DIR / 'pca.png'}")
