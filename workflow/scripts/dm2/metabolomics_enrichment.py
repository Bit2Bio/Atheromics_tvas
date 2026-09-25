#!/usr/bin/env python3
"""
Pathway over-representation analysis (ORA) — metabolomics DEM for dm2 (G1 vs G0).

Fisher exact test (one-sided, "greater") on super_pathway/sub_pathway
annotations, same logic as Atheromics' metabolomics_enrichment_ora.R
(run_pathway_enrichment), ported to Python. Pathway columns are already
joined into the dm2 contrast CSV (see metabolomics_limma.R), so no separate
chemical_metadata.csv load is needed here.

Standalone, one-off script (not wired into Snakemake) — run manually:
    python workflow/scripts/dm2/metabolomics_enrichment.py

Output:
    results/dm2/metabolomics_enrichment_super_pathway.csv
    results/dm2/metabolomics_enrichment_sub_pathway.csv
"""
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
from pathlib import Path

INPUT           = "results/metabolomics/limma/dm2/contrast_G1_vs_G0.csv"
OUT_DIR         = Path("results/dm2")
FDR_DEM         = 0.05   # threshold defining a "differentially expressed metabolite"
MIN_PATHWAY_SIZE = 3


def run_pathway_enrichment(metabolite_set, universe, chem_meta,
                            pathway_col="super_pathway", min_pathway_size=3):
    """Fisher exact ORA of `metabolite_set` within `universe`, per pathway."""
    chem_meta = chem_meta.copy()
    chem_meta["feature"] = chem_meta["feature"].astype(str)
    metabolite_set = {str(x) for x in metabolite_set}
    universe       = {str(x) for x in universe}

    annot = chem_meta[["feature", "plot_name", pathway_col]].dropna(subset=[pathway_col])
    annot = annot[annot[pathway_col] != ""]
    annot = annot[annot["feature"].isin(universe)]

    universe_annot = annot["feature"]
    set_annot      = metabolite_set & set(universe_annot)

    N = len(universe_annot)
    n = len(set_annot)
    if n == 0:
        return pd.DataFrame()

    pathway_counts = annot.groupby(pathway_col).size().reset_index(name="K")

    overlap = annot[annot["feature"].isin(set_annot)]
    overlap_counts = (
        overlap.groupby(pathway_col)["plot_name"]
        .agg(k="size", metabolites=lambda s: ",".join(s))
        .reset_index()
    )

    res = pathway_counts.merge(overlap_counts, on=pathway_col, how="left")
    res["k"] = res["k"].fillna(0).astype(int)
    res["metabolites"] = res["metabolites"].fillna("")
    res = res[res["K"] >= min_pathway_size].reset_index(drop=True)
    if len(res) == 0:
        return res

    p_values, expected, enrichment = [], [], []
    for _, row in res.iterrows():
        k, K = row["k"], row["K"]
        table = [[k, n - k], [K - k, N - K - (n - k)]]
        _, p = fisher_exact(table, alternative="greater")
        exp = n * (K / N)
        p_values.append(p)
        expected.append(exp)
        enrichment.append(k / exp if exp > 0 else np.nan)

    res["p_value"]   = p_values
    res["expected"]  = expected
    res["enrichment"] = enrichment
    res["FDR"] = multipletests(res["p_value"], method="fdr_bh")[1]
    res = res.sort_values("p_value").reset_index(drop=True)
    return res


if __name__ == "__main__":
    df = pd.read_csv(INPUT)
    chem_meta = df[["feature", "plot_name", "super_pathway", "sub_pathway"]]

    universe = df["feature"].unique()
    dem      = df.loc[df["adj.P.Val"] < FDR_DEM, "feature"].unique()
    print(f"Universe: {len(universe)}  |  DEM (FDR<{FDR_DEM}): {len(dem)}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    for pathway_col in ["super_pathway", "sub_pathway"]:
        res = run_pathway_enrichment(dem, universe, chem_meta,
                                      pathway_col=pathway_col,
                                      min_pathway_size=MIN_PATHWAY_SIZE)
        out_path = OUT_DIR / f"metabolomics_enrichment_{pathway_col}.csv"
        res.to_csv(out_path, index=False)
        n_sig = (res["FDR"] < 0.05).sum() if len(res) else 0
        print(f"{pathway_col}: {len(res)} pathways tested, {n_sig} FDR<0.05  ->  {out_path}")
