#!/usr/bin/env python3
"""
Clinical descriptive statistics — stratified by a grouping variable (Table 1).

For each clinical variable, computes per-group summary and overall p-value:
  - Continuous  : mean ± SD per group  | ANOVA omnibus p-value, Tukey HSD pairwise
  - Binary (0/1): n (%) of value=1     | Chi-square p-value
  - Multi-level : n (%) per level      | Chi-square p-value

Percentages are computed over non-missing observations of each variable
(not over the total group size).

Omnibus p-values (one per variable) are corrected for multiple testing with
Benjamini-Hochberg FDR, reported in the q_value column. Tukey HSD pairwise
p-values are left untouched since they are already multiple-comparison
adjusted.

The stratification variable and any variables listed in --exclude are omitted.

Output:
  - descriptive_table_stratified.csv   tidy table

Usage:
    python clinical_descriptive_stratified.py \
        --input       data/processed/clinical_tvas.csv \
        --guide       data/processed/clinical_guide.csv \
        --strat_var   plaque_type \
        --group_labels "0:Stable,1:Unstable" \
        --exclude     plaque_type \
        --out_dir     results/clinics/descriptive_stratified/
"""

import argparse
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input",        required=True)
    p.add_argument("--guide",        required=True)
    p.add_argument("--strat_var",    required=True)
    p.add_argument("--group_labels", default=None,
                   help="Label map e.g. '0:NGM,1:IGM,2:T2D'")
    p.add_argument("--exclude",      default="",
                   help="Comma-separated variables to exclude from rows")
    p.add_argument("--out_dir",      required=True)
    return p.parse_args()


def parse_labels(s):
    if not s:
        return None
    out = {}
    for tok in s.split(","):
        k, v = tok.split(":")
        out[float(k)] = v.strip()
    return out


def fmt_pval(p):
    if pd.isna(p):
        return ""
    if p < 0.001:
        return "<0.001"
    return f"{p:.3f}"


def continuous_stats(series, groups, group_vals):
    """ANOVA omnibus + Tukey HSD pairwise p-values across all groups at once.

    Tukey HSD adjusts jointly for the number of groups being compared, so it
    must be run once per variable with all groups together (not pair by pair).
    Returns (p_omnibus, {(i, j): p_pairwise}) with i, j indexing group_vals.
    """
    grp_data, valid_idx = [], []
    for i, g in enumerate(group_vals):
        s = series[groups == g].dropna()
        if len(s) >= 2:
            grp_data.append(s.values)
            valid_idx.append(i)

    if len(grp_data) < 2:
        return np.nan, {}

    _, p_omnibus = stats.f_oneway(*grp_data)

    pw = {}
    tukey = stats.tukey_hsd(*grp_data)
    for a, i in enumerate(valid_idx):
        for b, j in enumerate(valid_idx):
            if a < b:
                pw[(i, j)] = tukey.pvalue[a, b]

    return p_omnibus, pw


def pairwise_chisq(series, groups, g1, g2):
    """Chi-square between two groups for a categorical variable."""
    sub = series[groups.isin([g1, g2])]
    sub_g = groups[groups.isin([g1, g2])]
    ct = pd.crosstab(sub, sub_g)
    if ct.shape[0] < 2 or ct.shape[1] < 2:
        return np.nan
    _, p, _, _ = stats.chi2_contingency(ct.values)
    return p


def build_table(df, guide, strat_var, group_vals, label_map, exclude):
    var_types = dict(zip(guide["harmonized_name"], guide["type"]))
    skip = set([strat_var] + [v.strip() for v in exclude.split(",") if v.strip()])
    variables = [v for v in guide["harmonized_name"] if v in df.columns and v not in skip]

    group_labels = [label_map.get(g, str(int(g))) for g in group_vals] if label_map \
                   else [str(int(g)) for g in group_vals]
    groups = df[strat_var]
    group_ns = {g: (groups == g).sum() for g in group_vals}

    pair_idxs  = list(combinations(range(len(group_vals)), 2))
    pair_keys  = [f"p_{group_labels[i]}_vs_{group_labels[j]}" for i, j in pair_idxs]

    rows = []

    # Header row: N per group
    header_row = {"label": "N", "p_value": "", "is_header": True}
    for g, gl in zip(group_vals, group_labels):
        header_row[gl] = str(group_ns[g])
    for pk in pair_keys:
        header_row[pk] = ""
    rows.append(header_row)

    for v in variables:
        vtype = var_types.get(v, "categorical")
        series = df[v]
        n_missing = series.isna().sum()
        n_total   = series.notna().sum()

        if vtype == "continuous":
            row = {"label": v, "is_header": False}
            for g, gl in zip(group_vals, group_labels):
                s = series[groups == g].dropna()
                row[gl] = f"{s.mean():.1f} ± {s.std():.1f}" if len(s) > 0 else ""

            p_omnibus, pw_idx = continuous_stats(series, groups, group_vals)
            row["p_value"] = fmt_pval(p_omnibus)
            row["p_value_raw"] = p_omnibus
            row.update({pk: fmt_pval(pw_idx.get((i, j), np.nan))
                        for pk, (i, j) in zip(pair_keys, pair_idxs)})
            if n_missing > 0:
                row["label"] += f" (N={n_total})"
            rows.append(row)

        else:
            levels = sorted(series.dropna().unique())
            is_binary = (len(levels) == 2 and set(levels) <= {0.0, 1.0})

            # Chi-square on full contingency table
            ct = pd.crosstab(series, groups)
            if ct.shape[0] >= 2 and ct.shape[1] >= 2:
                _, p, _, _ = stats.chi2_contingency(ct.values)
            else:
                p = np.nan

            # Pairwise chi-square (same for all rows of a multi-level variable)
            pw = {pk: fmt_pval(pairwise_chisq(series, groups, group_vals[i], group_vals[j]))
                  for pk, (i, j) in zip(pair_keys, pair_idxs)}

            if is_binary:
                lv = 1.0
                row = {"label": f"{v}, n (%)", "is_header": False}
                for g, gl in zip(group_vals, group_labels):
                    sub_series = series[groups == g]
                    n_valid = sub_series.notna().sum()
                    n_pos   = (sub_series == lv).sum()
                    row[gl] = f"{n_pos} ({n_pos/n_valid*100:.1f}%)" if n_valid > 0 else ""
                row["p_value"] = fmt_pval(p)
                row["p_value_raw"] = p
                row.update(pw)
                if n_missing > 0:
                    row["label"] += f" (N={n_total})"
                rows.append(row)
            else:
                for i, lv in enumerate(levels):
                    try:
                        lv_str = str(int(lv)) if float(lv).is_integer() else str(lv)
                    except (ValueError, TypeError):
                        lv_str = str(lv)
                    label = v if i == 0 else f"  {lv_str}"
                    row   = {"label": label, "is_header": False}
                    for g, gl in zip(group_vals, group_labels):
                        sub_series = series[groups == g]
                        n_valid = sub_series.notna().sum()
                        n_lv    = (sub_series == lv).sum()
                        row[gl] = f"{n_lv} ({n_lv/n_valid*100:.1f}%)" if n_valid > 0 else ""
                    row["p_value"] = fmt_pval(p) if i == 0 else ""
                    if i == 0:
                        row["p_value_raw"] = p
                    row.update({pk: (pv if i == 0 else "") for pk, pv in pw.items()})
                    if i == 0 and n_missing > 0:
                        row["label"] += f" (N={n_total})"
                    rows.append(row)

    # Benjamini-Hochberg FDR correction on the omnibus p-values (one per variable).
    # Tukey HSD pairwise p-values are left as-is: already multiple-comparison adjusted.
    raw_idx, raw_pvals = [], []
    for idx, row in enumerate(rows):
        p_raw = row.get("p_value_raw", np.nan)
        if pd.notna(p_raw):
            raw_idx.append(idx)
            raw_pvals.append(p_raw)

    if raw_pvals:
        _, q_vals, _, _ = multipletests(raw_pvals, method="fdr_bh")
        for idx, q in zip(raw_idx, q_vals):
            rows[idx]["q_value"] = fmt_pval(q)

    for row in rows:
        row.pop("p_value_raw", None)
        row.setdefault("q_value", "")

    return rows, group_labels, pair_keys


def main():
    args     = parse_args()
    out_dir  = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df         = pd.read_csv(args.input, index_col=0)
    guide      = pd.read_csv(args.guide)
    label_map  = parse_labels(args.group_labels)
    group_vals = sorted(df[args.strat_var].dropna().unique())

    print(f"Stratification: {args.strat_var}")
    for g in group_vals:
        lbl = label_map.get(g, str(g)) if label_map else str(g)
        print(f"  {lbl} (n={(df[args.strat_var] == g).sum()})")

    rows, group_labels, pair_keys = build_table(
        df, guide, args.strat_var, group_vals, label_map, args.exclude
    )

    # CSV
    csv_path = out_dir / "descriptive_table_stratified.csv"
    pd.DataFrame(rows).drop(columns=["is_header"], errors="ignore").to_csv(csv_path, index=False)
    print(f"\nSaved: {csv_path}")


if __name__ == "__main__":
    main()
