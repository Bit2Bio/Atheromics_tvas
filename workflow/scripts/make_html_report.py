#!/usr/bin/env python3
"""
Prototype HTML report — clinical descriptive statistics only.

Standalone script (not wired into the Snakefile yet) — run manually:
    python workflow/scripts/make_html_report.py

Output:
    results/report/report_prototype.html

Fully self-contained: no external network calls (fonts/JS/CSS all inline),
so it opens correctly offline and never sends clinical data anywhere.
"""

import base64
import pandas as pd
from pathlib import Path

CONTINUOUS  = "results/clinics/descriptive/continuous.csv"
CATEGORICAL = "results/clinics/descriptive/categorical.csv"
STRATIFIED  = "results/clinics/descriptive_stratified/descriptive_table_stratified.csv"
QC_GRID     = "results/clinics/qc/continuous_distributions.png"
MET_BOXPLOT = "results/metabolomics/qc/sample_boxplot.png"
MET_PCA     = "results/metabolomics/qc/pca.png"
HEATMAP_METABOLOMICS   = "results/metabolomics/summary_heatmap.png"
HEATMAP_PROTEOMICS     = "results/proteomics/summary_heatmap.png"
HEATMAP_TRANSCRIPTOMICS = "results/transcriptomics/summary_heatmap.png"
OUT = Path("results/report/atheromics_tvas_report.html")


def img_to_data_uri(path: str) -> str:
    data = Path(path).read_bytes()
    return "data:image/png;base64," + base64.b64encode(data).decode("ascii")


def df_to_table(df: pd.DataFrame, table_id: str) -> str:
    thead = "".join(f"<th onclick=\"sortTable('{table_id}',{i})\">{c}</th>" for i, c in enumerate(df.columns))
    rows = []
    for _, r in df.iterrows():
        cells = "".join(f"<td>{'' if pd.isna(v) else v}</td>" for v in r)
        rows.append(f"<tr>{cells}</tr>")
    tbody = "".join(rows)
    return f'<table id="{table_id}"><thead><tr>{thead}</tr></thead><tbody>{tbody}</tbody></table>'


# --------------------------------------------------------------------------
# Load + format
# --------------------------------------------------------------------------
cont = pd.read_csv(CONTINUOUS)
cont[["mean", "sd", "median", "q25", "q75"]] = cont[["mean", "sd", "median", "q25", "q75"]].round(2)

cat = pd.read_csv(CATEGORICAL)
cat["pct"] = cat["pct"].round(1)
cat["category"] = cat["category"].apply(lambda x: str(int(x)) if pd.notna(x) and float(x).is_integer() else x)

strat = pd.read_csv(STRATIFIED)
pair_cols  = [c for c in strat.columns if "_vs_" in c]
group_cols = [c for c in strat.columns if c not in ("label", "p_value", "q_value") + tuple(pair_cols)]
strat = strat[["label"] + group_cols + ["p_value", "q_value"] + pair_cols]
strat = strat.rename(columns={"label": "Variable", "p_value": "p (omnibus)", "q_value": "q (FDR)"})
strat = strat.rename(columns={c: c.replace("p_", "p ").replace("_vs_", " vs ") for c in pair_cols})

# --------------------------------------------------------------------------
# Build HTML
# --------------------------------------------------------------------------
html = f"""<!DOCTYPE html>
<html lang="it">
<head>
<meta charset="utf-8">
<title>Atheromics TVAS — Clinical descriptive statistics</title>
<style>
  :root {{
    --bg: #f7f8fa; --card: #ffffff; --text: #1a1d23; --muted: #6b7280;
    --border: #e5e7eb; --accent: #2563eb; --accent-bg: #eff6ff;
    --stripe: #fafbfc;
  }}
  @media (prefers-color-scheme: dark) {{
    :root {{
      --bg: #0f1115; --card: #171a21; --text: #e6e8eb; --muted: #9aa3af;
      --border: #2a2e37; --accent: #5b9bff; --accent-bg: #14213a;
      --stripe: #1b1e26;
    }}
  }}
  * {{ box-sizing: border-box; }}
  body {{
    margin: 0; padding: 32px 16px; background: var(--bg); color: var(--text);
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
  }}
  .wrap {{ max-width: 1100px; margin: 0 auto; }}
  h1 {{ font-size: 1.6rem; margin-bottom: 4px; }}
  .subtitle {{ color: var(--muted); margin-bottom: 32px; font-size: 0.95rem; }}
  .card {{
    background: var(--card); border: 1px solid var(--border); border-radius: 12px;
    padding: 20px 24px; margin-bottom: 28px; box-shadow: 0 1px 3px rgba(0,0,0,0.04);
  }}
  h2 {{ font-size: 1.15rem; margin-top: 0; margin-bottom: 14px; }}
  table {{ width: 100%; border-collapse: collapse; font-size: 0.85rem; }}
  th, td {{ text-align: left; padding: 7px 10px; border-bottom: 1px solid var(--border); white-space: nowrap; }}
  th {{
    cursor: pointer; user-select: none; color: var(--muted); font-weight: 600;
    position: sticky; top: 0; background: var(--card);
  }}
  th:hover {{ color: var(--accent); }}
  th::after {{ content: " \\21C5"; font-size: 0.7em; opacity: 0.4; }}
  tbody tr:nth-child(even) {{ background: var(--stripe); }}
  tbody tr:hover {{ background: var(--accent-bg); }}
  .table-scroll {{ overflow-x: auto; max-height: 520px; overflow-y: auto; }}
  .meta {{ color: var(--muted); font-size: 0.8rem; margin-top: 8px; }}
</style>
</head>
<body>
<div class="wrap">
  <h1>Atheromics TVAS — Clinical descriptive statistics</h1>
  <div class="subtitle">Prototipo report HTML — clicca sull'intestazione di una colonna per ordinare</div>

  <div class="card">
    <h2>Variabili continue</h2>
    <div class="table-scroll">{df_to_table(cont, "tbl-cont")}</div>
    <div class="meta">{len(cont)} variabili</div>
  </div>

  <div class="card">
    <h2>Variabili categoriche</h2>
    <div class="table-scroll">{df_to_table(cat, "tbl-cat")}</div>
    <div class="meta">{cat['variable'].nunique()} variabili, {len(cat)} righe</div>
  </div>

  <div class="card">
    <h2>Statistiche cliniche stratificate</h2>
    <div class="table-scroll">{df_to_table(strat, "tbl-strat")}</div>
  </div>

  <div class="card">
    <h2>Variabili continue — raw vs log (QC)</h2>
    <img src="{img_to_data_uri(QC_GRID)}" style="width:100%; height:auto; border-radius:8px;">
    <div class="meta">Bordo blu = trasformazione usata in pipeline · linea rossa tratteggiata = outlier |z|&gt;3</div>
  </div>

  <div class="card">
    <h2>Metabolomica — Quality control</h2>
    <img src="{img_to_data_uri(MET_BOXPLOT)}" style="width:100%; height:auto; border-radius:8px;">
    <img src="{img_to_data_uri(MET_PCA)}" style="width:100%; height:auto; border-radius:8px; margin-top:16px;">
    <div class="meta">Boxplot per campione (normalizzazione) + PCA con ellisse di confidenza al 95% (outlier campione etichettati)</div>
  </div>

  <div class="card">
    <h2>Metabolomica — summary heatmap</h2>
    <img src="{img_to_data_uri(HEATMAP_METABOLOMICS)}" style="width:100%; height:auto; border-radius:8px;">
  </div>

  <div class="card">
    <h2>Proteomica — summary heatmap</h2>
    <img src="{img_to_data_uri(HEATMAP_PROTEOMICS)}" style="width:100%; height:auto; border-radius:8px;">
  </div>

  <div class="card">
    <h2>Trascrittomica — summary heatmap</h2>
    <img src="{img_to_data_uri(HEATMAP_TRANSCRIPTOMICS)}" style="width:100%; height:auto; border-radius:8px;">
  </div>
</div>

<script>
function sortTable(tableId, colIdx) {{
  const table = document.getElementById(tableId);
  const tbody = table.tBodies[0];
  const rows = Array.from(tbody.rows);
  const asc = table.dataset.sortCol == colIdx && table.dataset.sortDir !== 'asc';
  rows.sort((a, b) => {{
    let x = a.cells[colIdx].innerText, y = b.cells[colIdx].innerText;
    const fx = parseFloat(x), fy = parseFloat(y);
    if (!isNaN(fx) && !isNaN(fy)) {{ x = fx; y = fy; }}
    if (x < y) return asc ? -1 : 1;
    if (x > y) return asc ? 1 : -1;
    return 0;
  }});
  rows.forEach(r => tbody.appendChild(r));
  table.dataset.sortCol = colIdx;
  table.dataset.sortDir = asc ? 'asc' : 'desc';
}}
</script>
</body>
</html>
"""

OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(html, encoding="utf-8")
print(f"Saved: {OUT}")
