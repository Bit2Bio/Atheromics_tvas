# Atheromics_TVAS

Multi-omic analysis pipeline (metabolomics, proteomics) on the Atheromics TVAS cohort.

## Structure

```
Atheromics_TVAS/
├── data/
│   ├── raw/           ← raw data (not versioned)
│   └── processed/     ← preprocessed data (not versioned)
├── results/           ← analysis outputs (not versioned)
├── config/
│   └── config.yaml
└── workflow/
    ├── Snakefile
    └── scripts/
```

## Run

```bash
cd /home/lorenzo/Atheromics_TVAS
snakemake --cores 4 -s workflow/Snakefile
```
