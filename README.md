# Atheromics_v2

Multi-omic analysis pipeline (metabolomics, proteomics) on the Atheromics cohort.

## Structure

```
Atheromics_v2/
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
cd workflow
snakemake --cores 4
```
