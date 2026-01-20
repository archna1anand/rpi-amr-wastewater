# RPI Toolkit (Resistance Persistence Index)

A small, reproducible toolkit to compute the **Resistance Persistence Index (RPI)** for antimicrobial resistance (AMR) features in wastewater.

RPI is designed as a **surveillance prioritization** score that integrates:

- **Risk** (e.g., Zhang rank / pathogen association)
- **Persistence** across time/conditions (e.g., wet vs dry seasonal stability)
- **Mobility** (e.g., plasmid-associated annotation)

> This repo is intentionally lightweight so it’s easy to run, fork, and adapt.

---

## Quick start (recommended)

### Option A — Conda
```bash
conda env create -f environment.yml
conda activate rpi
```

### Option B — Pip
```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

---

## Run on the included sample data

```bash
python scripts/compute_rpi.py \
  --input data/czid_amr_sample.csv \
  --outdir outputs
```

This writes:
- `outputs/family_metrics_with_RPI.csv`
- `outputs/top25_genes_RPI.csv`
- `outputs/figure_RPI_vs_risk.png`

---

## Run on your full CZID export

Replace `--input` with your full CSV.

```bash
python scripts/compute_rpi.py --input /path/to/AA_Oct2024combined_amr_results.csv --outdir outputs
```

---

## Expected input columns

The script is flexible. It will try to infer column names, but it works best if your table includes:

- `sample` (or `sample_name`)
- `gene_name`
- `rpm`
- `season` (values like `Wet` / `Dry`)
- `species` (free text; used for pathogen keywords)

If you already have risk ranks (e.g., `zhang_rank`), the script will use them.

---

## What to cite

If you use this in a paper, include:
- Your manuscript citation (preferred)
- This repository (CITATION.cff)

---

## License

MIT — see `LICENSE`.
