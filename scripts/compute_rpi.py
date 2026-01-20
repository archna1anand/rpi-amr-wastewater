#!/usr/bin/env python
"""Compute Resistance Persistence Index (RPI) from a CZID-like AMR table.

This repository is meant to accompany the RPI manuscript and provides a
compact, reproducible reference implementation.

Expected (minimum) input columns (case-insensitive):
- gene_name
- sample_name
- season (wet/dry)
- rpm
- read_species (or any column containing a text annotation that may include
  'plasmid' and pathogen names)

Outputs:
- family_metrics_with_RPI.csv
- top25_genes_RPI.csv
- Figure_RPI_vs_risk.png

Notes
-----
The implementation follows the manuscript description:
- Risk score: Zhang rank mapped to numeric (I=4, II=3, III=2)
- Persistence: 1 - |RPM_wet - RPM_dry|/(RPM_wet + RPM_dry + eps)
- Mobility: 1 if plasmid-annotated, else 0 (family-level is max or mean)
- RPI at family level: risk_score + persistence + mobility

Because CZID outputs may not include explicit Zhang ranks, this script can
infer a simplified Zhang-like rank using pathogen-association and plasmid
annotations from a species/annotation string.
"""

from __future__ import annotations

import argparse
import os
from typing import Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ------------------------
# Helpers
# ------------------------

def _find_col(cols, candidates):
    """Return the first matching column name (case-insensitive)."""
    cols_lower = {c.lower(): c for c in cols}
    for cand in candidates:
        if cand.lower() in cols_lower:
            return cols_lower[cand.lower()]
    return None


def _normalize_season(s: str) -> str:
    if pd.isna(s):
        return "Unknown"
    s = str(s).strip().lower()
    if s.startswith("wet"):
        return "Wet"
    if s.startswith("dry"):
        return "Dry"
    return s.title()


def _mobility_from_text(x: str) -> int:
    """Heuristic: mobility=1 if annotation mentions plasmid."""
    if pd.isna(x):
        return 0
    return 1 if "plasmid" in str(x).lower() else 0


PATHOGEN_KEYWORDS = [
    # ESKAPE / common nosocomial
    "pseudomonas",
    "acinetobacter",
    "klebsiella",
    "enterobacter",
    "staphylococcus aureus",
    "enterococcus",
    # common enterics
    "escherichia",
    "salmonella",
    "shigella",
    "campylobacter",
    "vibrio",
]


def _pathogen_from_text(x: str) -> int:
    if pd.isna(x):
        return 0
    t = str(x).lower()
    return 1 if any(k in t for k in PATHOGEN_KEYWORDS) else 0


def infer_zhang_rank(df: pd.DataFrame, annotation_col: str) -> pd.Series:
    """Infer a simplified Zhang-like rank based on pathogen association and mobility.

    Rank I: pathogen-associated AND mobile
    Rank II: pathogen-associated OR mobile
    Rank III: neither

    Returns: 'I', 'II', 'III'
    """
    pathogen = df[annotation_col].apply(_pathogen_from_text)
    mobile = df[annotation_col].apply(_mobility_from_text)

    rank = np.where(
        (pathogen == 1) & (mobile == 1),
        "I",
        np.where((pathogen == 1) | (mobile == 1), "II", "III"),
    )
    return pd.Series(rank, index=df.index)


def zhang_score_from_rank(rank: pd.Series) -> pd.Series:
    mapping = {"I": 4, "II": 3, "III": 2, "IV": 1}
    return rank.map(mapping).fillna(2).astype(float)


# ------------------------
# Main computation
# ------------------------

def compute_family_metrics(
    df: pd.DataFrame,
    gene_col: str,
    rpm_col: str,
    season_col: str,
    annotation_col: str,
    family_col: Optional[str] = None,
    zhang_rank_col: Optional[str] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Compute family-level RPI and a top-25 gene table."""

    work = df.copy()
    work[season_col] = work[season_col].apply(_normalize_season)

    fam_col = family_col or ("gene_family" if "gene_family" in work.columns else gene_col)

    work["_family"] = work[fam_col].astype(str)
    work["_gene"] = work[gene_col].astype(str)

    # Rank and risk score
    if zhang_rank_col and zhang_rank_col in work.columns:
        rank_raw = work[zhang_rank_col].astype(str)
        rank_raw = rank_raw.str.replace("Rank", "", regex=False).str.strip()
        rank_raw = rank_raw.replace({"1": "I", "2": "II", "3": "III", "4": "IV"})
        work["_zhang_rank"] = rank_raw
    else:
        work["_zhang_rank"] = infer_zhang_rank(work, annotation_col)

    work["_risk_score"] = zhang_score_from_rank(work["_zhang_rank"])
    work["_mobility"] = work[annotation_col].apply(_mobility_from_text)

    # Persistence at family level (mean RPM in Wet vs Dry)
    pivot = (
        work.groupby(["_family", season_col])[rpm_col]
        .mean()
        .unstack(fill_value=0.0)
    )

    wet = pivot.get("Wet", pd.Series(0.0, index=pivot.index))
    dry = pivot.get("Dry", pd.Series(0.0, index=pivot.index))

    eps = 1e-9
    persistence = 1.0 - (np.abs(wet - dry) / (wet + dry + eps))
    persistence = persistence.clip(lower=0.0, upper=1.0)

    # Family risk score: max risk across observations within family
    fam_risk = work.groupby("_family")["_risk_score"].max()

    # Family mobility: mean mobility across observations (0..1)
    fam_mob = work.groupby("_family")["_mobility"].mean()

    family_metrics = pd.DataFrame(
        {
            "family": fam_risk.index,
            "risk_score": fam_risk.values,
            "persistence": persistence.reindex(fam_risk.index).values,
            "mobility": fam_mob.reindex(fam_risk.index).values,
            "wet_mean_rpm": wet.reindex(fam_risk.index).values,
            "dry_mean_rpm": dry.reindex(fam_risk.index).values,
        }
    )

    family_metrics["RPI"] = (
        family_metrics["risk_score"]
        + family_metrics["persistence"]
        + family_metrics["mobility"]
    )

    family_metrics = family_metrics.sort_values("RPI", ascending=False).reset_index(drop=True)

    # Build top 25 gene table from top families
    top_fams = set(family_metrics.head(25)["family"].astype(str))

    top_gene_table = (
        work[work["_family"].isin(top_fams)]
        .groupby(["_gene", "_family"])[rpm_col]
        .mean()
        .reset_index()
        .rename(columns={"_gene": "gene", "_family": "family", rpm_col: "mean_rpm"})
    )

    top_gene_table = top_gene_table.merge(
        family_metrics[["family", "RPI", "risk_score", "persistence", "mobility"]],
        on="family",
        how="left",
    )

    top_gene_table = top_gene_table.sort_values(["RPI", "mean_rpm"], ascending=[False, False]).head(25)

    return family_metrics, top_gene_table


def plot_rpi_vs_risk(family_metrics: pd.DataFrame, outpath: str):
    plt.figure(figsize=(6, 4))
    plt.scatter(family_metrics["risk_score"], family_metrics["RPI"], alpha=0.7)
    plt.xlabel("Risk score")
    plt.ylabel("RPI")
    plt.tight_layout()
    plt.savefig(outpath, dpi=300)
    plt.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="Input CZID AMR CSV")
    ap.add_argument("--outdir", default="outputs", help="Output directory")
    ap.add_argument("--gene-col", default=None, help="Column name for gene (default: auto)")
    ap.add_argument("--rpm-col", default=None, help="Column name for RPM (default: auto)")
    ap.add_argument("--season-col", default=None, help="Column name for season (default: auto)")
    ap.add_argument(
        "--annotation-col",
        default=None,
        help="Column name with read/contig species text (default: auto)",
    )
    ap.add_argument("--family-col", default=None, help="Optional column for ARG family")
    ap.add_argument("--zhang-rank-col", default=None, help="Optional column for Zhang rank")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    df = pd.read_csv(args.input)

    gene_col = args.gene_col or _find_col(df.columns, ["gene_name", "gene"]) or "gene_name"
    rpm_col = args.rpm_col or _find_col(df.columns, ["rpm", "dpm"]) or "rpm"
    season_col = args.season_col or _find_col(df.columns, ["season"]) or "season"
    annotation_col = args.annotation_col or _find_col(df.columns, ["read_species", "contig_species", "species"]) or "read_species"

    family_metrics, top25 = compute_family_metrics(
        df=df,
        gene_col=gene_col,
        rpm_col=rpm_col,
        season_col=season_col,
        annotation_col=annotation_col,
        family_col=args.family_col,
        zhang_rank_col=args.zhang_rank_col,
    )

    family_out = os.path.join(args.outdir, "family_metrics_with_RPI.csv")
    top25_out = os.path.join(args.outdir, "top25_genes_RPI.csv")
    plot_out = os.path.join(args.outdir, "Figure_RPI_vs_risk.png")

    family_metrics.to_csv(family_out, index=False)
    top25.to_csv(top25_out, index=False)
    plot_rpi_vs_risk(family_metrics, plot_out)

    print(f"Wrote: {family_out}")
    print(f"Wrote: {top25_out}")
    print(f"Wrote: {plot_out}")


if __name__ == "__main__":
    main()
