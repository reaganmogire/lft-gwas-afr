#!/usr/bin/env python3
"""
SuSiE fine-mapping (susieR::susie_rss) from GWAS summary statistics and an LD matrix computed by PLINK.

Inputs:
  - GWAS summary statistics (tab-delimited; can be .gz) containing at least:
      SNP, CHR, POS, Allele1, Allele2, BETA, SE, P.value
  - PLINK LD reference panel prefix (bfile; without .bed/.bim/.fam)
  - Target region (chr + lead position + window) or explicit start/end
Outputs:
  - Region-extracted summary statistics TSV
  - SNP list used for LD computation
  - LD matrix (r) from PLINK
  - SuSiE outputs: PIP per variant, credible set assignment, and a small summary JSON/TSV

This script is designed as a cohort-agnostic template; cohort-specific paths and inputs are supplied at runtime.
"""

import argparse
import gzip
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# rpy2 imports
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
import rpy2.robjects.numpy2ri as numpy2ri

numpy2ri.activate()


REQ_COLS = ["SNP", "CHR", "POS", "Allele1", "Allele2", "BETA", "SE", "P.value"]


def die(msg: str) -> None:
    print(f"[ERROR] {msg}", file=sys.stderr)
    sys.exit(1)


def read_sumstats(path: str) -> pd.DataFrame:
    # pandas can read .gz directly, but we keep explicit for clarity.
    df = pd.read_csv(path, sep="\t", dtype={"P.value": "string"})
    missing = [c for c in REQ_COLS if c not in df.columns]
    if missing:
        die(f"Missing required columns in sumstats: {missing}. Found: {list(df.columns)}")
    # Coerce types
    df["CHR"] = pd.to_numeric(df["CHR"], errors="coerce").astype("Int64")
    df["POS"] = pd.to_numeric(df["POS"], errors="coerce").astype("Int64")
    df["BETA"] = pd.to_numeric(df["BETA"], errors="coerce")
    df["SE"] = pd.to_numeric(df["SE"], errors="coerce")
    df["P.value"] = pd.to_numeric(df["P.value"], errors="coerce")
    df = df.dropna(subset=["CHR", "POS", "BETA", "SE", "P.value"])
    return df


def run_cmd(cmd, desc: str) -> None:
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        die(f"{desc} failed (exit={e.returncode}). Command: {' '.join(cmd)}")


def main():
    ap = argparse.ArgumentParser(
        description="Cohort-agnostic SuSiE RSS fine-mapping using PLINK LD + susieR.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("--sumstats", required=True, help="GWAS summary stats TSV (.gz ok)")
    ap.add_argument("--bfile", required=True, help="PLINK LD reference prefix (no .bed/.bim/.fam)")
    ap.add_argument("--chr", required=True, type=int, help="Chromosome")
    ap.add_argument("--pos", required=True, type=int, help="Lead position (bp)")
    ap.add_argument("--window", type=int, default=250000, help="Window around --pos (bp)")
    ap.add_argument("--n", required=True, type=int, help="Sample size for susie_rss")
    ap.add_argument("--plink", default="plink", help="Path to PLINK executable")
    ap.add_argument("--L", type=int, default=10, help="Maximum number of SuSiE effects (L)")
    ap.add_argument("--coverage", type=float, default=0.95, help="Credible set coverage")
    ap.add_argument("--min-abs-corr", type=float, default=0.5, help="susie_get_cs min_abs_corr")
    ap.add_argument("--outdir", default="results/susie", help="Output directory root")
    args = ap.parse_args()

    outdir = Path(args.outdir) / f"chr{args.chr}_{args.pos}"
    outdir.mkdir(parents=True, exist_ok=True)

    # 1) Load sumstats
    df = read_sumstats(args.sumstats)

    start = args.pos - args.window
    end = args.pos + args.window

    locus = df[(df["CHR"] == args.chr) & (df["POS"] > start) & (df["POS"] < end)].copy()
    if locus.empty:
        die(f"No variants found in region chr{args.chr}:{start}-{end}")

    # Ensure consistent ordering for downstream joins
    locus = locus.sort_values(["POS", "SNP"]).reset_index(drop=True)

    # Write region TSV + SNP list
    region_tsv = outdir / "region_sumstats.tsv"
    snplist_all = outdir / "snplist.all"
    locus.to_csv(region_tsv, sep="\t", index=False)
    locus["SNP"].to_csv(snplist_all, sep="\t", index=False, header=False)

    # 2) Compute LD matrix (r) using PLINK: --r square + --extract SNP list
    ld_prefix = outdir / "ld"
    ld_matrix_file = Path(str(ld_prefix) + ".ld")

    plink_cmd = [
        args.plink,
        "--bfile", args.bfile,
        "--keep-allele-order",
        "--r", "square",
        "--write-snplist",
        "--extract", str(snplist_all),
        "--out", str(ld_prefix),
    ]
    run_cmd(plink_cmd, "PLINK LD computation")

    # PLINK writes:
    #   ld.ld (matrix) and ld.snplist (filtered list of SNPs that survived)
    snplist_final = outdir / "ld.snplist"
    if not snplist_final.exists():
        die(f"Expected PLINK snplist not found: {snplist_final}")

    # 3) Align locus to PLINK SNP order
    kept = pd.read_csv(snplist_final, header=None, names=["SNP"])
    locus2 = locus[locus["SNP"].isin(kept["SNP"])].copy()
    locus2 = locus2.set_index("SNP").reindex(kept["SNP"]).reset_index()
    if locus2.empty:
        die("After PLINK filtering, no variants remain for SuSiE.")

    # 4) Read LD matrix
    if not ld_matrix_file.exists():
        die(f"Expected LD matrix not found: {ld_matrix_file}")
    R = pd.read_csv(ld_matrix_file, sep="\t", header=None).values
    if R.shape[0] != R.shape[1] or R.shape[0] != locus2.shape[0]:
        die(f"Dimension mismatch: LD={R.shape} vs variants={locus2.shape[0]}")

    # 5) Run SuSiE (susieR::susie_rss)
    susieR = importr("susieR")
    bhat = locus2["BETA"].to_numpy().reshape((-1, 1))
    shat = locus2["SE"].to_numpy().reshape((-1, 1))

    try:
        fit = susieR.susie_rss(
            bhat=bhat,
            shat=shat,
            R=R,
            L=args.L,
            n=args.n
        )
    except Exception as e:
        die(f"SuSiE susie_rss failed: {e}")

    converged = bool(fit.rx2("converged")[0])
    pip = np.array(susieR.susie_get_pip(fit))

    cs = susieR.susie_get_cs(
        fit,
        coverage=args.coverage,
        min_abs_corr=args.min_abs_corr,
        Xcorr=R
    )[0]

    # Convert credible sets (indices are 1-based in R)
    cs_assign = np.zeros(locus2.shape[0], dtype=int)
    try:
        # cs is a list of integer vectors (R indices, 1..m)
        for i in range(len(cs)):
            idx = np.array(cs[i], dtype=int) - 1
            cs_assign[idx] = i + 1
    except Exception:
        # If no CS returned, keep zeros
        pass

    # 6) Save results
    out_tsv = outdir / "susie_results.tsv"
    out_json = outdir / "susie_summary.json"

    locus2["PIP"] = pip
    locus2["CS"] = cs_assign
    locus2.to_csv(out_tsv, sep="\t", index=False)

    summary = {
        "chromosome": args.chr,
        "lead_position": args.pos,
        "window_bp": args.window,
        "region": f"chr{args.chr}:{start}-{end}",
        "n_variants_input": int(locus.shape[0]),
        "n_variants_final": int(locus2.shape[0]),
        "susie_L": args.L,
        "coverage": args.coverage,
        "min_abs_corr": args.min_abs_corr,
        "converged": converged,
        "n_credible_sets": int(np.max(cs_assign)) if locus2.shape[0] > 0 else 0,
    }
    out_json.write_text(json.dumps(summary, indent=2))

    print("[OK] SuSiE fine-mapping complete")
    print(f"[OK] Results TSV: {out_tsv}")
    print(f"[OK] Summary JSON: {out_json}")


if __name__ == "__main__":
    main()

