#!/usr/bin/env python3
"""
SuSiE fine-mapping (susieR::susie_rss) from GWAS summary statistics and an
LD matrix computed by PLINK.

Key design decisions (see inline comments for rationale):
  - GWASLab is used for allele harmonization against an Ensembl-style reference
    FASTA before LD computation.
  - PLINK --a1-allele ensures the LD r-matrix is signed consistently with the
    harmonized effect allele, avoiding sign mismatches that break SuSiE.
  - estimate_residual_variance=FALSE because the LD reference panel is
    out-of-sample relative to the summary statistics.
  - estimate_prior_variance=TRUE so that SuSiE prunes empty components
    internally; effective K is read from the number of returned credible sets.

Inputs:
  - GWAS summary statistics (tab-delimited; .gz ok) with columns:
      SNP, CHR, POS, Allele1, Allele2, BETA, SE, P.value
  - PLINK LD reference panel prefix (bfile)
  - Ensembl-style reference genome FASTA (chromosome names match plain integers)
  - Target region: chromosome + lead position + window

Outputs  (all written under --outdir/chr{CHR}_{POS}/):
  - region_sumstats.tsv          harmonized region summary statistics
  - snplist.all                  full SNP list sent to PLINK
  - ld.ld / ld_r2.ld             r and r² LD matrices
  - ld.snplist                   SNPs retained after PLINK QC
  - ld_heatmaps.png              r and r² heatmaps
  - susie_results.tsv            per-variant PIP and CS assignment
  - credible_set_variants.tsv    variants inside credible sets only
  - credible_set_plot.pdf        locus plot (-log10P and PIP panels)
  - credible_set_summary.pdf     bar charts of CS sizes and max PIPs
  - susie_summary.json           run metadata and key results
"""

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

import gwaslab as gl
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
import rpy2.robjects.numpy2ri as numpy2ri

numpy2ri.activate()

REQ_COLS = ["SNP", "CHR", "POS", "Allele1", "Allele2", "BETA", "SE", "P.value"]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def die(msg: str) -> None:
    print(f"[ERROR] {msg}", file=sys.stderr)
    sys.exit(1)


def run_cmd(cmd: list, desc: str) -> None:
    try:
        subprocess.run(cmd, check=True)
        print(f"[OK] {desc}")
    except subprocess.CalledProcessError as e:
        die(f"{desc} failed (exit={e.returncode}). Command: {' '.join(str(c) for c in cmd)}")


def read_sumstats(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype={"P.value": "string"})
    missing = [c for c in REQ_COLS if c not in df.columns]
    if missing:
        die(f"Missing required columns: {missing}. Found: {list(df.columns)}")
    df["CHR"]     = pd.to_numeric(df["CHR"],     errors="coerce").astype("Int64")
    df["POS"]     = pd.to_numeric(df["POS"],     errors="coerce").astype("Int64")
    df["BETA"]    = pd.to_numeric(df["BETA"],    errors="coerce")
    df["SE"]      = pd.to_numeric(df["SE"],      errors="coerce")
    df["P.value"] = pd.to_numeric(df["P.value"], errors="coerce")
    return df.dropna(subset=["CHR", "POS", "BETA", "SE", "P.value"])


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Cohort-agnostic SuSiE RSS fine-mapping: PLINK LD + GWASLab harmonization + susieR.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("--sumstats",    required=True,  help="GWAS summary stats TSV (.gz ok)")
    ap.add_argument("--bfile",       required=True,  help="PLINK LD reference prefix (no .bed/.bim/.fam)")
    ap.add_argument("--ref-fasta",   required=True,  help="Ensembl-style GRCh38 FASTA (plain integer chr names)")
    ap.add_argument("--chr",         required=True,  type=int, help="Chromosome")
    ap.add_argument("--pos",         required=True,  type=int, help="Lead position (bp, GRCh38)")
    ap.add_argument("--n",           required=True,  type=int, help="GWAS sample size for susie_rss")
    ap.add_argument("--window",      type=int,   default=250000, help="Window ± --pos (bp)")
    ap.add_argument("--build",       default="38",   help="Genome build passed to GWASLab (37 or 38)")
    ap.add_argument("--plink",       default="plink",help="Path to PLINK 1.9 executable")
    ap.add_argument("--L",           type=int,   default=10, help="Max SuSiE effects (L)")
    ap.add_argument("--coverage",    type=float, default=0.95, help="Credible set coverage")
    ap.add_argument("--min-abs-corr",type=float, default=0.5,  help="susie_get_cs min_abs_corr")
    ap.add_argument("--lambda-reg",  type=float, default=1e-4,
                    help="Diagonal regularization added to R (avoids non-PSD issues)")
    ap.add_argument("--outdir",      default="results/susie", help="Output directory root")
    args = ap.parse_args()

    outdir = Path(args.outdir) / f"chr{args.chr}_{args.pos}"
    outdir.mkdir(parents=True, exist_ok=True)

    start = args.pos - args.window
    end   = args.pos + args.window
    pfx   = f"chr{args.chr}_{args.pos}"

    # -----------------------------------------------------------------------
    # 1. Load sumstats and filter to region using GWASLab
    # -----------------------------------------------------------------------
    print(f"\n[1/7] Loading sumstats and filtering to chr{args.chr}:{start}-{end}")
    raw = read_sumstats(args.sumstats)

    mysumstats = gl.Sumstats(
        raw,
        snpid="SNP", chrom="CHR", pos="POS",
        ea="Allele1", nea="Allele2",
        beta="BETA", se="SE", p="P.value",
        build=args.build,
    )
    mysumstats.data["N"] = args.n
    mysumstats.basic_check()

    locus = mysumstats.filter_value(
        f"CHR=={args.chr} & POS>{start} & POS<{end}"
    )
    if locus.data.empty:
        die(f"No variants found in region chr{args.chr}:{start}-{end}")

    # -----------------------------------------------------------------------
    # 2. Harmonize alleles against Ensembl reference FASTA
    # -----------------------------------------------------------------------
    # GWASLab aligns NEA to the reference REF allele and flips EA/BETA where
    # necessary. Using an Ensembl-style FASTA (plain integer chromosome names)
    # ensures chromosome name matching with standard sumstats.
    print(f"\n[2/7] Harmonizing alleles against reference FASTA: {args.ref_fasta}")
    locus.harmonize(basic_check=False, ref_seq=args.ref_fasta)

    # Save harmonized region TSV + full SNP list
    region_tsv  = outdir / "region_sumstats.tsv"
    snplist_all = outdir / "snplist.all"
    a1_file     = outdir / f"{pfx}.a1alleles"

    locus.data.to_csv(region_tsv, sep="\t", index=False)
    locus.data["SNPID"].to_csv(snplist_all, sep="\t", index=False, header=False)

    # Write allele file for PLINK --a1-allele (col1=SNPID, col2=EA).
    # This aligns the LD r-matrix A1 with the harmonized effect allele so that
    # the sign of r(i,j) is consistent with BETA — critical for susie_rss.
    locus.data[["SNPID", "EA"]].to_csv(a1_file, sep="\t", index=False, header=False)
    print(f"[OK] Harmonized A1-allele file: {a1_file}")

    # -----------------------------------------------------------------------
    # 3. Compute LD matrices with PLINK
    # -----------------------------------------------------------------------
    print(f"\n[3/7] Computing LD matrices with PLINK")
    ld_prefix    = outdir / "ld"
    ld_r2_prefix = outdir / "ld_r2"

    for rtype, out_pfx in [("r", ld_prefix), ("r2", ld_r2_prefix)]:
        run_cmd([
            args.plink,
            "--bfile",      args.bfile,
            # --a1-allele: col 1 = variant ID, col 2 = effect allele.
            # Ensures R is signed consistently with harmonized BETA.
            "--a1-allele",  str(a1_file), "2", "1",
            f"--{rtype}",   "square",
            "--write-snplist",
            "--extract",    str(snplist_all),
            "--out",        str(out_pfx),
        ], f"PLINK --{rtype}")

    # -----------------------------------------------------------------------
    # 4. Align summary stats to PLINK SNP order and read LD
    # -----------------------------------------------------------------------
    print(f"\n[4/7] Aligning summary statistics to PLINK SNP order")
    snplist_final = outdir / "ld.snplist"
    if not snplist_final.exists():
        die(f"Expected PLINK snplist not found: {snplist_final}")

    region_df   = pd.read_csv(region_tsv, sep="\t")
    kept        = pd.read_csv(snplist_final, header=None, names=["SNPID"])
    df_filtered = region_df[region_df["SNPID"].isin(kept["SNPID"])].copy()
    df_filtered = df_filtered.set_index("SNPID").reindex(kept["SNPID"]).reset_index()

    if df_filtered.empty:
        die("After PLINK filtering, no variants remain for SuSiE.")

    ld_file    = Path(str(ld_prefix)    + ".ld")
    ld_r2_file = Path(str(ld_r2_prefix) + ".ld")
    if not ld_file.exists():
        die(f"LD matrix not found: {ld_file}")
    if not ld_r2_file.exists():
        die(f"LD r² matrix not found: {ld_r2_file}")

    R_df  = pd.read_csv(ld_file,    sep="\t", header=None).values
    R_df2 = pd.read_csv(ld_r2_file, sep="\t", header=None).values

    if R_df.shape[0] != df_filtered.shape[0]:
        die(f"Dimension mismatch: LD={R_df.shape} vs variants={df_filtered.shape[0]}")

    # Regularize: adds lambda to diagonal to ensure positive semi-definiteness.
    # Required when LD reference n is small relative to the number of variants.
    n_snps = R_df.shape[0]
    R_reg  = R_df + args.lambda_reg * np.eye(n_snps)
    print(f"[OK] LD matrix: {n_snps}x{n_snps}, regularized with lambda={args.lambda_reg}")

    # LD heatmaps
    fig, ax = plt.subplots(ncols=2, figsize=(20, 10))
    sns.heatmap(R_df,  cmap="Spectral", ax=ax[0]); ax[0].set_title("LD r matrix")
    sns.heatmap(R_df2, cmap="Spectral", ax=ax[1]); ax[1].set_title("LD r² matrix")
    fig.savefig(outdir / "ld_heatmaps.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    # -----------------------------------------------------------------------
    # 5. Probe SuSiE at L to determine effective K
    # -----------------------------------------------------------------------
    # When estimate_prior_variance=TRUE, SuSiE automatically zeros out
    # components not supported by the data, so the number of non-empty
    # credible sets returned is the effective K. A single run at L suffices;
    # ELBO-based K sweeps are redundant and produce identical ELBOs for all
    # L >= true_K.
    print(f"\n[5/7] Probe SuSiE at L={args.L} to determine effective K")
    susieR = importr("susieR")

    def run_susie(L: int) -> object:
        return susieR.susie_rss(
            bhat=df_filtered["BETA"].values.reshape((-1, 1)),
            shat=df_filtered["SE"].values.reshape((-1, 1)),
            R=R_reg,
            L=L,
            n=args.n,
            # estimate_residual_variance=FALSE: required when the LD reference
            # panel is a different (out-of-sample) dataset from the GWAS.
            # Setting TRUE causes the residual variance estimator to go negative
            # when sample sizes differ, crashing the model.
            estimate_residual_variance=False,
            # estimate_prior_variance=TRUE: SuSiE prunes empty components,
            # making the number of returned CSs the effective K directly.
            estimate_prior_variance=True,
            max_iter=1000,
        )

    try:
        fit_probe   = run_susie(args.L)
        cs_probe    = susieR.susie_get_cs(
            fit_probe, coverage=args.coverage,
            min_abs_corr=args.min_abs_corr, Xcorr=R_reg
        )[0]
        effective_k = len(cs_probe) if cs_probe else 1
        probe_conv  = bool(fit_probe.rx2("converged")[0])
        print(f"[OK] Probe converged: {probe_conv} | Effective K: {effective_k}")
    except Exception as e:
        print(f"[WARN] Probe failed ({e}), defaulting effective_k=1")
        effective_k = 1

    final_L = max(effective_k, args.L)
    print(f"[OK] Final SuSiE will run with L={final_L}")

    # -----------------------------------------------------------------------
    # 6. Final SuSiE run
    # -----------------------------------------------------------------------
    print(f"\n[6/7] Final SuSiE run (L={final_L})")
    try:
        fit = run_susie(final_L)
    except Exception as e:
        die(f"Final SuSiE run failed: {e}")

    converged     = bool(fit.rx2("converged")[0])
    pip_values    = np.array(susieR.susie_get_pip(fit))
    credible_sets = susieR.susie_get_cs(
        fit, coverage=args.coverage,
        min_abs_corr=args.min_abs_corr, Xcorr=R_reg
    )[0]

    print(ro.r.summary(fit))
    print(f"Converged: {converged}")
    if not converged:
        print("[WARN] Final model did not converge — interpret results with caution.")

    df_filtered = df_filtered.copy()
    df_filtered["pip"] = pip_values
    df_filtered["cs"]  = 0
    for i in range(len(credible_sets)):
        idx = np.array(credible_sets[i], dtype=int) - 1   # R is 1-indexed
        df_filtered.loc[idx, "cs"] = i + 1

    # -----------------------------------------------------------------------
    # 7. Save outputs and plots
    # -----------------------------------------------------------------------
    print(f"\n[7/7] Saving outputs to {outdir}")

    # Per-variant results
    df_filtered.to_csv(outdir / "susie_results.tsv", sep="\t", index=False)

    # Credible set variants only
    cs_out = df_filtered[df_filtered["cs"] > 0].sort_values(
        ["cs", "pip"], ascending=[True, False]
    )
    cs_out.to_csv(outdir / "credible_set_variants.tsv", sep="\t", index=False)

    # Locus plot (-log10P top, PIP bottom)
    df_filtered["MLOG10P"] = -np.log10(df_filtered["P"])
    lead_idx = int(df_filtered["P"].values.argmin())
    ld_color = R_df2[:, lead_idx]

    fig, axes = plt.subplots(nrows=2, sharex=True, figsize=(15, 7),
                              gridspec_kw={"height_ratios": (4, 1)})
    for ax, y_col, y_label in zip(axes,
                                   ["MLOG10P", "pip"],
                                   ["-log10(P)", "PIP"]):
        p = ax.scatter(df_filtered["POS"], df_filtered[y_col],
                       c=ld_color, cmap="viridis", s=20, zorder=2)
        ax.scatter(
            df_filtered.loc[df_filtered["cs"] == 1, "POS"],
            df_filtered.loc[df_filtered["cs"] == 1, y_col],
            marker="o", s=40, c="red", edgecolors="black",
            label="Credible set 1", zorder=3
        )
        if y_col == "MLOG10P":
            ax.scatter(
                df_filtered.loc[(df_filtered["CHR"] == args.chr) &
                                 (df_filtered["POS"] == args.pos), "POS"],
                df_filtered.loc[(df_filtered["CHR"] == args.chr) &
                                 (df_filtered["POS"] == args.pos), y_col],
                marker="*", s=120, c="gold", edgecolors="black",
                label=f"Lead SNP ({args.pos})", zorder=4
            )
            ax.set_title(
                f"Fine-mapping: chr{args.chr}:{start}-{end}  "
                f"(Effective K={effective_k}, L={final_L})"
            )
        plt.colorbar(p, ax=ax, label="r² with lead variant")
        ax.set_ylabel(y_label)
        ax.legend(fontsize=8)
        ax.set_xlim(start, end)
    axes[1].set_xlabel("Position (GRCh38)")
    plt.tight_layout()
    fig.savefig(outdir / "credible_set_plot.pdf", dpi=300)
    plt.close(fig)

    # CS summary bar charts
    cs_ids = sorted(df_filtered.loc[df_filtered["cs"] > 0, "cs"].unique())
    cs_summary_rows = []
    for cs_id in cs_ids:
        sub = df_filtered[df_filtered["cs"] == cs_id]
        cs_summary_rows.append({
            "CS": cs_id,
            "Size": len(sub),
            "Max_PIP": float(sub["pip"].max()),
        })
    cs_summary = pd.DataFrame(cs_summary_rows)

    if not cs_summary.empty:
        fig, axes = plt.subplots(ncols=2, figsize=(12, 4))
        axes[0].bar(cs_summary["CS"], cs_summary["Size"],  color="steelblue",  alpha=0.8)
        axes[0].set_xlabel("Credible Set"); axes[0].set_ylabel("Number of Variants")
        axes[0].set_title(f"Credible Set Sizes\nchr{args.chr}:{args.pos}")
        axes[1].bar(cs_summary["CS"], cs_summary["Max_PIP"], color="darkorange", alpha=0.8)
        axes[1].set_xlabel("Credible Set"); axes[1].set_ylabel("Max PIP")
        axes[1].set_ylim(0, 1)
        axes[1].set_title(f"Max PIP per Credible Set\n(K={effective_k}, L={final_L})")
        plt.tight_layout()
        fig.savefig(outdir / "credible_set_summary.pdf", dpi=300)
        plt.close(fig)

    # JSON summary
    summary = {
        "chromosome":          args.chr,
        "lead_position":       args.pos,
        "window_bp":           args.window,
        "region":              f"chr{args.chr}:{start}-{end}",
        "reference_fasta":     args.ref_fasta,
        "n_variants_input":    int(locus.data.shape[0]),
        "n_variants_final":    int(df_filtered.shape[0]),
        "susie_L":             final_L,
        "effective_K":         effective_k,
        "coverage":            args.coverage,
        "min_abs_corr":        args.min_abs_corr,
        "lambda_reg":          args.lambda_reg,
        "estimate_residual_variance": False,
        "estimate_prior_variance":    True,
        "converged":           converged,
        "n_credible_sets":     len(cs_ids),
        "sample_size":         args.n,
    }
    (outdir / "susie_summary.json").write_text(json.dumps(summary, indent=2))

    print(f"\n[DONE] Effective K={effective_k} | Converged={converged} | "
          f"Credible sets={len(cs_ids)} | Output: {outdir}")


if __name__ == "__main__":
    main()
