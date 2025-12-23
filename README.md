# Multi-cohort GWAS of Liver Enzyme Traits (ALP, ALT, AST, GGT)

This repository contains **minimal, cohort-agnostic template scripts** illustrating the main analytical steps used for genome-wide association and downstream statistical genetics analyses of circulating liver enzyme traits. The goal is to document **methods and execution patterns**; it is not a drop-in, end-to-end reproduction of any single cohort pipeline.

The repository accompanies the manuscript:

*Multi-cohort genome-wide association analyses reveal loci underlying circulating liver enzyme levels in African-ancestry populations.*

No individual-level data or cohort-specific configuration files are included.

---

## Repository scope

The repository documents the following analytical components:

1. Genome-wide association analyses (**SAIGE**)
2. Meta-analysis across cohorts (**METASOFT**, Han–Eskin random effects)
3. Conditional analyses (**GCTA–COJO**)
4. Statistical fine-mapping (**SuSiE RSS**)
5. Colocalisation with cis-eQTLs (**FastENLOC**, GTEx liver)
6. Secondary colocalisation using a signal-isolated liver eQTL resource (**COLOC / coloc.abf**)
7. Genetic correlation analyses (**LD score regression**)

Each component is represented by a single template script and/or an example command.

---

## End-to-end analysis outline 

**Step 1 — Run within-cohort GWAS (SAIGE).**  
Fit a null model per phenotype and cohort, then run association testing per chromosome.

**Step 2 — Harmonize cohort GWAS summary statistics.**  
Ensure consistent variant identifiers, alleles, effect direction, and column schema across cohorts.

**Step 3 — Meta-analyze across cohorts (METASOFT).**  
Perform Han–Eskin meta-analysis per phenotype using harmonized cohort-level summary statistics.

**Step 4 — Test independence of “novel” signals (GCTA–COJO).**  
Condition lead signals on previously reported variants within the locus using an ancestry-matched LD reference.

**Step 5 — Fine-map selected loci (SuSiE RSS).**  
Construct locus-specific LD matrices and estimate posterior inclusion probabilities (PIP) and credible sets.

**Step 6 — Colocalise GWAS loci with liver cis-eQTLs (FastENLOC).**  
Estimate regional colocalisation probabilities using liver eQTL resources (e.g., GTEx liver for primary analysis).

**Step 7 — Secondary colocalisation with signal-isolated liver eQTLs (COLOC).**  
Run coloc.abf per gene (and per independent eQTL signal where available) and generate PP4 summary plots/tables.

**Step 8 — Estimate genetic correlations (LDSC).**  
Munge summary statistics and compute rg using ancestry-matched LD scores.

---

## Step 1 — Genome-wide association analyses (SAIGE)

Genome-wide association analyses were performed within each cohort using **SAIGE**, which fits a generalized linear mixed model to account for population structure and sample relatedness. Quantitative phenotypes were analyzed with adjustment for cohort-specific covariates, and association testing was performed per chromosome.

A cohort-agnostic SAIGE template is provided in `gwas/run_saige.sh`. Users provide cohort-specific inputs (genotypes, phenotypes/covariates, sample IDs) and software paths at runtime.

### Example commands

```bash
# 1.1 Fit null model (once per phenotype per cohort)
export SAIGE_EXTDATA_DIR=<path_to_saige_extdata>
bash gwas/run_saige.sh --step1 \
  --plink <plink_prefix> \
  --pheno-file <pheno_cov.tsv> \
  --pheno-col <phenotype> \
  --covar-list <covariates_csv> \
  --out-prefix results/saige/<cohort>/<phenotype> \
  --trait-type quantitative \
  --sample-id-col IID \
  --nthreads 32

# 1.2 Association test (per chromosome)
bash gwas/run_saige.sh --step2 \
  --vcf <chr1.vcf.gz> \
  --vcf-index <chr1.vcf.gz.csi> \
  --chrom chr1 \
  --out-prefix results/saige/<cohort>/<phenotype> \
  --saige-out results/saige/<cohort>/<phenotype>.chr1.assoc.txt
```

---

## Step 2 — Meta-analysis across cohorts (METASOFT)

Meta-analyses were performed using **METASOFT** (Han–Eskin framework) on cohort-level GWAS summary statistics. For each phenotype, harmonized summary files were merged by variant identifier (allowing missingness across cohorts) and effect sizes with standard errors were meta-analyzed.

A cohort-agnostic template is provided in `meta_analysis/run_metasoft.sh`. Cohort-specific file paths and software locations are supplied at runtime; no cohort-specific logic is hard-coded.

---

## Step 3 — Conditional analyses (GCTA–COJO)

Conditional analyses were performed using **GCTA–COJO** to assess the independence of lead association signals from previously reported variants at the same locus. For each region, summary statistics were conditioned on known trait-associated SNPs using an LD reference panel matched to the study population.

### Example command

```bash
gcta64 \
  --bfile <LD_reference_prefix> \
  --cojo-file <summary_statistics> \
  --cojo-cond <known_snps_list.txt> \
  --out <output_prefix>
```

---

## Step 4 — Fine-mapping (SuSiE RSS)

Statistical fine-mapping was performed using **SuSiE** (susieR `susie_rss`) based on meta-analysis summary statistics. For each locus, variants within a fixed window around the lead position were extracted, and an LD correlation matrix was computed from a population-matched reference panel (e.g., PLINK). SuSiE posterior inclusion probabilities (PIP) and 95% credible sets were then inferred for each region.

A cohort-agnostic template implementation is provided in `finemapping/run_susie_rss.py`.

### Example command

```bash
python finemapping/run_susie_rss.py \
  --sumstats <sorted_sumstats.tsv.gz> \
  --bfile <LD_reference_plink_prefix> \
  --chr <chromosome> \
  --pos <lead_position_bp> \
  --window 250000 \
  --n <sample_size> \
  --plink <path_to_plink> \
  --outdir results/susie
```

---

## Step 5 — Colocalisation with liver cis-eQTLs (FastENLOC)

Colocalisation analyses were performed using **FastENLOC** to evaluate whether GWAS loci share a causal variant with liver cis-eQTL signals. GWAS summary statistics (effect sizes and standard errors) were integrated with annotated eQTL resources, and tissue-specific colocalisation probabilities were estimated for **liver**.

A cohort-agnostic template is provided in `colocalisation/run_fastenloc.sh`.

### Example command

```bash
bash colocalisation/run_fastenloc.sh \
  --fastenloc <path_to_fastenloc_binary> \
  --eqtl <annotated_eqtl.vcf.gz> \
  --gwas <gwas_beta_se.tsv.gz> \
  --total-variants <total_variant_count> \
  --tissue Liver \
  --out-prefix results/fastenloc/<trait>/<locus>/gtex_liver
```

---

## Step 6 — Secondary colocalisation (GWAS × signal-isolated liver eQTL; COLOC)

This repository provides a single script (`coloc_signal_isolated_eqtl.R`) that performs Bayesian colocalization using `coloc.abf` between a GWAS locus and a liver eQTL resource with **conditionally independent (signal-isolated)** eQTLs. The script outputs posterior probabilities **PP0–PP4** per gene (and per independent signal where available) and can optionally generate a **PP4 dot plot**.

### 6.1 Input 1 — GWAS summary statistics (TSV; optionally gzipped)

Required columns (tab-delimited; header required):

- `SNP`: variant identifier. **Must match the eQTL file exactly** (choose one convention and use it consistently).
- `CHR`: chromosome (integer)
- `POS`: base-pair position (integer)
- `BETA`: effect size (numeric)
- `SE`: standard error (numeric; >0)
- `P`: p-value (numeric; 0 < P ≤ 1)

Optional column:
- `MAF`: minor allele frequency (numeric; 0 < MAF < 1)

Notes:
- The script filters variants to the locus window defined by `--lead-chr`, `--lead-pos`, and `--window-bp` (default ±1 Mb).

### 6.2 Input 2 — Signal-isolated liver eQTL summary statistics (TSV; optionally gzipped)

Required columns:

- `SNP`: same variant identifier as GWAS
- `BETA`: eQTL effect size
- `SE`: eQTL standard error (>0)
- `P`: eQTL p-value (0 < P ≤ 1)
- `MAF`: minor allele frequency (0 < MAF < 1)
- `GeneSymbol`: target gene symbol

Optional columns:

- `SignalLead`: label for an independent eQTL signal (if present, colocalisation is run per **GeneSymbol × SignalLead**)
- `ENSG`: Ensembl gene ID
- `Gene_Biotype`: gene biotype annotation
- `CHR`, `POS`: if present, the eQTL file is also window-filtered to the locus region
- `N`: eQTL sample size; if absent, provide a constant via `--N-eqtl`

### 6.3 Output

A tab-delimited results table with columns:

`trait, locus, lead_chr, lead_pos, GeneSymbol, SignalLead, ENSG, Gene_Biotype, nsnps, PP0, PP1, PP2, PP3, PP4`

Optionally, a PP4 dot plot (all points labeled).

### 6.4 Run

```bash
Rscript coloc_signal_isolated_eqtl.R \
  --trait ALT \
  --locus ALT_rs738408 \
  --lead-chr 22 \
  --lead-pos 44324730 \
  --window-bp 1000000 \
  --gwas path/to/ALT.sumstats.tsv.gz \
  --eqtl path/to/liver_signal_isolated_eqtl.tsv.gz \
  --N-gwas 28989 \
  --out results/coloc/ALT_rs738408.coloc.tsv \
  --plot results/coloc/ALT_rs738408.PP4.png
```

Dependencies (R packages): `data.table`, `coloc`, `optparse`, `ggplot2`.

---

## Step 7 — LD score regression (genetic correlation)

Genetic correlations between traits were estimated using **LD score regression (LDSC)**. GWAS summary statistics were reformatted and munged using the LDSC pipeline, and pairwise genetic correlations were computed using ancestry-matched LD score reference panels.

A cohort-agnostic LDSC workflow template is provided in `ldsc/run_ldsc_rg.sh`.

### Example command

```bash
bash ldsc/run_ldsc_rg.sh \
  --gwas1 <trait1_gwas.tsv.gz> \
  --gwas2 <trait2_gwas.tsv.gz> \
  --n1 <sample_size_trait1> \
  --n2 <sample_size_trait2> \
  --ref-ld <ref_ld_prefix> \
  --w-ld <w_ld_prefix> \
  --out results/ldsc/trait1_vs_trait2_rg
```

---

## Notes on reproducibility

Scripts in this repository are provided as methodological templates. Exact cohort definitions, quality control procedures, covariate specifications, software versions, and parameter choices are described in the accompanying manuscript and supplementary materials.

---

## License

This repository is released under the **BSD 3-Clause License**, consistent with NIH open science and data-sharing requirements.

