# Multi-cohort GWAS of Liver Enzyme Traits

This repository provides template scripts used to perform genome-wide association
analyses and downstream statistical genetics analyses for circulating liver enzyme
traits. The code is intentionally minimalist and cohort-agnostic, designed to
document analytical methodology rather than reproduce cohort-specific pipelines.

The repository accompanies the manuscript:

*Multi-cohort genome-wide association analyses reveal loci underlying circulating
liver enzyme levels in African-ancestry populations.*

No individual-level data or cohort-specific configuration files are included.

---

## Repository scope

The repository documents the following analytical components:

- Genome-wide association analyses (SAIGE)
- Meta-analysis across cohorts (METASOFT)
- Conditional analyses (GCTA–COJO)
- Fine-mapping (SuSiE RSS)
- Colocalisation with cis-eQTLs (FastENLOC)
- Genetic correlation analyses (LD score regression)

Each component is represented by a single template script or example command.

---

## Genome-wide association analyses

Genome-wide association analyses were performed within each cohort using SAIGE,
which fits a generalized linear mixed model to account for population structure
and sample relatedness. Quantitative phenotypes were analyzed with adjustment for
cohort-specific covariates, and association testing was performed per chromosome.

A cohort-agnostic SAIGE template is provided in `gwas/run_saige.sh`. Users supply
cohort-specific genotype inputs, phenotype/covariate files, and SAIGE locations at
runtime.

### Example commands

```bash
# Step 1: fit null model (once per phenotype)
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

# Step 2: association test (per chromosome)
bash gwas/run_saige.sh --step2 \
  --vcf <chr1.vcf.gz> \
  --vcf-index <chr1.vcf.gz.csi> \
  --chrom chr1 \
  --out-prefix results/saige/<cohort>/<phenotype> \
  --saige-out results/saige/<cohort>/<phenotype>.chr1.assoc.txt
```

---

## Meta-analysis

Meta-analyses were performed using METASOFT (Han–Eskin framework) on cohort-level
GWAS summary statistics. For each phenotype, harmonized summary files were merged
by variant identifier, allowing missingness across cohorts, and effect sizes with
standard errors were meta-analyzed using METASOFT.

The script `meta_analysis/run_metasoft.sh` provides a cohort-agnostic template that
accepts an arbitrary number of cohorts and phenotypes. Cohort-specific file paths
and software locations are supplied at runtime; no cohort-specific logic is
hard-coded in the repository.

---

## Conditional analyses

Conditional analyses were performed using GCTA–COJO to assess the independence of
lead association signals from previously reported variants at the same locus.
For each region of interest, summary statistics were conditioned on known
trait-associated SNPs using a linkage disequilibrium reference panel matched to
the study population.

An example GCTA–COJO command illustrating the analytical approach is shown below:

```bash
gcta64 \
  --bfile <LD_reference_prefix> \
  --cojo-file <summary_statistics> \
  --cojo-cond <known_snps> \
  --out <output_prefix>
```

---

## Fine-mapping

Statistical fine-mapping was performed using SuSiE (susieR `susie_rss`) based on
meta-analysis summary statistics. For each locus, variants within a fixed window
around the lead position were extracted, and a linkage disequilibrium (LD)
correlation matrix was computed from a population-matched reference panel using
PLINK. SuSiE posterior inclusion probabilities (PIP) and 95% credible sets were
then inferred for each region.

A cohort-agnostic template implementation is provided in
`finemapping/run_susie_rss.py`.

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

## Colocalisation analysis

Colocalisation analyses were performed using FastENLOC to evaluate shared genetic
signals between GWAS loci and cis-eQTLs. GWAS summary statistics (effect sizes and
standard errors) were integrated with annotated eQTL VCFs, and colocalisation
probabilities were estimated in a tissue-specific manner.

A cohort-agnostic template is provided in `colocalisation/run_fastenloc.sh`.

### Example command

```bash
bash colocalisation/run_fastenloc.sh \
  --fastenloc <path_to_fastenloc_binary> \
  --eqtl <annotated_eqtl.vcf> \
  --gwas <gwas_beta_se.tsv> \
  --total-variants <total_variant_count> \
  --tissue Liver \
  --out-prefix eqtl_liver
```

---

## LD score regression

Genetic correlations between traits were estimated using LD score regression
(LDSC). GWAS summary statistics were reformatted and munged using the LDSC
pipeline, and pairwise genetic correlations were computed using ancestry-matched
LD score reference panels.

A cohort-agnostic LDSC workflow template is provided in `ldsc/run_ldsc_rg.sh`.

### Example command

```bash
bash ldsc/run_ldsc_rg.sh \
  --gwas1 <trait1_gwas.tsv.gz> \
  --gwas2 <trait2_gwas.tsv.gz> \
  --n1 <sample_size_trait1> \
  --n2 <sample_size_trait2> \
  --ref-ld UKBB.AFR \
  --w-ld UKBB.AFR \
  --out trait1_vs_trait2_rg
```

---

## Notes on reproducibility

Scripts in this repository are provided as methodological templates. Exact cohort
definitions, quality control procedures, covariate specifications, and software
versions are described in detail in the accompanying manuscript.

---

## License

This repository is released under the BSD 3-Clause License, consistent with NIH
open science and data-sharing requirements.

