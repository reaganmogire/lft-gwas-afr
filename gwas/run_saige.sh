#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Run SAIGE GWAS (cohort-agnostic template).

This script supports:
  --step1 : fit null GLMM model (recommended once per phenotype)
  --step2 : run single-variant association test (per chromosome)

Required inputs are supplied at runtime (no cohort-specific paths are hard-coded).

Usage (step1):
  bash gwas/run_saige.sh --step1 \
    --plink <plink_prefix> \
    --pheno-file <pheno_cov.tsv> \
    --pheno-col <phenotype> \
    --covar-list <covariates_csv> \
    --out-prefix <output_prefix> \
    --trait-type quantitative \
    --sample-id-col IID \
    --nthreads 32

Usage (step2):
  bash gwas/run_saige.sh --step2 \
    --vcf <chr.vcf.gz> \
    --vcf-index <chr.vcf.gz.csi> \
    --chrom chr1 \
    --out-prefix <same_output_prefix_from_step1> \
    --saige-out <assoc_output.txt> \
    --min-maf 0.01 \
    --min-mac 1
EOF
}

die(){ echo "[ERROR] $*" >&2; exit 1; }

STEP=""
PLINK=""
PHENO_FILE=""
PHENO_COL=""
COVAR_LIST=""
SAMPLE_ID_COL="IID"
TRAIT_TYPE="quantitative"
OUT_PREFIX=""
NTHREADS="16"
LOCO="TRUE"
RELATEDNESS_CUTOFF="0.0"
FEMALE_CODE="1"
MALE_CODE="0"
INV_NORM="FALSE"

VCF=""
VCF_INDEX=""
VCF_FIELD="GT"
CHROM=""
MIN_MAF="0.01"
MIN_MAC="1"
SAIGE_OUT=""

# SAIGE resolution: module (optional) + location of step scripts (required)
SAIGE_MODULE="${SAIGE_MODULE:-}"
SAIGE_EXTDATA_DIR="${SAIGE_EXTDATA_DIR:-}"

[[ $# -gt 0 ]] || { usage; exit 1; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --step1) STEP="1"; shift;;
    --step2) STEP="2"; shift;;

    --plink) PLINK="$2"; shift 2;;
    --pheno-file) PHENO_FILE="$2"; shift 2;;
    --pheno-col) PHENO_COL="$2"; shift 2;;
    --covar-list) COVAR_LIST="$2"; shift 2;;
    --sample-id-col) SAMPLE_ID_COL="$2"; shift 2;;
    --trait-type) TRAIT_TYPE="$2"; shift 2;;
    --out-prefix) OUT_PREFIX="$2"; shift 2;;
    --nthreads) NTHREADS="$2"; shift 2;;
    --loco) LOCO="$2"; shift 2;;
    --relatedness-cutoff) RELATEDNESS_CUTOFF="$2"; shift 2;;
    --female-code) FEMALE_CODE="$2"; shift 2;;
    --male-code) MALE_CODE="$2"; shift 2;;
    --inv-norm) INV_NORM="$2"; shift 2;;

    --vcf) VCF="$2"; shift 2;;
    --vcf-index) VCF_INDEX="$2"; shift 2;;
    --vcf-field) VCF_FIELD="$2"; shift 2;;
    --chrom) CHROM="$2"; shift 2;;
    --min-maf) MIN_MAF="$2"; shift 2;;
    --min-mac) MIN_MAC="$2"; shift 2;;
    --saige-out) SAIGE_OUT="$2"; shift 2;;

    -h|--help) usage; exit 0;;
    *) die "Unknown argument: $1";;
  esac
done

[[ -n "${STEP}" ]] || die "Specify --step1 or --step2"
[[ -n "${SAIGE_EXTDATA_DIR}" ]] || die "Set SAIGE_EXTDATA_DIR (env var) to directory containing step1_fitNULLGLMM.R and step2_SPAtests.R"

if [[ -n "${SAIGE_MODULE}" ]]; then
  module load "${SAIGE_MODULE}"
fi

cd "${SAIGE_EXTDATA_DIR}"

if [[ "${STEP}" == "1" ]]; then
  [[ -n "${PLINK}" ]] || die "--plink required for step1"
  [[ -f "${PHENO_FILE}" ]] || die "--pheno-file not found: ${PHENO_FILE}"
  [[ -n "${PHENO_COL}" ]] || die "--pheno-col required for step1"
  [[ -n "${COVAR_LIST}" ]] || die "--covar-list required for step1"
  [[ -n "${OUT_PREFIX}" ]] || die "--out-prefix required for step1"

  step1_fitNULLGLMM.R \
    --plinkFile="${PLINK}" \
    --phenoFile="${PHENO_FILE}" \
    --phenoCol="${PHENO_COL}" \
    --covarColList="${COVAR_LIST}" \
    --sampleIDColinphenoFile="${SAMPLE_ID_COL}" \
    --traitType="${TRAIT_TYPE}" \
    --outputPrefix="${OUT_PREFIX}" \
    --nThreads="${NTHREADS}" \
    --LOCO="${LOCO}" \
    --relatednessCutoff="${RELATEDNESS_CUTOFF}" \
    --FemaleCode="${FEMALE_CODE}" \
    --MaleCode="${MALE_CODE}" \
    --invNormalize="${INV_NORM}" \
    --IsOverwriteVarianceRatioFile=TRUE

  echo "[OK] SAIGE step1 complete: ${OUT_PREFIX}.rda"

elif [[ "${STEP}" == "2" ]]; then
  [[ -f "${VCF}" ]] || die "--vcf not found: ${VCF}"
  [[ -f "${VCF_INDEX}" ]] || die "--vcf-index not found: ${VCF_INDEX}"
  [[ -n "${CHROM}" ]] || die "--chrom required (e.g., chr1)"
  [[ -n "${OUT_PREFIX}" ]] || die "--out-prefix required (same used for step1)"
  [[ -n "${SAIGE_OUT}" ]] || die "--saige-out required (assoc output path)"

  MODEL="${OUT_PREFIX}.rda"
  VR="${OUT_PREFIX}.varianceRatio.txt"
  [[ -f "${MODEL}" ]] || die "Missing null model: ${MODEL} (run step1 first)"
  [[ -f "${VR}" ]] || die "Missing variance ratio: ${VR} (run step1 first)"

  step2_SPAtests.R \
    --vcfFile="${VCF}" \
    --vcfFileIndex="${VCF_INDEX}" \
    --vcfField="${VCF_FIELD}" \
    --chrom="${CHROM}" \
    --minMAF="${MIN_MAF}" \
    --minMAC="${MIN_MAC}" \
    --GMMATmodelFile="${MODEL}" \
    --varianceRatioFile="${VR}" \
    --SAIGEOutputFile="${SAIGE_OUT}" \
    --is_Firth_beta=TRUE \
    --pCutoffforFirth=0.05 \
    --LOCO="${LOCO}"

  echo "[OK] SAIGE step2 complete: ${SAIGE_OUT}"
fi

