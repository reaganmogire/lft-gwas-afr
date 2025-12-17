#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Run METASOFT meta-analysis for a given phenotype across an arbitrary set of cohorts.

Assumptions:
  - Each input file is a bgzipped (or gzipped) tab-delimited table with a header.
  - Column 1 is SNP (variant ID) and is identical across cohorts (or missing).
  - BETA is column 6 and SE is column 7 (as in: SNP CHR POS Allele1 Allele2 BETA SE P.value N).
  - The script will merge cohorts by SNP, allowing missingness, and will run METASOFT.

Usage:
  bash run_metasoft.sh \
    --phenotype ALB_int \
    --metasoft-jar /path/to/Metasoft.jar \
    --pvalue-table /path/to/HanEskinPvalueTable.txt \
    --outdir results/metasoft \
    --inputs \
      CARDIA:/path/to/sorted_ALB_int_data.tbl.gz \
      AOS:/path/to/sorted_ALB_int_data.tbl.gz \
      UKBB:/path/to/sorted_ALB_int_data.tbl.gz

Notes:
  - Input files are sorted internally on SNP (col1) after removing header.
  - Output:
      <outdir>/<phenotype>/<phenotype>_metasoft_input.txt
      <outdir>/<phenotype>/<phenotype>_metasoft_results.txt
      <outdir>/<phenotype>/<phenotype>_metasoft.log
EOF
}

die(){ echo "[ERROR] $*" >&2; exit 1; }

PHENO=""
METASOFT_JAR=""
PVAL_TABLE=""
OUTDIR="results/metasoft"
INPUT_SPECS=()

# -------- argument parsing --------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --phenotype)      PHENO="${2:-}"; shift 2;;
    --metasoft-jar)   METASOFT_JAR="${2:-}"; shift 2;;
    --pvalue-table)   PVAL_TABLE="${2:-}"; shift 2;;
    --outdir)         OUTDIR="${2:-}"; shift 2;;
    --inputs)         shift; while [[ $# -gt 0 && "$1" != --* ]]; do INPUT_SPECS+=("$1"); shift; done;;
    -h|--help)        usage; exit 0;;
    *) die "Unknown argument: $1 (use --help)";;
  esac
done

[[ -n "${PHENO}" ]] || die "Missing --phenotype"
[[ -n "${METASOFT_JAR}" ]] || die "Missing --metasoft-jar"
[[ -n "${PVAL_TABLE}" ]] || die "Missing --pvalue-table"
[[ -f "${METASOFT_JAR}" ]] || die "METASOFT jar not found: ${METASOFT_JAR}"
[[ -f "${PVAL_TABLE}" ]] || die "Pvalue table not found: ${PVAL_TABLE}"
[[ ${#INPUT_SPECS[@]} -ge 2 ]] || die "Provide >=2 cohorts via --inputs COHORT:file ..."

mkdir -p "${OUTDIR}/${PHENO}"
WORKDIR="$(mktemp -d -p "${OUTDIR}/${PHENO}" tmp.${PHENO}.XXXXXX)"
trap 'rm -rf "${WORKDIR}"' EXIT

echo "[INFO] phenotype=${PHENO}"
echo "[INFO] outdir=${OUTDIR}/${PHENO}"
echo "[INFO] workdir=${WORKDIR}"

# -------- preprocess: headerless + sorted for each cohort --------
COHORTS=()
FILES_SORTED=()

for spec in "${INPUT_SPECS[@]}"; do
  cohort="${spec%%:*}"
  file="${spec#*:}"

  [[ -n "${cohort}" && -n "${file}" ]] || die "Bad input spec: ${spec} (expected COHORT:/path/to/file.gz)"
  [[ -f "${file}" ]] || die "Input not found: ${file}"

  out="${WORKDIR}/${PHENO}_${cohort}.sorted.txt"
  zcat "${file}" | tail -n +2 | sort -k1,1 > "${out}"

  COHORTS+=("${cohort}")
  FILES_SORTED+=("${out}")

  echo "[INFO] prepared ${cohort}: ${out}"
done

# -------- iterative join: produce SNP, (BETA,SE)* across cohorts --------
# We only need columns 6 (BETA) and 7 (SE) from each cohort file.
# Layout of each cohort file: 1:SNP 2:CHR 3:POS 4:A1 5:A2 6:BETA 7:SE ...

# Start from first cohort: SNP, BETA1, SE1
MERGED="${WORKDIR}/merge_0.txt"
awk -F'\t' -v OFS='\t' '{print $1,$6,$7}' "${FILES_SORTED[0]}" > "${MERGED}"

for ((i=1; i<${#FILES_SORTED[@]}; i++)); do
  NEXT="${WORKDIR}/next_${i}.txt"
  awk -F'\t' -v OFS='\t' '{print $1,$6,$7}' "${FILES_SORTED[$i]}" > "${NEXT}"

  OUT="${WORKDIR}/merge_${i}.txt"
  join -a 1 -a 2 -e "NA" -t $'\t' -1 1 -2 1 -o auto "${MERGED}" "${NEXT}" > "${OUT}"
  MERGED="${OUT}"
done

# Final METASOFT input
META_INPUT="${OUTDIR}/${PHENO}/${PHENO}_metasoft_input.txt"
mv "${MERGED}" "${META_INPUT}"

echo "[INFO] METASOFT input written: ${META_INPUT}"
echo "[INFO] Cohorts: ${COHORTS[*]}"
echo "[INFO] Columns: SNP plus (BETA,SE) per cohort in the order above"

# -------- run METASOFT --------
META_OUT="${OUTDIR}/${PHENO}/${PHENO}_metasoft_results.txt"
META_LOG="${OUTDIR}/${PHENO}/${PHENO}_metasoft.log"

java -jar "${METASOFT_JAR}" \
  -input "${META_INPUT}" \
  -pvalue_table "${PVAL_TABLE}" \
  -output "${META_OUT}" \
  -log "${META_LOG}"

echo "[OK] METASOFT complete"
echo "[OK] Results: ${META_OUT}"
echo "[OK] Log:     ${META_LOG}"

