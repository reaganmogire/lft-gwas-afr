#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Run FastENLOC colocalisation analysis between GWAS and eQTL data.

Usage:
  bash run_fastenloc.sh \
    --fastenloc <path_to_fastenloc_binary> \
    --eqtl <eqtl_annotated_vcf> \
    --gwas <gwas_beta_se_tsv> \
    --total-variants <integer> \
    --tissue <tissue_name> \
    --out-prefix <output_prefix>

Inputs:
  --eqtl            Annotated eQTL VCF (e.g., GTEx v8 FastENLOC format)
  --gwas            GWAS summary statistics with beta and SE
  --total-variants  Total number of variants used by FastENLOC
  --tissue          Tissue label (for reporting)
  --out-prefix      Output prefix for FastENLOC results
EOF
}

die(){ echo "[ERROR] $*" >&2; exit 1; }

FASTENLOC=""
EQTL=""
GWAS=""
TOTAL=""
TISSUE=""
OUT=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --fastenloc)      FASTENLOC="$2"; shift 2;;
    --eqtl)           EQTL="$2"; shift 2;;
    --gwas)           GWAS="$2"; shift 2;;
    --total-variants) TOTAL="$2"; shift 2;;
    --tissue)         TISSUE="$2"; shift 2;;
    --out-prefix)     OUT="$2"; shift 2;;
    -h|--help)        usage; exit 0;;
    *) die "Unknown argument: $1";;
  esac
done

[[ -n "${FASTENLOC}" && -x "${FASTENLOC}" ]] || die "FastENLOC binary not found or not executable"
[[ -f "${EQTL}" ]] || die "Missing eQTL file: ${EQTL}"
[[ -f "${GWAS}" ]] || die "Missing GWAS summary file: ${GWAS}"
[[ -n "${TOTAL}" ]] || die "Missing --total-variants"
[[ -n "${TISSUE}" ]] || die "Missing --tissue"
[[ -n "${OUT}" ]] || die "Missing --out-prefix"

"${FASTENLOC}" \
  -eqtl "${EQTL}" \
  -sum "${GWAS}" \
  -total_variants "${TOTAL}" \
  -tissue "${TISSUE}" \
  -prefix "${OUT}"

echo "[OK] FastENLOC colocalisation complete"
echo "[OK] Output prefix: ${OUT}"

