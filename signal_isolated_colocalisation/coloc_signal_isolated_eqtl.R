#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(coloc)
  library(ggplot2)
})

# =========================
# CLI
# =========================
opt_list <- list(
  make_option("--trait", type="character", help="Trait name (e.g., ALP/ALT/AST/GGT)"),
  make_option("--locus", type="character", help="Locus label (e.g., ALT_rs738408)"),
  make_option("--lead-chr", type="integer", help="Lead chromosome (integer)"),
  make_option("--lead-pos", type="integer", help="Lead position (bp)"),
  make_option("--window-bp", type="integer", default=1000000, help="Window around lead (default: 1,000,000; i.e., ±1 Mb)"),

  make_option("--gwas", type="character", help="GWAS summary stats (TSV; may be .gz)"),
  make_option("--eqtl", type="character", help="Signal-isolated liver eQTL summary stats (TSV; may be .gz)"),

  make_option("--N-gwas", type="integer", help="GWAS sample size N"),
  make_option("--N-eqtl", type="integer", default=NA_integer_,
              help="Optional: eQTL sample size N (used if eqtl file has no N column)"),

  make_option("--out", type="character", help="Output coloc results TSV"),
  make_option("--plot", type="character", default=NA_character_,
              help="Optional: output plot path (.png); if omitted, no plot is written")
)

opt <- parse_args(OptionParser(option_list=opt_list))

stopifnot(!is.na(opt$trait), !is.na(opt$locus),
          !is.na(opt$`lead-chr`), !is.na(opt$`lead-pos`),
          !is.na(opt$gwas), !is.na(opt$eqtl),
          !is.na(opt$`N-gwas`), !is.na(opt$out))

msg <- function(...) cat(sprintf("[INFO] %s\n", sprintf(...)))

# =========================
# IO helpers
# =========================
read_tsv_auto <- function(path){
  if (grepl("\\.gz$", path, ignore.case=TRUE)) {
    fread(cmd = paste("zcat", shQuote(path)), sep="\t", header=TRUE, showProgress=FALSE)
  } else {
    fread(path, sep="\t", header=TRUE, showProgress=FALSE)
  }
}

require_cols <- function(dt, cols, label){
  miss <- setdiff(cols, names(dt))
  if(length(miss) > 0) stop(sprintf("%s missing required columns: %s", label, paste(miss, collapse=", ")))
}

coerce_numeric <- function(dt, cols){
  for (cc in cols) dt[, (cc) := as.numeric(get(cc))]
  dt
}

qc_common <- function(dt, label){
  # Basic QC on commonly used columns
  for (cc in intersect(c("BETA","SE","P","MAF"), names(dt))){
    dt <- dt[is.finite(get(cc))]
  }
  if("SE" %in% names(dt))  dt <- dt[SE > 0]
  if("P" %in% names(dt))   dt <- dt[P > 0 & P <= 1]
  if("MAF" %in% names(dt)) dt <- dt[MAF > 0 & MAF < 1]
  if(nrow(dt) == 0) stop(sprintf("%s has 0 variants after QC", label))
  dt
}

dedup_by_snp_minp <- function(dt, label){
  if(any(duplicated(dt$SNP))){
    msg("%s: duplicated SNP IDs detected; keeping smallest P per SNP", label)
    setorder(dt, SNP, P)
    dt <- dt[, .SD[1], by=SNP]
  }
  dt
}

# =========================
# Read + standardize GWAS
# =========================
msg("Reading GWAS: %s", opt$gwas)
gwas <- read_tsv_auto(opt$gwas)

# Minimal standardization to canonical names
if("P.value" %in% names(gwas)) setnames(gwas, "P.value", "P")
if("p" %in% names(gwas))       setnames(gwas, "p", "P")
if("beta" %in% names(gwas) && !"BETA" %in% names(gwas)) setnames(gwas, "beta", "BETA")
if("se" %in% names(gwas) && !"SE" %in% names(gwas))     setnames(gwas, "se", "SE")
if("Freq" %in% names(gwas) && !"MAF" %in% names(gwas))  setnames(gwas, "Freq", "MAF")
if("EAF" %in% names(gwas) && !"MAF" %in% names(gwas))   setnames(gwas, "EAF", "MAF")

require_cols(gwas, c("SNP","CHR","POS","BETA","SE","P"), "GWAS")
gwas <- coerce_numeric(gwas, c("CHR","POS","BETA","SE","P"))
if("MAF" %in% names(gwas)) gwas <- coerce_numeric(gwas, "MAF")

lo <- opt$`lead-pos` - opt$`window-bp`
hi <- opt$`lead-pos` + opt$`window-bp`
gwas <- gwas[CHR == opt$`lead-chr` & POS >= lo & POS <= hi]
if(nrow(gwas) == 0) stop("GWAS: 0 variants in window after filtering")

gwas <- qc_common(gwas, "GWAS")
gwas <- dedup_by_snp_minp(gwas, "GWAS")

# =========================
# Read + standardize eQTL
# =========================
msg("Reading eQTL: %s", opt$eqtl)
eqtl <- read_tsv_auto(opt$eqtl)

if("p" %in% names(eqtl))       setnames(eqtl, "p", "P")
if("p_value" %in% names(eqtl)) setnames(eqtl, "p_value", "P")
if("beta" %in% names(eqtl) && !"BETA" %in% names(eqtl)) setnames(eqtl, "beta", "BETA")
if("se" %in% names(eqtl) && !"SE" %in% names(eqtl))     setnames(eqtl, "se", "SE")
if("maf" %in% names(eqtl) && !"MAF" %in% names(eqtl))   setnames(eqtl, "maf", "MAF")

# Required eQTL cols
require_cols(eqtl, c("SNP","BETA","SE","P","MAF","GeneSymbol"), "eQTL")
eqtl <- coerce_numeric(eqtl, c("BETA","SE","P","MAF"))

# Optional locus filtering by position if present
if(all(c("CHR","POS") %in% names(eqtl))){
  eqtl <- coerce_numeric(eqtl, c("CHR","POS"))
  eqtl <- eqtl[CHR == opt$`lead-chr` & POS >= lo & POS <= hi]
}

eqtl <- qc_common(eqtl, "eQTL")
eqtl <- dedup_by_snp_minp(eqtl, "eQTL")

# Provide defaults for optional identifiers
if(!("ENSG" %in% names(eqtl))) eqtl[, ENSG := NA_character_]
if(!("Gene_Biotype" %in% names(eqtl))) eqtl[, Gene_Biotype := NA_character_]
if(!("SignalLead" %in% names(eqtl))) eqtl[, SignalLead := NA_character_]
if(!("N" %in% names(eqtl))) eqtl[, N := opt$`N-eqtl`]

# =========================
# Harmonize / overlap
# =========================
common <- intersect(gwas$SNP, eqtl$SNP)
if(length(common) < 100){
  stop(sprintf("Too few overlapping SNPs (n=%d). Ensure SNP IDs match exactly between GWAS and eQTL.", length(common)))
}

gwas2 <- gwas[SNP %in% common]
eqtl2 <- eqtl[SNP %in% common]

# Prepare GWAS dataset for coloc
d1 <- list(
  snp     = gwas2$SNP,
  beta    = gwas2$BETA,
  varbeta = (gwas2$SE)^2,
  N       = opt$`N-gwas`,
  type    = "quant"
)
if("MAF" %in% names(gwas2)) d1$MAF <- gwas2$MAF

# =========================
# Run coloc per gene (and per signal if provided)
# =========================
groups <- unique(eqtl2[, .(GeneSymbol, ENSG, Gene_Biotype, SignalLead)])
msg("Running coloc for %d gene/signal group(s); overlap=%d SNPs", nrow(groups), length(common))

res <- vector("list", nrow(groups))
keep <- 0L

for(i in seq_len(nrow(groups))){
  g <- groups[i]

  sub <- eqtl2[
    GeneSymbol == g$GeneSymbol &
      (is.na(g$ENSG) | ENSG == g$ENSG) &
      (is.na(g$SignalLead) | SignalLead == g$SignalLead)
  ]

  # Ensure same SNP order as d1
  sub <- sub[match(d1$snp, sub$SNP), ]
  if(any(is.na(sub$SNP))) next

  Neq <- unique(na.omit(sub$N))
  Neq <- if(length(Neq) >= 1) Neq[1] else NA_integer_

  d2 <- list(
    snp     = sub$SNP,
    beta    = sub$BETA,
    varbeta = (sub$SE)^2,
    MAF     = sub$MAF,
    N       = Neq,
    type    = "quant"
  )

  cc <- coloc.abf(dataset1 = d1, dataset2 = d2)

  s <- cc$summary
  keep <- keep + 1L
  res[[keep]] <- data.table(
    trait       = opt$trait,
    locus       = opt$locus,
    lead_chr    = opt$`lead-chr`,
    lead_pos    = opt$`lead-pos`,
    GeneSymbol  = g$GeneSymbol,
    SignalLead  = g$SignalLead,
    ENSG        = g$ENSG,
    Gene_Biotype = g$Gene_Biotype,
    nsnps       = length(common),
    PP0         = unname(s["PP.H0.abf"]),
    PP1         = unname(s["PP.H1.abf"]),
    PP2         = unname(s["PP.H2.abf"]),
    PP3         = unname(s["PP.H3.abf"]),
    PP4         = unname(s["PP.H4.abf"])
  )
}

out <- rbindlist(res[seq_len(keep)], fill=TRUE)
if(nrow(out) == 0) stop("No colocalization results produced (check gene grouping / overlap).")
setorder(out, -PP4, GeneSymbol)

# =========================
# Write results
# =========================
dir.create(dirname(opt$out), showWarnings=FALSE, recursive=TRUE)
fwrite(out, opt$out, sep="\t", quote=FALSE)
msg("Wrote results: %s", opt$out)

# =========================
# Optional plot: PP4 dotplot with all labels
# =========================
if(!is.na(opt$plot)){
  msg("Writing plot: %s", opt$plot)

  plot_dt <- copy(out)
  plot_dt[, label := ifelse(is.na(SignalLead) | SignalLead=="", GeneSymbol,
                            paste0(GeneSymbol, " (", SignalLead, ")"))]
  # order labels by PP4
  plot_dt[, label := factor(label, levels = plot_dt[order(PP4)]$label)]

  p <- ggplot(plot_dt, aes(x = PP4, y = label)) +
    geom_point() +
    geom_text(aes(label = label), hjust = -0.05, size = 3) +
    labs(x = "Colocalization posterior probability (PP4)", y = NULL,
         title = paste0(opt$trait, " — ", opt$locus)) +
    xlim(0, 1.05) +
    theme_bw() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(face="bold"))

  ggsave(opt$plot, p, width=10, height=7, dpi=300)
}

