#!/usr/bin/env Rscript
stopifnot(getRversion() >= "3.2.0")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  tidyverse,
  stringr,
  biomaRt
)

args <- commandArgs(trailingOnly = TRUE)
path <- list(
  raw_rasqual_input = args[1],
  annotated_output = args[2],
  run_mode = args[3],
  timepoint = args[4]
)

header = c(
  'gene_id',
  'rs_id',
  'chrom',
  'pos',
  'Ref_allele',
  'Alt_allele',
  'Allfreq',
  'HWE_ChiSq',
  'Imput_score',
  'L10_BH',
  'ChiSq',
  'effectSz', # (Pi)
  'Delta', # seq map error
  'Phi', # Reference allele mapping bias
  'Odisp', #overdispersion
  'SNPRegionID',
  'NoFeaturesNPs',
  'NotestedSNPs',
  'No. of iterations for null hypothesis',
  'No. of iterations for alternative hypothesis',
  'Random location of ties (tie lead SNP; only useful with -t option)',
  'LogLikH0',
  'Convergence',
  'Squared correlation between prior and posterior genotypes (fSNPs)',
  'Squared correlation between prior and posterior genotypes (rSNP)*')
#path <- list(raw_rasqual_input = 'rasqual_CeD_LD_output.txt')
df <- read_tsv(path$raw_rasqual_input, col_names = header) %>%
  dplyr::select(gene_id, rs_id, chrom, pos, ChiSq, effectSz, LogLikH0) %>%
  mutate(
    pval = pchisq(ChiSq, 1, lower.tail = FALSE),
    FDR = p.adjust(pval, method = 'BH'),
    ensembl_gene_id = str_extract(gene_id, '^ENSG\\d+'),
    run_mode = path$run_mode,
    timepoint = path$timepoint
  )


ensembl_version = 'feb2014.archive.ensembl.org'
ensembl_connect <- function(ensembl_version = 'feb2014.archive.ensembl.org') {
  useMart(
    'ENSEMBL_MART_ENSEMBL',
    host = ensembl_version,
    dataset = 'hsapiens_gene_ensembl')
}
mart = ensembl_connect()
gene_names <- getBM(attributes = c(
  'ensembl_gene_id', 'external_gene_id', 'gene_biotype'),
      filters = 'ensembl_gene_id',
      values = df$ensembl_gene_id,
      mart = mart)

df %>%
  left_join(gene_names) %>%
  dplyr::select(ensembl_gene_id, external_gene_id, rs_id, chrom, pos, pval, FDR, effectSz, everything()) %>%
  arrange(pval) %>%
  write_tsv(path$annotated_output)
