#!/usr/bin/env Rscript
stopifnot(getRversion() >= "3.2.0")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  tidyverse
)

args <- commandArgs(trailingOnly = TRUE)
path <- list(
  raw_rasqual_input = args[1],
  annotated_output = args[2]
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

df <- read_tsv(path$raw_rasqual_input, col_names = header) %>%
  dplyr::select(gene_id, rs_id, chrom, pos, ChiSq, effectSz, LogLikH0) %>%
  mutate(
    pval = pchisq(ChiSq, 1, lower.tail = FALSE),
    FDR = p.adjust(pval, method = 'BH')) %>%
  arrange(FDR)

df %>%
  write_tsv(path$annotated_output)
