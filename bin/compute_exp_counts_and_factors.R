#!/usr/bin/env Rscript
pacman::p_load(
  tidyverse,
  stringr,
  devtools
)
install_github("kauralasoo/rasqual/rasqualTools")
library(rasqualTools)

args <- commandArgs(trailingOnly = TRUE)
path <- list(
  sample_names = args[1],
  gene_gc = args[2],
  gene_counts = args[3],
  out = list(
    size_factors_bin = args[4],
    gene_counts_bin = args[5],
    gene_counts_txt = args[6]
  )
)
print(path)
# path <- list(
#   sample_names = 'out/timepoint/time_30/sample_genotype_map.list',
#   gene_gc = 'out/preprocessing/gene_gc_prct.tsv',
#   gene_counts = 'out/preprocessing/gene_counts.tsv')

sample_names <- read_tsv(path$sample_names,
                         col_names = c('rna_name', 'gt_name'))$rna_name
gene_gc <- read_tsv(path$gene_gc)
gene_counts <- read_tsv(path$gene_counts) %>%
  filter(!str_detect(gene_id, '_PAR_Y$')) %>%
  semi_join(gene_gc)

gene_gc <- gene_gc %>%
  inner_join(gene_counts) %>%
  dplyr::select(gene_id, percentage_gc_content)

count_matrix <- gene_counts %>%
  dplyr::select(one_of(c('gene_id', sample_names))) %>%
  as.data.frame() %>%
  column_to_rownames(var = 'gene_id')
size_factors = rasqualCalculateSampleOffsets(count_matrix, gene_gc, gc_correct = TRUE)

writeBin(as.double(c(t(count_matrix))), path$out$gene_counts_bin) # '/tmp/test.bin')
writeBin(as.double(c(t(size_factors))), path$out$size_factors_bin)
gene_counts %>% write_tsv(path$out$gene_counts_txt)
