#!/usr/bin/env Rscript
stopifnot(getRversion() >= "3.2.0")
pacman::p_load(
  assertthat,
  tidyverse,
  rtracklayer,
  GenomicFeatures,
  VariantAnnotation
  )

args <- commandArgs(trailingOnly = TRUE)
path <- list(
  gtf = args[1],
  output = args[2]
)

get_exon_unions <- function(uxons) {
  df <- as.data.frame(uxons) %>% as_tibble() %>%
    dplyr::rename(gene_id = group_name, chr = seqnames,
                  exon_start = start, exon_end = end)
  udf <- df %>%
    dplyr::select(gene_id, chr, strand) %>%
    distinct() %>%
    mutate(strand = ifelse(strand == '+', 1, -1))
  exon_df <- df %>% group_by(gene_id) %>% summarize(
      exon_starts = paste(exon_start, collapse = ','),
      exon_ends = paste(exon_start, collapse = ',')
    )
  inner_join(udf, exon_df, by = 'gene_id') %>%
    mutate(strand = as.integer(strand), chr=as.character(chr))
}

assert_that(file.exists(path$gtf))
txdb <- makeTxDbFromGFF(file = path$gtf)
uxons <- reduce(exonsBy(txdb, "gene"))
rasqual_df <- get_exon_unions(uxons)
rasqual_df %>% write_tsv(path$output)

