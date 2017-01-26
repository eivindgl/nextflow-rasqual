gene_counts = file('out/preprocessing/gene_counts.tsv')

process normalize_count_table {
  publishDir params.post_proc_dir, mode: 'copy'
  //module 'R/3.3.1-foss-2015b'

  input:
    file 'gene_counts.tsv' from gene_counts
  output:
    file 'normalized_gene_counts.tsv' into norm_gene_counts_ch
  '''
  #!/usr/bin/env Rscript
  if (!require("pacman")) install.packages("pacman")
  pacman::p_load(
    tidyverse,
    stringr,
    DESeq2
  )

  raw_counts <- read_tsv('gene_counts.tsv') %>% 
    dplyr::select(-ensembl_gene_id) %>% 
    as.data.frame() %>% 
    column_to_rownames('gene_id')

  meta <- tibble(sample_name = names(raw_counts)) 
  x <- meta$sample_name %>% 
    str_split_fixed('_', 2) %>% 
    as_tibble()
  colnames(x) <- c('timepoint', 'sampleID')
  meta <- meta %>% bind_cols(x)

  DESeqDataSetFromMatrix(raw_counts, colData = meta, design = ~ sampleID + timepoint) %>% 
    estimateSizeFactors() %>% 
    varianceStabilizingTransformation(blind = FALSE) %>%
    assay() %>% 
    as.data.frame() %>% 
    rownames_to_column(var = 'gene_id') %>% 
    write_tsv('normalized_gene_counts.tsv')
  '''
}
