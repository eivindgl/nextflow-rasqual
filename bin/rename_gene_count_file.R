#!/usr/bin/env Rscript
stopifnot(getRversion() >= "3.2.0")

if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  tidyverse,
  stringr
)
args <- commandArgs(trailingOnly = TRUE)

path <- list(
  count_table = args[1],
  sample_bam_map = args[2],
  output = args[3]
)
print("Script arguments are:")
print(path)

rename_columns_with_bad_names <- function(df, old_name, new_name) {
  # as.name quotes ugly names
  src_RID <- map(old_name, as.name)
  nm <- setNames(src_RID, new_name)
  df %>% rename_(.dots = nm)
}

sample_to_bamid <- read_tsv(path$sample_bam_map)

# First line is just the shell input command
df <- read_tsv(path$count_table, skip = 1)

# Maps (complex) column names to (simple) sample names
bam_name <- tibble(bam = colnames(df)[7:ncol(df)])
name_mapping <- right_join(bam_name, sample_to_bamid)


df <- df %>%
  # select proper but wrongly named subset
  select(gene_id = 1, one_of(name_mapping$bam)) %>%
  mutate(ensembl_gene_id = str_extract(gene_id, 'ENSG\\d+')) %>%
  select(gene_id, ensembl_gene_id, everything()) %>%
  rename_columns_with_bad_names(
    name_mapping$bam, name_mapping$sample_name)

df %>% write_tsv(path$output)
