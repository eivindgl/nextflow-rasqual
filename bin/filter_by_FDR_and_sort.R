#!/usr/bin/env Rscript
pacman::p_load(
  tidyverse
)

args <- commandArgs(trailingOnly = TRUE)

tsv_in = args[1]
tsv_out = args[2]
fdr_lim = as.numeric(args[3])

read_tsv(tsv_in) %>% 
  filter(FDR < fdr_lim) %>% 
  arrange(FDR) %>% 
  write_tsv(tsv_out)
