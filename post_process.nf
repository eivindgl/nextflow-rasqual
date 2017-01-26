CeD_LD_SNPs = file('input_data/CeD_LD_SNPs_Iris_maf-.001_r2-9.bed')
gene_counts = file('out/preprocessing/gene_counts.tsv')

Channel.fromPath('out/timepoint/rasqual_full_output/*_txt')
  .map { 
  (full_string, timepoint, run_mode) = (it.baseName =~ /time_([^_]+)_([^_]+).*/)[0]
  [run_mode, timepoint, it]
}
  .into { CeDSubInputCh }

Channel.fromPath('out/timepoint/rasqual_full_output/time_*_lead_SNPs.rasqual_txt')
  .map { 
  (full_string, timepoint, run_mode) = (it.baseName =~ /time_([^_]+)_([^_]+).*/)[0]
  [run_mode, timepoint, it]
}
  .into { leadSnpRawCh }

process subset_CeD_LD_SNPs_from_rasqual_raw {
  module 'Python/3.5.1-foss-2015b'
  input:
    file CeD_LD_SNPs
    set run_mode, timepoint, 'rasqual_output.txt' from CeDSubInputCh
  output:
    set run_mode, timepoint, 'rasqual_CeD_LD_output.txt' into CedSubCh
  """
  filter_rasqual_by_SNP_subset.py $CeD_LD_SNPs -i rasqual_output.txt -o rasqual_CeD_LD_output.txt
  """
}

process add_pvalue_FDR_and_annotations_genome_wide_for_lead_all {
  publishDir "${params.post_proc_dir}/lead_genome_wide_fdr", mode: 'copy'
  module 'R/3.3.1-foss-2015b'
  input:
    set run_mode, timepoint, file('rasqual_output.txt') from leadSnpRawCh
  output:
    file "${run_mode}_SNP_time${timepoint}_genome-wide.tsv" into lead_SNP_all_results
  """
  tsv="${run_mode}_SNP_time${timepoint}_genome-wide.tsv"
  compute_pvalue_and_FDR.R rasqual_output.txt \$tsv lead-gw all
  filter_by_FDR_and_sort.R \$tsv \$tsv 0.05
  """
}

process add_pvalue_FDR_and_annotations {
  publishDir "${params.post_proc_dir}/CeD_LD_subset", mode: 'copy'
  module 'R/3.3.1-foss-2015b'
  input:
    set run_mode, timepoint, 'rasqual_output.txt' from CedSubCh
  output:
    file "LD_subset_${run_mode}_${timepoint}.tsv" into CedSubResultsPre
  """
  compute_pvalue_and_FDR.R rasqual_output.txt "LD_subset_${run_mode}_${timepoint}.tsv" $run_mode $timepoint
  """
}
CedSubResultsPre.into { CedSubResults ; CedMinGenePre }

process recompute_FDR_with_best_SNP_per_gene_only{
  publishDir "${params.post_proc_dir}/CeD_LD_subset_1gene", mode: 'copy'
  module 'R/3.3.1-foss-2015b'
  input:
    file tsv from CedMinGenePre.filter{ it.baseName =~ /_all_/ }
  output:
    file("best_SNP_${tsv}") into bestLDSnpGeneCh
  """
  #!/usr/bin/env Rscript
  stopifnot(getRversion() >= "3.2.0")
  if (!require("pacman")) install.packages("pacman")
  pacman::p_load(
    tidyverse
  )
  df <- read_tsv("$tsv") %>% 
    group_by(ensembl_gene_id) %>% 
    slice(which.min(pval)) %>% 
    ungroup() %>% 
    mutate(
      FDR = p.adjust(pval, method = 'BH')
    ) %>% 
    arrange(pval) %>% 
    write_tsv("best_SNP_${tsv}")
  """
}

process merge_FDR_sig_genes {
  publishDir "${params.post_proc_dir}/CeD_LD_subset_1gene", mode: 'copy'
  module 'R/3.3.1-foss-2015b'
  input:
    file('tsv/*') from bestLDSnpGeneCh.toSortedList()
  output:
    file('best_SNP_Gene_merged.tsv') into MergedLDSigCh
  """
  #!/usr/bin/env Rscript
  stopifnot(getRversion() >= "3.2.0")
  if (!require("pacman")) install.packages("pacman")
  pacman::p_load(
    tidyverse,
    purrr
  )

  list.files("tsv", full.names = TRUE) %>% 
    map(~ read_tsv(.x, col_types = cols(
      .default = col_guess(),
      timepoint = 'c'
    ))) %>% 
    bind_rows() %>% 
    filter(FDR < 0.05) %>% 
    arrange(external_gene_id, run_mode, timepoint) %>% 
    write_tsv("best_SNP_Gene_merged.tsv")
  """
}
    
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
