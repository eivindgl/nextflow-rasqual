#!/usr/bin/env Rscript
pacman::p_load(
  tidyverse,
  forcats,
  stringr,
  VariantAnnotation
)

args <- commandArgs(trailingOnly = TRUE)
path <- list(
  norm_expr = args[1],
  eqtl_list = args[2],
  vcf = args[3],
  out_dir = args[4]
)
vst <- read_tsv(path$norm_expr) %>%
  gather(sample, vst, -gene_id)
eqtl <- read_tsv(path$eqtl_list)
vcf <- readVcf(path$vcf, genome = 'GRCh37')
snps <- unique(eqtl$rs_id)
gt <- geno(vcf)$GT %>%
  as.data.frame %>%
  rownames_to_column('SNP') %>%
  as_tibble %>%
  filter(SNP %in% snps) %>%
  gather(sample, genotype, -SNP) %>%
  mutate(genotype = fct_recode(genotype,
                               ref = '0|0',
                               heterozygous = '0|1',
                               heterozygous = '1|0',
                               alt = '1|1'))

gene_snp_map <- eqtl %>% 
  dplyr::select(external_gene_id, gene_id, SNP = rs_id)

df <- gene_snp_map %>%
  inner_join(vst) %>%
  inner_join(gt) %>% 
  dplyr::rename(hgnc = external_gene_id)

# print highly expressed genes
df %>%
  group_by(hgnc, gene_id) %>%
  summarise(n = n(), median = median(vst), mean = mean(vst)) %>%
  arrange(desc(median))

plot_eQTL <- function(df, gene) {
  df %>%
    filter(hgnc == gene) %>%
    ggplot(aes(genotype, vst, color = genotype)) +
    geom_boxplot() +
    geom_jitter(width = 0.2) +
    ggtitle(gene)
}

dir.create(path$out_dir, showWarnings = FALSE)
for (gene in unique(df$hgnc)) {
  p <- plot_eQTL(df, gene)
  outpath <- file.path(path$out_dir, paste(gene, '.png', sep = ''))
  ggsave(outpath, plot = p)
}
