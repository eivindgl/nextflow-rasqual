vcfFiles = Channel.fromPath("${params.vcfDir}/*.vcf")

/* Merges and sorts per_chrom vcf files into a single genome-wide vcf file.
 * The sorting is necessary for allele specific read counting.
 */
process unify_input_vcf {
  publishDir params.preprocessing_dir, mode: 'copy'
  input:
    file vcfs from vcfFiles.toList()
  output:
    file 'all.vcf' into vcf_stream 
  module 'picard'
  module 'VCFtools'
  """
  vcf-concat $vcfs > merged_unsorted.vcf
  java -jar \${EBROOTPICARD}/picard.jar SortVcf \
    I=merged_unsorted.vcf \
    O="all.vcf" \
    SEQUENCE_DICTIONARY="$params.refVcfDict"
  rm -f merged_unsorted.vcf
  """
}
vcf = vcf_stream.first()

process list_vcf_samples {
  module 'VCFtools'
  input:
    file vcf
  output:
    file 'sampleIDs.list' into sampleID_stream
  """
  vcf-query -l $vcf > sampleIDs.list
  """
}
sampleIDs = sampleID_stream.first()

process bind_genotype_to_rna {
  publishDir params.preprocessing_dir, mode: 'copy'
  input:
    file sampleIDs
    val timestring from params.timepoints.toList().join(' ')

  output:
    file 'sample_time_bam_map.tsv' into RnaWithGenotypeStream
    file 'bam_paths.txt' into bamPathStream

  """
  bind_genotype_and_rnaseq.p6 ${sampleIDs.name} ${params.bamDir} ${timestring} > sample_time_bam_map.tsv
  tail -n +2 sample_time_bam_map.tsv | cut -f 4 | tr '\n' ' ' > bam_paths.txt
  """
}
bamPaths = bamPathStream.first()
RnaWithGenotypeList = RnaWithGenotypeStream.first()
GtfFile = file(params.GTF)

process count_reads {
  publishDir params.preprocessing_dir, mode: 'copy'
  module 'Subread'
  executor 'slurm'
  memory '40 GB'
  cpus 15
  time '15 h'
  input:
    file bampaths from bamPaths
  output:
    file "gene_counts.txt" into rnaCounts
  """
  bam_paths=\$(cat $bampaths)
  featureCounts -T ${task.cpus} -s 1 -t exon -g gene_id -a "$params.GTF" \
    -o gene_counts.txt \$bam_paths
  """
}

process compute_exon_unions {
  publishDir params.preprocessing_dir, mode: 'copy'
  module 'R/3.3.1-foss-2015b'
  input:
    file GtfFile
  output:
    file 'exon_unions.tsv' into exonUnions
  """
  compute_exon_unions.R $GtfFile exon_unions.tsv
  """
}
exon_unions = exonUnions.first()



process rename_gene_counts {
  publishDir params.preprocessing_dir, mode: 'copy'
  module 'R/3.3.1-foss-2015b'
  input:
    file raw_count_table from rnaCounts
    file sample_bam_map from RnaWithGenotypeList
  output:
    file 'gene_counts.tsv' into rnaCountsFinal
  """
  rename_gene_count_file.R $raw_count_table $sample_bam_map gene_counts.tsv
  """
}
gene_counts = rnaCountsFinal.first()

process get_gene_gc_content_from_biomaRt {
  publishDir params.preprocessing_dir, mode: 'copy'
  module 'R/3.3.1-foss-2015b'
  input:
  output:
    file 'gene_gc_prct.tsv' into geneGcPercentageCh
  """
#!/usr/bin/env Rscript
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  biomaRt,
  tidyverse
)
ensembl_version = 'feb2014.archive.ensembl.org'
mart <- useMart(
    'ENSEMBL_MART_ENSEMBL',
    host = ensembl_version,
    dataset = 'hsapiens_gene_ensembl')
getBM(
      attributes = c('ensembl_gene_id', 'percentage_gc_content'), 
      mart = mart) %>%
  write_tsv('gene_gc_prct.tsv')
  """
}

process extract_SNP_coordinates {
  publishDir params.preprocessing_dir, mode: 'copy'
  module 'R/3.3.1-foss-2015b'
  input:
    file vcf
  output:
    file 'snp_coordinates.tsv' into SnpCoordinatesCh
  """
#!/usr/bin/env Rscript
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  tidyverse,
  VariantAnnotation
)
readVcf("$vcf", genome='GRCh37') %>%
  rowRanges() %>%
  as.data.frame() %>%
  rownames_to_column(var = 'snp_id') %>%
  dplyr::select(chr = seqnames, pos = start, snp_id) %>%
  write_tsv('snp_coordinates.tsv')
  """
}
snp_coord = SnpCoordinatesCh.first()

RnaWithGenotypeList
  .splitCsv(sep: '\t', header: true)
  .map { it.values() }
  .set { RnaGenoItems }

process count_ASE {
  publishDir params.preprocessing_ase, mode: 'copy'
  module 'GATK'
  executor 'slurm'
  cpus 2
  memory '6 GB'
  time '4 h'
  input:
    file vcf
    set sample_id, sample_name, timepoint, bam from RnaGenoItems
  output:
    set sample_id, sample_name, timepoint, "${sample_id}_ASE_count.csv" into AseCounts
    // This should be included to play nice, but I am in a hurry
    // -Djava.io.tmpdir="\$java_tmp" \
    """
    java \
     -XX:ParallelGCThreads=1  -Xmx4g \
     -jar \${EBROOTGATK-dummy}/GenomeAnalysisTK.jar \
     -R "$params.genome_ref" \
     -T ASEReadCounter \
     -o "${sample_id}_ASE_count.csv" \
     -I "$bam" \
     -sites "$vcf" \
     -L "$vcf" \
     -U ALLOW_N_CIGAR_READS \
     -dt NONE \
     --minMappingQuality 10
    """
}

// Creates a file with sample_name -> ASE_count_path, which is required for
// the mergeAseCounts scripts
AseCounts
  .map{"${it[0]}\t${it[3]}"}
  .collectFile(name: 'sampleID_to_ASECountPath.tsv', newLine: true)
  .set { ASEPathFile }

process merge_ASE_counts {
  publishDir params.preprocessing_dir, mode: 'copy'
  input:
    file sample_list from ASEPathFile
  output:
    file 'ASE_counts.tsv' into mergedASECounts

  """
  kauralasoo_merge_ASE_counts.py --sample_list $sample_list > ASE_counts.tsv
  """
}


/* Feature SNPs are within exons (used to estimate ASE)
 * Regulatory SNPs are outside exons but within a CIS-window of the gene
 */
process count_regulatory_and_feature_SNPs_per_gene_cis_window {
  publishDir params.preprocessing_dir, mode: 'copy'
  module 'R/3.3.1-foss-2015b'
  input:
    file exon_unions
    file snp_coord
  output:
    file 'gene_cis_snp_count.tsv' into geneCisSnpCountCh
  """
  #!/usr/bin/env Rscript
  if (!require("pacman")) install.packages("pacman")
  pacman::p_load(
    tidyverse,
    devtools
  )
  devtools::install_github("kauralasoo/rasqual/rasqualTools")
  library(rasqualTools)

  exon_unions <- read_tsv("$exon_unions", col_types = cols(
    .default = col_guess(), 
    exon_starts = 'c',
    exon_ends = 'c'))
  snp_coords  <- read_tsv("$snp_coord")
  countSnpsOverlapingExons(exon_unions, snp_coords, 
			   cis_window = ${params.snp_window_size}) %>%
    write_tsv('gene_cis_snp_count.tsv')
  """
}

