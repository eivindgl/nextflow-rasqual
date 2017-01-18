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

process count_ASE {
  module GATK
  cpus 2
  memory '6 GB'
  time '4 h'

    """
java \
   -XX:ParallelGCThreads=1  -Xmx4g \
   -Djava.io.tmpdir="\$java_tmp" \
   -jar \${EBROOTGATK-dummy}/GenomeAnalysisTK.jar \
   -R "$genome_ref" \
   -T ASEReadCounter \
   -o "$output_csv" \
   -I "$input_bam" \
   -sites "$input_vcf" \
   -L "$input_vcf" \
   -U ALLOW_N_CIGAR_READS \
   -dt NONE \
   --minMappingQuality 10


    """
}


