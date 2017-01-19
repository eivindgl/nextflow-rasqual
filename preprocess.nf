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

process rename_gene_counts {
  publishDir params.preprocessing_dir, mode: 'copy'
  module 'R/3.3.1-foss-2015b'
  input:
    file count_table_raw from rnaCounts
    file sample_bam_map from RnaWithGenotypeList
  output:
    file 'gene_counts.tsv' into rnaCountsFinal
  """
  rename_gene_count_file.R $count_table_raw $sample_bam_map gene_counts.tsv
  """
}

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


