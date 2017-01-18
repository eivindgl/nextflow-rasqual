params.vcfDir = 'input_data/vcf_by_chrom'
params.refVcfDict = '/apps/data/ftp.broadinstitute.org/bundle/2.8/b37/human_g1k_v37.dict'
params.bamDir = '/groups/umcg-pub/tmp04/projects/stimulated_gluten_specific_Tcell_clones_TCC23052016/pipelines/splice_junctions/results/sortedBam'
params.GTF = '/apps/data/ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.annotation.gtf'


vcfFiles = Channel.fromPath("${params.vcfDir}/*.vcf")
rnaBamFiles = Channel
		.fromPath("${params.bamDir}/*.bam")
		.map { file -> tuple(file.baseName, file) }

/* Merges and sorts per_chrom vcf files into a single genome-wide vcf file.
 * The sorting is necessary for allele specific read counting.
 */
process unify_input_vcf {
  input:
    file vcfs from vcfFiles.toList()
  output:
    file 'all.vcf' into vcf 
  module 'picard'
  module 'VCFtools'
  """
  vcf-concat $vcfs > merged_unsorted.vcf
  java -jar ${EBROOTPICARD}/picard.jar SortVcf \
    I=merged_unsorted.vcf \
    O="all.vcf" \
    SEQUENCE_DICTIONARY="$params.refVcfDict"
  rm -f merged_unsorted.vcf
  """
}

/* Counts reads per gene for all bam samples with regard to strand direction
 */
process count_reads {
  input:
    set raw_name, file(sample) from rnaBamFiles
  output:
    set raw_name, "counts.txt" into rnaCounts
  module 'Subread'
  executor 'slurm'
  memory '512 MB'
  cpus 1
  time '1h'

  """
  featureCounts -T ${task.cpus} -s 1 -t exon -g gene_id -a "$params.GTF" \
    -o counts.txt $sample
  """
}

process readVcfSampleNames {
  input:
    file vcf
  output:
    file 'vcf_samples.list' into vcfSamples
  module 'VCFtools'
  """
  vcf-query -l $vcf > vcf_samples.list
  """
}


//rnaCounts.subscribe { name, path -> println "RNA Got counts for: $name ($path)"}
vcf.subscribe { println "Got: $it" }
