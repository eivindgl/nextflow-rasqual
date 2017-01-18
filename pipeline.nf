params.bamDir = '/groups/umcg-pub/tmp04/projects/stimulated_gluten_specific_Tcell_clones_TCC23052016/pipelines/splice_junctions/results/sortedBam'


rnaBamFiles = Channel
		.fromPath("${params.bamDir}/*.bam")
		.map { file -> tuple(file.baseName, file) }

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
