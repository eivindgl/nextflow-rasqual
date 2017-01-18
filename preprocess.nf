
process list_vcf_samples {
  module 'VCFtools'
  input:
    file vcf from Channel.fromPath("$params.vcf")
  output:
    file 'sampleIDs.list' into sampleIDs
  """
  vcf-query -l $vcf > sampleIDs.list
  """
}

process bind_genotype_to_rna {
  publishDir params.preprocessing_dir// , mode: 'copy'
  input:
    file sampleIDs
    val timestring from params.timepoints.toList().join(' ')

  output:
    file 'sample_time_bam_map.tsv' into RnaWithGenotypeList

  """
  bind_genotype_and_rnaseq.p6 ${sampleIDs.name} ${params.bamDir} ${timestring} > sample_time_bam_map.tsv
  """

}

RnaWithGenotypeList.subscribe { println "got $it" }


