/* 
 * Input files are generated by preprocess.nf
 */
params.orig_vcf = 'out/preprocessing/all.vcf'
orig_vcf = file(params.orig_vcf)

params.gene_counts= 'out/preprocessing/gene_counts.tsv'
gene_counts = file(params.gene_counts)

params.ase_counts= 'out/preprocessing/ASE_counts.tsv'
ase_counts = file(params.ase_counts)

params.gene_gc_content = 'out/preprocessing/gene_gc_prct.tsv'
gene_gc = file(params.gene_gc_content)

params.gene_cis_snp_count = 'out/preprocessing/gene_cis_snp_count.tsv'
gene_cis_snp_count = file(params.gene_cis_snp_count)

params.ced_ld_snps = 'input_data/CeD_LD_SNPs_Iris_maf-.001_r2-9.bed'
CeD_LD_SNPs = file(params.ced_ld_snps)

process create_vcf_copy_with_timepoint_prefix {
  module 'VCFtools'
  module 'BCFtools'
  input:
    file orig_vcf from orig_vcf
    val timepoint from params.timepoints
  output:
    set timepoint, "time_${timepoint}.vcf" into timepointVcfCh
  shell:
  '''
  dst="time_!{timepoint}.vcf"
  mkdir -p tmp
  tmp_dir=$(mktemp -d -p tmp)
  tmp=$(mktemp -p tmp)
  src=$tmp_dir/orig.vcf.gz
  bgzip -c !{orig_vcf} > $src
  bcftools query -l $src | awk -v g=time!{timepoint}_ '{ print g $0 }'  > $tmp
  bcftools reheader -s $tmp $src > ${dst}.gz
  bcftools convert -O v -o $dst ${dst}.gz
  rm -fr $tmp_dir ${dst}.gz # deletes $tmp and $tmp_dir
  '''
}
timepointVcfCh.into { allTimepointsInputCh; timepointsVcfCh }

process compress_and_index_vcf_prior_to_merge {
  module 'BCFtools'
  input:
    set timepoint, file('tp.vcf') from allTimepointsInputCh
  output:
    file "time_${timepoint}.vcf.gz" into allTimepointsInputGzCh
    file "time_${timepoint}.vcf.gz.csi" into allTimepointsInputGzIndexCh
  """
  bcftools convert -O z -o time_${timepoint}.vcf.gz tp.vcf
  bcftools index time_${timepoint}.vcf.gz
  """
}

process merge_all_vcf_timepoints {
  //cache 'deep'
  module 'BCFtools'
  input:
    file('vcf/*') from allTimepointsInputGzCh.toSortedList()
    file('vcf/*') from allTimepointsInputGzIndexCh.toSortedList()
  output:
    set val('all'), 'all_timepoints.vcf' into mergedVcfCh
  """
  bcftools merge -o all_timepoints.vcf -O v vcf/*.vcf.gz
  """
}

timepointsVcfCh
  .mix(mergedVcfCh) 
  .set{vcfCh}

process extract_RNA_to_genotype_mapfile {
  publishDir "${params.timepoint_base_dir}/time_${timepoint}", mode: 'copy'
  module 'VCFtools'
  input:
    set val(timepoint), file(vcf) from vcfCh
    file gene_counts
  output:
    set timepoint, file(vcf), "sample_genotype_map.list" into svcfCh
  shell:
  '''
  # crude check that each genotype has an entry in RNA count file
  for gt in $(vcf-query -l !{vcf}) ; do
    if ! grep -q $gt <(head -1 !{gene_counts}) ; then
      echo "Genotype $gt not present in gene counts file: !{gene_counts}"
      exit 1
    fi
  done
  vcf-query -l !{vcf} | awk '{ print $1 "\t" $1 }' > "sample_genotype_map.list"
  '''
}

process combine_genotype_and_ASE_counts {
  publishDir "${params.timepoint_base_dir}/time_${timepoint}", mode: 'copy'
  module 'VCFtools'
  input:
  set val(timepoint), file(vcf), file(sample_map) from svcfCh
  file ase_counts
  output:
  set timepoint, file(vcf), file(sample_map), 'ASE_genotype.vcf.gz', 'ASE_genotype.vcf.gz.tbi' into csvcfCh
  """
  # piping into bgzip masks errors, so we make a tmp file
  add_ASE_to_vcf.py \
    --ASEcounts ${ase_counts} \
    --ASESampleGenotypeMap ${sample_map} \
    --VCFfile ${vcf} > tmp.vcf
  bgzip tmp.vcf -c > ASE_genotype.vcf.gz
  tabix -p vcf ASE_genotype.vcf.gz
  rm tmp.vcf
  """
}

process generate_timepoint_count_and_factor_tables {
  publishDir "${params.timepoint_base_dir}/time_${timepoint}", mode: 'copy'
  module 'R/3.3.1-foss-2015b'
  input:
    file gene_gc
    set val(timepoint), file(vcf), file(sample_map), file(ASE_vcf), file(ASE_vcf_idx) from csvcfCh
  output:
    set val(timepoint), file(sample_map), file(ASE_vcf), file(ASE_vcf_idx),
	file('gene_counts.bin'), file('size_factors.bin'), file('gene_counts.tsv') into rasqualInCh
  """
  compute_exp_counts_and_factors.R $sample_map $gene_gc $gene_counts size_factors.bin gene_counts.bin gene_counts.tsv
  """
}

process split_into_rasqual_batches {
  input:
    set val(timepoint), file(sample_map), file(ASE_vcf), file(ASE_vcf_idx),
	file(gene_counts_bin), file(size_factors_bin), file(gene_counts_tsv) from rasqualInCh
  output:
    set val(timepoint), file(sample_map), file(ASE_vcf), file(ASE_vcf_idx),
	file(gene_counts_bin), file(size_factors_bin), file('geneids.txt'), 
	file('geneid_batch_*') into runRasqualCh mode flatten
  """
  tail -n +2 $gene_counts_tsv | cut -f 1 > geneids.txt
  split -l $params.batch_size geneids.txt geneid_batch_

  """
}

runRasqualCh.into { runLeadCh ; runAllCh }

process run_rasqual_all_SNPs {
  module 'GSL'
  module 'VCFtools'
  executor 'slurm'
  memory '1 GB'
  cpus 1
  time '2 h'
  input:
    file gene_cis_snp_count
    set val(timepoint), file(sample_map), file(ASE_vcf), file(ASE_vcf_idx),
	file(gene_counts_bin), file(size_factors_bin), file(geneids), file(geneid_batch) from runAllCh
  output:
    set val(timepoint), file("rasqual_raw_time${timepoint}*.txt") into rasqRawAllCh
  """
  n=\$(wc -l < $sample_map)
  echo -ne "${geneid_batch}\t" | cat - <( paste -sd "," $geneid_batch) | \
    runRasqual.py \
    --readCounts $gene_counts_bin \
    --offsets $size_factors_bin \
    --n \$n \
    --vcf $ASE_vcf \
    --outprefix "rasqual_raw_time${timepoint}" \
    --geneids $geneids \
    --geneMetadata $gene_cis_snp_count \
    --execute True \
    --rasqualBin rasqual
  """
}

process run_rasqual_lead_SNPs {
  module 'GSL'
  module 'VCFtools'
  executor 'slurm'
  memory '1 GB'
  cpus 1
  time '2 h'
  input:
    file gene_cis_snp_count
    set val(timepoint), file(sample_map), file(ASE_vcf), file(ASE_vcf_idx),
	file(gene_counts_bin), file(size_factors_bin), file(geneids), file(geneid_batch) from runLeadCh
  output:
    set val(timepoint), file("rasqual_raw_time${timepoint}.*.txt") into rasqRawLeadCh
  """
  n=\$(wc -l < $sample_map)
  echo -ne "${geneid_batch}\t" | cat - <( paste -sd "," $geneid_batch) | \
    runRasqual.py \
    --readCounts $gene_counts_bin \
    --offsets $size_factors_bin \
    --n \$n \
    --vcf $ASE_vcf \
    --outprefix "rasqual_raw_time${timepoint}" \
    --geneids $geneids \
    --geneMetadata $gene_cis_snp_count \
    --execute True \
    --rasqualBin rasqual \
    --parameters '\\--lead-snp'
  """
}

rasqRawLeadCh
  .collectFile(
    storeDir: "${params.timepoint_base_dir}/rasqual_full_output") {
    [ 
      "time_${it[0]}_lead_SNPs.rasqual_txt",
      it[1]
    ]
  }
  .set { RawLeadFullCh }

rasqRawAllCh
  .collectFile(
    storeDir: "${params.timepoint_base_dir}/rasqual_full_output") {
    [ 
      "time_${it[0]}_all_SNPs.rasqual_txt",
      it[1]
    ]
  }
  .set { RawAllFullCh }

RawLeadFullCh.into { CeDSubLeadInputCh ; GenomeWideInputCh }

CeDSubLeadInputCh
  .mix(RawAllFullCh)
  .map { 
  (full_string, timepoint, run_mode) = (it.baseName =~ /time_([^_]+)_([^_]+).*/)[0]
  [run_mode, timepoint, it]
}.into { CeDSubInputCh }

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

process add_pvalue_FDR_and_annotations {
  publishDir "${params.timepoint_base_dir}/CeD_LD_subset", mode: 'copy'
  module 'R/3.3.1-foss-2015b'
  input:
    set run_mode, timepoint, 'rasqual_output.txt' from CedSubCh
  output:
    file "LD_subset_${run_mode}_${timepoint}.tsv" into CedSubResults
  """
  compute_pvalue_and_FDR.R rasqual_output.txt "LD_subset_${run_mode}_${timepoint}.tsv"
  """
}
