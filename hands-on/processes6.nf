#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.results = "results"

/*
 * Process 6A: process the VCF result to prepare variants file for allele
 * specific expression (ASE) analysis
 */
process "p6A_post_process_vcf" {
  tag "$sampleId"
  publishDir "${params.results}/$sampleId"

  input:
  tuple val(sampleId), path("final.vcf")
  tuple path(recoded_vcf), path(tab_index)

  output:
  tuple val(sampleId), path("final.vcf"), path("commonSNPs.diff.sites_in_files")

  script:
  """
  docker run -w \$(pwd) --volumes-from workspace cbcrg/callings-with-gatk:latest bash -c '
  grep -v "#" final.vcf | awk "\\\$7~/PASS/" | perl -ne "chomp(\\\$_); (\\\$dp)=\\\$_=~/DP\\\\=(\\\\d+)\\\\;/; if(\\\$dp>=8){print \\\$_.\\"\\\\n\\"};" > result.DP8.vcf
  cat result.DP8.vcf
  vcftools --vcf result.DP8.vcf --gzdiff $recoded_vcf --diff-site --out commonSNPs
  '
  """
}

/*
 * Process 6B: prepare the VCF for allele specific expression (ASE) and
 * generate a figure in R.
 */
process "p6B_prepare_vcf_for_ase" {
  tag "$sampleId"
  publishDir "${params.results}/$sampleId"

  input:
  tuple val(sampleId), path(vcf), path(snps)

  output:
  tuple val(sampleId), path("known_snps.vcf"), emit: known_snps
  path("AF.histogram.pdf"), emit: figure

  script:
  """
  docker run -w \$(pwd) --volumes-from workspace cbcrg/callings-with-gatk:latest bash -c '\

  awk "BEGIN{OFS=\\"\\\\t\\"} \\\$4~/B/{print \\\$1,\\\$2,\\\$3}" commonSNPs.diff.sites_in_files > test.bed

  vcftools --vcf final.vcf --bed test.bed --recode --keep-INFO-all --stdout > known_snps.vcf

  grep -v "#"  known_snps.vcf | awk -F "\\\\t" "{print \\\$10}" \
               | awk -F ":" "{print \\\$2}" | perl -ne "chomp(\\\$_); \
               @v=split(/\\\\,/,\\\$_); if(\\\$v[0]!=0 || \\\$v[1]!=0) \
               {print \\\$v[1]/(\\\$v[1]+\\\$v[0]).\\"\\\\n\\"; }" |awk "\\\$1!=1" \
               > AF.4R

  Rscript ${script_path}/gghist.R -i AF.4R -o AF.histogram.pdf
  '
  """
}

/*
 * Process 6C: calculate allele counts at a set of positions with GATK tools
 */
process "p6C_ASE_knownSNPs" {
  tag "$sampleId"
  publishDir "${params.results}/$sampleId", mode: "copy", overwrite: true

  input:
  path(genome)
  path(genome_index_sam)
  path(genome_picard_dict)
  tuple val(sampleId), path(vcf), path(bam), path(bai)

  output:
  path("ASE.tsv")

  script:
  """
  docker run -w \$(pwd) --volumes-from workspace cbcrg/callings-with-gatk:latest bash -c '\
  echo "${bam.join('\n')}" > bam.list

  java -jar $GATK -R $genome \
                  -T ASEReadCounter \
                  -o ASE.tsv \
                  -I bam.list \
                  -sites $vcf
  '
  """
}
