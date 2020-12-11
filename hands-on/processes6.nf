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

