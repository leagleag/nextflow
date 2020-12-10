#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.genome = "$baseDir/data/genome.fa"
params.variants   = "$baseDir/data/known_variants.vcf.gz"
params.blacklist  = "$baseDir/data/blacklist.bed"

/*
params.reads      = "$baseDir/data/reads/ENCSR000COQ1_{1,2}.fastq.gz"
params.results    = "results"
params.gatk       = "/opt/broad/GenomeAnalysisTK.jar"
reads_ch          = Channel.fromFilePairs(params.reads)
GATK              = params.gatk
*/

include {
  p1A_prepare_genome_samtools;
  p1B_prepare_genome_picard;
  p1C_prepare_star_genome_index;
  p1D_prepare_vcf_file;
  } from './processes1.nf'

workflow {
  p1A_prepare_genome_samtools(params.genome)
  p1B_prepare_genome_picard(params.genome)
  p1C_prepare_star_genome_index(params.genome)
  p1D_prepare_vcf_file(params.variants, params.blacklist)
}
