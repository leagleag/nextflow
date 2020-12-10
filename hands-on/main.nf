#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.genome = "$baseDir/data/genome.fa"
params.variants   = "$baseDir/data/known_variants.vcf.gz"
params.blacklist  = "$baseDir/data/blacklist.bed"
params.reads      = "$baseDir/data/reads/ENCSR000COQ1_{1,2}.fastq.gz"

/*
params.results    = "results"
*/

include {
  p1A_prepare_genome_samtools;
  p1B_prepare_genome_picard;
  p1C_prepare_star_genome_index;
  p1D_prepare_vcf_file;
  } from './processes1.nf'

include {
  p2_rnaseq_mapping_star;
} from './processes2.nf'

include {
  p3_rnaseq_gatk_splitNcigar;
} from './processes3.nf'

workflow {
  p1A_prepare_genome_samtools(params.genome)
  genome_index_sam = p1A_prepare_genome_samtools.out
  p1B_prepare_genome_picard(params.genome)
  genome_dict_picard = p1B_prepare_genome_picard.out
  p1C_prepare_star_genome_index(params.genome)
  genome_index_star = p1C_prepare_star_genome_index.out
  p1D_prepare_vcf_file(params.variants, params.blacklist)

  read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
  p2_rnaseq_mapping_star(params.genome, genome_index_star, read_pairs_ch)

  p3_rnaseq_gatk_splitNcigar(
    params.genome,
    genome_index_sam,
    genome_dict_picard,
    p2_rnaseq_mapping_star.out,
  )
}
