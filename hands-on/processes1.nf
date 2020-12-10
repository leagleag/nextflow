#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Process 1A: Create a FASTA genome index with samtools
 */
process "p1A_prepare_genome_samtools" {

  input:
    path genome

  output:
    path "${genome}.fai"

  script:
  """
  pwd=\$(pwd)
  docker run -w \$pwd --volumes-from workspace cbcrg/callings-with-gatk:latest bash -c 'samtools faidx ${genome}'
  """
}
