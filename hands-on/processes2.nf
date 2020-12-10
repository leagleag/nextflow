#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Process 2: Align RNA-Seq reads to the genome with STAR
 */
process "p2_rnaseq_mapping_star" {

  input:
  path(genome)
  path(genome_index)
  tuple val(replicateId), path(paired_reads)

  output:
  tuple val(replicateId), \
  path("Aligned.sortedByCoord.out.bam"), \
  path("Aligned.sortedByCoord.out.bam.bai")

  script:
  """
  docker run -w \$(pwd) --volumes-from workspace cbcrg/callings-with-gatk:latest bash -c '\
  # ngs-nf-dev Align reads to genome
  STAR --genomeDir $genome_index \
       --readFilesIn $paired_reads \
       --runThreadN ${task.cpus} \
       --readFilesCommand zcat \
       --outFilterType BySJout \
       --alignSJoverhangMin 8 \
       --alignSJDBoverhangMin 1 \
       --outFilterMismatchNmax 999

  # 2nd pass (improve alignments using table of splice junctions and create a new index)
  mkdir genomeDir
  STAR --runMode genomeGenerate \
       --genomeDir genomeDir \
       --genomeFastaFiles $genome \
       --sjdbFileChrStartEnd SJ.out.tab \
       --sjdbOverhang 75 \
       --runThreadN ${task.cpus}

  # Final read alignments
  STAR --genomeDir genomeDir \
       --readFilesIn $paired_reads \
       --runThreadN ${task.cpus} \
       --readFilesCommand zcat \
       --outFilterType BySJout \
       --alignSJoverhangMin 8 \
       --alignSJDBoverhangMin 1 \
       --outFilterMismatchNmax 999 \
       --outSAMtype BAM SortedByCoordinate \
       --outSAMattrRGline ID:$replicateId LB:library PL:illumina PU:machine SM:GM12878

  # Index the BAM file
  samtools index Aligned.sortedByCoord.out.bam
  '
  """
}
