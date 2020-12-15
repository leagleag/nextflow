#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.gatk = "/opt/broad/GenomeAnalysisTK.jar"
GATK = params.gatk

/*
 * Process 3: split reads on Ns in CIGAR string using GATK
 */
process "p3_rnaseq_gatk_splitNcigar" {
  tag "$replicateId"

  input:
  path(genome)
  path(genome_index_sam)
  path(genome_picard_dict)
  tuple val(replicateId), path(bam), path(bam_index)

  output:
  tuple val(replicateId), path("split.bam"), \
  path("split.bai")

  script:
  """
  docker run -w \$(pwd) --volumes-from workspace cbcrg/callings-with-gatk:latest bash -c '\

  # SplitNCigarReads and reassign mapping qualities
  java -jar $GATK -T SplitNCigarReads \
                  -R $genome -I $bam \
                  -o split.bam \
                  -rf ReassignOneMappingQuality \
                  -RMQF 255 -RMQT 60 \
                  -U ALLOW_N_CIGAR_READS \
                  --fix_misencoded_quality_scores\
  '
  """
}
