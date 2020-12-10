#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.gatk = "/opt/broad/GenomeAnalysisTK.jar"
GATK = params.gatk

/*
 * Process 4: uses GATK to detect systematic errors in the base quality scores, select unique alignments and then index the resulting bam file with samtools
 */
process "p4_rnaseq_gatk_recalibrate" {
  tag "$replicateId"

  input:
  path(genome)
  path(genome_index_sam)
  path(genome_picard_dict)
  tuple val(replicateId), path(bam), path(bam_index)
  tuple path(recoded_vcf), path(tab_index)

  output:
  tuple val(sampleId), path("${replicateId}.final.uniq.bam"), \
  path("${replicateId}.final.uniq.bam.bai")

  script:
  sampleId = replicateId.replaceAll(/[12]$/,'')
  """
  docker run -w \$(pwd) --volumes-from workspace cbcrg/callings-with-gatk:latest bash -c '\
  # Indel Realignment and Base Recalibration
  java -jar $GATK -T BaseRecalibrator \
                --default_platform illumina \
                -cov ReadGroupCovariate \
                -cov QualityScoreCovariate \
                -cov CycleCovariate \
                -knownSites ${recoded_vcf} \
                -cov ContextCovariate \
                -R ${genome} -I ${bam} \
                --downsampling_type NONE \
                -nct ${task.cpus} \
                -o final.rnaseq.grp

   java -jar $GATK -T PrintReads \
                -R ${genome} -I ${bam} \
                -BQSR final.rnaseq.grp \
                -nct ${task.cpus} \
                -o final.bam

    # Select only unique alignments, no multimaps
    (samtools view -H final.bam; samtools view final.bam | grep -w 'NH:i:1') \
    | samtools view -Sb -  > ${replicateId}.final.uniq.bam

    # Index BAM files
    samtools index ${replicateId}.final.uniq.bam
    '
    """
}
