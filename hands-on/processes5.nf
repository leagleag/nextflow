#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.gatk = "/opt/broad/GenomeAnalysisTK.jar"
GATK = params.gatk

/*
 * Process 5: variant calling of each sample using GATK
 */
process "p5_rnaseq_call_variants" {
  tag "$sampleId"

  input:
  path(genome)
  path(genome_index_sam)
  path(genome_picard_dict)
  tuple val(sampleId), path(bam), path(bai)

  output:
  tuple val(sampleId), path("final.vcf")

  script:
  """
  docker run -w \$(pwd) --volumes-from workspace cbcrg/callings-with-gatk:latest bash -c '\
  echo "${bam.join('\n')}" > bam.list

  # Variant calling
  java -jar $GATK -T HaplotypeCaller \
                  -R $genome -I bam.list \
                  -dontUseSoftClippedBases \
                  -stand_call_conf 20.0 \
                  -o output.gatk.vcf.gz

  # Variant filtering
  java -jar $GATK -T VariantFiltration \
                  -R $genome -V output.gatk.vcf.gz \
                  -window 35 -cluster 3 \
                  -filterName FS -filter "FS > 30.0" \
                  -filterName QD -filter "QD < 2.0" \
                  -o final.vcf
  '
  """
}
