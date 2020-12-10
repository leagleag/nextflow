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

/*
 * Process 1B: Create a FASTA genome sequence dictionary with Picard for GATK
 */
process "p1B_prepare_genome_picard" {

  input:
  path genome

  output:
  path "${genome.baseName}.dict"

  script:
  """
  docker run -w \$(pwd) --volumes-from workspace cbcrg/callings-with-gatk:latest bash -c 'PICARD=`which picard.jar`; java -jar \$PICARD CreateSequenceDictionary R= $genome O= ${genome.baseName}.dict'
  """
}

/*
 * Process 1C: Create the genome index file for STAR
 */
process "p1C_prepare_star_genome_index" {

  input:
  path genome

  output:
  path genome_dir

  script:
  """
  docker run -w \$(pwd) --volumes-from workspace cbcrg/callings-with-gatk:latest bash -c '
  mkdir genome_dir

  STAR --runMode genomeGenerate \
       --genomeDir genome_dir \
       --genomeFastaFiles ${genome} \
       --runThreadN ${task.cpus}
  '
  """
}

/*
 * Process 1D: Create a file containing the filtered and recoded set of variants
 */
process "p1D_prepare_vcf_file" {

  input:
  path variantsFile
  path blacklisted

  output:
  tuple path("${variantsFile.baseName}.filtered.recode.vcf.gz"), \
  path("${variantsFile.baseName}.filtered.recode.vcf.gz.tbi")

  script:
  """
  docker run -w \$(pwd) --volumes-from workspace cbcrg/callings-with-gatk:latest bash -c '
  vcftools --gzvcf $variantsFile -c \
           --exclude-bed ${blacklisted} \
           --recode | bgzip -c \
           > ${variantsFile.baseName}.filtered.recode.vcf.gz

  tabix ${variantsFile.baseName}.filtered.recode.vcf.gz
  '
  """
}
