<!-- MDTOC maxdepth:6 firsth1:1 numbering:1 flatten:0 bullets:0 updateOnSave:1 -->

1. [Summary](#summary)   
2. [Installation](#installation)   
3. [Sources](#sources)   
&emsp;3.1. [Material](#material)   
&emsp;3.2. [Products: nextflow tower](#products-nextflow-tower)   
4. [Concepts](#concepts)   
&emsp;4.1. [Containerisation](#containerisation)   
&emsp;4.2. [Nextflow and cloud deployment](#nextflow-and-cloud-deployment)   
&emsp;&emsp;4.2.1. [AWS batch and nextflow tower](#aws-batch-and-nextflow-tower)   
5. [Command line](#command-line)   
6. [Code](#code)   
&emsp;6.1. [Process](#process)   
&emsp;6.2. [Channels and processes](#channels-and-processes)   
&emsp;6.3. [Parameters](#parameters)   
&emsp;6.4. [Environment and profiles scopes](#environment-and-profiles-scopes)   
&emsp;6.5. [work directory and resume](#work-directory-and-resume)   
&emsp;6.6. [Error management](#error-management)   
7. [DSL2](#dsl2)   
&emsp;7.1. [Processes are modular](#processes-are-modular)   
&emsp;7.2. [Workflow (take/main/emit)](#workflow-takemainemit)   
&emsp;7.3. [Pipes](#pipes)   
&emsp;7.4. [Transition from dsl to dsl2](#transition-from-dsl-to-dsl2)   
&emsp;&emsp;7.4.1. [Some guidelines](#some-guidelines)   
&emsp;&emsp;7.4.2. [A workflow example](#a-workflow-example)   
8. [Bioinformatics pipeline example](#bioinformatics-pipeline-example)   
&emsp;8.1. [data description](#data-description)   
&emsp;&emsp;8.1.1. [genome assembly `genome.fa`](#genome-assembly-genomefa)   
&emsp;&emsp;8.1.2. [RNAseq reads](#rnaseq-reads)   
&emsp;&emsp;8.1.3. [variant file](#variant-file)   
&emsp;&emsp;8.1.4. [blacklisted regions](#blacklisted-regions)   
&emsp;8.2. [workflow](#workflow)   

<!-- /MDTOC -->

# Summary
- manage dependencies efficiently
- have portable pipelines
- reproduce results easily
- create pipeline with a test dataset and then scale
- separates application logic from deployment logic

# Installation
Using Ubuntu VirtualBox VM:
- enable copy/paste in VM settings then go ``Devices > insert Guest Addition CD image``; run software; restart; test copy/paste
- enable a shared directory pointing to ``/media``; add user to group ``sudo adduser $(whoami) vboxsf``; reboot vm for changes to take effect.
- Install java ``sudo apt install default-jre``; test ``java --version``
- Install docker ``sudo snap install docker``; test ``docker --version``
- Install [singularity](https://singularity.lbl.gov/install-linux); use neuro debian mirror to get package; test ``singularity --version``
- Install [conda](https://www.anaconda.com/products/individual#linux);
````
wget https://repo.anaconda.com/archive/Anaconda3-2020.11-Linux-x86_64.sh
md5sum Anaconda3-2020.11-Linux-x86_64.sh
bash Anaconda3-2020.11-Linux-x86_64.sh
conda init  # test
````
- Install graphviz ``sudo apt update -y; sudo apt install -y graphviz``
- Install [nextflow](https://seqera.io/training/#_environment_setup); ``curl get.nextflow.io | bash``; ``mv nextflow ~/bin``; add ~/bin to PATH; test ``nextflow info``
- For atom; install package ``language-nextflow``

# Sources
## Material
``wget -q -O- https://s3-eu-west-1.amazonaws.com/seqeralabs.com/public/nf-training.tar.gz | tar xvz``
The nf-core community builds nextflow pipelines.
## Products: nextflow tower
Helps monitoring pipelines. You can launch a pipeline from git specifying where to run it. You can create pipeline actions that can be used to trigger a pipeline using the API. They allow to launch the execution of the pipeline when a change is submitted to the associated code repository. Testing. Anything you can do is API driven.

# Concepts
Using nextflow, means using docker/singularity (safe/deploy) from docker and conda environments.

## Containerisation
NF uses docker and singularity (better for cluster and security). Singularity containers can be built using docker files.

## Nextflow and cloud deployment
To deploy nextflow on a cluster, you need to define the executor (e.g. slurm, Kubernetes, AWS batch, google lifesciences ...).

NF orchestrates workflow execution via AWS batch: nextflow submits job to EC2 VMs and requires a shared object storage to exchange data between VMs.
Multiple jobs can be assigned per VM. Similarly, this can be done with Google LifeSciences API but one VM is created per job/task.
Also, Kubernetes Cluster can be used, each task becomes a pod and a persistent volume claim is needed to exchange data between tasks.

On AWS, you can use FSX for Lustre for debugging and allow the system to run faster.

### AWS batch and nextflow tower
When running on AWS batch, Nextflow is linked to a job queue which lives in a computing environment that has a provisional model. AWS batch is managing server for us, and you can check the number of instances span off for the queue. EC2 spot requests.

Setting a cloud environment with nextflow tower requires:
 - platform
 - credential
 - region (like EU-west)
 - pipeline work directory (in s3://)
 - config mode (batch forge or manual)
 - provisional model: spot (defined limits of resources; interruptible Ec2 instances) or on-demand Ec2
 - FSx for Lustre: create a central filesystem, no longer have to use s3 ?

# Command line
On the command line, access using:
``-``: a nextflow command
``--``: a params from workflow and scripts

# Code
## Process
````
process task_name {
  directive1
  directive2

  input:
  file 'filename' from 'channel1'
  file 'filename2' from 'channel2'

  output:
  file 'filename3' into 'channel3'

  script:
  """
  some command line thingy
  other command number 2
  """
}
````

## Channels and processes
Channels link processes together. Channels are queues. Data is loaded into channels. Parallelism is implicit.
Queue channels are asynchronous, unidirectional and keep data order.

``Channel.fromFilePairs("*_{1,2}.fq")`` creates a tuple with a sample id and a list of files.

## Parameters
Parameters can be defined in ``nextflow.config`` file and in the nextflow scripts, and overridden using ``--<parameter line>`` on command line.

## Environment and profiles scopes
In the config file, use the ``env`` scope to inject the environment variables inside processes.
Closures can be used inside the config file, or the process itself. You can change process directives in the config file.
Use the special scope ``profiles`` to define profile and switch between them using ``-profile <profile_name>``.

## work directory and resume
The work directory stores ALL inputs and outputs from processes. It is necessary to have a policy regarding the removal of the files depending on how we want to use ``-resume`` option. So make sure to design your pipeline for ``--resume`` to work appropriately.

Things to care about for resume success:
 - how input files are modified
 - Some timestamps can be dependent on the filesystem, so it can change even if not modified
 - inappropriate use of global variables: race condition in global variables (use local variables for parallelism)
 - not deterministic input channels (remember that channels order is kept and results may depend on that order); solution: use a key to join channels for example.

## Error management
Use the ``errorStrategy`` directive, to specify the nextflow behavior when an error is encountered. You can also retry with backoff (sleep using a closure) or even allocate resources dynamically. Default is to kill all processes. You might want to wait for processes to finish.

# DSL2
Dsl2 is an improvement of implementation patterns allowing for pipeline modularisation and components (process/workflow) reuse.

## Processes are modular
- Processes are not written with channels anymore, meaning ``from`` and ``into`` syntax should not be used.
- modules of processes/workflows are defined in single files then included with ``include { PROCESS_NAME } from ./path/file``.

## Workflow (take/main/emit)
Channels connecting processes (and workflows) have to be defined in workflows.
It is preferable to name local variables for clarity.
It is preferable to name output variable for clarity, otherwise variables have to be accessed using ``output[<int>]``.
````
workflow {
  transcriptome = params.transcriptome
  read_pairs_ch = Channel.fromFilePairs( params.reads )

  INDEX(transcriptome)
}
````
Naming workflows allows for control over configuration, similar to processes.
````
workflow NAME {
  take:
    file
    channel

  main:
    CHANNEL_NAME(file, channel)
    CHANNEL_NAME2(CHANNEL_NAME.out, channel)

  emit:
    my_output_name = CHANNEL_NAME2.out
}
````

## Pipes
Also we can connect channels using pipes. It seems to be equivalent to the ``.`` operator.
````
workflow {
  Channel(params) | some_process_operator | some_process | view_operator { some_closure }
}
````

## Transition from dsl to dsl2
### Some guidelines
 - declare use of syntax of dsl2
 - Decouple processes by removing channel references in processes (``from`` and ``into``)
 - recombine processes together into ``workflow {}``
 - group processes into modules into file(s) then ``include`` each process by name
 - build a typical workflow structure; name workflows for easy reuse.
 - don't forget to name each emitted output.

### A workflow example
 ````
 workflow workflow_name {
   take:
   main:
    task1()
   emit:
   assign_name = task1.out
 }
 ````
Note:
 - task1 is named ``workflow_name:task1``, workflow name behave as scope
 - output is ``workflow_name.out``
 - use name assignment in ``emit:`` (``assign_name = task1.out``) if possible
 - workflows can be called by name using command line: ``nextflow run <filename.nf> -entry <workflow_name>``
 - workflow without names is ignored if imported into another one.

# Bioinformatics pipeline example
## data description
All data is available at this [git repo](https://github.com/nextflow-io/nf-hack18/tree/master/hands-on). Retrieve data for the [tutorial](https://seqera.io/training/handson/#_data_description).
### genome assembly `genome.fa`
Genome assembly hg19 from genbank - chr22 only [here](https://www.ncbi.nlm.nih.gov/nuccore/CM000684.1)

### RNAseq reads
Expect 6 bam files ranging 17-35gb per file; leading to 12 fastq files

 - Download the following files and check md5 sums.
| sample | isogenic replicate 1; alignment hg19| isogenic replicate 2; alignment hg19 |
|---|---|---|
| ENCSR000COQ (whole cell) |ENCFF664VNB|ENCFF790WER|
| ENCSR000CPO (nuclear) |ENCFF009TYW|ENCFF149PXU|
| ENCSR000COR (cytosol) |ENCFF705WOW|ENCFF367BKI|

 - Then we need to retrieve reads mapping chr 22 only into fastq; we use samtools. Get samtools with this docker [image](https://github.com/chrishah/samtools-docker/blob/master/Dockerfile).

 ``docker run -it -v /c/Users/User/Desktop/data:/test 288856f3b2ca /bin/bash``

 - create index
 ``samtools index ENCFF664VNB.bam``
 - extract reads for locus 22q11
 ``samtools view -h ENCFF664VNB.bam chr22:16000000-18000000 > ENCFF664VNB_chr22.bam``
 - output paired reads to separate files, discarding singletons into fastq
 ``samtools fastq -1 ENCSR000COQ1_1.fastq -2 ENCSR000COQ1_2.fastq -0 /dev/null -n ENCFF664VNB_chr22.bam``

### variant file
 - variants are coming from Illumina. Get the data and the index file ``tbi``.
 ``wget ftp://ussd-ftp.illumina.com/2017-1.0/hg19/small_variants/NA12878/NA12878.vcf.gz``
 - Then subset the variants to keep only those matching chr 22 using [vcftools](https://vcftools.github.io/man_latest.html#EXAMPLES).
 ``vcftools --gzvcf  NA12878.vcf.gz --chr chr22 --recode --out NA12878_chr22``
 Note: output file is NA12878_chr22.recode.vcf

### blacklisted regions
This is a bed file.
https://www.encodeproject.org/files/ENCFF001TDO/

## workflow
### run on win10 using gitbash
 ``winpty docker run -it --rm -v "C:\\Users\\User\\Desktop\\data\\nextflow:/test" nextflow/nextflow bash -c "nextflow run /test/hands-on/main.nf"``
