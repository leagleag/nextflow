<!-- MDTOC maxdepth:6 firsth1:1 numbering:0 flatten:0 bullets:1 updateOnSave:1 -->

- [Summary](#summary)   
- [DSL2](#dsl2)   
   - [Transition from dsl to dsl2.](#transition-from-dsl-to-dsl2)   

<!-- /MDTOC -->

# Summary


Problems
- manage dependencies efficiently
- have portable pipelines
- reproduce results can be hard
- create pipeline with a test dataset and then scale

Separates application logic from deployment logic


Data Flow programming??

nextflow script
nextflow runtime: deploy on cloud

````
process task_name {
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
Channels link processes together. Channels are queues. Data is loaded into channels.
Parallelism is implicit.

``Channel.fromFilePairs("*_{1,2}.fq")`` creates a tuple with a sample id and a list
of files.

Deployement
NF orchestrates workflow excution via AWS batch; nf submit job to EC2 VMs; this requires a shared object storage to exchange data between VMs. Multiple jobs can be assigned per VM.
Similarly, this can be done with Google LifeSciences API. One VM is created per job/task.
Also, Kubernetes Cluster can be used. Each task becomes a pod. Use a persistant volume claim to exchange data between tasks.

On AWS, you can use FSX for Lustre for debugging and allow the system to run faster.

Containerisation
NF use docker and singularity (better for cluster and security).

Nextflow has an execution report. And has syntax highlighting available for editors.

nf-core: community building pipelines ready to use using nf.

nextflow tower: help monitoring pipelines. You can launch a pipeline from git specifying where to run it. You can create pipeline actions that can be used to trigger a pipeline using the API. They allow to launch the execution of the pipeline when a change is submitted to the associated code repository. Testing. Anything you can do is API driven.

Run Ubuntu VirtualBox VM:
- enable copy/paste in VM settings then go ``Devices > insert Guest Addition CD image``; run software; restart; test copy/paste
- enable a shared directory pointing to ``/media``; add user to group ``sudo adduser $(whoami) vboxsf``; reboot for changes to take effect.
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
- For atom; intall package ``language-nextflow``

Material
``wget -q -O- https://s3-eu-west-1.amazonaws.com/seqeralabs.com/public/nf-training.tar.gz | tar xvz``

Access from command line:
``-``: nextflow command
``--``: params from workflow and scripts

Channels
 - queue channels: asynchronous; unidirectional; with data order kept.

Use path instead of file.
Do examples about combining inputs.

use docker images
use singularity from docker images
use conda environments
use nextflow config ``nextflow.config``. Parameters can be defined in ``nextflow.config`` and in the nextflow script, and using ``--<parameter line>`` on command line.
Use the ``env`` scope to inject the environment variables inside processes, closures can be used inside the config file, or the process itself. You can change process directives in the config file.
Use the special scope ``profiles`` to define profile and switch between them using ``-profile <profile_name>``.

note: tendency is to move to one container per process using dsl2.

Deployement
To deploy nextflow on a cluster, you need to define the executor (e.g. slurm, Kubernetes, AWS batch, google lifesciences ...).

AWS Batch
Nextflow is linked to a job queue which lives in a computing environment that has a provisional model (resources set on demand or on spot - fixed?). AWS batch is managing server for us, and you can check the number of instances span off for the queue. EC2 spot requests.

How to setup with nextflow tower.
1. Create compute environment
 - platform
 - credential
 - region (like EU-west)
 - pipeline work directory (in s3://)
 - config mode (batch forge or manual)
 - provisional model: spot (defined limits of resources; interruptible Ec2 instances) or on-demand Ec2
 - FSx for Lustre: create a central filesystem, no longer have to use s3 ?

The work directory stores all inputs and outputs from processes. It is necessary to have a policy regarding the removal of the files depending on how we want to use ``-resume`` option.

- Make sure to design your pipeline for ``--resume`` to work appropriately.
 - how input files are modified
 - timestamp can be dependent on the filesystem, so it can change even if not modified
 - race condition in global variables (use local variable for parallelism)
 - not deterministic input channels (remember that channels order is kept and results may depend on that order); solution: use a key to join channels for example.

Error management
Use the ``errorStrategy`` directive, to specify the nextflow behavior when an error is encountered. You can also retry with backoff (sleep using a closure) or even allocate resources dynamically.

# DSL2
 - pipeline modularisation
 - component reuse
 - fluent definition of recurrent implementation patterns

Process are linked with channels so it is hard to reuse.
Define ``nextflow.enable.dsl=2``
 - processes do not use channels
 - processes can be defines in single files then included in a workflow file ``include { PROCESS_NAME } from ./path/file``.

Workflow example:
```
workflow {
  transcriptome = params.transcriptome
  read_pairs_ch = Channel.fromFilePairs( params.reads )

  INDEX(trascriptome)
}

workflow NAME {
  take:
    file
    channel

  main:
    CHANNEL_NAME(file, channel)
    CHANNEL_NAME2(CHANNEL_NAME.out, channel)

  emit:
    CHANNEL_NAME2.out
}
```
DSL2 pipes example:
```
workflow {
  Channel(params) | process_operator | process | view_operator { some_closure }
}
```
## Transition from dsl to dsl2.
 - declare use of syntax of dsl2
 - remove channel references in processes (``from`` and ``into``)
 - recombine processes together into a ``workflow {}``
 - group processes into modules into file(s) then ``include`` each process by name
 - typical workflow structure; use a name to be able to reuse it.
 ```
 workflow workflow_name {
   take:
   main:
    task1()
   emit:
   assign_name = task1.out
 }
 ```
Note:
 - task1 is named ``workflow_name:task1``, workflow name behave as scope
 - output is ``workflow_name.out``
 - use name assignment in ``emit:`` (``assign_name = task1.out``) if possible
 - workflows can be called by name using command line: ``nextflow run <filename.nf> -entry <workflow_name>``
 - workflow without names is ignored if imported into another one.
