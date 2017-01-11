MANUAL for TOGGLEv0.3 using ONTHEFLY creation of PIPELINES
===========
# SUMMARY

* [PRE REQUISITES](#prerequisites)
* [LAUNCHING AN ANALYSIS](#launching)
* [RESULTING FOLDER](#finalFolder)
* [SENDING OPTIONS](#sendingOptions)
* [CREATING A PIPELINE](#creatingPipeline)
 *  [PROVIDING AN ORDER](#order)
 * [SAME SOFTWARE REPEATED MULTIPLE TIMES](#samesoftmultiple)
 * [GIVING A COMMON STEP TO ALL INDIVIDUALS](#commonstep)
 * [STARTING ONLY A COMMON TREATMENT](#commononly)
* [CLEANING STEPS](#cleaning)
* [COMPRESSING STEPS](#compressing)
* [SCHEDULER AND PARALLEL RUNS](#scheduling)
* [OUTPUT AND ERROR LOGS](#logs)
* [DRAWING THE PIPELINE](#drawing)
* [PROCESS_RADSTACKS STEPS](#processradstacks)


The TOGGLEv0.3 *onTheFly* version allows users to create their own customized pipelines.
You can modify not only the different options for each software but also define specific organization for your analysis.

You can therefore remove some steps compared to the previous version, starting from SAM/BAM or VCF files instead of FASTQ only, asking for individual treatments only, individual then common (such as mapping followed by common calling), or even only common treatment.

# <a name="prerequisites"></a>PRE REQUISITES

THE FILE NAMES MUST BE UNDER THE FORM **individual1_1.fastq**, **mapping1.sam**, **myVcf.vcf**

THE ONLY DOT ('.') ALLOWED IN NAMES IS FOR THE EXTENSION. The text before the dot will be used as individual names.

For *Fastq* files, only one underscore ('_') is authorized, before the direction of sequences ( _1, _2 or _R1, _R2), just before the extension (*.fastq*).

Please use only UTF-8 standard symbols, no weird characters or space, pipe, tilde or any type of commas.

# <a name="launching"></a>Launching an analysis
The current version is based on the **toggleGenerator.pl** script

````
$toggleGenerator.pl -d|--directory DIR -c|--config FILE -o|--outputdir DIR [-r|--reference FILE] [-k|--keyfile FILE] [-g|--gff FILE] [-nocheck|--nocheckFastq] [--help|-h]
````

| Required named arguments:       |                                                                                                                                |
| :------------------------------ | :----------------------------------------------------------------------------------------------------------------------------- |
| -d / --directory DIR:           | a folder with raw data to be treated (FASTQ, FASTQ.GZ, SAM, BAM, VCF)                                                          |
| -c / --config FILE:             | generally it is the *software.config.txt* file but it can be any text file structured as shown below.                          |
| -o / --outputdir DIR:           | the current version of TOGGLE will not modify the initial data folder but will create an output directory with all analyses in.|

| Optional named arguments:       |                                                                                                                                |
| :------------------------------ | :----------------------------------------------------------------------------------------------------------------------------- |
| -r / --reference FILE:          | the reference FASTA file to be used. (1)                                                                                           |
| -g / -gff FILE:                 | the GFF file to be used for some tools.                                                                                        |
| -k / --keyfile FILE:            | the keyfile use for demultiplexing step.                                                                                       |
| -nocheck / --nocheckFastq:      | by default toggle checks if fastq format is correct in every file. This option allows to skip this step.                       |
| -h / --help:                    | show help message and exit                                                                                                     |

(1): If no index exists it will be created accordingly to the pipeline requested index. If the index exist, they will not be re-created UNLESS the pipeline order (see below) expressively requests it (updating the index e.g.)

All the the paths (files and folders) can be provided as absolute (/home/mylogin/data/myRef.fasta) or relative (../data/myRef.fasta).

The *SNPdiscoveryPaired.config.txt* file is an example of how to customize your pipeline.

Any software configuration will start as follows:
 ````
 $First Software
 option1
 option2

 $Second Software
 option1
 option2
 ````

 # <a name="finalFolder"></a>What will you have in results ?

TOGGLE will generate an output folder containing different files and subfolders, as follows:

 ![TOGGLE Final Folder](/docs/images/toggleOutputFolder.png)

 The final results are contained in the **finalResults** folder.
 TOGGLE will also copy the *software config* file corresponding to the analysis, in order users can recover their options.
 The **output** folder contains all sub analyses, i.e. the individual analyses or intermediate data.

# <a name="sendingOptions"></a>Sending options
As for the previous version, you can address any option to any given software (as soon as the given option exists for the given software ^^) as follows:
````
$bwaAln
-e=1
-n=5
````

You can also write as follows
````
$BWA ALN
-e=1
-n=5
````

The software name is not case sensitive and the subprogram can be "glued" to the name (bwaALN is recognized, as well as bwa aln).

If your option has an equal symbol within, such as **-l hostname=MyNode**, you have to write the option as follows:

````
-l hostname=MyNode
````

**BE CAREFUL: for any file send as option (e.g. "-knownsite" option in GATK BaseRecalibrator), the path for this file must provided as an absolute one, and not relative. It must be thus noted as /my/absolute/path/FILE and not ../FILE, for instance.**

# <a name="creatingPipeline"></a>Creating a pipeline

The **toggleGenerator.pl** script will use the configuration file informations (generally named as *software.config.txt* file, but not mandatory) to create specific pipeline scripts in the output folder called **toggleBzz.pl** for individual treatments (individual mapping e.g.) and **toggleMultiple.pl** for global treatment (global SNP calling e.g.).

The order will be verified before creating the scripts, by checking if the output files from the step n-1 will be accepted as input for the step n.

![TOGGLE pipeline](/docs/images/TogglePipeline.png)

### <a name="order"></a>Providing an order
The order of a specific pipeline is provided in the *software.config.txt* as the software *order*

Thus, if you provide option such as:
````
$order
1=bwa aln
2=BWASAMPE
3=samtools view
````
You will obtain a pipeline taking fastq files (plain text of gzipped), aligning them upon a given reference using bwa aln/sampe, and finally cleaning the resulting sam and providing a BAM file (if samtools view option are as such).

Note that the way you write the name of the step is not really important, e.g. it is not case-sensitive, and can have space between the different words. A dedicated module will adjust everything.

You can start from any type of format recognized by TOGGLE at any block from TOGGLE. Thus, if you have a raw VCF file you can start from this one at step one and clean it (**even if your VCF contains informations for multiple individuals**):

````
$order
1=gatkVariantFiltrator
2=gatkSelectVariants
````
**BE CAREFUL**: You can comment any step using a '#' symbol before it (ex *#3=bwaAln*) to avoid to retype all numbers. However, such a comment will provoke anomalies in cleaning and compressing steps. The best approach is to provide a list of following numbers (e.g. 1, 2, 3, 4 and not 1, 3, 4). This is also true for steps higher than 1000.


### <a name="samesoftmultiple"></a>Same software repeated multiple times

If you want to adress multiple times the same step BUT with different configurations (e.g. a first samtools view to obtain a BAM then a second samtools view to clean this BAM using the -f option), you have to indicate as follows
````
$order
1=bwa aln
2=BWASAMPE
3=samtools view 1
4= samtools view 2
````

In the same time you have to provide the same informations in your configuration:
````
$samtoolsView 1
-Sb

$samtools View 2
-f 0x02
````
**BE CAREFUL**: in such multiple times repeated software, you have to put a space between the software name and the step (ex samtoolsview 1 and samtoolsView 2) in the order as well as in the options !!


### <a name="commonstep"></a>Giving a common step to all individuals (multiple entry files)

If you want for instance to map all individuals (multiple files) separately then perform a common SNP calling, please add a step number higher than 999.

````
$order
1= bwaAln
2=bwaSampe
3=picardTools SortSam
1000=gatkUnifiedGenotyper
````
This pipeline will map FASTQ then sort per coordinates the SAM for each individuals, and then will launch a common step for all as a global SNP calling. You can add after this calling other steps (1001=gatkVariantFiltrator for example).

**Warnings**: if you want to map and call SNP from only one file, you do not need to call a 1000+ step (here *gatkUnifiedGenotyper* would be then the 4th step).

### <a name="commononly"></a>Starting only as a common treatment
If you want only a global treatment **for multiple files** (you may have all your formatted BAM and you want to perform the calling and subsequent steps), you can create the following pipeline:
````
$order
1000=gatk UnifiedGenotyper
1001=gatk VariantFiltrator
1002=gatkSelectVariants
````

The pipeline will treat the data as a global pipeline in this case, without separating the individual files.

# <a name="cleaning"></a>Cleaning steps
In order to gain hard drive space, you may want to remove some steps from your complete analysis.

In this case, you can specify in the configuration file which step must be REMOVED along the pipeline (as soon as the step following them has been correctly completed), using the key *cleaner*.

As an example
````
$order
1=bwaAln
2=bwaSampe
3=samtools sort

$cleaner
1
2
````

There only the last step result (samtools sort) will be conserved. The folders and logs of the cleaned steps are conserved, but not the resulting files.

**NOTE**: the last step will never be cleaned. This is true for single analyses (steps lower than 999 - *toggleBzz.pl* script) as well as for common steps (higher than 1000 - *toggleMultiple.pl* script)

# <a name="compressing"></a>Compressing steps
In order to gain hard drive space but conserving data, you may want to compress some steps from your complete analysis.

In this case, you can specify in the configuration file which step must be COMPRESSED in tar.gz along the pipeline (as soon as the step following them has been correctly completed), using the key *compress*.

As an example
````
$order
1=bwaAln
2=bwaSampe
3=samtools sort

$compress
1
2
````

There only the last step result (samtools sort) will be conserved, the other being compressed in their respective tar.gz archive. The folders and logs of the compressed steps are conserved, but not the resulting files.

**NOTE**: as for the cleaner system, the last step will never be compressed. This is true for single analyses (steps lower than 999 - *toggleBzz.pl* script) as well as for common steps (higher than 1000 - *toggleMultiple.pl* script).

BE CAREFUL: CLEANING IS OF HIGHER ORDER THAN COMPRESS. If the same step is required to be cleaned AND compressed, it will be only cleaned!

# <a name="scheduling"></a>Scheduler and parallel runs
The current version is scheduler-aware (**SGE**, **MPRUN**, **Slurm** and **LSF** on v0.3+), and is able to decide by itself to launch on such a system.
The jobs will be launched in parallel however only if the *software.config* file contains informations for scheduling, *i.e.* the queue name, number of core/threads per jobs, etc...

As an example for **SGE**:
````
$sge
-q myqueue.q
-pe ompi 2
-b Y
-cwd
-V

````

If the *software.config* file contains those **SGE/Slurm/MPRUN/LSF** informations AND the machine is **SGE/Slurm/MPRUN/LSF** capable, the *toggleBzz.pl* and the *toggleMultiple.pl* scripts will be launched as parallel jobs.

Moreover, in parallel as in linear mode, an error in one individual will not stop the whole process. This individual will be marked as error in the error log, but the others will continue.

**NOTE**: If you need to provide specific ENVIRONMENT variables in your jobs (e.g. *export* or else), you can provide the $env key in your software.config file

````
$env
export PATH=$PATH:/my/other/PATH
export PERL5LIB=$PERL5LIB/my/new/path
````

# <a name="logs"></a>Output and Error Logs

TOGGLE will generate two main types of logs, the *.o* for normal output and the *.e* for the errors and warnings (these last ones are normally empty files).
Each level of TOGGLE will generate this pair of log:
* **GLOBAL_ANALYSIS_date.o/.e** logs represent the general output for the complete analysis (*toggleGenerator.pl* logs). They are located at the root of the output directory.
* **IndividualName_global_log.o/.e** logs represent the local output for sub analysis (*toggleBzz.pl* and *toggleMultiple.pl* logs). They are located in their respective subdirectories in the output folder.

# <a name="drawing"></a>Drawing the pipeline
If *Graphviz* is installed on the running server, the **toggleGenerator.pl** script will also generate a graphical view of the pipeline in the output folder.
If *Graphviz* is not installed, a .dot file is nevertheless created, allowing the user to create a graphical view on another machine having *Graphviz*

# <a name="processradstacks"></a>Process_Radstacks step
If process_radstacks from stacks is used in your pipeline, you have to put it as the fist step else TOGGLE return you an error. 

