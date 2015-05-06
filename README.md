TOGGLE : Toolbox for generic NGS analyses
=========================================

TOGGLE (TOolbox for Generic nGs anaLysEs) is a suite of 10 packages and more than 110 modules able to manage a large set of NGS softwares
and utilities to easily design pipelines able to handle hundreds of samples. Moreover, TOGGLE offers an easy way to manipulate the various
options of the different softwares through the pipelines in using a single basic configuration file, that can be changed for each assay without
having to change the code itself.

We present also the implementation of TOGGLE in a complete analysis pipeline designed for SNP discovery for large sets of NGS data, ready to use
in different environments (single machine to HPC clusters).

##  Contributing

* Licencied under CeCill-C (http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html) and GPLv3
* Intellectual property belongs to IRD, CIRAD, ADNid and SouthGreen development platform
* Written by Cecile Monat, Christine Tranchant, Ayite Kougbeadjo, Cedric Farcy, Mawusse Agbessi, Enrique Ortega-Abboud, Sébastien Ravel, Julie Orjuela-Bouniol, Souhila Amanzougarene, Gauthier Sarah, Marilyne Summo, and Francois Sabot

## REQUIREMENTS

#### Perl

* Data::Translate
* Data::Dumper
* Test::More
* Test::Deep
* Capture::Tiny


#### Bioinformatics software

* java 1.7
* fastQC v0.10.1
* cutadapt 1.2.1
* BWA 0.7.2
* gatk 3.3
* picardTools 1.63
* SAMtools 0.1.18

#### TOGGLE DIRECTORY

````
TOGGLE
|
|_ DATA-TEST
|
|_ DATA
    |———— arcadTest
            |———— arcad1_1.fastq
            |———— arcad1_2.fastq
            |———— arcad2_1.fastq
            |———— arcad3_1.fastq
            |———— arcad3_2.fastq
    |———— expectedData
    |———— iriginTest
            |———— irigin1_1.fastq
            |———— irigin1_2.fastq
            |———— irigin2_1.fastq
            |———— irigin3_1.fastq
            |———— irigin3_2.fastq  
    |____ referenceArcad.fasta
    |____ referenceIrigin.fasta
|
|_ Modules
|
            |———— bwa.pm
            |———— cufflinks.pm
            |———— cutadapt.pm
            |———— fastqUtils.pm
            |———— fastqc.pm
            |———— gatk.pm
            |———— localConfig.pm
            |———— pairing.pm
            |———— picardTools.pm
            |———— samTools.pm
            |____ toolbox.pm
            |____ tophat.pm
|
|_ TEST
            |———— adaptator.txt
            |————  adaptator.txt
            |———— all_tests.sh
            |———— bwa_test.t
            |———— cufflinks_test.t
            |———— cutadapt_test.t
            |———— fastqUtils_test.t
            |———— fastqc_test.t
            |———— gatk_test.t
            |———— localConfig_test.t
            |———— pairing_test.t
            |———— picardTools_test.t
            |———— pipelineTryExample.p
            |———— samTools_test.t
            |———— software.config.txt
            |———— toolbox_test.t
            |———— tophat_test.t
|_ LICENSE
|_ README.md
|_ adaptator.txt
|_ globalAnalysis.pl
|_ mergeAnaysis.pl
|_ pairAnalysis.pl
|_ singleAnalysis.pl
|_ software.config.txt
|_ software.config.txt.test
````
