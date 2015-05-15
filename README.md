<<<<<<< HEAD
TOGGLE : Toolbox for generic NGS analyses
=========================================
=======
TOGGLE : Toolbox for generic NGS analyses => version test: development for RnaSEQ 
===========
>>>>>>> dev

TOGGLE (TOolbox for Generic nGs anaLysEs) is a suite of 10 packages and more than 110 modules able to manage a large set of NGS softwares
and utilities to easily design pipelines able to handle hundreds of samples. Moreover, TOGGLE offers an easy way to manipulate the various
options of the different softwares through the pipelines in using a single basic configuration file, that can be changed for each assay without
having to change the code itself.

We present also the implementation of TOGGLE in a complete analysis pipeline designed for SNP discovery for large sets of NGS data, ready to use
in different environments (single machine to HPC clusters).

##  Contributing

* Licencied under CeCill-C (http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html) and GPLv3
<<<<<<< HEAD
* Intellectual property belongs to IRD, CIRAD, ADNid and SouthGreen development platform
* Written by Cecile Monat, Christine Tranchant, Ayite Kougbeadjo, Cedric Farcy, Mawusse Agbessi, Enrique Ortega-Abboud, Sébastien Ravel, Julie Orjuela-Bouniol, Souhila Amanzougarene, Gauthier Sarah, Marilyne Summo, and Francois Sabot
* Copyright 2014-2015

##INSTALLATION

see https://github.com/SouthGreenPlatform/TOGGLE/blob/master/INSTALL.md

=======
* Intellectual property belongs to IRD, CIRAD and SouthGreen developpement plateform
* Written by Cecile Monat, Christine Tranchant, Ayite Kougbeadjo, Cedric Farcy, Mawusse Agbessi, Marilyne Summo, and Francois Sabot
>>>>>>> dev

## REQUIREMENTS

#### Perl

* Data::Translate
* Data::Dumper
* Test::More
* Test::Deep
* Capture::Tiny


#### Bioinformatics software

<<<<<<< HEAD
* Java 1.7
* FastQC v0.11.1 or higher
* CutAdapt 1.2.1 or higher
* bwa 0.7.2 or higher
* GATK 3.3x or higher
* PicardTools 1.124 or higher
* SAMtools 0.1.18 or higher
=======
* java 1.7
* fastQC v0.10.1
* cutadapt 1.2.1
* BWA 0.7.2
* gatk 3.3
* picardTools 1.63
* SAMtools 0.1.18
>>>>>>> dev

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

