TOGGLE : Toolbox for generic NGS analyses
===========

TOGGLE (TOolbox for Generic nGs anaLysEs) is a suite of 10 packages and more than 110 modules able to manage a large set of NGS softwares
and utilities to easily design pipelines able to handle hundreds of samples. Moreover, TOGGLE offers an easy way to manipulate the various
options of the different softwares through the pipelines in using a single basic configuration file, that can be changed for each assay without
having to change the code itself.

We present also the implementation of TOGGLE in a complete analysis pipeline designed for SNP discovery for large sets of NGS data, ready to use
in different environments (single machine to HPC clusters).

##  Contributing sebastien

* Licencied under CeCill-C (http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html) and GPLv3 
* Intellectual property belongs to IRD, CIRAD and SouthGreen developpement plateform 
* Written by Cecile Monat, Christine Tranchant, Ayite Kougbeadjo, Cedric Farcy, Mawusse Agbessi, Marilyne Summo, and Francois Sabot

## REQUIREMENTS

#### Perl - Julie^^

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

#### INSTALL 

* Create the directory TOGGLE where you want to install TOGGLE and go into this directory

* Get the TOGGLE code Clone the git

````
$git clone https://github.com/SouthGreenPlatform/TOGGLE.git
````
* Add the Module path to the PERL5LIB environment variable

````
export PERL5LIB=$PERL5LIB:/pathToToggle/Modules
````
* In the same way add the TOGGLE directory to the PATH environment variable

````
export PATH=$PATH:/pathToToggle
````

Rq : you can add this to the ~/.bashrc to make it always available when you log-in.

* Run the test script

````
$cd TEST
$sh all_tests.sh
````

* Test the pipiline with the test data

````
$globalAnalysis.pl $PATH_INSTALL/DATA/arcardTest/ $PATH_INSTALL/software.config.txt $PATH_INSTALL/DATA/referenceArcad.fasta
````

* Check the good running 
> > * No error message
> > * BamDirectory has been well created into $PATH_INSTALL/DATA
> > * the data generated are good

````
tail $PATH_INSTALL/DATA/BamDirectory/GATKVARIANTFILTRATION.vcf

##contig=<ID=2299889,length=39753>
##contig=<ID=2299897,length=19555>
##reference=file:///home/tranchant/TOGGLE-19-02/DATA/referenceIrigin.fasta
#CHROM    POS    ID    REF    ALT    QUAL    FILTER    INFO    FORMAT    irigin1    irigin1Single    irigin2_1    irigin3    irigin3Single
2224477    996    .    TA    T    35.61    PASS    AC=2;AF=1.00;AN=2;DP=2;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;MQ0=0;QD=17.81;SOR=0.693    GT:AD:DP:GQ:PL    ./.    ./.    ./.    1/1:0,2:2:6:69,6,0    ./.
2248321    377    .    C    G    65.65    PASS    AC=2;AF=1.00;AN=2;DP=2;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=29.00;MQ0=0;QD=32.83;SOR=0.693    GT:AD:DP:GQ:PL    ./.    ./.    ./.    1/1:0,2:2:6:90,6,0    ./.
2248321    379    .    C    T    65.65    PASS    AC=2;AF=1.00;AN=2;DP=2;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=29.00;MQ0=0;QD=32.83;SOR=0.693    GT:AD:DP:GQ:PL    ./.    ./.    ./.    1/1:0,2:2:6:90,6,0    ./.
2281178    4213    .    G    A    65.65    PASS    AC=2;AF=1.00;AN=2;DP=2;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;MQ0=0;QD=32.83;SOR=0.693    GT:AD:DP:GQ:PL    ./.    ./.    ./.    1/1:0,2:2:6:90,6,0    ./.
2281178    4214    .    A    G    65.65    PASS    AC=2;AF=1.00;AN=2;DP=2;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;MQ0=0;QD=32.83;SOR=0.693    GT:AD:DP:GQ:PL    ./.    ./.    ./.    1/1:0,2:2:6:90,6,0    ./.
2290182    1013    .    A    G    45.65    PASS    AC=2;AF=1.00;AN=2;DP=2;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=29.00;MQ0=0;QD=22.83;SOR=0.693    GT:AD:DP:GQ:PL    ./.    ./.    ./.    1/1:0,2:2:6:70,6,0    ./.
````

##  Versions Notes - branche dev

Release 0.2, 14st of March, 2015

Second version release
