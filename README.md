TOGGLE : Toolbox for generic NGS analyses
===========

TOGGLE (TOolbox for Generic nGs anaLysEs) is a suite of 10 packages and more than 110 modules able to manage a large set of NGS softwares
and utilities to easily design pipelines able to handle hundreds of samples. Moreover, TOGGLE offers an easy way to manipulate the various
options of the different softwares through the pipelines in using a single basic configuration file, that can be changed for each assay without
having to change the code itself.

We present also the implementation of TOGGLE in a complete analysis pipeline designed for SNP discovery for large sets of NGS data, ready to use
in different environments (single machine to HPC clusters).

-------

## REQUIREMENTS

#### Perl

* Data::Translate
* Data::Dumper
* Test::More
* Test::Deep
* Capture::Tiny

#### Bioinformatics software

java 1.7
fastQC v0.10.1
cutadapt 1.2.1 
BWA 0.7.2 
gatk 3.3 
picardTools 1.63
SAMtools 0.1.18


## Versions Notes


Release 0.2, 14st of March, 2015

Second version release

## Contributing

Licencied under CeCill-C (http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html) and GPLv3 
Intellectual property belongs to IRD, CIRAD and SouthGreen developpement plateform 
Written by Cecile Monat, Christine Tranchant, Ayite Kougbeadjo, Cedric Farcy, Mawusse Agbessi, Marilyne Summo, and Francois Sabot
