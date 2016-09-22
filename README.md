TOGGLE : Toolbox for generic NGS analyses
===========

![TOGGLE Logo](toggleLogo.png)

Dear Biologist, have you ever dream of using the whole power of those numerous NGS tools that your bioinformatician colleagues use through this awful list of command line ?

Dear Bioinformatician, have you ever guess how to design really fastly a new NGS pipeline without having to retype again dozens of code lines to readapt your scripts or starting from scratch ?

**So, be Happy ! TOGGLE is for you !!**

TOGGLE (TOolbox for Generic nGs anaLysEs) is a suite of 19 packages and more than 110 modules able to manage a large set of NGS softwares
and utilities to easily design pipelines able to handle hundreds of samples. Moreover, TOGGLE offers an easy way to manipulate the various
options of the different softwares through the pipelines in using a single basic configuration file, that can be changed for each assay without
having to change the code itself.

Users can also create their own pipeline through an easy and user-friendly approach. The pipelines can start from Fastq (plain or gzipped), SAM, BAM or VCF (plain or gzipped) files, with parallel and global analyses. Samples pipelines are provided for SNP discovery and RNAseq counts.

The system is able to detect parallel/scheduling launching and to manage large amount of samples on large cluster machines.


##  Contributing

* Licencied under CeCill-C (http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html) and GPLv3
* Intellectual property belongs to IRD, CIRAD, ADNid and SouthGreen development platform
* Written by Cecile Monat, Christine Tranchant, Ayite Kougbeadjo, Cedric Farcy, Mawusse Agbessi, Enrique Ortega-Abboud, Sébastien Ravel, Julie Orjuela-Bouniol, Souhila Amanzougarene, Gauthier Sarah, Marilyne Summo, and Francois Sabot
* Copyright 2014-2016

## Contact

For bug tracking purpose you can use the GitHub or questions about TOGGLE, you can contact the mainteners using the following email addresses:

* christine.tranchant_at_ird.fr
* francois.sabot_at_ird.fr

##  Citation
**TOGGLE: Toolbox for generic NGS analyses**. Cécile Monat, Christine Tranchant-Dubreuil, Ayité Kougbeadjo, Cédric Farcy, Enrique
Ortega-Abboud, Souhila Amanzougarene, Sébastien Ravel, Mawussé Agbessi, Julie Orjuela-Bouniol, Maryline Summo and François Sabot.

[*BMC Bioinformatics* 2015, 16:374  doi:10.1186/s12859-015-0795-6][paperLink]

##  INSTALLATION

[Follow the INSTALLATION instructions][installLink]

## MANUAL

[You can find a detailed MANUAL here][manualLink]

## KNOWN ISSUES

[You can find detailed known issues][knownIssues]

## Release Notes

[Current Release Notes][releaseLink]

## REQUIREMENTS

#### Perl


* [Data::Translate](http://search.cpan.org/~davieira/Data_Translate-0.3/Translate.pm)
* [Data::Dumper](http://search.cpan.org/~smueller/Data-Dumper-2.154/Dumper.pm)
* [Test::More](http://search.cpan.org/~exodist/Test-Simple-1.001014/lib/Test/More.pm)
* [Test::Deep](http://search.cpan.org/~rjbs/Test-Deep-0.119/lib/Test/Deep.pm)
* [Capture::Tiny](http://search.cpan.org/~dagolden/Capture-Tiny-0.30/lib/Capture/Tiny.pm)
* [List::Compare](http://search.cpan.org/~jkeenan/List-Compare-0.53/lib/List/Compare.pm)
* [Switch](https://metacpan.org/pod/Switch)


#### Bioinformatics software (minimal version)

* [Perl 5.16](https://www.perl.org/)
* [java 1.7](https://www.java.com/fr/)
* [BWA 0.7.2](http://bio-bwa.sourceforge.net/)
* [SAMtools 0.1.18](http://samtools.sourceforge.net/)
* [picardTools 1.63](http://broadinstitute.github.io/picard/)
* [gatk 3.3](https://www.broadinstitute.org/gatk/)
* [fastQC v0.10.1](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [cutadapt 1.2.1](https://pypi.python.org/pypi/cutadapt)
* [FastxToolkit](http://hannonlab.cshl.edu/fastx_toolkit/)
* [Tophat](https://ccb.jhu.edu/software/tophat/index.shtml)
* [Snpeff](http://snpeff.sourceforge.net/)

#### Bioinformatics tools included

##### BWA (http://bio-bwa.sourceforge.net/)

- bwaAln
- bwaSampe
- bwaSamse
- bwaIndex
- bwaMem

##### SamTools (http://samtools.sourceforge.net/)

- samToolsFaidx
- samToolsIndex
- samToolsView
- samToolsSort
- mergeHeader
- samToolsMerge
- samToolsIdxstats
- samToolsDepth
- samToolsFlagstat
- samToolsMpileUp

##### PicardTools (http://broadinstitute.github.io/picard/)

- picardToolsMarkDuplicates
- picardToolsCreateSequenceDictionary
- picardToolsSortSam
- picardToolsAddOrReplaceReadGroup
- picardToolsValidateSamFile
- picardToolsCleanSam
- picardToolsSamFormatConverter


##### GATK (https://www.broadinstitute.org/gatk/)

- gatkBaseRecalibrator
- gatkRealignerTargetCreator
- gatkIndelRealigner
- gatkHaplotypeCaller
- gatkSelectVariants
- gatkVariantFiltration
- gatkReadBackedPhasing
- gatkUnifiedGenotyper
- gatkBaseRecalibrator
- gatkPrintReads

##### Fastqc (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

- fastqc

##### FastxToolkit (http://hannonlab.cshl.edu/fastx_toolkit/)

- fastxTrimmer

##### Tophat (https://ccb.jhu.edu/software/tophat/index.shtml)

- bowtiebuild
- bowtie2build
- tophat2

##### Snpeff (http://snpeff.sourceforge.net/)

- snpeffAnnotation

##### Cutadapt (https://pypi.python.org/pypi/cutadapt)

- cutadapt

#### OPTIONAL
- Graphviz v2.xx (http://www.graphviz.org/)



[paperLink]:http://www.biomedcentral.com/1471-2105/16/374
[installLink]:https://github.com/SouthGreenPlatform/TOGGLE/blob/master/INSTALL.md
[manualLink]:https://github.com/SouthGreenPlatform/TOGGLE/blob/master/MANUAL.md
[knownIssues]:https://github.com/SouthGreenPlatform/TOGGLE-DEV/blob/master/KnownIssues.md
[releaseLink]:https://github.com/SouthGreenPlatform/TOGGLE/blob/master/ReleaseNotes.md

