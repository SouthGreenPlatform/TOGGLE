TOGGLE : Toolbox for generic NGS analyses
===========

![TOGGLE Logo](docs/images/toggleLogo.png)

Dear Biologist, have you ever dream of using the whole power of those numerous NGS tools that your bioinformatician colleagues use through this awful list of command line ?

Dear Bioinformatician, have you ever guess how to design really fastly a new NGS pipeline without having to retype again dozens of code lines to readapt your scripts or starting from scratch ?

**So, be Happy ! TOGGLE is for you !!**

TOGGLE (TOolbox for Generic nGs anaLysEs) is a suite of 19 packages and more than 110 modules able to manage a large set of NGS softwares
and utilities to easily design pipelines able to handle hundreds of samples. Moreover, TOGGLE offers an easy way to manipulate the various
options of the different softwares through the pipelines in using a single basic configuration file, that can be changed for each assay without
having to change the code itself.

Users can also create their own pipeline through an easy and user-friendly approach. The pipelines can start from Fastq (plain or gzipped), SAM, BAM or VCF (plain or gzipped) files, with parallel and global analyses. Samples pipelines are provided for SNP discovery and RNAseq counts.

The system is able to detect parallel/scheduling launching and to manage large amount of samples on large cluster machines.

## Development version
You are currently on the production version, *i.e.* the last stable one. A more advanced version, the dev version is available on the [TOGGLE-DEV gitHub](https://github.com/SouthGreenPlatform/TOGGLE-DEV)


##  Contributing

* Licencied under CeCill-C (http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html) and GPLv3
* Intellectual property belongs to IRD, CIRAD, ADNid and SouthGreen development platform
* Written by Cecile Monat, Christine Tranchant, Ayite Kougbeadjo, Cedric Farcy, Mawusse Agbessi, Enrique Ortega-Abboud, Sébastien Ravel, Julie Orjuela-Bouniol, Souhila Amanzougarene, Gauthier Sarah, Marilyne Summo, and Francois Sabot
* Copyright 2014-2016

## Contact

For bug tracking purpose you can use the GitHub or questions about TOGGLE, you can contact the mainteners using the following email addresses:

* christine.tranchant@ird.fr
* francois.sabot@ird.fr

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

* [Capture::Tiny](http://search.cpan.org/~dagolden/Capture-Tiny-0.30/lib/Capture/Tiny.pm)
* [Data::Translate](http://search.cpan.org/~davieira/Data_Translate-0.3/Translate.pm)
* [Data::Dumper](http://search.cpan.org/~smueller/Data-Dumper-2.154/Dumper.pm)
* [Getopt::ArgParse](http://search.cpan.org/~mytram/Getopt-ArgParse-1.0.2/lib/Getopt/ArgParse.pm)
* [List::Compare](http://search.cpan.org/~jkeenan/List-Compare-0.53/lib/List/Compare.pm)
* [Switch](https://metacpan.org/pod/Switch)
* [Test::More](http://search.cpan.org/~exodist/Test-Simple-1.001014/lib/Test/More.pm)
* [Test::Deep](http://search.cpan.org/~rjbs/Test-Deep-0.119/lib/Test/Deep.pm)

#### Bioinformatics software (minimal version)

<table border="1" cellpadding="5" cellspacing="1" >
<thead>
<tr bgcolor="#68AFDF">
	<th> Category        </th>
	<th> Software        </th>
	<th> Minimal version </th>
	<th> Tools included  </th>
</tr>
</thead>
<tbody>
<tr bgcolor="#B9F0FB" >
	<td align="center" rowspan="2" > <b>System<b>   </td>
	<td> <a href="https://www.java.com">Java</a></td>
	<td> 1.7 non open JDK </td>
	<td></td>
</tr>
<tr bgcolor="#B9F0FB">
	<td> <a href="https://www.perl.org/">Perl</a></td>
	<td> 5.16 </td>
	<td>  </td>
</tr>
<tr bgcolor="#54D1C0">
	<td align="center" rowspan="3"> <b> Mapping <b> </td>
	<td> <a href="http://bowtie-bio.sourceforge.net/bowtie2/index.shtml">Bowtie2</a></td>
	<td> 2.2.9 </td>
	<td> bowtie2build </td>
</tr>
<tr bgcolor="#54D1C0">
	<td> <a href="http://bio-bwa.sourceforge.net/">BWA</a></td>
	<td> 0.7.2 </td>
	<td> bwaAln </br> bwaSampe </br> bwaSamse </br> bwaMem </br> bwaIndex</td>
</tr>
<tr bgcolor="#54D1C0">
	<td> <a href="https://ccb.jhu.edu/software/tophat/index.shtml">Tophat</a></td>
	<td> 2.1.1 </td>
	<td> bowtiebuild </br> tophat2 </td>
</tr>
<tr bgcolor="#B9F0FB">
	<td align="center" rowspan="6"> <b> Fastq/BAM/SAM </br> tools <b> </td>
	<td> <a href="https://www.broadinstitute.org/gatk/">GATK</a></td>
	<td> 3.3 </td>
	<td> gatkBaseRecalibrator </br> gatkRealignerTargetCreator </br> gatkIndelRealigner </br> gatkHaplotypeCaller </br> gatkSelectVariants </br> gatkVariantFiltration </br> gatkReadBackedPhasing </br> gatkUnifiedGenotyper </br> gatkBaseRecalibrator </br> gatkPrintReads
 </td>
</tr>
<tr bgcolor="#B9F0FB">
	<td> <a href="http://samtools.sourceforge.net/">SAMtools</a></td>
	<td> 0.1.18 </td>
	<td> samToolsFaidx </br> samToolsIndex </br> samToolsView </br> samToolsSort </br> mergeHeader </br> samToolsMerge </br> samToolsIdxstats </br> samToolsDepth </br> samToolsFlagstat </br> samToolsMpileUp
</td>
</tr>
<tr bgcolor="#B9F0FB">
	<td> <a href="http://broadinstitute.github.io/picard/">PicardTools</a></td>
	<td> 1.63 </td>
	<td> picardToolsMarkDuplicates </br> picardToolsCreateSequenceDictionary </br> picardToolsSortSam </br> picardToolsAddOrReplaceReadGroup </br> picardToolsValidateSamFile </br> picardToolsCleanSam </br> picardToolsSamFormatConverter
 </td>
</tr>
<tr bgcolor="#B9F0FB">
	<td> <a href="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/">FastQC</a></td>
	<td> 0.10.1 </td>
	<td> fastqc </td>
</tr>
<tr bgcolor="#B9F0FB">
	<td> <a href="https://pypi.python.org/pypi/cutadapt">Cutadapt</a></td>
	<td> 1.2.1 </td>
	<td> cutadapt </td>
</tr>
<tr bgcolor="#B9F0FB">
	<td> <a href="http://hannonlab.cshl.edu/fastx_toolkit/">FastxToolkit</a></td>
	<td> 0.0.13 </td>
	<td> fastxTrimmer </td>
</tr>
<tr bgcolor="#54D1C0">
	<td align="center"> <b> Demultiplexing <b> </td>
	<td> <a href="http://catchenlab.life.illinois.edu/stacks/">Stacks</a></td>
	<td> 1.43 </td>
	<td> process_radtags </td>
</tr>
<tr bgcolor="#B9F0FB">
	<td align="center" > <b> VCF <b> </td>
	<td> <a href="http://snpeff.sourceforge.net/">Snpeff</a></td>
	<td> 4.2 </td>
	<td> snpeffAnnotation </td>
</tr>
<tr bgcolor="#54D1C0">
	<td align="center" rowspan="2"> <b> RNA Seq <b> </td>
	<td> <a href="http://www-huber.embl.de/HTSeq/doc/count.html">HTSeq-Count</a></td>
	<td>  </td>
	<td>  </td>
</tr>
<tr bgcolor="#54D1C0">
	<td> <a href="https://github.com/trinityrnaseq/trinityrnaseq/wiki">Trinity</a></td>
	<td> 2.2.0 </td>
	<td>  </td>
</tr>
<tr bgcolor="#B9F0FB">
	<td align="center" rowspan="2"> <b> Assembly <b> </td>
	<td> <a href="https://sourceforge.net/projects/tgicl/files/">TGICL</a></td>
	<td> 2.2.0 </td>
	<td>  </td>
</tr>
<tr bgcolor="#B9F0FB">
	<td> <a href="http://gmt.genome.wustl.edu/packages/pindel/">Pindel</a></td>
	<td> 0.2.4 </td>
	<td>  </td>
</tr>
<tr bgcolor="#54D1C0">
	<td align="center" rowspan="1"> <b> Structural Variants <b> </td>
	<td> <a href="http://breakdancer.sourceforge.net/">Breakdancer</a></td>
	<td> 1.4.5 </td>
	<td>  </td>
</tr>
<tr bgcolor="#B9F0FB">
	<td align="center" rowspan="1"> <b> Optional <b> </td>
	<td> <a href="http://www.graphviz.org/">Graphviz</a></td>
	<td> v2.xx </td>
	<td>  </td>
</tr>
</tbody>
</table>



[paperLink]:http://www.biomedcentral.com/1471-2105/16/374
[installLink]:https://github.com/SouthGreenPlatform/TOGGLE/blob/master/docs/INSTALL.md
[manualLink]:https://github.com/SouthGreenPlatform/TOGGLE/blob/master/docs/MANUAL.md
[knownIssues]:https://github.com/SouthGreenPlatform/TOGGLE-DEV/blob/master/docs/KnownIssues.md
[releaseLink]:https://github.com/SouthGreenPlatform/TOGGLE/blob/master/docs/ReleaseNotes.md
