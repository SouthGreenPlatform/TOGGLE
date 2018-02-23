package versionSofts;

###################################################################################################################################
#
# Copyright 2014-2018 IRD-CIRAD-INRA-ADNid
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/> or
# write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston,
# MA 02110-1301, USA.
#
# You should have received a copy of the CeCILL-C license with this program.
#If not see <http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt>
#
# Intellectual property belongs to IRD, CIRAD and South Green developpement plateform for all versions also for ADNid for v2 and v3 and INRA for v3
# Version 1 written by Cecile Monat, Ayite Kougbeadjo, Christine Tranchant, Cedric Farcy, Mawusse Agbessi, Maryline Summo, and Francois Sabot
# Version 2 written by Cecile Monat, Christine Tranchant, Cedric Farcy, Enrique Ortega-Abboud, Julie Orjuela-Bouniol, Sebastien Ravel, Souhila Amanzougarene, and Francois Sabot
# Version 3 written by Cecile Monat, Christine Tranchant, Laura Helou, Abdoulaye Diallo, Julie Orjuela-Bouniol, Sebastien Ravel, Gautier Sarah, and Francois Sabot
#
###################################################################################################################################

use strict;
use warnings;
use localConfig;
use toolbox;
use Data::Dumper;
use Switch;

sub bowtieVersion
{   #We works with the STDIN output
	my $version = `$bowtie --version | grep "bowtie version"` or die toolbox::exportLog("ERROR: versionSoft::bowtieVersion : Can not grep bowtie version.\nPlease check your bowtie installation.\n", 0);
	chomp($version);
	return $version;
}
sub bowtie2Version
{   #We works with the STDIN output
	my $version = `$bowtie2 --version | grep "bowtie2-align-s version"` or die toolbox::exportLog("ERROR: versionSoft::bowtie2Version : Can not grep bowtie2 version.\nPlease check your bowtie installation.\n", 0);
	chomp($version);
	return $version;
}
sub bowtieBuildVersion
{   #We works with the STDOUT output
	my $version = `$bowtieBuild --version 2>&1 | grep "bowtie-build version"` or die toolbox::exportLog("ERROR: versionSoft::bowtieBuildVersion : Can not grep bowtieBuild version.\nPlease check your bowtie installation.\n", 0);
	chomp($version);
	return $version;
}

sub bowtie2BuildVersion
{   #We works with the STDOUT output
	my $version = `$bowtie2Build --version 2>&1 | grep "bowtie2-build version"` or die toolbox::exportLog("ERROR: versionSoft::bowtie2BuildVersion : Can not grep bowtie2Build version.\nPlease check your bowtie2 installation.\n", 0);
	chomp($version);
	return $version;
}

sub bwaVersion
{    #We works with the STDERR output
	my $version = `$bwa 2>&1 | grep "Version"` or die toolbox::exportLog("ERROR: versionSoft::bwaVersion : Can not grep bwa version\nPlease check your bwa installation.\n", 0);
	chomp($version);
	return $version;
}

sub cufflinksVersion
{   #We works with the STDERR output
	my $version = `$cufflinks/cufflinks 2>&1 | grep "cufflinks v"`or die toolbox::exportLog("ERROR: versionSoft::cufflinksVersion : Can not grep cufflinks version\nPlease check your cufflinks installation.\n", 0);
	chomp($version);
	return $version;
}

sub cutadaptVersion
{   #We works with the STDERR output
	my $version = `$cutadapt 2>&1 | grep "cutadapt version"` or die toolbox::exportLog("ERROR: versionSoft::cutadaptVersion : Can not grep cutadapt version\nPlease check your cutadapt installation.\n", 0);
	chomp($version);
	return $version;
}

sub atroposVersion
{   #We works with the STDERR output
	my $version = `$atropos 2>&1 | grep "Atropos version"` or die toolbox::exportLog("ERROR: versionSoft::atroposVersion : Can not grep Atropos version\nPlease check your Atropos installation.\n", 0);
	chomp($version);
	return $version;
}

sub fastqcVersion
{
	my $version = `$fastqc -v ` or die toolbox::exportLog("ERROR: versionSoft::fastqcVersion : Can not grep fastqc version\nPlease check your fastqc installation.\n", 0);
	chomp($version);
	return $version;
}

sub fastxToolkitVersion
{   #We works with the STDOUT output
	my $version = `$fastxTrimmer -h | grep "FASTX Toolkit"` or die toolbox::exportLog("ERROR: versionSoft::fastqxToolkitVersion : Can not grep fastxToolkit version\nPlease check your fastxToolkit installation.\n", 0);
	chomp($version);
	return $version;
}

sub gatkVersion
{
	my $version = `$GATK -version` or die toolbox::exportLog("ERROR: versionSoft::gatkVersion : Can not grep gatk version\nPlease check your GATK installation.\n", 0); #We works with the STDOUT output
	chomp($version);
	return $version;
}

sub htseqcountVersion
{    #We works with the STDOUT output
	my $version = `$htseqcount -h | grep "version" | cut -d"," -f 2,2` or die toolbox::exportLog("ERROR: versionSoft::htseqcountVersion : Can not grep htseqcount version\nPlease check your HTseq-Count installation.\n", 0);
	chomp($version);
	return $version;
}

sub javaVersion
{   #We works with the STDERR output
	my $version = `$java -version 2>&1 | grep "version"`or die toolbox::exportLog("ERROR: versionSoft::javaVersion : Can not grep java version\nPlease check your java installation.\n", 0);
	chomp($version);
	return $version;
}

sub picardToolsVersion
{
	my $version = `$picard CheckFingerprint --version 2>&1` or die toolbox::exportLog("ERROR: versionSoft::picardToolsVersion : Can not grep picardTools version\nPlease check your picardtools installation.\n", 0); #We works with the STDOUT output
	chomp($version);
	return $version;
}

sub stacksVersion
{
	my $version = `$stacks -v 2>&1` or die toolbox::exportLog("ERROR: versionSoft::stacksVersion : Can not grep stacks version\nPlease check your stacks installation.\n", 0); #We works with the STDOUT output
	chomp($version);
	return $version;
}

sub samtoolsVersion
{   #We works with the STDERR output
	my $version = `$samtools 2>&1 | grep "Version"` or die toolbox::exportLog("ERROR: versionSoft::samtoolsVersion : Can not grep samtools version\nPlease check your samtools installation.\n", 0);
	chomp($version);
	return $version;
}

sub snpeffVersion
{
	my $version = `$snpEff -version 2>&1` or die toolbox::exportLog("ERROR: versionSoft::snpeffVersion : Can not grep snpeff version\nPlease check your snpeff installation.\n", 0); #We works with the STDOUT output
	chomp($version);
	return $version;
}

sub tgiclVersion
{
	return "No version available";
}

sub tophatVersion
{
	my $version = `$tophat2 -v` or die toolbox::exportLog("ERROR: versionSoft::tophatVersion : Can not grep tophat version\nPlease check your tophat installation.\n", 0); #We works with the STDOUT output
	chomp($version);
	return $version;
}

sub trinityVersion
{   #We works with the STDOUT output
	my $version = `$trinity --version | grep "Trinity version"` or die toolbox::exportLog("ERROR: versionSoft::trinityVersion : Can not grep trinity version\nPlease check your trinity installation.\n", 0);
	chomp($version);
	return $version;
}

sub bamutilsVersion
{   #We works with the STDOUT output
	my $version = `$bamutils | tail -1` or die toolbox::exportLog("ERROR: versionSoft::bamutilsVersion : Can not grep bamutils version\nPlease check your bamutils/ngsutils installation.\n", 0);
	chomp($version);
	return $version;
}
sub plinkVersion
{   #We works with the STDOUT output
	my $version = `$plink -version` or die toolbox::exportLog("ERROR: versionSoft::plinkVersion : Can not grep plink version\nPlease check your plink installation.\n", 0);
	chomp($version);
	return $version;
}

sub fastmeVersion
{   #We works with the STDOUT output
	my $version = `$fastme --version` or die toolbox::exportLog("ERROR: versionSoft::fastmeVersion : Can not grep fastme version\nPlease check your fastme installation.\n", 0);
	chomp($version);
	return $version;
}

sub readseqVersion
{   #We works with the STDOUT output
	my $version = `$readseqjar -h | head -1` or die toolbox::exportLog("ERROR: versionSoft::readseqVersion : Can not grep readseq version\nPlease check your fastme installation.\n", 0);
	chomp($version);
	return $version;
}

sub cracVersion
{ #We works with the STDOUT output
	my $version = `$crac -version | grep -m 1 "version"` or die toolbox::exportLog("ERROR: versionSoft::cracVersion : Can not grep CRAC version\nPlease check your CRAC installation.\n", 0);
	chomp($version);
	return $version;
}

sub bedToolsVersion
{ #We works with the STDOUT output
	my $version = `$bedtools --version` or die toolbox::exportLog("ERROR: versionSoft::bedToolsVersion : Can not grep BEDtools version\nPlease check your BEDtools installation.\n", 0);
	chomp($version);
	return $version;

}

sub abyssVersion
{ #We works with the STDOUT output
	my $version = `$abyss --version | grep 'GNU Make'` or die toolbox::exportLog("ERROR: versionSoft::abyssVersion : Can not grep Abyss version\nPlease check your Abyss installation.\n", 0);
	chomp($version);
	return $version;

}

#sub transAbyssVersion
#{ #We works with the STDOUT output
#	my $version = `$transAbyss --version 2>&1` or die toolbox::exportLog("ERROR: versionSoft::transAbyssVersion : Can not grep transAbyss version\nPlease check your transAbyss installation.\n", 0);
#	chomp($version);
#	return $version;
#
#}

sub breakDancerVersion
{   #We works with the STDERR output
	my $version = `$breakDancer 2>&1 | grep Version` or die toolbox::exportLog("ERROR: versionSoft::breakDancerVersion : Can not grep breakDancer version\nPlease check your breakDancer installation.\n", 0);
	chomp($version);
	return $version;
}

sub pindelVersion
{
	my $version = `$pindel | grep -m 1 version` or die toolbox::exportLog("ERROR: versionSoft::pindelVersion : Can not grep breakDancer version\nPlease check your pindel installation.\n", 0);
	chomp($version);
	return $version;
}

sub fastqStatsVersion
{
	my $version =`$fastqStats -h | grep Version` or die toolboox::exportLog("ERROR: versionSoft::fastqStatsVersion: Cannot grep fastq-stats version\nPlease check your fastq-stats/ea-utils installation.\n",0);
	chomp $version;
	return $version;
}

sub writeLogVersion
{
	my ($fileConf, $version, $reportDir,$report) = @_; #recovery $report boolean value: set to 1 if report is requested. $reportDir is the path were software.txt is generated
	$report=0 if not defined $report; # by default $report does not generate sofware.txt

	my %softPathVersion = ("toggle"	=> $version);
	my %softPath = ("toggle"	=> $toggle);

	toolbox::checkFile($fileConf);
	my $configInfo=toolbox::readFileConf($fileConf);
	my $hashOrder=toolbox::extractHashSoft($configInfo,"order");	#Picking up the options for the order of the pipeline

	for my $softOrder ( values %{ $hashOrder } )
	{
		#DEBUG: print $softOrder." DANS LA BOUCLE\n";
		switch (1)
		{
			#FOR bwa.pm
			case ($softOrder =~ m/^bwa.*/i){$softPathVersion{"bwa"}= bwaVersion if not defined $softPathVersion{"bwa"};
											$softPath{"bwa"}= $bwa if not defined $softPath{"bwa"};
											}

			#FOR samTools.pm
			case ($softOrder =~ m/^samtools.*/i){$softPathVersion{"samtools"}= samtoolsVersion if not defined $softPathVersion{"samtools"};
												 $softPath{"samtools"}= $samtools if not defined $softPath{"samtools"};
												}

			#FOR picardTools.pm
			case ($softOrder =~ m/^picard.*/i){ $softPathVersion{"java"}= javaVersion if not defined $softPathVersion{"java"};
												$softPath{"java"}= $java if not defined $softPath{"java"};
												$softPathVersion{"picard"}= picardToolsVersion if not defined $softPathVersion{"picard"};
												$softPath{"picard"}= $picard if not defined $softPath{"picard"};
											   }

			#FOR gatk.pm
			case ($softOrder =~ m/^gatk.*/i){$softPathVersion{"java"}= javaVersion if not defined $softPathVersion{"java"};
											 $softPathVersion{"GATK"}= gatkVersion if not defined $softPathVersion{"GATK"};
											 $softPath{"java"}= $java if not defined $softPath{"java"};
											 $softPath{"GATK"}= $GATK if not defined $softPath{"GATK"};
											 }

			#FOR fastqc
			case ($softOrder =~ m/^fastqc/i){$softPathVersion{"fastqc"}= fastqcVersion if not defined $softPathVersion{"fastqc"};
											 $softPath{"fastqc"}= $fastqc if not defined $softPath{"fastqc"};
											 }

			#FOR fastxToolkit
			case ($softOrder =~ m/^fastx.*/i){$softPathVersion{"fastxTrimmer"}= fastxToolkitVersion if not defined $softPathVersion{"fastxTrimmer"};
											  $softPath{"fastxTrimmer"}= $fastxTrimmer if not defined $softPath{"fastxTrimmer"};
											  }

			#FOR tophat.pm
			case ($softOrder =~ m/^bowtie2.*/i){$softPathVersion{"bowtie2Build"}= bowtie2BuildVersion if not defined $softPathVersion{"bowtie2Build"};
												$softPath{"bowtie2Build"}= $bowtie2Build if not defined $softPath{"bowtie2Build"};
												$softPathVersion{"bowtie2"}= bowtie2Version if not defined $softPathVersion{"bowtie2"};
												$softPath{"bowtie2"}= $bowtie2 if not defined $softPath{"bowtie2"};
												}
			case ($softOrder =~ m/^bowtie$/i){$softPathVersion{"bowtieBuild"}= bowtieBuildVersion if not defined $softPathVersion{"bowtieBuild"};
											 $softPath{"bowtieBuild"}= $bowtieBuild if not defined $softPath{"bowtieBuild"};
											 $softPathVersion{"bowtie"}= bowtieVersion if not defined $softPathVersion{"bowtie"};
											 $softPath{"bowtie"}= $bowtie if not defined $softPath{"bowtie"};
											 }

			case ($softOrder =~ m/^tophat.*/i){$softPathVersion{"tophat2"}= tophatVersion if not defined $softPathVersion{"tophat2"};
											   $softPathVersion{"bowtieBuild"}= bowtieBuildVersion if not defined $softPathVersion{"bowtieBuild"};
											   $softPathVersion{"bowtie2Build"}= bowtie2BuildVersion if not defined $softPathVersion{"bowtie2Build"};
											   $softPath{"tophat2"}= $tophat2 if not defined $softPath{"tophat2"};
											   $softPath{"bowtieBuild"}= $bowtieBuild if not defined $softPath{"bowtieBuild"};
											   $softPath{"bowtie2Build"}= $bowtie2Build if not defined $softPath{"bowtie2Build"};
											   }

			#FOR cufflinks.pm
			case ($softOrder =~ m/^cufflinks.*/i){$softPathVersion{"cufflinks"}= cufflinksVersion if not defined $softPathVersion{"cufflinks"};
												  $softPath{"cufflinks"}= $cufflinks if not defined $softPath{"cufflinks"};
												  }
			case ($softOrder =~ m/^cuffdiff.*/i){$softPathVersion{"cuffdiff"}= cufflinksVersion if not defined $softPathVersion{"cuffdiff"};
												 $softPath{"cuffdiff"}= $cufflinks if not defined $softPath{"cuffdiff"};
												 }
			case ($softOrder =~ m/^cuffmerge.*/i){$softPathVersion{"cuffmerge"}= cufflinksVersion if not defined $softPathVersion{"cuffmerge"};
												  $softPath{"cuffmerge"}= $cufflinks if not defined $softPath{"cuffmerge"};
												  }

			#FOR HTSeq
			case ($softOrder =~ m/^htseq.*/i){$softPathVersion{"htseqcount"}= htseqcountVersion if not defined $softPathVersion{"htseqcount"};
											  $softPath{"htseqcount"}= $htseqcount if not defined $softPath{"htseqcount"};
											  }

			#FOR snpEff.pm
			case ($softOrder =~ m/^snp.*/i){$softPathVersion{"java"}= javaVersion if not defined $softPathVersion{"java"};
											$softPathVersion{"snpEff"}= snpeffVersion if not defined $softPathVersion{"snpEff"};
											$softPath{"java"}= $java if not defined $softPath{"java"};
											$softPath{"snpEff"}= $snpEff if not defined $softPath{"snpEff"};
											}

			#FOR processRadtags.pm
			case ($softOrder =~ m/process.*/i){$softPathVersion{"stacks"}= stacksVersion if not defined $softPathVersion{"stacks"};
											   $softPath{"stacks"}= $stacks if not defined $softPath{"stacks"};
											   }

			#FOR cutadapt functions
			case ($softOrder =~ m/^cutadapt/i){$softPathVersion{"cutadapt"}= cutadaptVersion if not defined $softPathVersion{"cutadapt"};
											   $softPath{"cutadapt"}= $cutadapt if not defined $softPath{"cutadapt"};
											  }
			#FOR atropos functions
			case ($softOrder =~ m/^atropos/i){$softPathVersion{"atropos"}= atroposVersion if not defined $softPathVersion{"atropos"};
											   $softPath{"atropos"}= $atropos if not defined $softPath{"atropos"};
											  }

			#FOR TGICL
			case ($softOrder =~ m/^tgicl/i){$softPathVersion{"tgicl"}= tgiclVersion if not defined $softPathVersion{"tgicl"};
											$softPath{"tgicl"}= $tgicl if not defined $softPath{"tgicl"};
											}

			#FOR trinity
			case ($softOrder =~ m/^trinity/i){$softPathVersion{"trinity"}= trinityVersion if not defined $softPathVersion{"trinity"};
											  $softPath{"trinity"}= $trinity if not defined $softPath{"trinity"};
											  }
			#FOR bamutils
			case ($softOrder =~ m/^bamutils.*/i){$softPathVersion{"bamutils"}= bamutilsVersion if not defined $softPathVersion{"bamutils"};
											  $softPath{"bamutils"}= $bamutils if not defined $softPath{"bamutils"};
											  }
			#FOR plink
			case ($softOrder =~ m/^plink.*/i){$softPathVersion{"plink"}= plinkVersion if not defined $softPathVersion{"plink"};
											  $softPath{"plink"}= $plink if not defined $softPath{"plink"};
											  }
			#FOR fastme
			case ($softOrder =~ m/^fastme.*/i){$softPathVersion{"fastme"}= fastmeVersion if not defined $softPathVersion{"fastme"};
											  $softPath{"fastme"}= $fastme if not defined $softPath{"fastme"};
											  }

			#FOR readseq
			case ($softOrder =~ m/^readseq.*/i){$softPathVersion{"readseq"}= readseqVersion if not defined $softPathVersion{"readseq"};
											  $softPath{"readseq"}= $readseqjar if not defined $softPath{"readseq"};
											  }

			#FOR SNIPLAY
			case ($softOrder =~ m/^sniplay.*/i){$softPathVersion{"sniplay"}= "v1.0" if not defined $softPathVersion{"sniplay"};
												$softPath{"sniplay"}= "sniplay" if not defined $softPath{"sniplay"};
											   }

			#For format checking
			case($softOrder =~ m/^check/i){next;}

			#FOR crac.pm
			case ($softOrder =~ m/^crac.*/i){$softPathVersion{"crac"}= cracVersion if not defined $softPathVersion{"crac"};
											$softPath{"crac"}= $crac if not defined $softPath{"crac"};
											}
			#FOR duplicationDetector
			case ($softOrder =~ m/^duplicationDetector/i){$softPathVersion{"duplicationDetector"}= "v1.0" if not defined $softPathVersion{"duplicationDetector"};
											$softPath{"duplicationDetector"}= $duplicationDetector if not defined $softPath{"duplicationDetector"};}

			#FOR checkFormatFasta
			case ($softOrder =~ m/^checkFormatFasta/i){$softPathVersion{"checkFormatFasta"}= "v1.0" if not defined $softPathVersion{"checkFormatFasta"};
											$softPath{"checkFormatFasta"}= "checkFormatFasta" if not defined $softPath{"checkFormatFasta"};}
			#FOR checkFormatFastq
			case ($softOrder =~ m/^checkFormatFastq/i){$softPathVersion{"checkFormatFastq"}= "v1.0" if not defined $softPathVersion{"checkFormatFasta"};
											$softPath{"checkFormatFastq"}= "checkFormatFastq" if not defined $softPath{"checkFormatFastq"};}
			#FOR checkFormatVcf
			case ($softOrder =~ m/^checkFormatVcf/i){$softPathVersion{"checkFormatVcf"}= "v1.0" if not defined $softPathVersion{"checkFormatVcf"};
											$softPath{"checkFormatVcf"}= "checkFormatVcf" if not defined $softPath{"checkFormatVcf"};}
			#FOR checkFormatSamOrBam
			case ($softOrder =~ m/^checkFormatSamOrBam/i){$softPathVersion{"checkFormatSamOrBam"}= "v1.0" if not defined $softPathVersion{"checkFormatSamOrBam"};
											$softPath{"checkFormatSamOrBam"}= "checkFormatSamOrBam" if not defined $softPath{"checkFormatSamOrBam"};}
			#FOR checkEncodeByASCIIcontrol
			case ($softOrder =~ m/^checkEncodeByASCIIcontrol/i){$softPathVersion{"checkEncodeByASCIIcontrol"}= "v1.0" if not defined $softPathVersion{"checkEncodeByASCIIcontrol"};
																 $softPath{"checkEncodeByASCIIcontrol"}= "checkEncodeByASCIIcontrol" if not defined $softPath{"checkEncodeByASCIIcontrol"};}

			#FOR checkFormatGff
			case ($softOrder =~ m/^checkFormatGff/i){$softPathVersion{"checkFormatGff"}= "v1.0" if not defined $softPathVersion{"checkFormatGff"};
											$softPath{"checkFormatGff"}= "checkFormatGff" if not defined $softPath{"checkFormatGff"};}

			#FOR checkFormatBed
			case ($softOrder =~ m/^checkFormatBed/i){$softPathVersion{"checkFormatBed"}= "v1.0" if not defined $softPathVersion{"checkFormatBed"};
											$softPath{"checkFormatBed"}= "checkFormatBed" if not defined $softPath{"checkFormatBed"};}

			#FOR BEDtools
			case ($softOrder =~ m/^bedtools/i){$softPathVersion{"bedtools"}= bedToolsVersion if not defined $softPathVersion{"bedtools"};
												$softPath{"bedtools"}= $bedtools if not defined $softPath{"bedtools"};}
			case ($softOrder =~ m/.*bed$/i){$softPathVersion{"bedtools"}= bedToolsVersion if not defined $softPathVersion{"bedtools"};
												$softPath{"bedtools"}= $bedtools if not defined $softPath{"bedtools"};}

			#FOR generic command
			case ($softOrder =~ m/^generic/i){$softPathVersion{"generic"}= "v0.1" if not defined $softPathVersion{"generic"};
												$softPath{"generic"}= "" if not defined $softPath{"generic"};}

			#FOR Abyss
			case ($softOrder =~ m/^abyss/i){$softPathVersion{"abyss"}= abyssVersion if not defined $softPathVersion{"abyss"};
												$softPath{"abyss"}= $abyss if not defined $softPath{"abyss"};}

			#FOR transAbyss
			#case ($softOrder =~ m/^transAbyss/i){$softPathVersion{"transAbyss"}= transAbyssVersion if not defined $softPathVersion{"transAbyss"};
												#$softPath{"transAbyss"}= $abyss if not defined $softPath{"transAbyss"};}

			#For breakDancer
			case ($softOrder =~ m/^bam2cfg/i){$softPathVersion{"breakDancer"}= breakDancerVersion if not defined $softPathVersion{"breakDancer"};
												$softPath{"bam2cfg"}= $bam2cfg if not defined $softPath{"bam2cfg"};}
			case ($softOrder =~ m/^breakDancer/i){$softPathVersion{"breakDancer"}= breakDancerVersion if not defined $softPathVersion{"breakDancer"};
												$softPath{"breakDancer"}= $breakDancer if not defined $softPath{"breakDancer"};}

			#For Pindel
			case ($softOrder =~ m/^pindel/i){$softPathVersion{"pindel"}= pindelVersion if not defined $softPathVersion{"pindel"};
												$softPath{"pindel"}= $pindel if not defined $softPath{"pindel"};}

			#For ea-Utils
			case ($softOrder =~ m/^fastqStats/i){$softPathVersion{"fastqStats"}=fastqStatsVersion if not defined $softPathVersion{"fastqStats"};
												$softPath{"fastqStats"}=$fastqStats if not defined $softPath{"fastqStats"};}

			else
			{
				toolbox::exportLog("ERROR VERSIONSOFTS: $0 : the $softOrder function or software is unknown to TOGGLE, cannot continue",0);
			}; # Name unknown to TOGGLE, must stop
		}
	}
	## DEBUG print Dumper(%softPathVersion);

	open (my $fhConfig, "<", "$toggle/modules/localConfig.pm");
	open (my $fhSoft, ">", "$reportDir/software.txt") if $report;
	#open (my $fhSoft, ">", "$reportDir/software.tex") if $report;
    #print $fhSoft "\\begin{itemize}\n";

	while (my $line = <$fhConfig>)
	{
		no strict "vars";
		chomp $line;
		chop $line; #Remove the last character, ie ";"
		next unless $line =~ m/^our \$/;
		my ($soft,$value) = split /=/, $line;
		$soft =~ s/our| |\$//g;
		if (defined $softPathVersion{$soft})
		{
			toolbox::exportLog(uc($soft)." : $softPath{$soft} : $softPathVersion{$soft}",1);
			#print $fhSoft "\\item"  .			uc($soft) ." : \n \\begin{verbatim}\n$softPathVersion{$soft} \n \\end{verbatim}" if $report;
			print $fhSoft uc($soft) ." : ".$softPathVersion{$soft}." (*\@{\\cite{$soft}}\@*)\n" if $report;
		}
	}

	#print $fhSoft "\\end{itemize}\n";
	close $fhConfig;
	close $fhSoft if $report;
}
1;



=head1 NAME

    Package I<versionSofts>

=head1 SYNOPSIS

        use versionSofts;

        versionSofts::writeLogVersion ();

=head1 DESCRIPTION

    Package Version Softs

=head2 FUNCTIONS

=head3 versionSofts::writeLogVersion

This module return soft version

=head1 AUTHORS

Intellectual property belongs to IRD, CIRAD and South Green developpement plateform
Written by Cecile Monat, Ayite Kougbeadjo, Marilyne Summo, Cedric Farcy, Mawusse Agbessi, Christine Tranchant and Francois Sabot

=head1 SEE ALSO

L<http://www.southgreen.fr/>

=cut
