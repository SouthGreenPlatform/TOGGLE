package softwareManagement;

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

###################################################################################################################################
#
# This package will change the names of softwares to be coherent through the toggle hashes.
#
###################################################################################################################################

use strict;
use warnings;
use Data::Dumper;
use Switch;

use localConfig;

sub softwareNomenclature # Will rewrite the correct name in the hash of configuration
{
    my ($hash) = @_;

    foreach my $currentSoft (keys %{$hash})
    {
        my $correctName;
        if ($currentSoft eq "order") # We aredealing with the order hash...
        {
            #Specific treatment
            my $hashOrder = $$hash{$currentSoft};
            foreach my $step (keys %{$hashOrder})
            {
                $$hashOrder{$step}=correctName($$hashOrder{$step}); # will change in order accordingly
            }

            next;
        }
        ##DEBUG print "---------------$currentSoft-->";
        $correctName=correctName($currentSoft);
        ##DEBUG print "$correctName--------\n";
        if ($currentSoft ne $correctName) # the name has changed
        {
            ##DEBUG print Dumper($hash);
            $$hash{$correctName}=\%{$$hash{$currentSoft}};
            ##DEBUG print Dumper($hash);
            delete $$hash{$currentSoft};
            ##DEBUG print Dumper($hash);
        }

    }
    return $hash;
}

sub correctName
{
    my ($name)=@_;
    my $correctedName="NA";
    my $order;
    ## DEBUG toolbox::exportLog("++++++++++++++$name\n",1);
    my @list = split /\s/,$name;
    $order = pop @list if ($list[-1] =~ m/^\d+/); # This is for a repetition of the same step
    switch (1)
    {
        #FOR cleaner
        case ($name =~ m/^cleaner/i){$correctedName="cleaner";} #Correction for cleaner step

        #FOR compressor
        case ($name =~ m/^compress/i){$correctedName="compress";} #correction for compressor step

        #FOR merge
        case ($name =~ m/^merge/i){$correctedName="merge";} #correction for merge step

        #FOR env
        case ($name =~ m/^env/i){$correctedName="env";} #correction for merge step

        #FOR SGE
        case ($name =~ m/^sge/i){$correctedName="sge";} #Correction for sge configuration

        #FOR SLURM
        case ($name =~ m/^slurm/i){$correctedName="slurm";} #Correction for slurm configuration

        #FOR MPRUN
        case ($name =~ m/^mprun/i){$correctedName="mprun";} #Correction for mprun configuration

        #FOR LSF
        case ($name =~ m/^lsf/i){$correctedName="lsf";} #Correction for lsf configuration

        #FOR SCP
        case ($name =~ m/^scp/i or $name =~ m/^rsync/i or $name =~ m/^transfer/i){$correctedName="scp";} #Correction for scp transfer

        #NEW SOFT ADDED AUTOMATICALLY
        
		#FOR nanoplot
		case ($name =~ m/^nanoplot[\s|\.|\-| \/|\|\|]*/i){$correctedName="nanoplot";} #Correction for nanoplot

        #FOR hisat2
        case ($name =~ m/^hisat2[\s|\.|\-| \/|\|\|]*/i){$correctedName="hisat2";} #Correction for hisat2
	
        #FOR stringtie
        case ($name =~ m/^stringtie[\s|\.|\-| \/|\|\|]*/i){$correctedName="stringtie";} #Correction for stringtie
    
    
        #FOR bwa.pm
        case ($name =~ m/^bwa[\s|\.|\-| \/|\\|\|]*aln/i){$correctedName="bwaAln"; } #Correction for bwaAln
        case ($name =~ m/^bwa[\s|\.|\-| \/|\\|\|]*sampe/i){$correctedName="bwaSampe"} # Correction for bwaSampe
        case ($name =~ m/^bwa[\s|\.|\-| \/|\\|\|]*samse/i){$correctedName="bwaSamse"} # Correction for bwaSamse
        case ($name =~ m/^bwa[\s|\.|\-| \/|\\|\|]*index/i){$correctedName="bwaIndex"} # Correction for bwaIndex
        case ($name =~ m/^bwa[\s|\.|\-| \/|\\|\|]*mem/i){$correctedName="bwaMem"} # Correction for bwaMem
        case ($name =~ m/^bwa[\s|\.|\-| \/|\\|\|]*sw/i){$correctedName="bwaSw"} # Correction for bwaSw


        #FOR samTools.pm
        case ($name =~ m/^samtools[\s|\.|\-| \/|\\|\|]*faidx/i){$correctedName="samToolsFaidx"} # Correction for samToolsFaidx
        case ($name =~ m/^samtools[\s|\.|\-| \/|\\|\|]*index/i){$correctedName="samToolsIndex"} # Correction for samToolsIndex
        case ($name =~ m/^samtools[\s|\.|\-| \/|\\|\|]*view/i){$correctedName="samToolsView"} # Correction for samToolsView
        case ($name =~ m/^samtools[\s|\.|\-| \/|\\|\|]*sort/i){$correctedName="samToolsSort"} # Correction for samToolsSort
        case ($name =~ m/^merge[\s|\.|\-| \/|\\|\|]*header/i){$correctedName="mergeHeader"} # Correction for mergeHeader
        case ($name =~ m/^samtools[\s|\.|\-| \/|\\|\|]*merge/i){$correctedName="samToolsMerge"} # Correction for samToolsMerge
        case ($name =~ m/^samtools[\s|\.|\-| \/|\\|\|]*idxstats/i){$correctedName="samToolsIdxstats"} # Correction for samToolsIdxstats
        case ($name =~ m/^samtools[\s|\.|\-| \/|\\|\|]*depth/i){$correctedName="samToolsDepth"} # Correction for samToolsDepth
        case ($name =~ m/^samtools[\s|\.|\-| \/|\\|\|]*flagstat/i){$correctedName="samToolsFlagstat"} # Correction for samToolsFlagstat
        case ($name =~ m/^samtools[\s|\.|\-| \/|\\|\|]*mpileup/i){$correctedName="samToolsMpileUp"} # Correction for samToolsMpileUp

        #FOR picardTools.pm
        case ($name =~ m/^picardtools[\s|\.|\-| \/|\\|\|]*mark[\s|\.|\-| \/|\\|\|]*duplicates/i){$correctedName="picardToolsMarkDuplicates"} # Correction for picardToolsMarkDuplicates
        case ($name =~ m/^picardtools[\s|\.|\-| \/|\\|\|]*create[\s|\.|\-| \/|\\|\|]*sequence[\s|\.|\-| \/|\\|\|]*dictionary/i){$correctedName="picardToolsCreateSequenceDictionary"} # Correction for picardToolsCreateSequenceDictionary
        case ($name =~ m/^picardtools[\s|\.|\-| \/|\\|\|]*sort[\s|\.|\-| \/|\\|\|]*sam/i){$correctedName="picardToolsSortSam"} # Correction for picardToolsSortSam
        case ($name =~ m/^picardtools[\s|\.|\-| \/|\\|\|]*validate[\s|\.|\-| \/|\\|\|]*sam[\s|\.|\-| \/|\\|\|]*file/i){$correctedName="picardToolsValidateSamFile"} # Correction for picardToolsValidateSamFile
        case ($name =~ m/^picardtools[\s|\.|\-| \/|\\|\|]*clean[\s|\.|\-| \/|\\|\|]*sam/i){$correctedName="picardToolsCleanSam"} # Correction for picardToolsCleanSam
        case ($name =~ m/^picardtools[\s|\.|\-| \/|\\|\|]*sam[\s|\.|\-| \/|\\|\|]*format[\s|\.|\-| \/|\\|\|]*converter/i){$correctedName="picardToolsSamFormatConverter"} # Correction for picardToolsSamFormatConverter
        case ($name =~ m/^picardtools[\s|\.|\-| \/|\\|\|]*add[\s|\.|\-| \/|\\|\|]*or[\s|\.|\-| \/|\\|\|]*replace[\s|\.|\-| \/|\\|\|]*read[\s|\.|\-| \/|\\|\|]*groups/i){$correctedName="picardToolsAddOrReplaceReadGroups"} # Correction for picardToolsAddOrReplaceReadGroups


        #FOR gatk.pm
        case ($name =~ m/^gatk[\s|\.|\-| \/|\\|\|]*base[\s|\.|\-| \/|\\|\|]*recalibrator/i){$correctedName="gatkBaseRecalibrator"} # Correction for gatkBaseRecalibrator
        case ($name =~ m/^gatk[\s|\.|\-| \/|\\|\|]*print[\s|\.|\-| \/|\\|\|]*reads/i){$correctedName="gatkPrintReads"} # Correction for gatkPrintReads
        case ($name =~ m/^gatk[\s|\.|\-| \/|\\|\|]*realigner[\s|\.|\-| \/|\\|\|]*target[\s|\.|\-| \/|\\|\|]*creator/i){$correctedName="gatkRealignerTargetCreator"} # Correction for gatkRealignerTargetCreator
        case ($name =~ m/^gatk[\s|\.|\-| \/|\\|\|]*indel[\s|\.|\-| \/|\\|\|]*realigner/i){$correctedName="gatkIndelRealigner"} # Correction for gatkIndelRealigner
        case ($name =~ m/^gatk[\s|\.|\-| \/|\\|\|]*haplotype[\s|\.|\-| \/|\\|\|]*caller/i){$correctedName="gatkHaplotypeCaller"} # Correction for gatkHaplotypeCaller
        case ($name =~ m/^gatk[\s|\.|\-| \/|\\|\|]*select[\s|\.|\-| \/|\\|\|]*variants/i){$correctedName="gatkSelectVariants"} # Correction for gatkSelectVariants
        case ($name =~ m/^gatk[\s|\.|\-| \/|\\|\|]*variant[\s|\.|\-| \/|\\|\|]*filtration/i){$correctedName="gatkVariantFiltration"} # Correction for gatkVariantFiltration
        case ($name =~ m/^gatk[\s|\.|\-| \/|\\|\|]*unified[\s|\.|\-| \/|\\|\|]*genotyper/i){$correctedName="gatkUnifiedGenotyper"} # Correction for gatkUnifiedGenotyper
        case ($name =~ m/^gatk[\s|\.|\-| \/|\\|\|]*read[\s|\.|\-| \/|\\|\|]*backed[\s|\.|\-| \/|\\|\|]*phasing/i){$correctedName="gatkReadBackedPhasing"} # Correction for gatkReadBackedPhasing

        #FOR fastqc
        case ($name =~ m/^fastqc/i){$correctedName="fastqc"} # Correction for fastqc

        #FOR fastqUtils.pm
        case ($name =~ m/^check[\s|\.|\-| \/|\\|\|]*encode[\s|\.|\-| \/|\\|\|]*by[\s|\.|\-| \/|\\|\|]*ascii[\s|\.|\-| \/|\\|\|]*control/i){$correctedName="checkEncodeByASCIIcontrol"} # Correction for checkEncodeByASCIIcontrol
        case ($name =~ m/^change[\s|\.|\-| \/|\\|\|]*encode/i){$correctedName="changeEncode"} # Correction for changeEncode

        #FOR fastxToolkit
        case ($name =~ m/^fastx[\s|\.|\-| \/|\\|\|]*trimmer/i){$correctedName="fastxTrimmer"} # Correction for fastxTrimmer

        #FOR tophat.pm
        case ($name =~ m/^tophat[\s|\.|\-| \/|\\|\|]*2/i){$correctedName="tophat2"; } #Correction for tophat2

        #FOR cufflinks.pm

        #FOR HTSeq.pm
        case ($name =~ m/^htseq[\s|\.|\-| \/|\\|\|]*count/i){$correctedName="htseqCount"; } #Correction for htseq-count

        #FOR snpEff.pm
        case ($name =~ m/^snp[\s|\.|\-| \/|\\|\|]*Eff[\s|\.|\-| \/|\\|\|]*annotation/i){$correctedName="snpEffAnnotation"} # Correction for snpEffAnnotation

		#FOR processRadtags.pm
        case ($name =~ m/process[\s|\.|\-| \/|\\|\|]*Radtags/i){$correctedName="processRadtags"} # Correction for processRadtags

        #FOR cutadapt functions
        case ($name =~ m/^cutadapt/i){$correctedName="cutadapt"} # Correction for cutadapt step

        #FOR atropos functions
        case ($name =~ m/^atropos/i){$correctedName="atropos"} # Correction for atropos step

        #FOR TGICL
        case ($name =~ m/^tgicl/i){$correctedName="tgicl"}

        #FOR trinity
        case ($name =~ m/^trinity/i){$correctedName="trinity"}  # Correction for Trinity step

        case ($name =~ m/^bamutils[\s|\.|\-| \/|\\|\|]*/i){$name =~ s/bamutils[\s|\.|\-| \/|\\|\|]*//gi; $correctedName="bamutils".$name;}  # Correction for bamutils from ngsutils tools

		#FOR checkFormat
		case ($name =~ m/^check[\s|\.|\-| \/|\\|\|]*format[\s|\.|\-| \/|\\|\|]*fasta/i){$correctedName="checkFormatFasta"}  # Correction for checkFormatFasta step
		case ($name =~ m/^check[\s|\.|\-| \/|\\|\|]*format[\s|\.|\-| \/|\\|\|]*fastq/i){$correctedName="checkFormatFastq"}  # Correction for checkFormatFastq step
		case ($name =~ m/^check[\s|\.|\-| \/|\\|\|]*format[\s|\.|\-| \/|\\|\|]*vcf/i){$correctedName="checkFormatVcf"}  # Correction for checkFormatVcf step
		case ($name =~ m/^check[\s|\.|\-| \/|\\|\|]*format[\s|\.|\-| \/|\\|\|]*sam[\s|\.|\-| \/|\\|\|]*or[\s|\.|\-| \/|\\|\|]*bam/i){$correctedName="checkFormatSamOrBam"}  # Correction for checkSamOrBam step
        case ($name =~ m/^check[\s|\.|\-| \/|\\|\|]*format[\s|\.|\-| \/|\\|\|]*sam/i){$correctedName="checkFormatSamOrBam"}  # Correction for checkSamOrBam step
        case ($name =~ m/^check[\s|\.|\-| \/|\\|\|]*format[\s|\.|\-| \/|\\|\|]*bam/i){$correctedName="checkFormatSamOrBam"}  # Correction for checkSamOrBam step
        case ($name =~ m/^check[\s|\.|\-| \/|\\|\|]*format[\s|\.|\-| \/|\\|\|]*gff/i){$correctedName="checkFormatGff"}  # Correction for checkFormatGff step
		case ($name =~ m/^check[\s|\.|\-| \/|\\|\|]*format[\s|\.|\-| \/|\\|\|]*bed/i){$correctedName="checkFormatBed"}  # Correction for checkFormatBed step

        #FOR BOWTIE SUITE
        case ($name =~ m/^bowtie$/i){$correctedName="bowtie"}#Correction for bowtie
        case ($name =~ m/^bowtie2$/i){$correctedName="bowtie2"}#Correction for bowtie2
        case ($name =~ m/^bowtie[\s|\.|\-| \/|\\|\|]*build/i){$correctedName="bowtieBuild"; } #Correction for bowtiebuild
        case ($name =~ m/^bowtie2[\s|\.|\-| \/|\\|\|]*build/i){$correctedName="bowtie2Build"; } #Correction for bowtie2build

        #FOR CRAC SUITE
        case ($name =~ m/^crac[\s|\.|\-| \/|\\|\|]*index/i){$correctedName="cracIndex"; } #Correction for cracIndex
        case ($name =~ m/^crac/i){$correctedName="crac"; } #Correction for crac

        #FOR DUPLICATIONDETECTOR
        case ($name =~ m/^duplication[\s|\.|\-| \/|\\|\|]*detector/i){$correctedName="duplicationDetector";} #Correction for duplicationDetector


        #FOR BEDTOOLS
        case ($name =~ m/^bed[\s|\.|\-| \/|\\|\|]*tools[\s|\.|\-| \/|\\|\|]*intersect/i){$correctedName="bedToolsIntersectBed";} #Correction for intersectBed version 1
        case ($name =~ m/^intersect[\s|\.|\-| \/|\\|\|]*bed/i){$correctedName="bedToolsIntersectBed";} #Correction for intersectBed version 2
        case ($name =~ m/^bed[\s|\.|\-| \/|\\|\|]*tools[\s|\.|\-| \/|\\|\|]*window/i){$correctedName="bedToolsWindowBed";} #Correction for windowBed version 1
        case ($name =~ m/^window[\s|\.|\-| \/|\\|\|]*bed/i){$correctedName="bedToolsWindowBed";} #Correction for windowBed version 2
        case ($name =~ m/^bed[\s|\.|\-| \/|\\|\|]*tools[\s|\.|\-| \/|\\|\|]*generic/i){$correctedName="bedToolsGeneric";} #Correction for bedtools generic

        #FOR GENERIC COMMAND SYSTEM
        case ($name =~ m/^generic/i){$correctedName="generic";} #Correction for generic command

        #FOR abyss functions
        #case ($name =~ m/^trans[\s|\.|\-| \/|\\|\|]*abyss/i){$correctedName="transAbyss"} # Correction for transAbyss step
        case ($name =~ m/^abyss/i){$correctedName="abyssSimple";} # Correction for abyss step

        #FOR pindel
        case ($name =~ m/^pindel/i){$correctedName="pindel";}

        #FOR breakDancer
        case ($name =~ m/^break[\s|\.|\-| \/|\\|\|]*dancer/i){$correctedName="breakDancer";}
        case ($name =~ m/^bam[\s|\.|\-| \/|\\|\|]*2[\s|\.|\-| \/|\\|\|]*cfg/i){$correctedName="bam2cfg";}
        case ($name =~ m/^bam[\s|\.|\-| \/|\\|\|]*2[\s|\.|\-| \/|\\|\|]*cfg[\s|\.|\-| \/|\\|\|]*pl/i){$correctedName="bam2cfg";}

        #FOR pindel
        case ($name =~ m/^pindel/i){$correctedName="pindel";}

        # SNIPLAY
        case ($name =~ m/^plink[\s|\.|\-| \/|\\|\|]*vcf[\s|\.|\-| \/|\\|\|]*2[\s|\.|\-| \/|\\|\|]*ped/i){$correctedName="plinkVcf2Ped";}
        case ($name =~ m/^sniplay[\s|\.|\-| \/|\\|\|]*ped[\s|\.|\-| \/|\\|\|]*2[\s|\.|\-| \/|\\|\|]*fasta/i){$correctedName="sniplayPed2fasta";}
        case ($name =~ m/^fastme/i){$correctedName="fastme";}
		case ($name =~ m/^readseq/i){$correctedName="readseq";}

        #ea-Utils
        case ($name =~ m/^fastq[\s|\.|\-| \/|\\|\|]*stats/i){$correctedName="fastqStats";}

        else
        {
            toolbox::exportLog("ERROR NAMING CONVENTION : $0 : the $name function or software is unknown to TOGGLE, cannot continue",0);
        }; # Name unknown to TOGGLE, must stop
    }
    $correctedName .= " ".$order if ($order);
    ##DEBUG toolbox::exportLog("$correctedName\n",1);
    return $correctedName;
}



#######################################################
## HASH for tools
#######################################################
sub returnSoftInfos
{
	my %softInfos = (
        
    #INFOS FOR NEW TOOLS

	'nanoplot'=>{		'IN' => 'fastq',
				'OUT'=>'NA',
				'cmdVersion' => "$nanoplot -v"},
    
	'stringtie'=>{		'IN' => 'bam,gtf',
				'OUT'=>'gtf',
				'cmdVersion' => "$stringtie -v 2>&1 | grep 'StringTie v' "},
    
	'hisat2'=>{		'IN' => 'fastq',
				'OUT'=>'sam',
				'cmdVersion' => "$hisat2 -v 2>&1 | grep 'HISAT2 version'"},

	'abyss' =>{		'IN' => 'fasta,fastq,sam,bam',
				'OUT' => 'fasta',
				'cmdVersion' => "$abyss --version | grep 'GNU Make' " },
	
	'atropos' =>{		'IN' => 'fastq',
				'OUT' => 'fastq',
				'cmdVersion' => "$atropos  2>&1 | grep 'Atropos version' " },

	'bam2cfg' =>{		'IN' => 'bam',
				'OUT' => 'NA',
				'cmdVersion' => "$breakDancer 2>&1 | grep Version' " },

	'bamutilsTool' =>{	'IN' => 'sam,bam',
				'OUT' => 'bam,fastq,fasta,bed',
				'cmdVersion' => "$bamutils | tail -1" },
	
	'bedTools' =>{		'cmdVersion' => "$bedtools --version " },
	
	'bedToolsGeneric' =>{	'IN' => 'bam,bed,vcf,gff',
				'OUT' => 'bed'},

	'bedToolsIntersect' =>{	'IN' => 'bam,bed,vcf,gff',
				'OUT' => 'bed'},

	'bedToolsWindow' =>{	'IN' => 'bam,bed,vcf,gff',
				'OUT' => 'bed' },

	'bowtieBuild' =>{	'IN' => 'fasta',
				'OUT' => 'NA',
				'MANDATORY' => 'reference',
				'cmdVersion' => "$bowtieBuild --version 2>&1 | grep 'bowtie-build version'" },

	'bowtie2Build' =>{	'IN' => 'fasta',
				'OUT' => 'NA',
				'MANDATORY' => 'reference',
				'cmdVersion' =>  "$bowtie2Build --version 2>&1 | grep 'bowtie2-build version'"},

	'bowtie' =>{		'IN' => 'fastq',
				'OUT' => 'sam',
				'MANDATORY' => 'reference',
				'cmdVersion' => "$bowtie --version | grep 'bowtie version' " },

	'bowtie2' =>{		'IN' => 'fastq',
				'OUT' => 'sam',
				'MANDATORY' => 'reference',
				'cmdVersion' => "$bowtie2 --version | grep 'bowtie2-align-s version'" },

	'breakDancer' =>{	'IN' => 'bam',
				'OUT' => 'NA',
				'cmdVersion' => "$breakDancer 2>&1 | grep 'Version'" },

	'bwa' =>{ 		'cmdVersion' => "$bwa 2>&1 | grep 'Version'" },

	'bwaAln' =>{		'IN' => 'fastq',
				'OUT' => 'sai',
				'MANDATORY' => 'reference'},

	'bwaIndex' =>{		'IN' => 'fasta',
				'OUT' => 'NA',
				'MANDATORY' => 'reference'},

	'bwaMem' =>{		'IN' => 'fastq,fasta',
				'OUT' => 'sam,bam',
				'MANDATORY' => 'reference'},

	'bwaSampe' =>{		'IN' => 'fastq,sai',
				'OUT' => 'sam',
				'MANDATORY' => 'reference'},

	'bwaSamse' =>{		'IN' => 'fastq,sai',
				'OUT' => 'sam',
				'MANDATORY' => 'reference'},

	'bwaSw' =>{		'IN' => 'fastq,fasta',
				'OUT' => 'sam',
				'MANDATORY' => 'reference'},

	'checkEncodeByASCIIcontrol' =>{'IN' => 'fastq',
						'OUT' => 'NA',
						'cmdVersion' => "echo 'v1.0'" },

	'checkFormatBed' =>{	'IN' => 'bed',
				'OUT' => 'NA'},

	'checkFormatFasta' =>{	'IN' => 'fasta',
				'OUT' => 'NA'},

	'checkFormatFastq' =>{	'IN' => 'fastq',
				'OUT' => 'NA'},


	'checkFormatGff' =>{	'IN' => 'gff,gtf',
				'OUT' => 'NA'},

	'checkFormatVcf' =>{	'IN' => 'vcf',
				'OUT' => 'NA'},

	'checkFormatSamOrBam' =>{'IN' => 'sam,bam',
				'OUT' => 'NA'},

	'checkFormatVcf' =>{	'IN' => 'vcf',
				'OUT' => 'NA'},

	'crac' =>{		'IN' => 'fastq',
				'OUT' => 'sam',
				'MANDATORY' => 'reference',
				'cmdVersion' => "$crac -version | grep -m 1 'version'" },

	'cracIndex' =>{		'IN' => 'fasta',
				'OUT' => 'NA',
				'MANDATORY' => 'reference'},

	'cutadapt' =>{		'IN' => 'fastq',
				'OUT' => 'fastq',
				'cmdVersion' => "$cutadapt 2>&1 | grep 'cutadapt version'" },

	'duplicationDetector' =>{'IN' => 'vcf',
				'OUT' => 'csv,bed',
				'cmdVersion' => "echo 'v1.0'"  },

	'fastme' =>{		'IN' => 'phy,phylip',
				'OUT' => 'nwk,newik,nk',
				'cmdVersion' => "$fastme --version"  },

	'fastqc' =>{		'IN' => 'fastq,bam',
				'OUT' => 'NA',
				'cmdVersion' => "$fastqc -v"  },

	'fastxTrimmer' =>{	'IN' => 'fastq',
				'OUT' => 'fastq',
				'cmdVersion' => "$fastxTrimmer -h | grep 'FASTX Toolkit'"  },

	'gatk' =>{		'cmdVersion' => "$GATK -version"  },

	'gatkBaseRecalibrator' =>{'IN' => 'bam',
				'OUT' => 'NA',
				'MANDATORY' => 'reference'},

	'gatkHaplotypeCaller' =>{'IN' => 'bam,vcf',
				'OUT' => 'vcf',	
				'MANDATORY' => 'reference'},

	'gatkIndelRealigner' =>{'IN' => 'intervals,bam,fasta',
				'OUT' => 'bam',
				'MANDATORY' => 'reference'},

	'gatkPrintReads' =>{	'IN' => 'bam',
				'OUT' => 'bam',
				'MANDATORY' => 'reference'},

	'gatkRealignerTargetCreator' =>{'IN' => 'bam,fasta',
					'OUT' => 'intervals',
					'MANDATORY' => 'reference'},

	'gatkSelectVariants' =>{'IN' => 'vcf',
				'OUT' => 'vcf',
				'MANDATORY' => 'reference'},

	'gatkUnifiedGenotyper' =>{'IN' => 'bam,vcf',
				'OUT' => 'vcf',
				'MANDATORY' => 'reference'},

	'gatkVariantFiltration' =>{	'IN' => 'vcf',
					'OUT' => 'vcf',
					'MANDATORY' => 'reference'},

	##################
	## ADD * for all input/output
	'generic' =>{		'IN' => '*',
				'OUT' => '*',
				'cmdVersion' => "echo 'v1.0'"  },
	##################

	'htseqCount' =>{	'IN' => 'sam,bam',
				'OUT' => 'NA',
				'MANDATORY' => 'gff',
				'cmdVersion' => "$htseqcount -h | grep 'version' | cut -d',' -f 2,2" },

	'java' => {		'cmdVersion' => "$java -version 2>&1 | grep 'version'" },

	'picardTools' =>{	'cmdVersion' => "$picard CheckFingerprint --version 2>&1" },

	'picardToolsAddOrReplaceReadGroups' =>{	'IN' => 'sam,bam',
						'OUT' => 'sam,bam'},

	'picardToolsCleanSam' =>{		'IN' => 'sam,bam',
						'OUT' => 'sam,bam'},

	'picardToolsCreateSequenceDictionary' =>{	'IN' => 'fasta',
							'OUT' => 'NA',
							'MANDATORY' => 'reference'},

	'picardToolsMarkDuplicates' =>{		'IN' => 'bam',
						'OUT' => 'bam'},

	'picardToolsSamFormatConverter' =>{	'IN' => 'sam,bam',
						'OUT' => 'sam,bam'},

	'picardToolsSortSam' =>{		'IN' => 'sam,bam',
						'OUT' => 'sam,bam'},

	'picardToolsValidateSamFile' =>{	'IN' => 'sam,bam',
						'OUT' => 'NA'},

	'pindel' =>{				'IN' => 'bam',
						'OUT' => 'NA',
						'cmdVersion' => "$pindel | grep -m 1 version" },

	'plinkVcf2Ped' =>{			'IN' => 'vcf',
						'OUT' => 'ped',
						'cmdVersion' => "$plink -version" },

	'stacks' =>{			'IN' => 'fastq',
						'OUT' => 'fastq',
						'MANDATORY' => 'keyfile',
						'cmdVersion' => "$stacks -v 2>&1" },

	'readseq' =>{				'IN' => 'fasta',
						'OUT' => 'phylip,phy',
						'cmdVersion' => "$readseqjar -h | head -1" },
						
	'samtools' =>{				'cmdVersion' => "$samtools 2>&1 | grep 'Version'" },
	
	'samToolsDepth' =>{			'IN' => 'bam',
						'OUT' => 'NA'},

	'samToolsFaidx' =>{			'IN' => 'fasta',
						'OUT' => 'NA',
						'MANDATORY' => 'reference'},

	'samToolsFlagstat' =>{			'IN' => 'bam',
						'OUT' => 'NA'},

	'samToolsIdxstats' =>{			'IN' => 'bam',
						'OUT' => 'NA'},

	'samToolsIndex' =>{			'IN' => 'bam',
						'OUT' => 'NA'},

	'samToolsMerge' =>{			'IN' => 'bam',
						'OUT' => 'bam'},

	'samToolsMpileUp' =>{			'IN' => 'bam',
						'OUT' => 'mpileup'},

	'samToolsSort' =>{			'IN' => 'bam',
						'OUT' => 'bam'},

	'samToolsView' =>{			'IN' => 'sam,bam',
						'OUT' => 'sam,bam'},

	'sniplayPed2fasta' =>{			'IN' => 'ped',
						'OUT' => 'fasta',
						'cmdVersion' => "echo 'v3'" },

	'snpEffAnnotation' =>{			'IN' => 'vcf',
						'OUT' => 'vcf',
						'MANDATORY' => 'vcf',
						'cmdVersion' => "$snpEff -version 2>&1" },

	'tgicl' =>{				'IN' => 'fasta',
						'OUT' => 'fasta',
						'cmdVersion' => "echo 'No version available" },

	'tophat2' =>{				'IN' => 'fastq',
						'OUT' => 'bam',
						'MANDATORY' => 'reference',
						'cmdVersion' => "$tophat2 -v" },

	'trinity' =>{				'IN' => 'fastq',
						'OUT' => 'fasta',
						'cmdVersion' => "$trinity --version | grep 'Trinity version'" },

	'fastqStats' =>{			'IN' => 'fastq',
						'OUT' => 'NA',
						'cmdVersion' => "$fastqStats -h | grep 'Version'" }

	);
# 	print Dumper( \%softInfos );
	return \%softInfos;
}


sub writeLogVersion
{
	my ($hashOrder, $toggleVersion, $reportDir,$report) = @_; #recovery $report boolean value: set to 1 if report is requested. $reportDir is the path were software.txt is generated
	$report=0 if not defined $report; # by default $report does not generate sofware.txt

	my %softInfos = %{returnSoftInfos()};

	my %softPathVersion = ();
	my %softPath = ();

	foreach my $softOrder ( values %{ $hashOrder } )
	{
		#DEBUG:print $softOrder." DANS LA BOUCLE\n";

		switch (1)
		{

            #LOG INFOS FOR NEW TOOLS
            
            #FOR stringtie
            case ($softOrder =~ m/^stringtie.*/i){$softPathVersion{"stringtie"}= `$softInfos{$softOrder}{"cmdVersion"}` if not defined $softPathVersion{"stringtie"};
            		                        $softPath{"stringtie"}= $stringtie if not defined $softPath{"stringtie"};
											}

			#FOR hisat2.pm
			case ($softOrder =~ m/^hisat2.*/i){$softPathVersion{"hisat2"}= `$softInfos{$softOrder}{"cmdVersion"}` if not defined $softPathVersion{"hisat2"};
			                               $softPath{"hisat2"}= $hisat2 if not defined $softPath{"hisat2"};
											}
            
			#FOR bwa.pm
			case ($softOrder =~ m/^bwa.*/i){$softPathVersion{"bwa"}= `$softInfos{"bwa"}{"cmdVersion"}` if not defined $softPathVersion{"bwa"};
											$softPath{"bwa"}= $bwa if not defined $softPath{"bwa"};
											}

			#FOR samTools.pm
			case ($softOrder =~ m/^samtools.*/i){$softPathVersion{"samtools"}= `$softInfos{"samtools"}{"cmdVersion"}` if not defined $softPathVersion{"samtools"};
												 $softPath{"samtools"}= $samtools if not defined $softPath{"samtools"};
												}

			#FOR picardTools.pm
			case ($softOrder =~ m/^picard.*/i){ $softPathVersion{"java"}= `$softInfos{"java"}{"cmdVersion"}` if not defined $softPathVersion{"java"};
												$softPath{"java"}= $java if not defined $softPath{"java"};
												$softPathVersion{"picard"}= `$softInfos{"picardTools"}{"cmdVersion"}` if not defined $softPathVersion{"picard"};
												$softPath{"picard"}= $picard if not defined $softPath{"picard"};
											   }

			#FOR gatk.pm
			case ($softOrder =~ m/^gatk.*/i){$softPathVersion{"java"}= `$softInfos{"java"}{"cmdVersion"}` if not defined $softPathVersion{"java"};
											 $softPathVersion{"GATK"}= `$softInfos{"gatk"}{"cmdVersion"}` if not defined $softPathVersion{"GATK"};
											 $softPath{"java"}= $java if not defined $softPath{"java"};
											 $softPath{"GATK"}= $GATK if not defined $softPath{"GATK"};
											 }

			#FOR fastqc
			case ($softOrder =~ m/^fastqc/i){$softPathVersion{"fastqc"}= `$softInfos{"fastqc"}{"cmdVersion"}` if not defined $softPathVersion{"fastqc"};
											 $softPath{"fastqc"}= $fastqc if not defined $softPath{"fastqc"};
											 }

			#FOR fastxToolkit
			case ($softOrder =~ m/^fastx.*/i){$softPathVersion{"fastxTrimmer"}= `$softInfos{"fastxTrimmer"}{"cmdVersion"}` if not defined $softPathVersion{"fastxTrimmer"};
											  $softPath{"fastxTrimmer"}= $fastxTrimmer if not defined $softPath{"fastxTrimmer"};
											  }

			#FOR tophat.pm
			case ($softOrder =~ m/^bowtie2.*/i){$softPathVersion{"bowtie2Build"}= `$softInfos{"bowtie2Build"}{"cmdVersion"}` if not defined $softPathVersion{"bowtie2Build"};
												$softPath{"bowtie2Build"}= $bowtie2Build if not defined $softPath{"bowtie2Build"};
												$softPathVersion{"bowtie2"}= `$softInfos{"bowtie2"}{"cmdVersion"}` if not defined $softPathVersion{"bowtie2"};
												$softPath{"bowtie2"}= $bowtie2 if not defined $softPath{"bowtie2"};
												}
			case ($softOrder =~ m/^bowtie$/i){$softPathVersion{"bowtieBuild"}= `$softInfos{"bowtieBuild"}{"cmdVersion"}` if not defined $softPathVersion{"bowtieBuild"};
											 $softPath{"bowtieBuild"}= $bowtieBuild if not defined $softPath{"bowtieBuild"};
											 $softPathVersion{"bowtie"}= `$softInfos{"bowtie"}{"cmdVersion"}` if not defined $softPathVersion{"bowtie"};
											 $softPath{"bowtie"}= $bowtie if not defined $softPath{"bowtie"};
											 }

			case ($softOrder =~ m/^tophat.*/i){$softPathVersion{"tophat2"}= `$softInfos{"tophat2"}{"cmdVersion"}` if not defined $softPathVersion{"tophat2"};
											   $softPathVersion{"bowtieBuild"}= `$softInfos{"bowtieBuild"}{"cmdVersion"}` if not defined $softPathVersion{"bowtieBuild"};
											   $softPathVersion{"bowtie2Build"}= `$softInfos{"bowtie2Build"}{"cmdVersion"}` if not defined $softPathVersion{"bowtie2Build"};
											   $softPath{"tophat2"}= $tophat2 if not defined $softPath{"tophat2"};
											   $softPath{"bowtieBuild"}= $bowtieBuild if not defined $softPath{"bowtieBuild"};
											   $softPath{"bowtie2Build"}= $bowtie2Build if not defined $softPath{"bowtie2Build"};
											   }

			#FOR cufflinks.pm
			case ($softOrder =~ m/^cufflinks.*/i){$softPathVersion{"cufflinks"}= `$softInfos{'cufflinks'}{"cmdVersion"}` if not defined $softPathVersion{"cufflinks"};
												  $softPath{"cufflinks"}= $cufflinks if not defined $softPath{"cufflinks"};
												  }
			case ($softOrder =~ m/^cuffdiff.*/i){$softPathVersion{"cuffdiff"}= `$softInfos{'cufflinks'}{"cmdVersion"}` if not defined $softPathVersion{"cuffdiff"};
												 $softPath{"cuffdiff"}= $cufflinks if not defined $softPath{"cuffdiff"};
												 }
			case ($softOrder =~ m/^cuffmerge.*/i){$softPathVersion{"cuffmerge"}= `$softInfos{'cufflinks'}{"cmdVersion"}` if not defined $softPathVersion{"cuffmerge"};
												  $softPath{"cuffmerge"}= $cufflinks if not defined $softPath{"cuffmerge"};
												  }

			#FOR HTSeq
			case ($softOrder =~ m/^htseq.*/i){$softPathVersion{"htseqcount"}= `$softInfos{'htseqCount'}{"cmdVersion"}` if not defined $softPathVersion{"htseqcount"};
											  $softPath{"htseqcount"}= $htseqcount if not defined $softPath{"htseqcount"};
											  }

			#FOR snpEff.pm
			case ($softOrder =~ m/^snp.*/i){$softPathVersion{"java"}= `$softInfos{'java'}{"cmdVersion"}` if not defined $softPathVersion{"java"};
											$softPathVersion{"snpEff"}= `$softInfos{'snpEff'}{"cmdVersion"}` if not defined $softPathVersion{"snpEff"};
											$softPath{"java"}= $java if not defined $softPath{"java"};
											$softPath{"snpEff"}= $snpEff if not defined $softPath{"snpEff"};
											}

			#FOR processRadtags.pm
			case ($softOrder =~ m/process.*/i){$softPathVersion{"stacks"}= `$softInfos{'stacks'}{"cmdVersion"}` if not defined $softPathVersion{"stacks"};
											   $softPath{"stacks"}= $stacks if not defined $softPath{"stacks"};
											   }
            

			#FOR cutadapt functions
			case ($softOrder =~ m/^cutadapt/i){$softPathVersion{"cutadapt"}= `$softInfos{'cutadapt'}{"cmdVersion"}` if not defined $softPathVersion{"cutadapt"};
											   $softPath{"cutadapt"}= $cutadapt if not defined $softPath{"cutadapt"};
											  }
			#FOR atropos functions
			case ($softOrder =~ m/^atropos/i){$softPathVersion{"atropos"}= `$softInfos{'atropos'}{"cmdVersion"}` if not defined $softPathVersion{"atropos"};
											   $softPath{"atropos"}= $atropos if not defined $softPath{"atropos"};
											  }

			#FOR TGICL
			case ($softOrder =~ m/^tgicl/i){$softPathVersion{"tgicl"}= `$softInfos{'tgicl'}{"cmdVersion"}` if not defined $softPathVersion{"tgicl"};
											$softPath{"tgicl"}= $tgicl if not defined $softPath{"tgicl"};
											}

			#FOR trinity
			case ($softOrder =~ m/^trinity/i){$softPathVersion{"trinity"}= `$softInfos{'trinity'}{"cmdVersion"}` if not defined $softPathVersion{"trinity"};
											  $softPath{"trinity"}= $trinity if not defined $softPath{"trinity"};
											  }
			#FOR bamutils
			case ($softOrder =~ m/^bamutils.*/i){$softPathVersion{"bamutils"}= `$softInfos{'bamutils'}{"cmdVersion"}` if not defined $softPathVersion{"bamutils"};
											  $softPath{"bamutils"}= $bamutils if not defined $softPath{"bamutils"};
											  }
			#FOR plink
			case ($softOrder =~ m/^plink.*/i){$softPathVersion{"plink"}= `$softInfos{'plink'}{"cmdVersion"}` if not defined $softPathVersion{"plink"};
											  $softPath{"plink"}= $plink if not defined $softPath{"plink"};
											  }
			#FOR fastme
			case ($softOrder =~ m/^fastme.*/i){$softPathVersion{"fastme"}= `$softInfos{'fastme'}{"cmdVersion"}` if not defined $softPathVersion{"fastme"};
											  $softPath{"fastme"}= $fastme if not defined $softPath{"fastme"};
											  }

			#FOR readseq
			case ($softOrder =~ m/^readseq.*/i){$softPathVersion{"readseq"}= `$softInfos{'readseq'}{"cmdVersion"}` if not defined $softPathVersion{"readseq"};
											  $softPath{"readseq"}= $readseqjar if not defined $softPath{"readseq"};
											  }

			#FOR SNIPLAY
			case ($softOrder =~ m/^sniplay.*/i){$softPathVersion{"sniplay"}= "v1.0" if not defined $softPathVersion{"sniplay"};
												$softPath{"sniplay"}= "sniplay" if not defined $softPath{"sniplay"};
											   }

			#For format checking
			case($softOrder =~ m/^check/i){next;}

			#FOR crac.pm
			case ($softOrder =~ m/^crac.*/i){$softPathVersion{"crac"}= `$softInfos{'crac'}{"cmdVersion"}` if not defined $softPathVersion{"crac"};
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
			case ($softOrder =~ m/^bedtools/i){$softPathVersion{"bedtools"}= `$softInfos{'bedTools'}{"cmdVersion"}` if not defined $softPathVersion{"bedtools"};
												$softPath{"bedtools"}= $bedtools if not defined $softPath{"bedtools"};}
			case ($softOrder =~ m/.*bed$/i){$softPathVersion{"bedtools"}= `$softInfos{'bedTools'}{"cmdVersion"}` if not defined $softPathVersion{"bedtools"};
												$softPath{"bedtools"}= $bedtools if not defined $softPath{"bedtools"};}

			#FOR generic command
			case ($softOrder =~ m/^generic/i){$softPathVersion{"generic"}= "v0.1" if not defined $softPathVersion{"generic"};
												$softPath{"generic"}= "" if not defined $softPath{"generic"};}

			#FOR Abyss
			case ($softOrder =~ m/^abyss/i){$softPathVersion{"abyss"}= `$softInfos{'abyss'}{"cmdVersion"}` if not defined $softPathVersion{"abyss"};
												$softPath{"abyss"}= $abyss if not defined $softPath{"abyss"};}

			#FOR transAbyss
			#case ($softOrder =~ m/^transAbyss/i){$softPathVersion{"transAbyss"}= transAbyssVersion if not defined $softPathVersion{"transAbyss"};
												#$softPath{"transAbyss"}= $abyss if not defined $softPath{"transAbyss"};}

			#For breakDancer
			case ($softOrder =~ m/^bam2cfg/i){$softPathVersion{"breakDancer"}= `$softInfos{'breakDancer'}{"cmdVersion"}` if not defined $softPathVersion{"breakDancer"};
												$softPath{"bam2cfg"}= $bam2cfg if not defined $softPath{"bam2cfg"};}
			case ($softOrder =~ m/^breakDancer/i){$softPathVersion{"breakDancer"}= `$softInfos{'breakDancer'}{"cmdVersion"}` if not defined $softPathVersion{"breakDancer"};
												$softPath{"breakDancer"}= $breakDancer if not defined $softPath{"breakDancer"};}

			#For Pindel
			case ($softOrder =~ m/^pindel/i){$softPathVersion{"pindel"}= `$softInfos{'pindel'}{"cmdVersion"}` if not defined $softPathVersion{"pindel"};
												$softPath{"pindel"}= $pindel if not defined $softPath{"pindel"};}

			#For ea-Utils
			case ($softOrder =~ m/^fastqStats/i){$softPathVersion{"fastqStats"}=`$softInfos{'fastqStats'}{"cmdVersion"}` if not defined $softPathVersion{"fastqStats"};
												$softPath{"fastqStats"}=$fastqStats if not defined $softPath{"fastqStats"};}

			else
			{
				toolbox::exportLog("ERROR SOFTWARE Managenemnt: $0 : the $softOrder function or software is unknown to TOGGLE, cannot continue",0);
			}; # Name unknown to TOGGLE, must stop
		}
	}
	## DEBUG print Dumper(%softPathVersion);

	open (my $fhSoft, ">", "$reportDir/software.txt") if $report;

	toolbox::exportLog("TOGGLE : $toggle : $toggleVersion\n",1);
	print $fhSoft "TOGGLE : $toggle : $toggleVersion (*\@{\\cite{toggle}}\@*)\n" if $report;

	foreach my $soft ( sort( keys %softPath ) )
	{
		my $version = $softPathVersion{$soft};
		toolbox::exportLog(uc($soft)." : $softPath{$soft} : $version",1);
		print $fhSoft uc($soft) ." : $version (*\@{\\cite{$soft}}\@*)\n" if $report;
	}

	close $fhSoft if $report;
}



1;

