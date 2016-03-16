package namingConvention;

###################################################################################################################################
#
# Copyright 2014-2015 IRD-CIRAD-INRA-ADNid
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
# Version 3 written by Cecile Monat, Christine Tranchant, Cedric Farcy, Maryline Summo, Julie Orjuela-Bouniol, Sebastien Ravel, Gautier Sarah, and Francois Sabot
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
	case ($name =~ m/cleaner/i){$correctedName="cleaner";} #Correction for cleaner step
	
	#FOR compressor
	case ($name =~ m/compress/i){$correctedName="compress";} #correction for compressor step
	
	#FOR SGE
	case ($name =~ m/sge/i){$correctedName="sge";} #Correction for sge configuration
	
	#FOR SLURM
	case ($name =~ m/slurm/i){$correctedName="slurm";} #Correction for slurm configuration
	
	#FOR MPRUN
	case ($name =~ m/mprun/i){$correctedName="mprun";} #Correction for mprun configuration
	
        #FOR bwa.pm
        case ($name =~ m/bwa[\s|\.|\-| \/|\\|\|]*aln/i){$correctedName="bwaAln"; } #Correction for bwaAln
        case ($name =~ m/bwa[\s|\.|\-| \/|\\|\|]*sampe/i){$correctedName="bwaSampe"} # Correction for bwaSampe
        case ($name =~ m/bwa[\s|\.|\-| \/|\\|\|]*samse/i){$correctedName="bwaSamse"} # Correction for bwaSamse
        case ($name =~ m/bwa[\s|\.|\-| \/|\\|\|]*index/i){$correctedName="bwaIndex"} # Correction for bwaIndex
        case ($name =~ m/bwa[\s|\.|\-| \/|\\|\|]*mem/i){$correctedName="bwaMem"} # Correction for bwaMem

        
        #FOR samTools.pm
        case ($name =~ m/samtools[\s|\.|\-| \/|\\|\|]*faidx/i){$correctedName="samToolsFaidx"} # Correction for samToolsFaidx
        case ($name =~ m/samtools[\s|\.|\-| \/|\\|\|]*index/i){$correctedName="samToolsIndex"} # Correction for samToolsIndex
        case ($name =~ m/samtools[\s|\.|\-| \/|\\|\|]*view/i){$correctedName="samToolsView"} # Correction for samToolsView
        case ($name =~ m/samtools[\s|\.|\-| \/|\\|\|]*sort/i){$correctedName="samToolsSort"} # Correction for samToolsSort
        case ($name =~ m/merge[\s|\.|\-| \/|\\|\|]*header/i){$correctedName="mergeHeader"} # Correction for mergeHeader
        case ($name =~ m/samtools[\s|\.|\-| \/|\\|\|]*merge/i){$correctedName="samToolsMerge"} # Correction for samToolsMerge
        case ($name =~ m/samtools[\s|\.|\-| \/|\\|\|]*idxstats/i){$correctedName="samToolsIdxstats"} # Correction for samToolsIdxstats
        case ($name =~ m/samtools[\s|\.|\-| \/|\\|\|]*depth/i){$correctedName="samToolsDepth"} # Correction for samToolsDepth
        case ($name =~ m/samtools[\s|\.|\-| \/|\\|\|]*flagstat/i){$correctedName="samToolsFlagstat"} # Correction for samToolsFlagstat
        case ($name =~ m/samtools[\s|\.|\-| \/|\\|\|]*mpileup/i){$correctedName="samToolsMpileUp"} # Correction for samToolsMpileUp

        #FOR picardTools.pm
        case ($name =~ m/picardtools[\s|\.|\-| \/|\\|\|]*mark[\s|\.|\-| \/|\\|\|]*duplicates/i){$correctedName="picardToolsMarkDuplicates"} # Correction for picardToolsMarkDuplicates
        case ($name =~ m/picardtools[\s|\.|\-| \/|\\|\|]*create[\s|\.|\-| \/|\\|\|]*sequence[\s|\.|\-| \/|\\|\|]*dictionary/i){$correctedName="picardToolsCreateSequenceDictionary"} # Correction for picardToolsCreateSequenceDictionary
        case ($name =~ m/picardtools[\s|\.|\-| \/|\\|\|]*sort[\s|\.|\-| \/|\\|\|]*sam/i){$correctedName="picardToolsSortSam"} # Correction for picardToolsSortSam
	case ($name =~ m/picardtools[\s|\.|\-| \/|\\|\|]*validate[\s|\.|\-| \/|\\|\|]*sam[\s|\.|\-| \/|\\|\|]*file/i){$correctedName="picardToolsValidateSamFile"} # Correction for picardToolsValidateSamFile
	case ($name =~ m/picardtools[\s|\.|\-| \/|\\|\|]*clean[\s|\.|\-| \/|\\|\|]*sam/i){$correctedName="picardToolsCleanSam"} # Correction for picardToolsCleanSam
	case ($name =~ m/picardtools[\s|\.|\-| \/|\\|\|]*sam[\s|\.|\-| \/|\\|\|]*format[\s|\.|\-| \/|\\|\|]*converter/i){$correctedName="picardToolsSamFormatConverter"} # Correction for picardToolsSamFormatConverter
	case ($name =~ m/picardtools[\s|\.|\-| \/|\\|\|]*add[\s|\.|\-| \/|\\|\|]*or[\s|\.|\-| \/|\\|\|]*replace[\s|\.|\-| \/|\\|\|]*group/i){$correctedName="picardToolsAddOrReplaceGroup"} # Correction for picardToolsAddOrReplaceGroup

        
        #FOR gatk.pm
        case ($name =~ m/gatk[\s|\.|\-| \/|\\|\|]*base[\s|\.|\-| \/|\\|\|]*recalibrator/i){$correctedName="gatkBaseRecalibrator"} # Correction for gatkBaseRecalibrator
        case ($name =~ m/gatk[\s|\.|\-| \/|\\|\|]*print[\s|\.|\-| \/|\\|\|]*reads/i){$correctedName="gatkPrintReads"} # Correction for gatkPrintReads
        case ($name =~ m/gatk[\s|\.|\-| \/|\\|\|]*realigner[\s|\.|\-| \/|\\|\|]*target[\s|\.|\-| \/|\\|\|]*creator/i){$correctedName="gatkRealignerTargetCreator"} # Correction for gatkRealignerTargetCreator
        case ($name =~ m/gatk[\s|\.|\-| \/|\\|\|]*indel[\s|\.|\-| \/|\\|\|]*realigner/i){$correctedName="gatkIndelRealigner"} # Correction for gatkIndelRealigner
        case ($name =~ m/gatk[\s|\.|\-| \/|\\|\|]*haplotype[\s|\.|\-| \/|\\|\|]*caller/i){$correctedName="gatkHaplotypeCaller"} # Correction for gatkHaplotypeCaller
        case ($name =~ m/gatk[\s|\.|\-| \/|\\|\|]*select[\s|\.|\-| \/|\\|\|]*variants/i){$correctedName="gatkSelectVariants"} # Correction for gatkSelectVariants
        case ($name =~ m/gatk[\s|\.|\-| \/|\\|\|]*variant[\s|\.|\-| \/|\\|\|]*filtration/i){$correctedName="gatkVariantFiltration"} # Correction for gatkVariantFiltration
        case ($name =~ m/gatk[\s|\.|\-| \/|\\|\|]*unified[\s|\.|\-| \/|\\|\|]*genotyper/i){$correctedName="gatkUnifiedGenotyper"} # Correction for gatkUnifiedGenotyper
        case ($name =~ m/gatk[\s|\.|\-| \/|\\|\|]*read[\s|\.|\-| \/|\\|\|]*backed[\s|\.|\-| \/|\\|\|]*phasing/i){$correctedName="gatkReadBackedPhasing"} # Correction for gatkReadBackedPhasing
        
        #FOR fastqc
        case ($name =~ m/fastqc/i){$correctedName="fastqc"} # Correction for fastqc
        
        #FOR fastqUtils.pm
        case ($name =~ m/check[\s|\.|\-| \/|\\|\|]*encode[\s|\.|\-| \/|\\|\|]*by[\s|\.|\-| \/|\\|\|]*ascii[\s|\.|\-| \/|\\|\|]*control/i){$correctedName="checkEncodeByASCIIcontrol"} # Correction for checkEncodeByASCIIcontrol
        case ($name =~ m/change[\s|\.|\-| \/|\\|\|]*encode/i){$correctedName="changeEncode"} # Correction for changeEncode
        
        #FOR fastxToolkit
        case ($name =~ m/fastx[\s|\.|\-| \/|\\|\|]*trimmer/i){$correctedName="fastxTrimmer"} # Correction for fastxTrimmer

        #FOR tophat.pm
        case ($name =~ m/bowtie[\s|\.|\-| \/|\\|\|]*build/i){$correctedName="bowtieBuild"; } #Correction for bowtiebuild
	case ($name =~ m/bowtie2[\s|\.|\-| \/|\\|\|]*build/i){$correctedName="bowtie2Build"; } #Correction for bowtie2build
	case ($name =~ m/tophat[\s|\.|\-| \/|\\|\|]*2/i){$correctedName="tophat2"; } #Correction for tophat2
	
        #FOR cufflinks.pm
        
        #FOR HTSeq.pm
        case ($name =~ m/htseq[\s|\.|\-| \/|\\|\|]*count/i){$correctedName="htseqCount"; } #Correction for htseq-count
	
        #FOR snpeff.pm
        case ($name =~ m/snpeff[\s|\.|\-| \/|\\|\|]*annotation/i){$correctedName="snpeffAnnotation"} # Correction for snpeffAnnotation

        
        #FOR cutadapt functions
        case ($name =~ m/cutadapt/i){$correctedName="cutadapt"} # Correction for cutadapt step
     
        
        
        else {toolbox::exportLog("ERROR : $0 : the $name function or software is unknown to TOGGLE, cannot continue",0);}; # Name unknown to TOGGLE, must stop
    }
    $correctedName .= " ".$order if ($order);
    ##DEBUG toolbox::exportLog("$correctedName\n",1);
    return $correctedName;
}

1;

=head1 NAME

    Package I<namingConvention> 

=head1 SYNOPSIS

	use namingConvention;
    
	namingConvention::softwareNomenclature ($hashSoftware);
    
	namingConvention::correctName ($softwareName);

=head1 DESCRIPTION

    Package namingConvention will reformat the softwares/steps names in the software.config.txt to be able to manage them correctly in the $configInfo hash.

=head2 FUNCTIONS

=head3 namingConvention::softwareNomenclature

This module will check and modify the $configInfo hash.
It takes a single argument: the configInfo hash (from toolbox::readFileConf).

=head3 namingConvention::correctName

This module will correct the name in order to be usable in TOGGLE
It takes a single argument: the software name.



=head1 AUTHORS

Intellectual property belongs to IRD, CIRAD and South Green developpement plateform for all versions also for ADNid for v2 and v3 and INRA for v3
Written by Cecile Monat, Christine Tranchant, Cedric Farcy, Maryline Summo, Julie Orjuela-Bouniol, Sebastien Ravel, Gautier Sarah, and Francois Sabot

=head1 SEE ALSO

L<http://www.southgreen.fr/>

=cut
