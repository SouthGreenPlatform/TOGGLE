package fileConfigurator;

###################################################################################################################################
#
# Copyright 2014-2017 IRD-CIRAD-INRA-ADNid
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
use Data::Dumper;

use lib qw(.);
use localConfig;
use toolbox;

our %testParams=    (
        #BWA soft
        bwaAln => ["-n 5"],
        bwaSampe => ["-a 500"],
        bwaSamse => [""],
        bwaMem => [""],
        bwaIndex => [""],
        #TopHat soft
        tophat2 => ["-i=30","-I=20000","-a=8","-m=1","--no-coverage-search","-g=10","--bowtie-n","--library-type=fr-unstranded","--microexon-search"],
        #GATK soft
        gatkRealignerTargetCreator => [""],
        gatkIndelRealigner => [""],
        gatkHaplotypeCaller => ["-rf BadCigar"],
        gatkUnifiedGenotyper => ["-rf BadCigar"],
        gatkVariantFiltration => ["--filterName 'FILTER-DP' --filterExpression 'DP<10 || DP>600' --filterName 'LowQual' --filterExpression 'QUAL<30'"],
        gatkSelectVariants => ["-selectType=SNP"],
        gatkBaseRecalibrator => ["-knownSites=$toggle/data/expectedData/GATKVARIANTFILTRATION.vcf"],
        gatkReadBackedPhasing => [""],
        gatkPrintReads => [""],
        #PicardTools soft
        picardToolsSortSam => ["SORT_ORDER=coordinate","VALIDATION_STRINGENCY=SILENT","CREATE_INDEX=TRUE"],
        picardToolsValidateSamFile => [""],
        picardToolsAddOrReplaceReadGroups => ["ID=Test","LB=Irigin","PL=Illumina","SM=glaberrima","VALIDATION_STRINGENCY=SILENT","PU=unit1"],
        picardToolsMarkDuplicates => ["VALIDATION_STRINGENCY=SILENT","CREATE_INDEX=TRUE","REMOVE_DUPLICATES=TRUE"],
        picardToolsCreateSequenceDictionary => [""],
        picardToolsCleanSam => ["VALIDATION_STRINGENCY=SILENT"],
        picardToolsSamFormatCOnverter => ["VALIDATION_STRINGENCY=SILENT"],
        #SAMtools
        samToolsView => ["-h","-b","-f=0x02"],
        samToolsSort => [""],
        samToolsFaidx => ["TRINITY_DN75358_c0_g1_i1:1-100"],
        samToolsIndex => [""],
        mergeHeader => [""],
        samToolsMerge => [""],
        samToolsIdxStats => [""],
        samToolsDepth => [""],
        samToolsFlagstat => [""],
        samToolsMpileUp => [""],
        #Fastqc
        fastqc => [""],
        #Cutadapt
        cutadapt => ["-O=10","-m=35","-q=20,20","--overlap=7","-b GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG -B GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG","-b GTTCGTCTTCTGCCGTATGCTCTAGCACTACACTGACCTCAAGTCTGCACACGAGAAGGCTAG -B GTTCGTCTTCTGCCGTATGCTCTAGCACTACACTGACCTCAAGTCTGCACACGAGAAGGCTAG"],
        #FastXToolKit
        fastxTrimmer => ["-f 8", "-Q33"],
        #SnpEFF
        snpEffAnn => [""],
        #Stacks
        processRadtags => ["-i fastq","-e 'apeKI'","--retain_header"],
        #TGICL
        tgicl => ["-c 6","-p 90","-l 20"],
        #Trinity
        trinity => ["--seqType fq","--max_memory 20G","--full_cleanup"],
        #HTSeqCount
        htseqcount => ["-r=name","-s=no","-t=mRNA","-m=union","-i=ID"],
        #Pindel
        pindel => ["-c ALL","-x 4","-l","-k","-s","-M 10","-A 20","--MIN_DD_MAP_DISTANCE 2000"],
        #BreakDancer
        breakdancer => [""],
        #SGE
        sge => ["-q bioinfo.q","-b Y"],
        #checkFormat
        checkFormatFasta => [""],
        checkFormatFastq => [""],
        checkFormatVcf => [""],
        checkFormatSamOrBam => [""],
        #fastqUtils
        checkEncodeByASCIIcontrol => [""]
        );

sub softParams
{
    #This sub gathers all the soft parameters for tests and return the text to add
    die("ERROR: fileConfigurator::softParams : should provide only one argument\n",0) if (@_ > 1);
    my ($softName) = @_ ;


    ##DEBUG    print Dumper (\%testParams);

    if (exists $testParams{$softName})
    {
        ##DEBUG    print @{$testParams{$softName}},"\n";
        my $returnValue = join ("\n",@{$testParams{$softName}});
        return $returnValue ;
    }
    else
    {
        die("ERROR: fileConfigurator::softParams : unknown software $softName\n",0);
    }
}

sub createFileConf
{
    #This subprogram will generate a config file for soft tests.
    die("ERROR: fileConfigurator::createFileConf : should provide two arguments, order of softs and output file name\n",0) if (@_ != 2);
    my ($listOrder,$outputName) = @_ ;

    #Order initialization
    my $i = 1;
    open (my $fhOut, ">", $outputName) or die("ERROR: fileConfigurator::createFileConf : Cannot create output file $outputName\n: $!\n",0);


    my $orderConfig = "\n\n\$order\n";

    while (@{$listOrder})
    {
        my $currentSoft = shift @{$listOrder};
        ##DEBUG print $currentSoft,"\n";
        if ($currentSoft eq "1000")
        {
                $i = 1000;
                next;
        }
        $orderConfig .= $i."=".$currentSoft."\n";
        my $currentParams = softParams($currentSoft);
        my $lineOut = "\n\$".$currentSoft."\n".$currentParams."\n";
        print $fhOut $lineOut;
        $i++; #increasing the order level for next turn
    }
    print $fhOut $orderConfig;
    close $fhOut;
    return 1;
}

sub generateAllConf
{
    #Will generate a file conf gathering ALL the configs...
}

=head1 NOM

package I<fileConfigurator>

=head1 SYNOPSIS

	use fileConfigurator;

	fileConfigurator::softParams($softName);

	fileConfigurator::createFileConf($listOrder,$outputName);

	fileConfigurator::generateAllConf($outputName);


=head1 DESCRIPTION

This module is a set of internal functions related to tests

1;
