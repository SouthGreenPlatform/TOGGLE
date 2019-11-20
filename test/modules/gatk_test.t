#!/usr/bin/perl

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

######################################################################################################################################
######################################################################################################################################
## COMMON MODULE TEST HEADER
######################################################################################################################################
######################################################################################################################################

use strict;
use warnings;
use Data::Dumper;

use Test::More 'no_plan'; #Number of tests, to modify if new tests implemented. Can be changed as 'no_plan' instead of tests=>11 .
use Test::Deep;

# Load localConfig if primary test is successful
use_ok('localConfig') or exit;
use localConfig;


########################################
# Extract automatically tool name and sub name list
########################################
my ($toolName,$tmp) = split /_/ , $0;
my $subFile=$toggle."/modules/".$toolName.".pm";
my @sub = `grep "^sub" $subFile`or die ("ERROR: $0 : Cannot extract automatically sub name list by grep command \n$!\n");


########################################
#Automatically module test with use_ok and can_ok
########################################

use_ok($toolName) or exit;
eval "use $toolName";

foreach my $subName (@sub)
{
    chomp ($subName);
    $subName =~ s/sub //;
    can_ok($toolName,$subName);
}

#########################################
#Preparing test directory
#########################################
my $testDir="$toggle/dataTest/$toolName"."TestModule";
my $cmd="rm -Rf $testDir ; mkdir -p $testDir";
system($cmd) and die ("ERROR: $0 : Cannot execute the test directory $testDir ($toolName) with the following cmd $cmd\n$!\n");
chdir $testDir or die ("ERROR: $0 : Cannot go into the test directory $testDir ($toolName) with the chdir cmd \n$!\n");


#########################################
#Creating log file
#########################################
my $logFile=$toolName."_log.o";
my $errorFile=$toolName."_log.e";
system("touch $testDir/$logFile $testDir/$errorFile") and die "\nERROR: $0 : cannot create the log files $logFile and $errorFile: $!\nExiting...\n";

######################################################################################################################################
######################################################################################################################################

##########################################
##### gatk::gatkRealignerTargetCreator
##########################################
my $bankData="$toggle/data/Bank/";
my $bamData="$toggle/data/testData/samBam/";
my $vcfData="$toggle/data/testData/vcf/singleVCF/";

# input file
my $bamIn="$bamData/oneBam/RC3-SAMTOOLSVIEW.bam";
# index BAM
system ("samtools index $bamIn") and die ("ERROR: $0 : Cannot index file $bamIn\n$!\n");

my $fastaRef="$bankData/referenceIrigin.fasta";
my $fastaRefFai="$bankData/referenceIrigin.fasta.fai";
my $fastaRefDict="$bankData/referenceIrigin.dict";

# output file
my $intervalsFile="RC3.GATKREALIGNERTARGETCREATOR.intervals";

# execution test
is(gatk::gatkRealignerTargetCreator($fastaRef, $bamIn, $intervalsFile),1, 'gatk::gatkRealignerTargetCreator');

# expected output test
my $observedOutput = `ls`;
my @observedOutput = split /\n/,$observedOutput;
my @expectedOutput = ('gatk_log.e','gatk_log.o','RC3.GATKREALIGNERTARGETCREATOR.intervals');

is_deeply(\@observedOutput,\@expectedOutput,'gatk::gatkRealignerTargetCreator - output list');

# expected content test
my $expectedMD5sum="ad6d03974c93118e47e2149c7a5f916e";      # structure of the ref file
my $observedMD5sum=`md5sum $intervalsFile`;       # structure of the test file
my @withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'gatk::gatkRealignerTargetCreator - output content');


##########################################
#### Test for gatkIndelRealigner
##########################################

# output File
my $bamOut="RC3.GATKINDELREALIGNER.bam";

# execution test
is(gatk::gatkIndelRealigner($fastaRef, $bamIn, $intervalsFile, $bamOut),1, 'gatk::gatkIndelRealigner');

# expected output test
$observedOutput = `ls`;
@observedOutput = split /\n/,$observedOutput;
@expectedOutput = ('gatk_log.e','gatk_log.o','RC3.GATKINDELREALIGNER.bai','RC3.GATKINDELREALIGNER.bam','RC3.GATKREALIGNERTARGETCREATOR.intervals');

is_deeply(\@observedOutput,\@expectedOutput,'gatk::gatkIndelRealigner - output list');


# expected output content
$observedOutput = `samtools view $bamOut | wc -l `;
chomp $observedOutput;
my $expectedOutput = 1998;
is($observedOutput,$expectedOutput, 'gatk::gatkIndelRealigner - output content 1');


$observedOutput = `samtools view $bamOut |grep "H2:C381HACXX:5:1101:1433:2214"`;
chomp $observedOutput;
$expectedOutput = "H2:C381HACXX:5:1101:1433:2214	99	2187676	10	29	101M	=	199	269	AGCCCGAAGACCCGCAGTGCGAGGATTTCGAGGATCAAGCTCAAGATCTCGAGCAAAGCAAGTCACCTTTGATCATCTTGCACCTATAATTTAAATCTAAG	CCCFFFFFHHHHHJJJJGIJJIIJIJJJJJIJIIJJGJIJJJJJIHGHHHFFFFDEEEDDDDCDEDDDDDDDDDECDDDDDDDDDDDDDEEECDCDDDDDD	X0:i:1	X1:i:0	MC:Z:6S80M15S	MD:Z:38T0G16G44	RG:Z:RC3	XG:i:0	AM:i:29	NM:i:3	SM:i:29	XM:i:3	XO:i:0	MQ:i:29	XT:A:U
H2:C381HACXX:5:1101:1433:2214	147	2187676	199	29	6S80M15S	=	10	-269	ATCCTAAATTGCTGCAAATACCCTCCGTGAATTATTGAACACTTAAACCTCCTTTGTCGACCGTTGTGCTTCGATGCACGGGCCTTCGGACACGCGCATCA	DCDDDDDDDDDDEDEEDC>2BDDDDDDDEDEDEEEEDDCCCDDDCA<DDBCCDDDDFFHCHIJJJJJIJJJJJJJJIHJJJIJJJJJJHHHHHFFFFFCCC	MC:Z:101M	MD:Z:11T8T0G21C36	RG:Z:RC3	XG:i:0	AM:i:29	NM:i:4	SM:i:29	XM:i:4	XO:i:0	MQ:i:29	XT:A:M";
is($observedOutput,$expectedOutput, 'gatk::gatkIndelRealigner - output content 2');

# RM index bai
system ("rm $bamIn.bai") and die ("ERROR: $0 : Cannot remove index file $bamIn.bai\n$!\n");


##########################################
#####Test for gatk Unified Genotyper
##########################################

undef($bamIn);
push(@{$bamIn},$bamOut);
my $vcfOut="RC3.GATKUNIFIEDGENOTYPER.vcf";

# execution test
is(gatk::gatkUnifiedGenotyper($fastaRef,$bamIn, $vcfOut),1,'gatk::gatkUnifiedGenotyper');


# expected output test
$observedOutput = `ls`;
@observedOutput = split /\n/,$observedOutput;
@expectedOutput = ('gatk_log.e','gatk_log.o','RC3.GATKINDELREALIGNER.bai','RC3.GATKINDELREALIGNER.bam','RC3.GATKREALIGNERTARGETCREATOR.intervals','RC3.GATKUNIFIEDGENOTYPER.vcf','RC3.GATKUNIFIEDGENOTYPER.vcf.idx');

is_deeply(\@observedOutput,\@expectedOutput,'gatk::gatkUnifiedGenotyper - output list');

# expected output content
$observedOutput=`tail -n 1 $vcfOut`;
chomp $observedOutput;
$expectedOutput="2299572	2417	.	C	T	15.65	.	AC=2;AF=1.00;AN=2;DP=1;Dels=0.00;ExcessHet=3.0103;FS=0.000;HaplotypeScore=0.0000;MLEAC=2;MLEAF=1.00;MQ=60.00;MQ0=0;QD=15.65;SOR=1.609	GT:AD:DP:GQ:PL	1/1:0,1:1:3:42,3,0";

is($observedOutput,$expectedOutput,'gatk::gatkUnifiedGenotyper - output content');


##########################################
##### Test for gatkBaseRecalibrator
##########################################

# execution test
$bamIn=$bamOut;
my $controlVCF=$vcfData."/GATKVARIANTFILTRATION.vcf";
my $tableReport="recal_data.table";

my %optionsHachees = (
			"-knownSites" => $controlVCF,
			);        # Hash containing informations
my $optionsHachees = \%optionsHachees;
is(gatk::gatkBaseRecalibrator($fastaRef,$bamIn,$tableReport, $optionsHachees),1, 'gatk::gatkBaseRecalibrator');

# expected output test
$observedOutput = `ls`;
@observedOutput = split /\n/,$observedOutput;
@expectedOutput = ('gatk_log.e','gatk_log.o','RC3.GATKINDELREALIGNER.bai','RC3.GATKINDELREALIGNER.bam','RC3.GATKREALIGNERTARGETCREATOR.intervals','RC3.GATKUNIFIEDGENOTYPER.vcf','RC3.GATKUNIFIEDGENOTYPER.vcf.idx','recal_data.table');

is_deeply(\@observedOutput,\@expectedOutput,'gatk::gatkBaseRecalibrator - output list');

# expected output content
$expectedOutput="6273";     # nb of line from the ref file
$observedOutput=`wc -l $tableReport`;     # structure of the test file
@withoutName =split (" ", $observedOutput);     # to separate the structure and the name of the test file
$observedOutput = $withoutName[0];      #just to have the wc -l result
is($observedOutput,$expectedOutput,'gatk::gatkBaseRecalibrator - output content');


##########################################
##### Test for gatkPrintReads
##########################################

# execution test
$bamOut=$bamIn;
$bamOut =~ s/\.bam/\.RECALIBRATED\.bam/;

is(gatk::gatkPrintReads($fastaRef,$bamIn,$bamOut,$tableReport),1, 'gatk::gatkPrintReads');


# expected output test
$observedOutput = `ls`;
@observedOutput = split /\n/,$observedOutput;
@expectedOutput = ('gatk_log.e','gatk_log.o','RC3.GATKINDELREALIGNER.bai','RC3.GATKINDELREALIGNER.bam','RC3.GATKINDELREALIGNER.RECALIBRATED.bai','RC3.GATKINDELREALIGNER.RECALIBRATED.bam','RC3.GATKREALIGNERTARGETCREATOR.intervals','RC3.GATKUNIFIEDGENOTYPER.vcf','RC3.GATKUNIFIEDGENOTYPER.vcf.idx','recal_data.table');

is_deeply(\@observedOutput,\@expectedOutput,'gatk::gatkPrintReads - output list');


# expected output content
#For number of lines
$observedOutput = `samtools view $bamOut | wc -l `;
chomp $observedOutput;
$expectedOutput = 1998;
is($observedOutput,$expectedOutput,'gatk::gatkPrintReads - output content 1');

#For correction of quality
$observedOutput = `samtools view $bamOut |grep "H2:C381HACXX:5:1101:1433:2214"`;
chomp $observedOutput;
$expectedOutput = "H2:C381HACXX:5:1101:1433:2214	99	2187676	10	29	101M	=	199	269	AGCCCGAAGACCCGCAGTGCGAGGATTTCGAGGATCAAGCTCAAGATCTCGAGCAAAGCAAGTCACCTTTGATCATCTTGCACCTATAATTTAAATCTAAG	88888888888888888888888888888888888888888888888888888878887777878777777777887777777777777888878777777	X0:i:1	X1:i:0	MC:Z:6S80M15S	BD:Z:LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL	MD:Z:38T0G16G44	RG:Z:RC3	XG:i:0	BI:Z:LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL	AM:i:29	NM:i:3	SM:i:29	XM:i:3	XO:i:0	MQ:i:29	XT:A:U
H2:C381HACXX:5:1101:1433:2214	147	2187676	199	29	6S80M15S	=	10	-269	ATCCTAAATTGCTGCAAATACCCTCCGTGAATTATTGAACACTTAAACCTCCTTTGTCGACCGTTGTGCTTCGATGCACGGGCCTTCGGACACGCGCATCA	78777777777787887888777777778787888877888777888777887777888887888888888888888888888888888888888888888	MC:Z:101M	BD:Z:LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL	MD:Z:11T8T0G21C36	RG:Z:RC3	XG:i:0	BI:Z:LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL	AM:i:29	NM:i:4	SM:i:29	XM:i:4	XO:i:0	MQ:i:29	XT:A:M";
is($observedOutput,$expectedOutput, 'gatk::gatkPrintReads - output content 2');


#################################################################################################
######Test for gatk Haplotype caller
##########################################
$bamIn=$bamData."/twoBamsIrigin/irigin1-PICARDTOOLSMARKDUPLICATES.bam";
my $bamSingle=$bamData."/twoBamsIrigin/irigin3-PICARDTOOLSMARKDUPLICATES.bam";

# index BAM
system ("samtools index $bamIn") and die ("ERROR: $0 : Cannot index file $bamIn\n$!\n");
system ("samtools index $bamSingle") and die ("ERROR: $0 : Cannot index file $bamSingle\n$!\n");


# execution test
my @bamsToCall=($bamIn,$bamSingle);
my $vcfCalled="GATKHAPLOTYPECALLER.vcf";
is(gatk::gatkHaplotypeCaller($fastaRef, $vcfCalled, \@bamsToCall) ,1,"gatk::gatkHaplotypeCaller");

# expected output test
$observedOutput = `ls`;
@observedOutput = split /\n/,$observedOutput;
@expectedOutput = ('GATKHAPLOTYPECALLER.vcf','GATKHAPLOTYPECALLER.vcf.idx','gatk_log.e','gatk_log.o','RC3.GATKINDELREALIGNER.bai','RC3.GATKINDELREALIGNER.bam','RC3.GATKINDELREALIGNER.RECALIBRATED.bai','RC3.GATKINDELREALIGNER.RECALIBRATED.bam','RC3.GATKREALIGNERTARGETCREATOR.intervals','RC3.GATKUNIFIEDGENOTYPER.vcf','RC3.GATKUNIFIEDGENOTYPER.vcf.idx','recal_data.table');

is_deeply(\@observedOutput,\@expectedOutput,'gatk::gatkHaplotypeCaller - output list');

# expected output structure test
my $expectedSNPLines="2224477	996	.	TA	T	34.13	.	AC=2;AF=1.00;AN=2;DP=2;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=17.07;SOR=0.693	GT:AD:DP:GQ:PL	./.	1/1:0,2:2:6:69,6,0
2248321	377	.	C	G	64.17	.	AC=2;AF=1.00;AN=2;DP=2;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=29.00;QD=32.08;SOR=0.693	GT:AD:DP:GQ:PL	./.	1/1:0,2:2:6:90,6,0
2248321	379	.	C	T	64.17	.	AC=2;AF=1.00;AN=2;DP=2;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=29.00;QD=32.08;SOR=0.693	GT:AD:DP:GQ:PL	./.	1/1:0,2:2:6:90,6,0
2281178	4213	.	G	A	64.17	.	AC=2;AF=1.00;AN=2;DP=2;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=32.08;SOR=0.693	GT:AD:DP:GQ:PL	./.	1/1:0,2:2:6:90,6,0
2281178	4214	.	A	G	64.17	.	AC=2;AF=1.00;AN=2;DP=2;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=32.08;SOR=0.693	GT:AD:DP:GQ:PL	./.	1/1:0,2:2:6:90,6,0
2290182	1013	.	A	G	44.17	.	AC=2;AF=1.00;AN=2;DP=2;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=29.00;QD=22.09;SOR=0.693	GT:AD:DP:GQ:PL	./.	1/1:0,2:2:6:70,6,0
2291957	21946	.	A	G	20	.	AC=2;AF=1.00;AN=2;DP=1;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=29.00;QD=20.00;SOR=1.609	GT:AD:DP:GQ:PL	./.	1/1:0,1:1:3:45,3,0
2291957	21948	.	A	G	20	.	AC=2;AF=1.00;AN=2;DP=1;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=29.00;QD=20.00;SOR=1.609	GT:AD:DP:GQ:PL	./.	1/1:0,1:1:3:45,3,0
2291957	21957	.	A	C	20	.	AC=2;AF=1.00;AN=2;DP=1;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=29.00;QD=20.00;SOR=1.609	GT:AD:DP:GQ:PL	./.	1/1:0,1:1:3:45,3,0
2291957	21958	.	T	C	20	.	AC=2;AF=1.00;AN=2;DP=0;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=NaN;SOR=0.693	GT:AD:DP:GQ:PL	./.	1/1:0,0:0:3:45,3,0";     # structure of the ref file
my $observedSNPLines=`grep -v "#" $vcfCalled`;      # structure of the test file
chomp $observedSNPLines;
is($observedSNPLines,$expectedSNPLines,'gatk::gatkHaplotypeCaller - structure of file');



##########################################
######Test for gatkVariant filtrator
##########################################

#execution Test
my $variantFiltered="GATKVARIANTFILTRATION.vcf";
is(gatk::gatkVariantFiltration($fastaRef, $variantFiltered, $vcfCalled), 1 , 'gatk::gatkVariantFiltration');

# expected output test
$observedOutput = `ls`;
@observedOutput = split /\n/,$observedOutput;
@expectedOutput = ('GATKHAPLOTYPECALLER.vcf','GATKHAPLOTYPECALLER.vcf.idx','gatk_log.e','gatk_log.o','GATKVARIANTFILTRATION.vcf','GATKVARIANTFILTRATION.vcf.idx','RC3.GATKINDELREALIGNER.bai','RC3.GATKINDELREALIGNER.bam','RC3.GATKINDELREALIGNER.RECALIBRATED.bai','RC3.GATKINDELREALIGNER.RECALIBRATED.bam','RC3.GATKREALIGNERTARGETCREATOR.intervals','RC3.GATKUNIFIEDGENOTYPER.vcf','RC3.GATKUNIFIEDGENOTYPER.vcf.idx','recal_data.table');

is_deeply(\@observedOutput,\@expectedOutput,'gatk::gatkVariantFiltration - output list');

# expected output structure test
$expectedSNPLines="2224477	996	.	TA	T	34.13	PASS	AC=2;AF=1.00;AN=2;DP=2;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=17.07;SOR=0.693	GT:AD:DP:GQ:PL	./.	1/1:0,2:2:6:69,6,0
2248321	377	.	C	G	64.17	PASS	AC=2;AF=1.00;AN=2;DP=2;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=29.00;QD=32.08;SOR=0.693	GT:AD:DP:GQ:PL	./.	1/1:0,2:2:6:90,6,0
2248321	379	.	C	T	64.17	PASS	AC=2;AF=1.00;AN=2;DP=2;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=29.00;QD=32.08;SOR=0.693	GT:AD:DP:GQ:PL	./.	1/1:0,2:2:6:90,6,0
2281178	4213	.	G	A	64.17	PASS	AC=2;AF=1.00;AN=2;DP=2;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=32.08;SOR=0.693	GT:AD:DP:GQ:PL	./.	1/1:0,2:2:6:90,6,0
2281178	4214	.	A	G	64.17	PASS	AC=2;AF=1.00;AN=2;DP=2;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=32.08;SOR=0.693	GT:AD:DP:GQ:PL	./.	1/1:0,2:2:6:90,6,0
2290182	1013	.	A	G	44.17	PASS	AC=2;AF=1.00;AN=2;DP=2;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=29.00;QD=22.09;SOR=0.693	GT:AD:DP:GQ:PL	./.	1/1:0,2:2:6:70,6,0
2291957	21946	.	A	G	20	PASS	AC=2;AF=1.00;AN=2;DP=1;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=29.00;QD=20.00;SOR=1.609	GT:AD:DP:GQ:PL	./.	1/1:0,1:1:3:45,3,0
2291957	21948	.	A	G	20	PASS	AC=2;AF=1.00;AN=2;DP=1;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=29.00;QD=20.00;SOR=1.609	GT:AD:DP:GQ:PL	./.	1/1:0,1:1:3:45,3,0
2291957	21957	.	A	C	20	PASS	AC=2;AF=1.00;AN=2;DP=1;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=29.00;QD=20.00;SOR=1.609	GT:AD:DP:GQ:PL	./.	1/1:0,1:1:3:45,3,0
2291957	21958	.	T	C	20	PASS	AC=2;AF=1.00;AN=2;DP=0;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=NaN;SOR=0.693	GT:AD:DP:GQ:PL	./.	1/1:0,0:0:3:45,3,0
";     # structure of the ref file
$observedSNPLines=`grep -v "#" $variantFiltered`;      # structure of the test file
is($observedSNPLines,$expectedSNPLines,'gatk::gatkVariantFiltration - structure of file');


#################################################################################################
######Test for gatk Select Variant
#################################################################################################


#execution Test
my $vcfVariantsSelected="GATKSELECTVARIANT.vcf";
is(gatk::gatkSelectVariants($fastaRef, $variantFiltered, $vcfVariantsSelected),1, 'gatk::gatkSelectVariants');

# expected output test
$observedOutput = `ls`;
@observedOutput = split /\n/,$observedOutput;
@expectedOutput = ('GATKHAPLOTYPECALLER.vcf','GATKHAPLOTYPECALLER.vcf.idx','gatk_log.e','gatk_log.o','GATKSELECTVARIANT.vcf','GATKSELECTVARIANT.vcf.idx','GATKVARIANTFILTRATION.vcf','GATKVARIANTFILTRATION.vcf.idx','RC3.GATKINDELREALIGNER.bai','RC3.GATKINDELREALIGNER.bam','RC3.GATKINDELREALIGNER.RECALIBRATED.bai','RC3.GATKINDELREALIGNER.RECALIBRATED.bam','RC3.GATKREALIGNERTARGETCREATOR.intervals','RC3.GATKUNIFIEDGENOTYPER.vcf','RC3.GATKUNIFIEDGENOTYPER.vcf.idx','recal_data.table');

is_deeply(\@observedOutput,\@expectedOutput,'gatk::gatkSelectVariants - output list');

# expected output structure test
$expectedSNPLines="2224477	996	.	TA	T	34.13	PASS	AC=2;AF=1.00;AN=2;DP=2;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=17.07;SOR=0.693	GT:AD:DP:GQ:PL	./.	1/1:0,2:2:6:69,6,0
2248321	377	.	C	G	64.17	PASS	AC=2;AF=1.00;AN=2;DP=2;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=29.00;QD=32.08;SOR=0.693	GT:AD:DP:GQ:PL	./.	1/1:0,2:2:6:90,6,0
2248321	379	.	C	T	64.17	PASS	AC=2;AF=1.00;AN=2;DP=2;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=29.00;QD=32.08;SOR=0.693	GT:AD:DP:GQ:PL	./.	1/1:0,2:2:6:90,6,0
2281178	4213	.	G	A	64.17	PASS	AC=2;AF=1.00;AN=2;DP=2;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=32.08;SOR=0.693	GT:AD:DP:GQ:PL	./.	1/1:0,2:2:6:90,6,0
2281178	4214	.	A	G	64.17	PASS	AC=2;AF=1.00;AN=2;DP=2;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=32.08;SOR=0.693	GT:AD:DP:GQ:PL	./.	1/1:0,2:2:6:90,6,0
2290182	1013	.	A	G	44.17	PASS	AC=2;AF=1.00;AN=2;DP=2;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=29.00;QD=22.09;SOR=0.693	GT:AD:DP:GQ:PL	./.	1/1:0,2:2:6:70,6,0
2291957	21946	.	A	G	20	PASS	AC=2;AF=1.00;AN=2;DP=1;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=29.00;QD=20.00;SOR=1.609	GT:AD:DP:GQ:PL	./.	1/1:0,1:1:3:45,3,0
2291957	21948	.	A	G	20	PASS	AC=2;AF=1.00;AN=2;DP=1;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=29.00;QD=20.00;SOR=1.609	GT:AD:DP:GQ:PL	./.	1/1:0,1:1:3:45,3,0
2291957	21957	.	A	C	20	PASS	AC=2;AF=1.00;AN=2;DP=1;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=29.00;QD=20.00;SOR=1.609	GT:AD:DP:GQ:PL	./.	1/1:0,1:1:3:45,3,0
2291957	21958	.	T	C	20	PASS	AC=2;AF=1.00;AN=2;DP=0;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=NaN;SOR=0.693	GT:AD:DP:GQ:PL	./.	1/1:0,0:0:3:45,3,0
";     # structure of the ref file
$observedSNPLines=`grep -v "#" $vcfVariantsSelected`;      # structure of the test file
is($observedSNPLines,$expectedSNPLines,'gatk::gatkSelectVariants - structure of file');


#################################################################################################
#######Test for gatk Read Backed Phasing
#################################################################################################

#execution Test

my $vcfFileOut="GATKPHASED.vcf";
is(gatk::gatkReadBackedPhasing($fastaRef, $bamIn, $variantFiltered, $vcfFileOut),1, 'gatk::gatkReadBackedPhasing');

# expected output test
$observedOutput = `ls`;
@observedOutput = split /\n/,$observedOutput;
@expectedOutput = ('GATKHAPLOTYPECALLER.vcf','GATKHAPLOTYPECALLER.vcf.idx','gatk_log.e','gatk_log.o','GATKPHASED.vcf','GATKPHASED.vcf.idx','GATKSELECTVARIANT.vcf','GATKSELECTVARIANT.vcf.idx','GATKVARIANTFILTRATION.vcf','GATKVARIANTFILTRATION.vcf.idx','RC3.GATKINDELREALIGNER.bai','RC3.GATKINDELREALIGNER.bam','RC3.GATKINDELREALIGNER.RECALIBRATED.bai','RC3.GATKINDELREALIGNER.RECALIBRATED.bam','RC3.GATKREALIGNERTARGETCREATOR.intervals','RC3.GATKUNIFIEDGENOTYPER.vcf','RC3.GATKUNIFIEDGENOTYPER.vcf.idx','recal_data.table');

is_deeply(\@observedOutput,\@expectedOutput,'gatk::gatkReadBackedPhasing - output list');

# expected output structure test
#$expectedSNPLines="2224477	996	.	TA	T	32.71	PASS	AC=2;AF=1.00;AN=2;DP=2;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=16.35;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,2:2:6:69,6,0
#2248321	377	.	C	G	62.74	PASS	AC=2;AF=1.00;AN=2;DP=2;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=29.00;QD=31.37;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,2:2:6:90,6,0
#2248321	379	.	C	T	62.74	PASS	AC=2;AF=1.00;AN=2;DP=2;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=29.00;QD=31.37;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,2:2:6:90,6,0
#2281178	4213	.	G	A	62.74	PASS	AC=2;AF=1.00;AN=2;DP=2;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=31.37;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,2:2:6:90,6,0
#2281178	4214	.	A	G	62.74	PASS	AC=2;AF=1.00;AN=2;DP=2;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=31.37;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,2:2:6:90,6,0
#2290182	1013	.	A	G	42.74	PASS	AC=2;AF=1.00;AN=2;DP=2;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=29.00;QD=21.37;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,2:2:6:70,6,0
#";     # structure of the ref file
$expectedSNPLines="";
$observedSNPLines=`grep -v "#" $vcfFileOut`;      # structure of the test file
is($observedSNPLines,$expectedSNPLines,'gatk::gatkReadBackedPhasing - structure of file');

# rm index
system ("rm $bamIn.bai") and die ("ERROR: $0 : Cannot remove index file $bamIn.bai\n$!\n");
system ("rm $bamSingle.bai") and die ("ERROR: $0 : Cannot remove index file $bamSingle.bai\n$!\n");
#rm idx file on
my $dataOneVcf = "$toggle/data/testData/vcf/singleVCF";
my $cleaningCmd="rm -f $dataOneVcf/GATKVARIANTFILTRATION.vcf.idx";
system ($cleaningCmd) and die ("ERROR: $0 : Cannot remove GATKVARIANTFILTRATION.vcf.idx with the command $cleaningCmd \n$!\n");



exit;
__END__
