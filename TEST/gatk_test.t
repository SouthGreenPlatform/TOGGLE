#!/usr/bin/perl

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

#Will test if GATK works correctly
use strict;
use warnings;
use Test::More 'no_plan'; #Number of tests, to modify if new tests implemented. Can be changed as 'no_plan' instead of tests=>11 .
use Test::Deep;
use Data::Dumper;
use lib qw(../Modules/);


########################################
#Test of the use of gatk modules
########################################
use_ok('toolbox') or exit;
use_ok('gatk') or exit;

can_ok('gatk','gatkBaseRecalibrator');
can_ok('gatk','gatkRealignerTargetCreator');
can_ok('gatk','gatkIndelRealigner');
can_ok('gatk','gatkHaplotypeCaller');
can_ok('gatk','gatkSelectVariants');
can_ok('gatk','gatkVariantFiltration');
can_ok('gatk','gatkUnifiedGenotyper');
can_ok('gatk','gatkReadBackedPhasing');
use toolbox;
use gatk;


#######################################
#Creating the IndividuSoft.txt file
#######################################
my $creatingCommand="echo \"gatk\nTEST\" > individuSoft.txt";
system($creatingCommand) and die ("ERROR: $0: Cannot create the individuSoft.txt file with the command $creatingCommand \n$!\n");

#######################################
#Cleaning the logs for the test
#######################################
my $cleaningCommand="rm -Rf gatk_TEST_log.*";
system($cleaningCommand) and die ("ERROR: $0: Cannot clean the previous log files for this test with the command $cleaningCommand \n$!\n");

########################################
#initialisation and setting configs
########################################
my $testingDir="../DATA-TEST/gatkTestDir";
my $creatingDirCom="rm -Rf $testingDir ; mkdir -p $testingDir";                                    #Allows to have a working directory for the tests
system($creatingDirCom) and die ("ERROR: $0 : Cannot execute the command $creatingDirCom\n$!\n");

my $OriginalFastaRef="../DATA/expectedData/Reference.fasta";
my $fastaRef="$testingDir/Reference.fasta";
my $refCopyCom="cp $OriginalFastaRef $fastaRef";
system($refCopyCom) and die ("ERROR: $0 : Cannot copy the Reference $OriginalFastaRef with the command $refCopyCom\n$!\n");     #Now we have a ref to be tested

my $OriginalFastaRefFai="../DATA/expectedData/Reference.fasta.fai";
my $fastaRefFai="$testingDir/Reference.fasta.fai";
my $refCopyComFai="cp $OriginalFastaRefFai $fastaRefFai";
system($refCopyComFai) and die ("ERROR: $0 : Cannot copy the Reference $OriginalFastaRefFai with the command $refCopyComFai\n$!\n");     #Now we have a ref to be tested

my $OriginalFastaRefDict="../DATA/expectedData/Reference.dict";
my $fastaRefDict="$testingDir/Reference.dict";
my $refCopyComDict="cp $OriginalFastaRefDict $fastaRefDict";
system($refCopyComDict) and die ("ERROR: $0 : Cannot copy the Reference $OriginalFastaRefDict with the command $refCopyComDict\n$!\n");     #Now we have a ref to be tested

my $originalBamFile="../DATA/expectedData/RC3.SAMTOOLSVIEW.bam";
my $bamFile="$testingDir/RC3.SAMTOOLSVIEW.bam";
my $bamFileCopyCom="cp $originalBamFile $bamFile";
system($bamFileCopyCom) and die("ERROR: $0 : Cannot copy the initial bam file $originalBamFile with the command $bamFileCopyCom\n$!\n");

my $originalBamFileIndex="../DATA/expectedData/RC3.SAMTOOLSVIEW.bam.bai";
my $bamFileIndex="$testingDir/RC3.SAMTOOLSVIEW.bam.bai";
my $bamFileCopyComIndex="cp $originalBamFileIndex $bamFileIndex";
system($bamFileCopyComIndex) and die("ERROR: $0 : Cannot copy the initial bam file index $originalBamFileIndex with the command $bamFileCopyComIndex\n$!\n");

toolbox::readFileConf("software.config.txt");





################################################################################################
#####Test for gatkRealignerTargetCreator
my %optionsHachees = ("-T" => "RealignerTargetCreator","-nt" => "4");        # Hash containing informations
my $optionHachees = \%optionsHachees;                           # Ref of the hash
my $bamToRealigne=$bamFile;
my $intervalsFile="../DATA-TEST/gatkTestDir/RC3.GATKREALIGNERTARGETCREATOR.intervals";
is(gatk::gatkRealignerTargetCreator($fastaRef, $bamToRealigne, $intervalsFile, $optionHachees),1, 'Test for gatk::gatkRealignerTargetCreator');

#####Checking the correct structure for the output file using md5sum
my $expectedMD5sum="ad6d03974c93118e47e2149c7a5f916e";      # structure of the ref file
my $observedMD5sum=`md5sum ../DATA-TEST/gatkTestDir/RC3.GATKREALIGNERTARGETCREATOR.intervals`;       # structure of the test file
my @withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'Test for gatk::gatkRealignerTargetCreator structure of file');


################################################################################################
####Test for gatkIndelRealigner
%optionsHachees = ("-T" => "IndelRealigner");        # Hash containing informations
$optionHachees = \%optionsHachees;                           # Ref of the hash
my $bamRealigned="$testingDir/RC3.GATKINDELREALIGNER.bam";

is(gatk::gatkIndelRealigner($fastaRef, $bamToRealigne, $intervalsFile, $bamRealigned, $optionHachees),1, 'Test for gatk::gatkIndelRealigner');

#####Checking the correct structure for the output file using md5sum
$expectedMD5sum="3296b45af3890583e257447cb8d94f2d";     # structure of the ref file
$observedMD5sum=`md5sum ../DATA-TEST/gatkTestDir/RC3.GATKINDELREALIGNER.bam`;     # structure of the test file
@withoutName =split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];      #just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'Test fot gatk::gatkIndelRealigner structure of file');


################################################################################################
#####Test for gatk Unified Genotyper
#my %optionsHachees = ("-T" => "IndelRealigner");        # Hash containing informations
#my $optionHachees = \%optionsHachees;                           # Ref of the hash
#my $rawvcf="$testingDir/raw.vcf";
#is(gatk::gatkUnifiedGenotyper($fastaRef, $bamFile, $rawvcf, $optionsHachees),1,'gatkUnifiedGenotyper::Running');
#
####Test for correct file value of bwa sampe
#my $grepResult=`grep -c "^LOC" ../DATA-TEST/gatkTestDir/raw.vcf`;
#chomp $grepResult;
#is($grepResult,76,'Test for the result of gatk::gatkUnifiedGenotyper');

################################################################################################
#####Test for gatkBaseRecalibrator
#$optionsHachees=$configInfos->{'GATK gatkBaseRecalibrator'};
#my $bamToRecalibrate="$testingDir/RC3_bamRealigned.bam";
#my $tableReport="$testingDir/recal_data.table";
#my $vcfSnpKnownFile="$testingDir/raw.vcf";
#is(gatk::gatkBaseRecalibrator($fastaRef, $bamFile,  $vcfSnpKnownFile, $tableReport, $optionsHachees),1, 'gatk BaseRecalibrator::Running');
#
#####Checking the correct structure for the output file using md5sum
#$expectedMD5sum="c625e8f8a10cebbc9eda06a22e762eda";     # structure of the ref file
#$observedMD5sum=`md5sum ../DATA-TEST/gatkTestDir/recal_data.table`;     # structure of the test file
#@withoutName =split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
#$observedMD5sum = $withoutName[0];      #just to have the md5sum result
#is($observedMD5sum,$expectedMD5sum,'gatkBaseRecalibrator output\'s content');


#################################################################################################
######Test for gatk Haplotype caller

my $originalBam="../DATA/expectedData/RC3.PICARDTOOLSMARKDUPLICATES.bam";
my $bam="$testingDir/RC3.PICARDTOOLSMARKDUPLICATES.bam";
my $bamCopyCom="cp $originalBam $bam";
system($bamCopyCom) and die("ERROR: $0 : Cannot copy the initial bam file $originalBam with the command $bamCopyCom\n$!\n");

my $originalBamIndex="../DATA/expectedData/RC3.PICARDTOOLSMARKDUPLICATES.bam.bai";
my $bamIndex="$testingDir/RC3.PICARDTOOLSMARKDUPLICATES.bam.bai";
my $bamIndexCopyCom="cp $originalBamIndex $bamIndex";
system($bamIndexCopyCom) and die("ERROR: $0 : Cannot copy the initial bam file $originalBamIndex with the command $bamIndexCopyCom\n$!\n");

my $originalBamSingle="../DATA/expectedData/RC3Single.PICARDTOOLSMARKDUPLICATES.bam";
my $bamSingle="$testingDir/RC3Single.PICARDTOOLSMARKDUPLICATES.bam";
my $bamSingleCopyCom="cp $originalBamSingle $bamSingle";
system($bamSingleCopyCom) and die("ERROR: $0 : Cannot copy the initial bam file $originalBamSingle with the command $bamSingleCopyCom\n$!\n");

my $originalBamSingleIndex="../DATA/expectedData/RC3Single.PICARDTOOLSMARKDUPLICATES.bam.bai";
my $bamSingleIndex="$testingDir/RC3Single.PICARDTOOLSMARKDUPLICATES.bam.bai";
my $bamSingleIndexCopyCom="cp $originalBamSingleIndex $bamSingleIndex";
system($bamSingleIndexCopyCom) and die("ERROR: $0 : Cannot copy the initial bam file $originalBamSingleIndex with the command $bamSingleIndexCopyCom\n$!\n");

%optionsHachees = ("-T" => "HaplotypeCaller");        # Hash containing informations
$optionHachees = \%optionsHachees;                           # Ref of the hash
my @bamsToCall=($bam,$bamSingle);
my $vcfCalled="$testingDir/GATKHAPLOTYPECALLER.vcf";
is(gatk::gatkHaplotypeCaller($fastaRef, $vcfCalled, \@bamsToCall, $optionHachees) ,1,"Test for gatk::gatkHaplotypeCaller");

######Checking the correct structure for the output file using diff
my $expectedSNPLines="2224477	996	.	TA	T	32.71	.	AC=2;AF=1.00;AN=2;DP=2;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;MQ0=0;QD=16.35;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,2:2:6:69,6,0
2248321	377	.	C	G	62.74	.	AC=2;AF=1.00;AN=2;DP=2;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=29.00;MQ0=0;QD=31.37;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,2:2:6:90,6,0
2248321	379	.	C	T	62.74	.	AC=2;AF=1.00;AN=2;DP=2;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=29.00;MQ0=0;QD=31.37;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,2:2:6:90,6,0
2281178	4213	.	G	A	62.74	.	AC=2;AF=1.00;AN=2;DP=2;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;MQ0=0;QD=31.37;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,2:2:6:90,6,0
2281178	4214	.	A	G	62.74	.	AC=2;AF=1.00;AN=2;DP=2;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;MQ0=0;QD=31.37;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,2:2:6:90,6,0
2290182	1013	.	A	G	42.74	.	AC=2;AF=1.00;AN=2;DP=2;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=29.00;MQ0=0;QD=21.37;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,2:2:6:70,6,0
";     # structure of the ref file
my $observedSNPLines=`grep -v "#" $vcfCalled`;      # structure of the test file
is($observedSNPLines,$expectedSNPLines,'Test for gatk::gatkHaplotypeCaller structure of file');


#################################################################################################
######Test for gatkVariant filtrator
%optionsHachees = ("-T" => "VariantFiltration");        # Hash containing informations
$optionHachees = \%optionsHachees;                           # Ref of the hash
my $vcfToFilter="$testingDir/GATKHAPLOTYPECALLER.vcf";
my $variantFiltered="$testingDir/GATKVARIANTFILTRATION.vcf";
is(gatk::gatkVariantFiltration($fastaRef, $variantFiltered, $vcfToFilter, $optionHachees), 1 , 'Test for gatk::gatkVariantFiltration');
##
######Checking the correct structure for the output file using md5sum
$expectedSNPLines="2224477	996	.	TA	T	32.71	PASS	AC=2;AF=1.00;AN=2;DP=2;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;MQ0=0;QD=16.35;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,2:2:6:69,6,0
2248321	377	.	C	G	62.74	PASS	AC=2;AF=1.00;AN=2;DP=2;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=29.00;MQ0=0;QD=31.37;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,2:2:6:90,6,0
2248321	379	.	C	T	62.74	PASS	AC=2;AF=1.00;AN=2;DP=2;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=29.00;MQ0=0;QD=31.37;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,2:2:6:90,6,0
2281178	4213	.	G	A	62.74	PASS	AC=2;AF=1.00;AN=2;DP=2;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;MQ0=0;QD=31.37;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,2:2:6:90,6,0
2281178	4214	.	A	G	62.74	PASS	AC=2;AF=1.00;AN=2;DP=2;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;MQ0=0;QD=31.37;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,2:2:6:90,6,0
2290182	1013	.	A	G	42.74	PASS	AC=2;AF=1.00;AN=2;DP=2;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=29.00;MQ0=0;QD=21.37;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,2:2:6:70,6,0
";     # structure of the ref file
$observedSNPLines=`grep -v "#" $variantFiltered`;      # structure of the test file
is($observedSNPLines,$expectedSNPLines,'Test for gatk::gatkVariantFiltration structure of file');



#################################################################################################
######Test for gatk Select Variant
# my ($refFastaFileIn, $vcfSnpKnownFile, $vcfVariantsSelected, $optionsHachees);
%optionsHachees = ("-T" => "SelectVariants");        # Hash containing informations
$optionHachees = \%optionsHachees;                           # Ref of the hash
my $vcfIn="../DATA-TEST/gatkTestDir/GATKVARIANTFILTRATION.vcf";
my $vcfVariantsSelected="$testingDir/GATKSELECTVARIANTS.vcf";
is(gatk::gatkSelectVariants($fastaRef, $vcfIn, $vcfVariantsSelected, $optionHachees),1, 'Test for gatk::gatkSelectVariants');
##
######Checking the correct structure for the output file using md5sum
$expectedSNPLines="2224477	996	.	TA	T	32.71	PASS	AC=2;AF=1.00;AN=2;DP=2;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;MQ0=0;QD=16.35;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,2:2:6:69,6,0
2248321	377	.	C	G	62.74	PASS	AC=2;AF=1.00;AN=2;DP=2;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=29.00;MQ0=0;QD=31.37;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,2:2:6:90,6,0
2248321	379	.	C	T	62.74	PASS	AC=2;AF=1.00;AN=2;DP=2;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=29.00;MQ0=0;QD=31.37;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,2:2:6:90,6,0
2281178	4213	.	G	A	62.74	PASS	AC=2;AF=1.00;AN=2;DP=2;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;MQ0=0;QD=31.37;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,2:2:6:90,6,0
2281178	4214	.	A	G	62.74	PASS	AC=2;AF=1.00;AN=2;DP=2;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;MQ0=0;QD=31.37;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,2:2:6:90,6,0
2290182	1013	.	A	G	42.74	PASS	AC=2;AF=1.00;AN=2;DP=2;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=29.00;MQ0=0;QD=21.37;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,2:2:6:70,6,0
";     # structure of the ref file
$observedSNPLines=`grep -v "#" $vcfVariantsSelected`;      # structure of the test file
is($observedSNPLines,$expectedSNPLines,'Test for gatk::gatkVariantFiltration structure of file');


#################################################################################################
#######Test for gatk Read Backed Phasing
#$optionsHachees=$configInfos->{'GATK ReadBackedPhasing'};
#my $vcfFileOut="$testingDir/phased.vcf";
#is(gatk::gatkReadBackedPhasing($fastaRef, $bamFile, $variantFiltered, $vcfFileOut, $optionsHachees),1, 'Gatk BackedPhasing::Running');
###
#######Checking the correct structure for the output file using md5sum
#$expectedMD5sum="e528b5da85a26663b22f1fb7d9fb04cb";     # structure of the ref file
#$observedMD5sum=`md5sum ../DATA-TEST/gatkTestDir/phased.vcf`;       # structure of the test file
#@withoutName =split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
#$observedMD5sum = $withoutName[0];      #just to have the md5sum result
#is($observedMD5sum,$expectedMD5sum,'gatk readBackedPhasing output\'scontent');

exit;
