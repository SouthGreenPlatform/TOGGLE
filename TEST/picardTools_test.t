#!/usr/bin/perl -w

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

#Will test if picardsTools module works correctly
use strict;
use warnings;

use Test::More 'no_plan'; 
use Test::Deep;
use Data::Dumper;
use FindBin qw($Bin);
use lib qw(../Modules/);

########################################
#use of picardTools modules ok
########################################
use_ok('toolbox') or exit;
use_ok('picardTools') or exit;
can_ok( 'picardTools','picardToolsMarkDuplicates');
can_ok( 'picardTools','picardToolsCreateSequenceDictionary');
can_ok( 'picardTools','picardToolsSortSam');
can_ok( 'picardTools','picardToolsAddOrReplaceReadGroups');
can_ok( 'picardTools','picardToolsCleanSam');
can_ok( 'picardTools','picardToolsSamFormatConverter');
can_ok( 'picardTools','picardToolsValidateSamFile');

use toolbox;
use picardTools;

#########################################
#Remove files and directory created by previous test
#########################################
my $testingDir="../DATA-TEST/picardtoolsTestDir";
my $cleaningCmd="rm -Rf $testingDir"; 
system ($cleaningCmd) and die ("ERROR: $0 : Cannot remove the previous test directory with the command $cleaningCmd \n$!\n");

my $expectedData="../../DATA/expectedData/";

########################################
#Creation of test directory
########################################
my $makeDirCmd = "mkdir $testingDir";
system ($makeDirCmd) and die ("ERROR: $0 : Cannot create the new directory with the command $makeDirCmd\n$!\n");
chdir $testingDir or die ("ERROR: $0 : Cannot go into the new directory with the command \"chdir $testingDir\"\n$!\n");

#######################################
#Creating the IndividuSoft.txt file
#######################################
my $creatingCmd="echo \"picardtools\nTEST\" > individuSoft.txt";
system($creatingCmd) and die ("ERROR: $0 : Cannot create the individuSoft.txt file with the command $creatingCmd\n$!\n");

#######################################
#Cleaning the logs for the test
#######################################
$cleaningCmd="rm -Rf picardtools_TEST_log.*";
system($cleaningCmd) and die ("ERROR: $0 : Cannot remove the previous log files with the command $cleaningCmd \n$!\n");

##########################################
#picardToolsCreateSequenceDictionary test
##########################################

#Input file
my $refFile = "Reference.fasta";  
my $originalRefFile = $expectedData."/".$refFile;    
my $cpCmd = "cp $originalRefFile ."; # command to copy the original Ref fasta file into the test directory
system ($cpCmd) and die ("ERROR: $0 : Cannot copy the file $originalRefFile in the test directory with the command $cpCmd\n$!\n"); 

#Output file
my $refFileDict = "Reference.dict";

#execution test
is(picardTools::picardToolsCreateSequenceDictionary($refFile,$refFileDict),1,'picardTools::picardToolsCreateSequenceDictionary');

# expected output test
my $observedOutput = `ls`;
my @observedOutput = split /\n/,$observedOutput;
my @expectedOutput = ('individuSoft.txt','picardtools_TEST_log.e','picardtools_TEST_log.o','Reference.dict','Reference.fasta');

is_deeply(\@observedOutput,\@expectedOutput,'picardTools::picardToolsCreateSequenceDictionary - output list');

# expected content test
my $expectedLastLine="\@SQ	SN:2299897	";  
my $observedLastLine=`tail -n 1 $refFileDict`; 
my @withoutName = split ("LN:", $observedLastLine); 
$observedLastLine = $withoutName[0];       # just to have the md5sum result
is($observedLastLine,$expectedLastLine,'picardTools::picardToolsCreateSequenceDictionary - output structure');



##########################################
#picardToolsSortSam test
##########################################
#input data
## Input files test for single analysis
my $samFile = $expectedData."RC3.BWASAMPE.sam";

#output data
my $bamFileOut = "RC3.PICARDTOOLSSORT.bam";

my %optionsRef = ("SORT_ORDER" => "coordinate","VALIDATION_STRINGENCY" => "SILENT");   
my $optionsHachees = \%optionsRef;

#execution test
is(picardTools::picardToolsSortSam($samFile,$bamFileOut,$optionsHachees),1,'picardTools::picardToolsSortSam');  

# expected output test
$observedOutput = `ls`;
@observedOutput = split /\n/,$observedOutput;
@expectedOutput = ('individuSoft.txt','picardtools_TEST_log.e','picardtools_TEST_log.o','RC3.PICARDTOOLSSORT.bam','Reference.dict','Reference.fasta');

is_deeply(\@observedOutput,\@expectedOutput,'picardTools::picardToolsSortSam - output list');

# expected content test
#my $expectedMD5sum="df7c7657d50f6dbf93a5ba5b6900b981";      # structure of the ref file
#my $observedMD5sum=`md5sum $bamFileOut`;       # structure of the test file
#@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
#$observedMD5sum = $withoutName[0];       # just to have the md5sum result
#is($observedMD5sum,$expectedMD5sum,'picardTools::picardToolsSortSam - output structure');

#test end of file because change in version 130
my $expectedEndLine="H2:C381HACXX:5:1101:9881:2219	141	*	0	0	*	*	0	0	CTCTTAGATCTTCTTTCTCCAATCTTGGATTAGGGAAGAAGGAGATATTCGCGACTCCTGGTGGTTTCATTATGGGGCAGCTCATGATCTTCATATCGATC	=;?DDAFBF>?<,<EF\@CIH:EHG4,<3+<293AF;:;?DBE8B\@9B<BF;\@D';@\@CDA).71'56??6(.6632;3>;8:(:@@(5:3(>\@:>BC?<?<	RG:Z:RC3
";
my $observedEndLine=`samtools view RC3.PICARDTOOLSSORT.bam | tail -1`  or die ("ERROR: $0 : Cannot execute: samtools view RC3.PICARDTOOLSSORT.bam | tail -1  \n$!\n");
is($observedEndLine,$expectedEndLine,'picardTools::picardToolsSortSam - output endFile');


###########################################
##picardToolsMarkDuplicates test
###########################################
#input file
my $bamFile = $expectedData."RC3.GATKINDELREALIGNER.bam"; 

#output files
$bamFileOut = "RC3.PICARDTOOLSMARKDUPLICATES.bam";
my $duplicatesFileOut = "RC3.PICARDTOOLSMARKDUPLICATES.bamDuplicates";

%optionsRef = ("VALIDATION_STRINGENCY" => "SILENT");        # Hash containing informations
$optionsHachees = \%optionsRef;                           # Ref of the hash

#execution test
is(picardTools::picardToolsMarkDuplicates($bamFile, $bamFileOut, $duplicatesFileOut, $optionsHachees),1,'picardTools::picardToolsMarkDuplicates');
# expected output test
$observedOutput = `ls`;
@observedOutput = split /\n/,$observedOutput;
@expectedOutput = ('individuSoft.txt','picardtools_TEST_log.e','picardtools_TEST_log.o','RC3.PICARDTOOLSMARKDUPLICATES.bam','RC3.PICARDTOOLSMARKDUPLICATES.bamDuplicates','RC3.PICARDTOOLSSORT.bam','Reference.dict','Reference.fasta');

is_deeply(\@observedOutput,\@expectedOutput,'picardTools::picardToolsMarkDuplicates - output list');

# expected content test
my $expectedLastLines="H2:C381HACXX:5:1101:9881:2219	77	*	0	0	*	*	0	0	AGTCCATGATATAACCAAATTGGATGGATCTTCCACCCGTTTAGCTAAGAAAGAATAGATGCAGAGGTGGATAATAGATCGATATGAAGATCATGAGCTGC	?7=DBB;=DF>C?CG<?FFEIIF3E<?EE\@FBFF<EGB6):D?4<9?D309??\@BF<B)8BBF).6;=CEF<E?A7?>@)?7==?A:AA>ABA5,:>>A:A	PG:Z:MarkDuplicates	RG:Z:RC3
H2:C381HACXX:5:1101:9881:2219	141	*	0	0	*	*	0	0	CTCTTAGATCTTCTTTCTCCAATCTTGGATTAGGGAAGAAGGAGATATTCGCGACTCCTGGTGGTTTCATTATGGGGCAGCTCATGATCTTCATATCGATC	=;?DDAFBF>?<,<EF\@CIH:EHG4,<3+<293AF;:;?DBE8B\@9B<BF;\@D';@\@CDA).71'56??6(.6632;3>;8:(:@@(5:3(>\@:>BC?<?<	PG:Z:MarkDuplicates	RG:Z:RC3";      
my $observedLastLines=`samtools view $bamFileOut |tail -n 2`;       
chomp $observedLastLines;       
is($observedLastLines,$expectedLastLines,'picardTools::picardToolsMarkDuplicates - output structure');


###########################################
# picardToolsCleanSam test
###########################################
#output files
$bamFileOut = "RC3.PICARDTOOLSCLEANSAM.bam";

%optionsRef = ("VALIDATION_STRINGENCY" => "SILENT");        # Hash containing informations
$optionsHachees = \%optionsRef;                           # Ref of the hash

#execution test
is(picardTools::picardToolsCleanSam($bamFile, $bamFileOut,$optionsHachees),1,'picardTools::picardToolsCleanSam');

# expected output test
$observedOutput = `ls`;
@observedOutput = split /\n/,$observedOutput;
@expectedOutput = ('individuSoft.txt','picardtools_TEST_log.e','picardtools_TEST_log.o','RC3.PICARDTOOLSCLEANSAM.bam','RC3.PICARDTOOLSMARKDUPLICATES.bam','RC3.PICARDTOOLSMARKDUPLICATES.bamDuplicates','RC3.PICARDTOOLSSORT.bam',,'Reference.dict','Reference.fasta');

is_deeply(\@observedOutput,\@expectedOutput,'picardTools::picardToolsCleanSam - output list');

# expected content test
$expectedLastLines="H2:C381HACXX:5:1101:9881:2219	77	*	0	0	*	*	0	0	AGTCCATGATATAACCAAATTGGATGGATCTTCCACCCGTTTAGCTAAGAAAGAATAGATGCAGAGGTGGATAATAGATCGATATGAAGATCATGAGCTGC	?7=DBB;=DF>C?CG<?FFEIIF3E<?EE\@FBFF<EGB6):D?4<9?D309??\@BF<B)8BBF).6;=CEF<E?A7?>@)?7==?A:AA>ABA5,:>>A:A	RG:Z:RC3
H2:C381HACXX:5:1101:9881:2219	141	*	0	0	*	*	0	0	CTCTTAGATCTTCTTTCTCCAATCTTGGATTAGGGAAGAAGGAGATATTCGCGACTCCTGGTGGTTTCATTATGGGGCAGCTCATGATCTTCATATCGATC	=;?DDAFBF>?<,<EF\@CIH:EHG4,<3+<293AF;:;?DBE8B\@9B<BF;\@D';@\@CDA).71'56??6(.6632;3>;8:(:@@(5:3(>\@:>BC?<?<	RG:Z:RC3";      
$observedLastLines=`samtools view $bamFileOut |tail -n 2`;       
chomp $observedLastLines;       
is($observedLastLines,$expectedLastLines,'picardTools::picardToolsCleanSam - output structure');


###########################################
# picardToolsSamFormatConverter test
###########################################
#output files
my $samFileOut = "RC3.PICARDTOOLSSAMFORMATCONVERTER.sam";

%optionsRef = ("VALIDATION_STRINGENCY" => "SILENT");        # Hash containing informations
$optionsHachees = \%optionsRef;                           # Ref of the hash

#execution test
is(picardTools::picardToolsSamFormatConverter($bamFile, $samFileOut,$optionsHachees),1,'picardTools::picardToolsSamFormatConverter');

# expected output test
$observedOutput = `ls`;
@observedOutput = split /\n/,$observedOutput;
@expectedOutput = ('individuSoft.txt','picardtools_TEST_log.e','picardtools_TEST_log.o','RC3.PICARDTOOLSCLEANSAM.bam','RC3.PICARDTOOLSMARKDUPLICATES.bam','RC3.PICARDTOOLSMARKDUPLICATES.bamDuplicates','RC3.PICARDTOOLSSAMFORMATCONVERTER.sam','RC3.PICARDTOOLSSORT.bam',,'Reference.dict','Reference.fasta');

is_deeply(\@observedOutput,\@expectedOutput,'picardTools::picardToolsSamFormatConverter - output list');

# expected content test
$expectedLastLines="H2:C381HACXX:5:1101:9881:2219	77	*	0	0	*	*	0	0	AGTCCATGATATAACCAAATTGGATGGATCTTCCACCCGTTTAGCTAAGAAAGAATAGATGCAGAGGTGGATAATAGATCGATATGAAGATCATGAGCTGC	?7=DBB;=DF>C?CG<?FFEIIF3E<?EE\@FBFF<EGB6):D?4<9?D309??\@BF<B)8BBF).6;=CEF<E?A7?>@)?7==?A:AA>ABA5,:>>A:A	RG:Z:RC3
H2:C381HACXX:5:1101:9881:2219	141	*	0	0	*	*	0	0	CTCTTAGATCTTCTTTCTCCAATCTTGGATTAGGGAAGAAGGAGATATTCGCGACTCCTGGTGGTTTCATTATGGGGCAGCTCATGATCTTCATATCGATC	=;?DDAFBF>?<,<EF\@CIH:EHG4,<3+<293AF;:;?DBE8B\@9B<BF;\@D';@\@CDA).71'56??6(.6632;3>;8:(:@@(5:3(>\@:>BC?<?<	RG:Z:RC3";      
$observedLastLines=`tail -n 2 $samFileOut`;       
chomp $observedLastLines;       
is($observedLastLines,$expectedLastLines,'picardTools::picardToolsSamFormatConverter - output structure');

###########################################
# picardToolsAddOrReplaceReadGroups test
###########################################
#output files
$bamFileOut = "RC3.PICARDTOOLSADDORREPLACEREADGROUPS.bam";

%optionsRef = ("ID" => "Test","LB" => "Irigin","PL" => "Illumina","SM" => "glaberrima","VALIDATION_STRINGENCY" => "SILENT","PU" => "unit1");        # Hash containing informations
$optionsHachees = \%optionsRef;                           # Ref of the hash

#execution test
is(picardTools::picardToolsAddOrReplaceReadGroups($bamFile, $bamFileOut,$optionsHachees),1,'picardTools::picardToolsAddOrReplaceReadGroups');

# expected output test
$observedOutput = `ls`;
@observedOutput = split /\n/,$observedOutput;
@expectedOutput = ('individuSoft.txt','picardtools_TEST_log.e','picardtools_TEST_log.o','RC3.PICARDTOOLSADDORREPLACEREADGROUPS.bam','RC3.PICARDTOOLSCLEANSAM.bam','RC3.PICARDTOOLSMARKDUPLICATES.bam','RC3.PICARDTOOLSMARKDUPLICATES.bamDuplicates','RC3.PICARDTOOLSSAMFORMATCONVERTER.sam','RC3.PICARDTOOLSSORT.bam',,'Reference.dict','Reference.fasta');

is_deeply(\@observedOutput,\@expectedOutput,'picardTools::picardToolsAddOrReplaceReadGroups - output list');

# expected content test
$expectedLastLine="\@RG	ID:Test	LB:Irigin	PL:Illumina	SM:glaberrima	PU:unit1";
$observedOutput=`samtools view -H RC3.PICARDTOOLSADDORREPLACEREADGROUPS.bam| grep \@RG`; # We pick up only the position field
chomp $observedOutput;
is($observedOutput,$expectedLastLine,'picardTools::picardToolsAddOrReplaceReadGroups - output structure');


###########################################
# picardToolsValidateSamFile test

#THE SOFT WILL STOP JOB AS SOON AS THE SAM IS NOT COMPLETELY PERFECT... CANNOT BE TESTED
###########################################
#input file
#$bamFile = "RC3.PICARDTOOLSCLEANSAM.bam";
#
##output files
#my $infoFileOut = "RC3.PICARDTOOLSVALIDATESAMFILE.infos";
#
#%optionsRef = ("VALIDATION_STRINGENCY" => "SILENT");        # Hash containing informations
#$optionsHachees = \%optionsRef;                           # Ref of the hash
#
#TODO:{
##execution test
#    is(picardTools::picardToolsValidateSamFile($bamFile, $infoFileOut,$optionsHachees),1,'picardTools::picardToolsValidateSamFile');
#}
#
## expected output test
#$observedOutput = `ls`;
#@observedOutput = split /\n/,$observedOutput;
#@expectedOutput = ('individuSoft.txt','picardtools_TEST_log.e','picardtools_TEST_log.o','RC3.PICARDTOOLSCLEANSAM.bam','RC3.PICARDTOOLSMARKDUPLICATES.bam','RC3.PICARDTOOLSMARKDUPLICATES.bamDuplicates','RC3.PICARDTOOLSSORT.bam','RC3.PICARDTOOLSVALIDATESAMFILE.infos','Reference.dict','Reference.fasta');
#
#is_deeply(\@observedOutput,\@expectedOutput,'picardTools::picardToolsValidateSamFile - output list');
#
## expected content test
#$expectedMD5sum="23180154325eaf0bedb97488846a3592";      # structure of the ref file
#$observedMD5sum=`md5sum $infoFileOut`;       # structure of the test file
#@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
#$observedMD5sum = $withoutName[0];       # just to have the md5sum result
#is($observedMD5sum,$expectedMD5sum,'picardTools::picardToolsValidateSamFile - output structure');

exit;
__END__


#PicardToolsSamFormatConverter
