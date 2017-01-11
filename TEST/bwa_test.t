#!/usr/bin/perl

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

#Will test if bwa works correctly
use strict;
use warnings;
use Test::More 'no_plan'; #Number of tests, to modify if new tests implemented. Can be changed as 'no_plan' instead of tests=>11 .
use Test::Deep;
use Data::Dumper;
use lib qw(../Modules/);


########################################
#use of bwa module ok
########################################
use_ok('toolbox') or exit;
use_ok('bwa') or exit;
can_ok( 'bwa','bwaIndex');
can_ok('bwa','bwaAln');
can_ok('bwa','bwaSampe');
can_ok('bwa','bwaSamse');
can_ok('bwa','bwaMem');

use toolbox;
use bwa;

my $expectedData="../../DATA/expectedData/";
#my $configInfos = toolbox::readFileConf("software.config.txt");

#########################################
#Remove files and directory created by previous test
#########################################
my $testingDir="../DATA-TEST/bwaTestDir";
my $creatingDirCom="rm -Rf $testingDir ; mkdir -p $testingDir";                                    #Allows to have a working directory for the tests
system($creatingDirCom) and die ("ERROR: $0 : Cannot execute the command $creatingDirCom\n$!\n");

chdir $testingDir or die ("ERROR: $0 : Cannot go into the new directory with the command \"chdir $testingDir\"\n$!\n");


#######################################
#Creating the IndividuSoft.txt file
#######################################
my $creatingCommand="echo \"bwa\nTEST\" > individuSoft.txt";
system($creatingCommand) and die ("ERROR: $0: Cannot create the individuSoft.txt file with the command $creatingCommand \n$!\n");


#######################################
#Cleaning the logs for the test
#######################################
my $cleaningCommand="rm -Rf bwa_TEST_log.*";
system($cleaningCommand) and die ("ERROR: $0: Cannot clean the previous log files for this test with the command $cleaningCommand \n$!\n");


##########################################
##### bwa::index
##########################################


# input file
my $fastaRef="Reference.fasta";

my $originalFastaRef=$expectedData."/Reference.fasta";
my $copyCmd= "cp $originalFastaRef $fastaRef";           # command to copy the original fasta file into the test directory
system ($copyCmd) and die ("ERROR: $0 : Cannot link the file $originalFastaRef in the test directory with the command $copyCmd\n$!\n");    # RUN the copy command

# output file
my $fastaRefBWT="Reference.fasta.bwt";
my $fastaRefPAC="Reference.fasta.pac";
my $fastaRefANN="Reference.fasta.ann";
my $fastaRefAMB="Reference.fasta.amb";
my $fastaRefSA="Reference.fasta.sa";

# execution test
my %optionsHachees = (
			"-a" => "is",
			);        # Hash containing informations
my $optionsHachees = \%optionsHachees;  

is(bwa::bwaIndex($fastaRef,$optionsHachees),1,'bwa::bwaIndex - running');

# expected output test
#Check if files created
my @expectedOutput = ("bwa_TEST_log.e","bwa_TEST_log.o","individuSoft.txt","Reference.fasta","Reference.fasta.amb","Reference.fasta.ann","Reference.fasta.bwt","Reference.fasta.pac","Reference.fasta.sa");
my $observedOutput = `ls`;
my @observedOutput = split /\n/,$observedOutput;
is_deeply(\@observedOutput,\@expectedOutput,'bwa::bwaIndex - Filetree created');


# expected content test $fastaRefBWT
my $expectedMD5sum = "e4fdc0af9540ee8365e7e324fc5c0cc3";                                            # structure of the ref file for checking
my $observedMD5sum = `md5sum $fastaRefBWT`;                                                        # structure of the test file for checking
my @withoutName = split (" ", $observedMD5sum);                                                     # to separate the structure and the name of file
$observedMD5sum = $withoutName[0];     										                        # just to have the md5sum result
is($observedMD5sum, $expectedMD5sum, "bwa::index - output content file BWT");               # TEST IF THE STRUCTURE OF THE FILE OUT IS GOOD

# expected content test $fastaRefPAC
$expectedMD5sum = "03454e7242900c436d9d7126f492e4d5";                                            # structure of the ref file for checking
$observedMD5sum = `md5sum $fastaRefPAC`;                                                        # structure of the test file for checking
@withoutName = split (" ", $observedMD5sum);                                                     # to separate the structure and the name of file
$observedMD5sum = $withoutName[0];     										                        # just to have the md5sum result
is($observedMD5sum, $expectedMD5sum, "bwa::index - output content file PAC");               # TEST IF THE STRUCTURE OF THE FILE OUT IS GOOD

# expected content test $fastaRefANN
$expectedMD5sum = "a51b6a5152f51b13833a40fe609474ea";                                            # structure of the ref file for checking
$observedMD5sum = `md5sum $fastaRefANN`;                                                        # structure of the test file for checking
@withoutName = split (" ", $observedMD5sum);                                                     # to separate the structure and the name of file
$observedMD5sum = $withoutName[0];     										                        # just to have the md5sum result
is($observedMD5sum, $expectedMD5sum, "bwa::index - output content file ANN");               # TEST IF THE STRUCTURE OF THE FILE OUT IS GOOD

# expected content test $fastaRefAMB
$expectedMD5sum = "b86728bb71903f8641530e61e9687b59";                                            # structure of the ref file for checking
$observedMD5sum = `md5sum $fastaRefAMB`;                                                        # structure of the test file for checking
@withoutName = split (" ", $observedMD5sum);                                                     # to separate the structure and the name of file
$observedMD5sum = $withoutName[0];     										                        # just to have the md5sum result
is($observedMD5sum, $expectedMD5sum, "bwa::index - output content file AMB");               # TEST IF THE STRUCTURE OF THE FILE OUT IS GOOD

# expected content test $fastaRefSA
$expectedMD5sum = "9243bf066de0cc18aa0d3813f174cae8";                                            # structure of the ref file for checking
$observedMD5sum = `md5sum $fastaRefSA`;                                                        # structure of the test file for checking
@withoutName = split (" ", $observedMD5sum);                                                     # to separate the structure and the name of file
$observedMD5sum = $withoutName[0];     										                        # just to have the md5sum result
is($observedMD5sum, $expectedMD5sum, "bwa::index - output content file SA");               # TEST IF THE STRUCTURE OF THE FILE OUT IS GOOD


##########################################
##### bwa::bwaAln
##########################################

# input file
my $forwardFastq=$expectedData."RC3_1.REPAIRING.fastq";
my $reverseFastq=$expectedData."RC3_2.REPAIRING.fastq";

# output file
my $forwardSaiFileIn="RC3_1.BWAALN.sai";
my $reverseSaiFileIn="RC3_2.BWAALN.sai";

# execution test
%optionsHachees = (
			"-n" => 5,
			"-o" => 1,
			);        # Hash containing informations
$optionsHachees = \%optionsHachees;

is (bwa::bwaAln($fastaRef,$forwardFastq,$forwardSaiFileIn,$optionsHachees),'1',"bwa::bwaAln - Test for bwa Aln running for forward");
is (bwa::bwaAln($fastaRef,$reverseFastq,$reverseSaiFileIn,$optionsHachees),'1',"bwa::bwaAln - Test for bwa Aln running for reverse");

# expected output test
#Check if files created
@expectedOutput = ("bwa_TEST_log.e","bwa_TEST_log.o","individuSoft.txt","RC3_1.BWAALN.sai","RC3_2.BWAALN.sai","Reference.fasta","Reference.fasta.amb","Reference.fasta.ann","Reference.fasta.bwt","Reference.fasta.pac","Reference.fasta.sa");
$observedOutput = `ls`;
@observedOutput = split /\n/,$observedOutput;
is_deeply(\@observedOutput,\@expectedOutput,'bwa::aln - Files created');

# expected content test $forwardSaiFileIn
#$expectedMD5sum = "8b6a2da9c90bd105f8b55bf3867e7f64";                                            # structure of the ref file for checking
#$observedMD5sum = `md5sum $forwardSaiFileIn`;                                                        # structure of the test file for checking
#@withoutName = split (" ", $observedMD5sum);                                                     # to separate the structure and the name of file
#$observedMD5sum = $withoutName[0];     										                        # just to have the md5sum result
ok((-s $forwardSaiFileIn)>0, "bwa::aln - output content file sai forward");               # TEST IF THE STRUCTURE OF THE FILE OUT IS GOOD

# expected content test $reverseSaiFileIn
ok((-s $reverseSaiFileIn)>0, "bwa::aln - output content file sai forward");               # TEST IF THE STRUCTURE OF THE FILE OUT IS GOOD

##########################################
##### bwa::sampe
##########################################

# output file
my $samFileOut="RC3.BWASAMPE.sam";
my $readGroupLine="RC3";

# execution test
is(bwa::bwaSampe($samFileOut,$fastaRef,$forwardSaiFileIn,$reverseSaiFileIn,$forwardFastq,$reverseFastq,$readGroupLine),'1',"bwa::sampe - Test for bwa sampe running");

# expected output test
#Check if files created
@expectedOutput = ("bwa_TEST_log.e","bwa_TEST_log.o","individuSoft.txt","RC3_1.BWAALN.sai","RC3_2.BWAALN.sai","RC3.BWASAMPE.sam","Reference.fasta","Reference.fasta.amb","Reference.fasta.ann","Reference.fasta.bwt","Reference.fasta.pac","Reference.fasta.sa");
$observedOutput = `ls`;
@observedOutput = split /\n/,$observedOutput;
is_deeply(\@observedOutput,\@expectedOutput,'bwa::sampe - Files created');


# expected content test $samFileOut
my $expectedLineNumber = "2947 $samFileOut";                                            # structure of the ref file for checking
my $observedLineNumber = `wc -l $samFileOut`;                                                        # structure of the test file for checking
chomp $observedLineNumber;                                                     # to separate the structure and the name of file
is($observedLineNumber, $expectedLineNumber, "bwa::sampe - output content file sam");               # TEST IF THE STRUCTURE OF THE FILE OUT IS GOOD

###Test for correct file value of bwa sampe
#GREP command result
my $grepResult=`grep -c "XT:A:U" $samFileOut`;
chomp $grepResult;
is($grepResult,1699,'bwa::sampe - output grep in file sam');


##########################################
##### bwa::samse
##########################################

# input file
my $fastqFile=$expectedData."RC3.REPAIRING.fastq";

# output files
my $singleSaiFileIn="RC3.BWAALN.sai";
my $samseFileOut="RC3.BWASAMSE.sam";

is (bwa::bwaAln($fastaRef,$fastqFile,$singleSaiFileIn,$optionsHachees),'1',"bwa::aln - Test for bwa Aln running for single");
is (bwa::bwaSamse($samseFileOut,$fastaRef,$singleSaiFileIn,$fastqFile,$readGroupLine),'1',"bwa::samse - Test for bwa samse running");

# expected output test
#Check if files created
@expectedOutput = ("bwa_TEST_log.e","bwa_TEST_log.o","individuSoft.txt","RC3_1.BWAALN.sai","RC3_2.BWAALN.sai","RC3.BWAALN.sai","RC3.BWASAMPE.sam","RC3.BWASAMSE.sam","Reference.fasta","Reference.fasta.amb","Reference.fasta.ann","Reference.fasta.bwt","Reference.fasta.pac","Reference.fasta.sa");
$observedOutput = `ls`;
@observedOutput = split /\n/,$observedOutput;
is_deeply(\@observedOutput,\@expectedOutput,'bwa::samse - Files created');


# expected content test $singleSaiFileIn
ok((-s $singleSaiFileIn)>0, "bwa::aln - output content file sai forward");               # TEST IF THE STRUCTURE OF THE FILE OUT IS GOOD

# expected content test $samseFileOut
$expectedLineNumber = "954 $samseFileOut";                                            # structure of the ref file for checking
$observedLineNumber = `wc -l $samseFileOut`;                                                        # structure of the test file for checking
chomp $observedLineNumber;     										                        # just to have the md5sum result
is($observedLineNumber, $expectedLineNumber, "bwa::samse - output content file sam");               # TEST IF THE STRUCTURE OF THE FILE OUT IS GOOD

####Test for correct file value of bwa samse - 
$grepResult= `grep -c "XT:A:U" $samseFileOut`;
chomp $grepResult;
is($grepResult,1,'bwa::samse - output grep in file sam');

##########################################
##### bwa::mem single
##########################################

#output files
$samFileOut="RC3.BWAMEM.sam";

##Running test
is (bwa::bwaMem($samFileOut,$fastaRef,$forwardFastq,undef,$readGroupLine),'1',"bwa::bwaMem - Test for bwa mem running single");


###Verify if output are correct for mem single
@expectedOutput = ("bwa_TEST_log.e","bwa_TEST_log.o","individuSoft.txt","RC3_1.BWAALN.sai","RC3_2.BWAALN.sai","RC3.BWAALN.sai","RC3.BWAMEM.sam","RC3.BWASAMPE.sam","RC3.BWASAMSE.sam","Reference.fasta","Reference.fasta.amb","Reference.fasta.ann","Reference.fasta.bwt","Reference.fasta.pac","Reference.fasta.sa");
$observedOutput = `ls`;
@observedOutput = split /\n/,$observedOutput;
is_deeply(\@observedOutput,\@expectedOutput,'bwa::mem Single - Files created');

##Output value test
$grepResult= `grep -c 'XT:A:U' $samFileOut`;
chomp $grepResult;
is($grepResult,0,'bwa::mem - Test for the result of bwa mem single');

##########################################
##### bwa::mem single
##########################################

#output files
$samFileOut="RC3.BWAMEMPaired.sam";

##Running test
is (bwa::bwaMem($samFileOut,$fastaRef,$forwardFastq,$reverseFastq,$readGroupLine),'1',"bwa::mem - Test for bwa mem running paired");


###Verify if output are correct for mem single
@expectedOutput = ("bwa_TEST_log.e","bwa_TEST_log.o","individuSoft.txt","RC3_1.BWAALN.sai","RC3_2.BWAALN.sai","RC3.BWAALN.sai","RC3.BWAMEMPaired.sam","RC3.BWAMEM.sam","RC3.BWASAMPE.sam","RC3.BWASAMSE.sam","Reference.fasta","Reference.fasta.amb","Reference.fasta.ann","Reference.fasta.bwt","Reference.fasta.pac","Reference.fasta.sa");
$observedOutput = `ls`;
@observedOutput = split /\n/,$observedOutput;
is_deeply(\@observedOutput,\@expectedOutput,'bwa::mem Paired - Files created');

##Output value test
$grepResult= `grep -c 'XT:A:U' $samFileOut`;
chomp $grepResult;
is($grepResult,0,'bwa::mem - Test for the result of bwa mem paired');

exit;


