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
### Test for bowtie::bowtieBuild
##########################################
my $bankData="$toggle/data/Bank/";
my $fastqData="$toggle/data/testData/fastq/pairedTwoIndividusIrigin/";

# input file

my $fastaRefIni=$bankData."/referenceIrigin.fasta";
my $fastaRef="referenceIrigin.fasta";

#copy fasta reference into test directory where the index will be created
my $copyCommand="cp $fastaRefIni ./$fastaRef";
system ($copyCommand) and die "ERROR: $0: Cannot copy the $fastaRefIni file with the command $copyCommand \n$!\n";

# execution test
is(bowtie::bowtieBuild($fastaRef),$fastaRef,'bowtie::bowtieBuild');

# expected output test
my $observedOutput = `ls`;
my @observedOutput = split /\n/,$observedOutput;
my @expectedOutput = ('bowtie_log.e','bowtie_log.o','referenceIrigin.fasta','referenceIrigin.fasta.1.ebwt','referenceIrigin.fasta.2.ebwt','referenceIrigin.fasta.3.ebwt','referenceIrigin.fasta.4.ebwt','referenceIrigin.fasta.rev.1.ebwt','referenceIrigin.fasta.rev.2.ebwt');
is_deeply(\@observedOutput,\@expectedOutput,'bowtie::bowtieBuild - output list');

# expected output content
my $expectedMD5sum="de1ef57892bd9f508fb466521bd5a5b6";
my $observedMD5sum=`md5sum $fastaRef.1.ebwt`;# structure of the test file
my @withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'bowtie::bowtieBuild - output content 1.ebwt');

$expectedMD5sum="5fe542df841de8685b4ee1c694b52f64";
$observedMD5sum=`md5sum $fastaRef.2.ebwt`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'bowtie::bowtieBuild - output content 2.ebwt');

$expectedMD5sum="dc12cca8433dfb22df23bc78bc6aeef6";
$observedMD5sum=`md5sum $fastaRef.3.ebwt`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'bowtie::bowtieBuild - output content 3.ebwt');

$expectedMD5sum="3d11892beee30c866ee5e2a06bbbc3d8";
$observedMD5sum=`md5sum $fastaRef.4.ebwt`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'bowtie::bowtieBuild - output content 4.ebwt');

$expectedMD5sum="cdf0694f4adfc7c5773f59c234081e98";
$observedMD5sum=`md5sum $fastaRef.rev.1.ebwt`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'bowtie::bowtieBuild - output content rev.1.ebwt');

$expectedMD5sum="f55fc9bd3bc5298fb0946289db6cff66";
$observedMD5sum=`md5sum $fastaRef.rev.2.ebwt`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'bowtie::bowtieBuild - output content rev.2.ebwt');


###############################################################################################
##bowtie::bowtie2Build
###############################################################################################

my %optionsHachees = ();                # Hash containing informations
my $optionHachees = \%optionsHachees;   # Ref of the hash

# execution test
is(bowtie::bowtie2Build($fastaRef,$optionHachees),$fastaRef, 'bowtie::bowtie2Build');

# expected output test
$observedOutput = `ls`;
@observedOutput = split /\n/,$observedOutput;
@expectedOutput = ('bowtie_log.e','bowtie_log.o','referenceIrigin.fasta','referenceIrigin.fasta.1.bt2','referenceIrigin.fasta.1.ebwt','referenceIrigin.fasta.2.bt2','referenceIrigin.fasta.2.ebwt','referenceIrigin.fasta.3.bt2','referenceIrigin.fasta.3.ebwt','referenceIrigin.fasta.4.bt2','referenceIrigin.fasta.4.ebwt','referenceIrigin.fasta.rev.1.bt2','referenceIrigin.fasta.rev.1.ebwt','referenceIrigin.fasta.rev.2.bt2','referenceIrigin.fasta.rev.2.ebwt');
##print Dumper(\@observedOutput);
is_deeply(\@observedOutput,\@expectedOutput,'bowtie::bowtie2Build - output list');

# expected output content
###Checking the correct structure for the output file using md5sum
$expectedMD5sum="b7a6d65d4bcefe2332dcdc8e9c0cb9c1";
$observedMD5sum=`md5sum $fastaRef.1.bt2`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'bowtie::bowtie2Build - output content 1.bt2');

$expectedMD5sum="481b0055258e98825bb4a8c52c3e90c0";
$observedMD5sum=`md5sum $fastaRef.2.bt2`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'bowtie::bowtie2Build - output content 2.bt2');

$expectedMD5sum="dc12cca8433dfb22df23bc78bc6aeef6";
$observedMD5sum=`md5sum $fastaRef.3.bt2`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'bowtie::bowtie2Build - output content 3.bt2');

$expectedMD5sum="3d11892beee30c866ee5e2a06bbbc3d8";
$observedMD5sum=`md5sum $fastaRef.4.bt2`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'bowtie::bowtie2Build - output content 4.bt2');

$expectedMD5sum="c5c13aab02f5bf0d2701b8de21df32ec";
$observedMD5sum=`md5sum $fastaRef.rev.1.bt2`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'bowtie::bowtie2Build - output content rev.1.bt2');

$expectedMD5sum="f53340fee1bdbd14d0da74565975c29d";
$observedMD5sum=`md5sum $fastaRef.rev.2.bt2`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'bowtie::bowtie2Build - output content rev.2.bt2');

##########################################
##### bowtie::bowtie
##########################################

# input file
my $forwardFastq=$fastqData."irigin1_1.fastq";
my $reverseFastq=$fastqData."irigin1_2.fastq";

# output file
my $samFileOut="irigin.BOWTIE.sam";
my $readGroupLine="irigin";

# execution test
is(bowtie::bowtie($samFileOut,$readGroupLine,$fastaRef,$forwardFastq,$reverseFastq),'1',"bowtie::bwotie - Test for bowtie running");

# expected output test
#Check if files created
@expectedOutput = ('bowtie_log.e','bowtie_log.o','irigin.BOWTIE.sam','referenceIrigin.fasta','referenceIrigin.fasta.1.bt2','referenceIrigin.fasta.1.ebwt','referenceIrigin.fasta.2.bt2','referenceIrigin.fasta.2.ebwt','referenceIrigin.fasta.3.bt2','referenceIrigin.fasta.3.ebwt','referenceIrigin.fasta.4.bt2','referenceIrigin.fasta.4.ebwt','referenceIrigin.fasta.rev.1.bt2','referenceIrigin.fasta.rev.1.ebwt','referenceIrigin.fasta.rev.2.bt2','referenceIrigin.fasta.rev.2.ebwt');
$observedOutput = `ls`;
@observedOutput = split /\n/,$observedOutput;
is_deeply(\@observedOutput,\@expectedOutput,'bowtie::bowtie - Files created');


# expected content test $samFileOut
my $expectedLineNumber = "8 $samFileOut";                        # structure of the ref file for checking
my $observedLineNumber = `wc -l $samFileOut`;                       # structure of the test file for checking
chomp $observedLineNumber;                                          # to separate the structure and the name of file
is($observedLineNumber, $expectedLineNumber, "bowtie::bowtie - output content file sam");               # TEST IF THE STRUCTURE OF THE FILE OUT IS GOOD

###Test for correct file value of bwa sampe
#GREP command result
my $grepResult=`grep -c ">" $samFileOut`;
chomp $grepResult;
is($grepResult,4,'bowtie::bowtie - output grep in file sam');

##########################################
##### bowtie::bowtie2
##########################################

# output file
$samFileOut="irigin.BOWTIE2.sam";


# execution test
is(bowtie::bowtie($samFileOut,$readGroupLine,$fastaRef,$forwardFastq,$reverseFastq),'1',"bowtie::bowtie2 - Test for bowtie2 running");

# expected output test
#Check if files created
@expectedOutput = ('bowtie_log.e','bowtie_log.o','irigin.BOWTIE2.sam','irigin.BOWTIE.sam','referenceIrigin.fasta','referenceIrigin.fasta.1.bt2','referenceIrigin.fasta.1.ebwt','referenceIrigin.fasta.2.bt2','referenceIrigin.fasta.2.ebwt','referenceIrigin.fasta.3.bt2','referenceIrigin.fasta.3.ebwt','referenceIrigin.fasta.4.bt2','referenceIrigin.fasta.4.ebwt','referenceIrigin.fasta.rev.1.bt2','referenceIrigin.fasta.rev.1.ebwt','referenceIrigin.fasta.rev.2.bt2','referenceIrigin.fasta.rev.2.ebwt');
$observedOutput = `ls`;
@observedOutput = split /\n/,$observedOutput;
is_deeply(\@observedOutput,\@expectedOutput,'bowtie::bowtie2 - Files created');


# expected content test $samFileOut
$expectedLineNumber = "8 $samFileOut";                        # structure of the ref file for checking
$observedLineNumber = `wc -l $samFileOut`;                       # structure of the test file for checking
chomp $observedLineNumber;                                          # to separate the structure and the name of file
is($observedLineNumber, $expectedLineNumber, "bowtie::bowtie2 - output content file sam");               # TEST IF THE STRUCTURE OF THE FILE OUT IS GOOD

###Test for correct file value of bwa sampe
#GREP command result
$grepResult=`grep -c ">" $samFileOut`;
chomp $grepResult;
is($grepResult,4,'bowtie::bowtie2 - output grep in file sam');
