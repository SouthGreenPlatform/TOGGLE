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
### Test for hisat2::hisat2Build
##########################################
my $bankData="$toggle/data/Bank/";
my $fastqData="$toggle/data/testData/rnaseq/pairedOneIndividu/";

# input file

my $fastaRefIni=$bankData."/referenceRnaseq.fa";
my $fastaRef="referenceRnaseq.fa";

#copy fasta reference into test directory where the index will be created
my $copyCommand="cp $fastaRefIni ./$fastaRef";
system ($copyCommand) and die "ERROR: $0: Cannot copy the $fastaRefIni file with the command $copyCommand \n$!\n";

# execution test
is(hisat2::hisat2Build($fastaRef),$fastaRef,'hisat2::hisat2-build');

# expected output test
my $observedOutput = `ls`;
my @observedOutput = split /\n/,$observedOutput;
my @expectedOutput = ('hisat2_log.e','hisat2_log.o','referenceRnaseq.fa','referenceRnaseq.fa.1.ht2','referenceRnaseq.fa.2.ht2','referenceRnaseq.fa.3.ht2','referenceRnaseq.fa.4.ht2','referenceRnaseq.fa.5.ht2','referenceRnaseq.fa.6.ht2','referenceRnaseq.fa.7.ht2','referenceRnaseq.fa.8.ht2');
is_deeply(\@observedOutput,\@expectedOutput,'hisat2::hisat2Build - output list');

# expected output content
my $expectedMD5sum="54c495bbba42cd61c1610442971982c0";
my $observedMD5sum=`md5sum $fastaRef.1.ht2`;# structure of the test file
my @withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'hisat2::hisat2Build - output content 1.ht2');

$expectedMD5sum="4ad45a523ecaef6310bb6f7f608eb311";
$observedMD5sum=`md5sum $fastaRef.2.ht2`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'hisat2::hisat2Build - output content 2.ht2');

$expectedMD5sum="8aa5c56a0ba0b0ab7e9e7f3fb7ee4a76";
$observedMD5sum=`md5sum $fastaRef.3.ht2`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'hisat2::hisat2Build - output content 3.ht2');

$expectedMD5sum="a4ebbf39ff457e410253b571ee79088d";
$observedMD5sum=`md5sum $fastaRef.4.ht2`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'hisat2::hisat2Build - output content 4.ht2');

##########################################
##### hisat2::hisat2
##########################################
# input file
my $forwardFastq=$fastqData."RNASeq_1.fastq";
my $reverseFastq=$fastqData."RNASeq_2.fastq";

# output file
my $samFileOut="RNASeq.HISAT2.sam";
my $readGroupLine="RNASeq";

# execution test
is(hisat2::hisat2($samFileOut,$readGroupLine,$fastaRef,$forwardFastq,$reverseFastq),'1',"hisat2::hisat2 - Test for hisat2 running");

# expected output test
#Check if files created
@expectedOutput = ('hisat2_log.e','hisat2_log.o','referenceRnaseq.fa','referenceRnaseq.fa.1.ht2','referenceRnaseq.fa.2.ht2','referenceRnaseq.fa.3.ht2','referenceRnaseq.fa.4.ht2','referenceRnaseq.fa.5.ht2','referenceRnaseq.fa.6.ht2','referenceRnaseq.fa.7.ht2','referenceRnaseq.fa.8.ht2','RNASeq.HISAT2.sam',);
$observedOutput = `ls`;
@observedOutput = split /\n/,$observedOutput;
is_deeply(\@observedOutput,\@expectedOutput,'hisat2::hisat2 - Files created');


# expected content test $samFileOut
my $expectedLineNumber = "2020 $samFileOut";                        # structure of the ref file for checking
my $observedLineNumber = `wc -l $samFileOut`;                       # structure of the test file for checking
chomp $observedLineNumber;                                          # to separate the structure and the name of file
is($observedLineNumber, $expectedLineNumber, "hisat2::hisat2 - output content file sam");               # TEST IF THE STRUCTURE OF THE FILE OUT IS GOOD

###Test for correct file value of bwa sampe
#GREP command result
my $grepResult=`grep -c ">" $samFileOut`;
chomp $grepResult;
is($grepResult,0,'hisat2::hisat2 - output grep in file sam');
