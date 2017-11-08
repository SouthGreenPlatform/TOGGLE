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
use lib qw(../../modules/);


########################################
#use of bwa module ok
########################################
use_ok('localConfig') or exit;
use_ok('crac') or exit;
can_ok('crac','cracIndex');
can_ok('crac','crac');

use localConfig;
use crac;

my $bankData="$toggle/data/Bank/";
my $fastqData="$toggle/data/testData/fastq/pairedTwoIndividusIrigin/";

#########################################
#Remove files and directory created by previous test
#########################################
my $testingDir="$toggle/dataTest/cracTestDir";
my $creatingDirCom="rm -Rf $testingDir ; mkdir -p $testingDir";                                    #Allows to have a working directory for the tests
system($creatingDirCom) and die ("ERROR: $0 : Cannot execute the command $creatingDirCom\n$!\n");

chdir $testingDir or die ("ERROR: $0 : Cannot go into the new directory with the command \"chdir $testingDir\"\n$!\n");


#######################################
#Creating the IndividuSoft.txt file
#######################################
my $creatingCommand="echo \"crac\nTEST\" > individuSoft.txt";
system($creatingCommand) and die ("ERROR: $0: Cannot create the individuSoft.txt file with the command $creatingCommand \n$!\n");


#######################################
#Cleaning the logs for the test
#######################################
my $cleaningCommand="rm -Rf crac_TEST_log.*";
system($cleaningCommand) and die ("ERROR: $0: Cannot clean the previous log files for this test with the command $cleaningCommand \n$!\n");


##########################################
##### crac::index
##########################################


# input file
my $fastaRef="referenceIrigin.fasta";

my $originalFastaRef=$bankData."/referenceIrigin.fasta";
my $copyCmd= "cp $originalFastaRef $fastaRef";           # command to copy the original fasta file into the test directory
system ($copyCmd) and die ("ERROR: $0 : Cannot link the file $originalFastaRef in the test directory with the command $copyCmd\n$!\n");    # RUN the copy command

# execution test
my %optionsHachees = ();
my $optionsHachees = \%optionsHachees;

is(crac::cracIndex($fastaRef,$optionsHachees),1,'crac::cracIndex - running');

# expected output test
#Check if files created
my @expectedOutput = ("crac_TEST_log.e","crac_TEST_log.o","individuSoft.txt","referenceIrigin.fasta","referenceIrigin.fasta.CRAC.index.conf","referenceIrigin.fasta.CRAC.index.ssa");
my $observedOutput = `ls`;
my @observedOutput = split /\n/,$observedOutput;
is_deeply(\@observedOutput,\@expectedOutput,'crac::cracIndex - Filetree created');

# expected content test $fastaRefBWT
my $expectedMD5sum = "2827db89a8f57cc6298c18b6eb6f9ac8";                                            # structure of the ref file for checkin
my $observedMD5sum = `md5sum referenceIrigin.fasta.CRAC.index.conf`;                                                        # structure of the test file for checking
my @withoutName = split (" ", $observedMD5sum);                                                     # to separate the structure and the name of file
$observedMD5sum = $withoutName[0];     										                        # just to have the md5sum result
is($observedMD5sum, $expectedMD5sum, "crac::cracIndex - output content file");               # TEST IF THE STRUCTURE OF THE FILE OUT IS GOOD


##########################################
##### crac::crac
##########################################

# input file
my $forwardFastq=$fastqData."irigin1_1.fastq";
my $reverseFastq=$fastqData."irigin1_2.fastq";
my $cracIndex="referenceIrigin.fasta.CRAC.index";
# execution test
%optionsHachees = (
			"-k" => 22
			);        # Hash containing informations
$optionsHachees = \%optionsHachees;

# output file
my $samFileOut="irigin.CRAC.sam";

# execution test
is(crac::crac($samFileOut,$cracIndex,$forwardFastq,$reverseFastq,$optionsHachees),'1',"crac::crac - Test for crac running");

# expected output test
#Check if files created
@expectedOutput = ("crac_TEST_log.e","crac_TEST_log.o","individuSoft.txt","irigin.CRAC.sam","referenceIrigin.fasta","referenceIrigin.fasta.CRAC.index.conf","referenceIrigin.fasta.CRAC.index.ssa");
$observedOutput = `ls ./`;
@observedOutput = split /\n/,$observedOutput;
is_deeply(\@observedOutput,\@expectedOutput,'crac::crac - Files created');

# expected content test $samFileOut
my $expectedLineNumber = "2954 $samFileOut";                                            # structure of the ref file for checking
my $observedLineNumber = `wc -l $samFileOut`;                                                        # structure of the test file for checking
chomp $observedLineNumber;                                                     # to separate the structure and the name of file
is($observedLineNumber, $expectedLineNumber, "crac::crac - output content file sam");               # TEST IF THE STRUCTURE OF THE FILE OUT IS GOOD

###Test for correct file value of crac crac
#GREP command result
my $grepResult=`grep -c "MQ:i:254" $samFileOut`;
chomp $grepResult;
is($grepResult,341,'crac::crac - output grep in file sam');

exit;
