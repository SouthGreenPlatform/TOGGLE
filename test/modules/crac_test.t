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
##### crac::index
##########################################
my $bankData="$toggle/data/Bank/";
my $fastqData="$toggle/data/testData/fastq/pairedTwoIndividusIrigin/";


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
my @expectedOutput = ("crac_log.e","crac_log.o","referenceIrigin.fasta","referenceIrigin.fasta.CRAC.index.conf","referenceIrigin.fasta.CRAC.index.ssa");
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
@expectedOutput = ("crac_log.e","crac_log.o","irigin.CRAC.sam","referenceIrigin.fasta","referenceIrigin.fasta.CRAC.index.conf","referenceIrigin.fasta.CRAC.index.ssa");
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
