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

########################################
##### fastxToolkit::fastxTrimmer with fastq
########################################
my $fastqFile="$toggle/data/testData/fastq/pairedOneIndividuArcad/arcad1_1.fastq";

# output file
my $fastqFileOut = "RC3_1.FASTXTRIMMER.fastq";  

# execution test
my %optionsHachees = ("-f" => "8", "-Q33" => "");                             # Hash containing informations
my $optionHachees = \%optionsHachees;                           # Ref of the hash

is(fastxToolkit::fastxTrimmer($fastqFile, $fastqFileOut, $optionHachees),1, 'fastxToolkit::fastxTrimmer');


# expected output test
my $observedOutput = `ls`;
my @observedOutput = split /\n/,$observedOutput;
my @expectedOutput=('fastxToolkit_log.e','fastxToolkit_log.o','RC3_1.FASTXTRIMMER.fastq');

is_deeply(\@observedOutput,\@expectedOutput,'fastxToolkit::fastxTrimmer - output list');

###Test for correct file value of fastq trimmed file using a md5sum file control
# expected content test
my $expectedMD5sum = "4550ec29721933d7dbd437952c535453";                                            # structure of the ref file for checking
my $observedMD5sum = `md5sum $fastqFileOut`;                                                        # structure of the test file for checking
my @withoutName = split (" ", $observedMD5sum);                                                     # to separate the structure and the name of file
$observedMD5sum = $withoutName[0];                                                     #just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'fastxToolkit::fastxTrimmer - output content');
##############################

########################################
##### fastxToolkit::fastxTrimmer with gz
########################################

# input file
$fastqFile="$toggle/data/testData/fastq/pairedTwoIndividusGzippedIrigin/irigin1_1.fastq.gz";

# output file
$fastqFileOut = "RC3_1.FASTXTRIMMER.fastq.gz";  

# execution test
%optionsHachees = ("-f" => "8", "-Q33" => "");                             # Hash containing informations
$optionHachees = \%optionsHachees;                           # Ref of the hash

is(fastxToolkit::fastxTrimmer($fastqFile, $fastqFileOut, $optionHachees),1, 'fastxToolkit::fastxTrimmer gzip');


# expected output test
$observedOutput = `ls`;
@observedOutput = split /\n/,$observedOutput;
@expectedOutput=('fastxToolkit_log.e','fastxToolkit_log.o','RC3_1.FASTXTRIMMER.fastq','RC3_1.FASTXTRIMMER.fastq.gz');

is_deeply(\@observedOutput,\@expectedOutput,'fastxToolkit::fastxTrimmer - output list gzip');

###Test for correct file value of fastq trimmed file using a md5sum file control
# expected content test
$expectedMD5sum = "\@H3:C39R6ACXX:3:1101:1215:1877/1";                                            # structure of the ref file for checking
$observedMD5sum = `zcat $fastqFileOut | head -n1`;                                                        # structure of the test file for checking

chomp $observedMD5sum;                                                    #just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'fastxToolkit::fastxTrimmer - output content gzip');
##############################



exit;