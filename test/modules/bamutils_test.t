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

#Will test if samTools module work correctly works correctly
use strict;
use warnings;

use Test::More 'no_plan'; #Number of tests, to modify if new tests implemented. Can be changed as 'no_plan' instead of tests=>11 .
use Test::Deep;
use Data::Dumper;
use lib qw(../../modules/);

########################################
#use of bamuttils modules ok
########################################
use_ok('localConfig') or exit;
use_ok('bamutils') or exit;

can_ok('bamutils','bamutilsTool');

use localConfig;
use bamutils;

#########################################
#Remove files and directory created by previous test
#########################################

my $testingDir="$toggle/dataTest/bamutilsTestDir";
my $cleaningCmd="rm -Rf $testingDir"; 
system ($cleaningCmd) and die ("ERROR: $0 : Cannot remove the previous test directory with the command $cleaningCmd \n$!\n");

########################################
#Creation of test directory
########################################
my $makeDirCmd = "mkdir $testingDir";
system ($makeDirCmd) and die ("ERROR: $0 : Cannot create the new directory with the command $makeDirCmd\n$!\n");
chdir $testingDir or die ("ERROR: $0 : Cannot go into the new directory with the command \"chdir $testingDir\"\n$!\n");

#######################################
#Creating the IndividuSoft.txt file
#######################################
my $creatingCmd="echo \"bamutils\nTEST\" > individuSoft.txt";
system($creatingCmd) and die ("ERROR: $0 : Cannot create the individuSoft.txt file with the command $creatingCmd\n$!\n");

#######################################
#Cleaning the logs for the test
#######################################
$cleaningCmd="rm -Rf bamutils_TEST_log.*";
system($cleaningCmd) and die ("ERROR: $0 : Cannot remove the previous log files with the command $cleaningCmd \n$!\n");




################################################################################################
##bamutils tool
################################################################################################
# tools for bamutils bamutilsFilter
my $toolName = "bamutilsFilter";

# input file
my $bamFileIn = "$toggle/data/testData/samBam/oneSam/RC3-SAMTOOLSVIEW.sam";

#Output file
my $bamFileOut="RC3.$toolName.bam";

my %optionsRef = ("-gte" => 'MAPQ 30');
my $optionsHachees = \%optionsRef; 

#execution test
is(bamutils::bamutilsTool($toolName, $bamFileIn, $bamFileOut, $optionsHachees),1,'bamutils::bamutilsTool - bamutilsFilter');

# expected output test
my $observedOutput = `ls`;
my @observedOutput = split /\n/,$observedOutput;
my @expectedOutput = ('bamutils_TEST_log.e','bamutils_TEST_log.o','individuSoft.txt','RC3.bamutilsFilter.bam');

is_deeply(\@observedOutput,\@expectedOutput,'bamutils::bamutilsTool - bamutilsFilter - output list');

# expected output structure
my $expectedMD5sum = "667002537d89d3a1b79d4628085fdbb6";
my $observedMD5sum=`md5sum $bamFileOut`;# structure of the test file
my @withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'bamutils::bamutilsTool - bamutilsFilter - output structure');



# tools for bamutils bamutilstobed
 $toolName = "bamutilstobed";
 
# input file
$bamFileIn = "$toggle/data/testData/samBam/oneBamUnsorted/unsorted.bam";

#Output file
$bamFileOut="RC3.$toolName.bed";

%optionsRef = ();
$optionsHachees = \%optionsRef; 

#execution test
is(bamutils::bamutilsTool($toolName, $bamFileIn, $bamFileOut, $optionsHachees),1,'bamutils::bamutilsTool - bamutilstobed');

# expected output test
$observedOutput = `ls`;
@observedOutput = split /\n/,$observedOutput;
@expectedOutput = ('bamutils_TEST_log.e','bamutils_TEST_log.o','individuSoft.txt','RC3.bamutilsFilter.bam','RC3.bamutilstobed.bed');

is_deeply(\@observedOutput,\@expectedOutput,'bamutils::bamutilsTool - bamutilstobed - output list');

# expected output structure
$expectedMD5sum = "1cdaccadf9fa80dd8822924685ec1291";
$observedMD5sum=`md5sum $bamFileOut`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'bamutils::bamutilsTool - bamutilstobed - output structure');

exit;
