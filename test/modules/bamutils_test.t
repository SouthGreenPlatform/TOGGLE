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


################################################################################################
##bamutils tool
################################################################################################
# tools for bamutils bamutilsFilter
my $toolName = "bamutilsFilter";

#########################################
#Creating log file
#########################################
my $logFile=$toolName."_log.o";
my $errorFile=$toolName."_log.e";
system("touch $testDir/$logFile $testDir/$errorFile") and die "\nERROR: $0 : cannot create the log files $logFile and $errorFile: $!\nExiting...\n";

######################################################################################################################################
######################################################################################################################################


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
my @expectedOutput = ('bamutilsFilter_log.e','bamutilsFilter_log.o','RC3.bamutilsFilter.bam');

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
@expectedOutput = ('bamutilsFilter_log.e','bamutilsFilter_log.o','RC3.bamutilsFilter.bam','RC3.bamutilstobed.bed');

is_deeply(\@observedOutput,\@expectedOutput,'bamutils::bamutilsTool - bamutilstobed - output list');

# expected output structure
$expectedMD5sum = "1cdaccadf9fa80dd8822924685ec1291";
$observedMD5sum=`md5sum $bamFileOut`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'bamutils::bamutilsTool - bamutilstobed - output structure');

exit;
