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
#Input files
########################################

#Correct files for verifying the outputs
my $CorrectPHRED33file="$toggle/data/testData/fastq/singleIndividualPHRED33/RC1_1.SANGER.fastq";
my $CorrectPHRED64file="$toggle/data/testData/fastq/singleIndividualPHRED64/RC1_1.ILLUMINA.fastq";

#########################################
# changeEncode
#########################################

#generate Test data file
my $fastqFileOut33 = "RC1_1.SANGER.fastq";
my $fastqFileOut64 = "RC1_1.ILLUMINA.fastq";

my $formatInit = 64;
my $formatOut = 33;

#execution test 64 to 33
is (fastqUtils::changeEncode($CorrectPHRED64file,$fastqFileOut33,$formatInit,$formatOut),'1','fastqUtils::changeEncode 64 to 33');

# expected output test
my $observedOutput = `ls`;
my @observedOutput = split /\n/,$observedOutput;
my @expectedOutput = ('fastqUtils_log.e','fastqUtils_log.o','RC1_1.SANGER.fastq');

is_deeply(\@observedOutput,\@expectedOutput,'fastqUtils::changeEncode 64 to 33 - output list');

# expected content test
my $diffResult = `diff $CorrectPHRED33file $fastqFileOut33`;
chomp $diffResult;
is ($diffResult,'','fastqUtils::changeEncode 64 to 33 - output structure');


$formatInit = 33;
$formatOut = 64;
#execution test 33 to 64
is (fastqUtils::changeEncode($fastqFileOut33,$fastqFileOut64,$formatInit,$formatOut),'1','fastqUtils::changeEncode 33 to 64');

# expected output test
$observedOutput = `ls`;
@observedOutput = split /\n/,$observedOutput;
@expectedOutput = ('fastqUtils_log.e','fastqUtils_log.o','RC1_1.ILLUMINA.fastq','RC1_1.SANGER.fastq');

is_deeply(\@observedOutput,\@expectedOutput,'fastqUtils::changeEncode 33 to 64 - output list');

# expected content test
my $diffResult2 = `diff $CorrectPHRED64file $fastqFileOut64`;
chomp $diffResult2;
is ($diffResult2,'','fastqUtils::changeEncode 33 to 64 - output content');




#########################################
# convertLinePHRED64ToPHRED33
#########################################
#input/output data
my $PHRED64Data="dddd`ddddd`dd\\dd^bbdaa`b\\dddd`abdbTbaYabBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB";
my $PHRED33Data="EEEEAEEEEEAEE=EE?CCEBBAC=EEEEABCEC5CB:BC####################################";

#execution test
is (fastqUtils::convertLinePHRED64ToPHRED33($PHRED64Data),$PHRED33Data,'fastqUtils::convertLinePHRED64ToPHRED33'); # Will verify that the re-encoded quality from 64 to 33 if correct

#########################################
# convertLinePHRED64ToPHRED33
#####################################

#For PHRED33 to PHRED64
is (fastqUtils::convertLinePHRED33ToPHRED64($PHRED33Data),$PHRED64Data,'fastqUtils::convertLinePHRED33ToPHRED64'); # Will verify that the re-encoded quality from 33 to 64 if correct


#########################################
# fastqUtils::checkEncodeByASCIIcontrol
#########################################

# Will verify that the file is a correct PHRED33 encoded one
is (fastqUtils::checkEncodeByASCIIcontrol($fastqFileOut33),'1','fastqUtils::checkEncode - positive test');

# Will verify that the file is not a correct PHRED33 encoded one
is (fastqUtils::checkEncodeByASCIIcontrol($fastqFileOut64),'0','fastqUtils::checkEncode - negative test');

exit;
__END__
