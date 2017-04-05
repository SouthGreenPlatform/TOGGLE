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

#Will test if fastq_utils works correctly
use strict;
use warnings;
use Data::Translate;#To convert ASCII to decimal
use Data::Dumper;
use Test::More 'no_plan'; #Number of tests, to modify if new tests implemented. Can be changed as 'no_plan' instead of tests=>11 .
use Test::Deep;
use lib qw(../../modules);
use toolbox;

########################################
#use of fastq_utils module ok
########################################

use_ok('localConfig') or exit;
use_ok('fastqUtils') or exit;
can_ok('fastqUtils','changeEncode');
can_ok('fastqUtils','checkEncodeByASCIIcontrol');
can_ok('fastqUtils','convertLinePHRED33ToPHRED64');
can_ok('fastqUtils','convertLinePHRED64ToPHRED33');

use localConfig;
use fastqUtils;

########################################
#initialisation and setting configs
########################################
my $testingDir="$toggle/dataTest/fastqUtilsTestDir";
my $creatingDirCom="rm -Rf $testingDir ; mkdir -p $testingDir";                                    #Allows to have a working directory for the tests
system($creatingDirCom) and die ("ERROR: $0 : Cannot execute the command $creatingDirCom\n$!\n");
chdir $testingDir or die("ERROR: $0 : Cannot change dir to $testingDir\n$!\n");

#######################################
#Creating the IndividuSoft.txt file
#######################################
my $creatingCommand="echo \"fastqUtils\nTEST\" > individuSoft.txt";
system($creatingCommand) and die ("ERROR: $0: Cannot create the individuSoft.txt file with the command $creatingCommand \n$!\n");

#######################################
#Cleaning the logs for the test
#######################################
my $cleaningCommand="rm -Rf fastqUtils_TEST_log.*";
system($cleaningCommand) and die ("ERROR: $0: Cannot clean the previous log files for this test with the command $cleaningCommand \n$!\n");

;


########################################
#Input files
########################################
my $expectedData="$toggle/data/expectedData/";

#Correct files for verifying the outputs
my $CorrectPHRED33file=$expectedData."RC1_1.fastq";
my $CorrectPHRED64file=$expectedData."RC1_1.ILLUMINA.fastq";


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
my @expectedOutput = ('fastqUtils_TEST_log.e','fastqUtils_TEST_log.o','individuSoft.txt','RC1_1.SANGER.fastq');

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
@expectedOutput = ('fastqUtils_TEST_log.e','fastqUtils_TEST_log.o','individuSoft.txt','RC1_1.ILLUMINA.fastq','RC1_1.SANGER.fastq');

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
