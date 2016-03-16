#!/usr/bin/perl

###################################################################################################################################
#
# Copyright 2014-2015 IRD-CIRAD-INRA-ADNid
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
# Version 3 written by Cecile Monat, Christine Tranchant, Cedric Farcy, Maryline Summo, Julie Orjuela-Bouniol, Sebastien Ravel, Gautier Sarah, and Francois Sabot
#
###################################################################################################################################

#Will test if fastq_utils works correctly
use strict;
use warnings;
use Data::Translate;#To convert ASCII to decimal
use Data::Dumper;
use Test::More 'no_plan'; #Number of tests, to modify if new tests implemented. Can be changed as 'no_plan' instead of tests=>11 .
use lib qw(../Modules);
use toolbox;

########################################
#use of fastq_utils module ok
########################################

use_ok('fastqUtils');
can_ok('fastqUtils','changeEncode');
can_ok('fastqUtils','checkEncodeByASCIIcontrol');
#can_ok('fastqUtils','checkNumberByWC');
can_ok('fastqUtils','convertLinePHRED33ToPHRED64');
can_ok('fastqUtils','convertLinePHRED64ToPHRED33');

use fastqUtils;

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

########################################
#initialisation and setting configs
########################################
my $testingDir="../DATA-TEST/fastqUtilsTestDir";
my $creatingDirCom="rm -Rf $testingDir ; mkdir -p $testingDir";                                    #Allows to have a working directory for the tests
system($creatingDirCom) and die ("ERROR: $0 : Cannot execute the command $creatingDirCom\n$!\n");


########################################
#Input files
########################################
my $originalFastaRef="../DATA/expectedData/Reference.fasta";
my $fastaRef="$testingDir/Reference.fasta";
my $refCopyCom="cp $originalFastaRef $fastaRef";
system($refCopyCom) and die ("ERROR: $0 : Cannot copy the Reference for test with the command $refCopyCom \n$!\n");

my $originalFastqFile1="../DATA/expectedData/RC3_1.fastq";
my $originalFastqFile2="../DATA/expectedData/RC3_2.fastq";
my $fastqFile1="$testingDir/RC3_1.fastq";
my $fastqFile2="$testingDir/RC3_2.fastq";
my $seqCopyCom="cp $originalFastqFile1 $fastqFile1";
system($seqCopyCom) and die ("ERROR: $0 : Cannot copy the file $originalFastqFile1 for test with the command $seqCopyCom \n$!\n");    #The sequences are copied for testing
$seqCopyCom="cp $originalFastqFile2 $fastqFile2";
system($seqCopyCom) and die ("ERROR: $0 : Cannot copy the file $originalFastqFile2 for test with the command $seqCopyCom \n$!\n");    #The sequences are copied for testing



#########################################
#Sequence count test
#########################################
#checkNumberByWC test with a fastq file
#my $count = fastqUtils::checkNumberByWC($fastqFile1);
#is($count,'1000',"Test for checkNumberByWC... $fastqFile1");

#$count = fastqUtils::checkNumberByWC($fastqFile2);
#is($count,'1000',"Test for checkNumberByWC... $fastqFile2");


#########################################
# Encode conversion test
#########################################

#########################################
# General mode, ie for a complete file, using changeEncode
#########################################
#Correct files for verifying the outputs
my $CorrectPHRED33file="../DATA/expectedData/RC1_1.fastq";
my $CorrectPHRED64file="../DATA/expectedData/RC1_1.ILLUMINA.fastq";
#generate Test data file
my $fastqFileOut33 = "$testingDir/RC1_1.SANGER.fastq";
my $formatInit = 64;
my $formatOut = 33;
is (fastqUtils::changeEncode($CorrectPHRED64file,$fastqFileOut33,$formatInit,$formatOut),'1','Test for changeEncode 64 to 33 running'); # Will verify that the changeEncode will run correctly from 64 to 33 on a file
#Testing the output content

my $diffResult = `diff $CorrectPHRED33file $fastqFileOut33`;
chomp $diffResult;
is ($diffResult,'','Test for changeEncode 64 to 33 output content');


my $fastqFileOut64 = "$testingDir/RC1_1.ILLUMINA.fastq";
$formatInit = 33;
$formatOut = 64;
#Testing the running for 33 to 64
is (fastqUtils::changeEncode($fastqFileOut33,$fastqFileOut64,$formatInit,$formatOut),'1','Test for changeEncode 64 to 33 running'); # Will verify that the changeEncode will run correctly from 33 to 64 on a file
#Testing the output content
my $diffResult2 = `diff $CorrectPHRED64file $fastqFileOut64`;
chomp $diffResult2;
is ($diffResult2,'','Test for changeEncode 33 to 64 output content');



#########################################
# Line by Line mode, using the convertLinePHRED64ToPHRED33 and convertLinePHRED33ToPHRED64
#########################################
#Generic values for encoding in 64 and 33
my $PHRED64Data="dddd`ddddd`dd\\dd^bbdaa`b\\dddd`abdbTbaYabBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB";
my $PHRED33Data="EEEEAEEEEEAEE=EE?CCEBBAC=EEEEABCEC5CB:BC####################################";

#For PHRED64 to PHRED33
is (fastqUtils::convertLinePHRED64ToPHRED33($PHRED64Data),$PHRED33Data,'Test for changing encoding in PHRED64 to PHRED33'); # Will verify that the re-encoded quality from 64 to 33 if correct

#For PHRED33 to PHRED64
is (fastqUtils::convertLinePHRED33ToPHRED64($PHRED33Data),$PHRED64Data,'Test for changing encoding in PHRED33 to PHRED64'); # Will verify that the re-encoded quality from 33 to 64 if correct

#########################################
# Encode check test general
#########################################

#checkEncodeByASCIIcontrol test with a fastq file
is (fastqUtils::checkEncodeByASCIIcontrol($fastqFileOut33),'1','Test for checkEncode on a correct file'); # Will verify that the file is a correct PHRED33 encoded one
is (fastqUtils::checkEncodeByASCIIcontrol($fastqFileOut64),'0','Test for checkEncode on a non correct file');# Will verify that the file is not a correct PHRED33 encoded one

#########################################
# Cleaning
#########################################
$cleaningCommand="rm -Rf $testingDir/* individuSoft.txt fastqUtils_TEST_log.* && rm -Rf $testingDir";
system ($cleaningCommand) and die ("ERROR: $0: Cannot erase the testing files with the command $cleaningCommand \n$!\n");

