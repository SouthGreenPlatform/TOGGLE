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

#Will test if pairing.pm works correctly
use strict;
use warnings;
use Test::More 'no_plan'; #Number of tests, to modify if new tests implemented. Can be changed as 'no_plan' instead of tests=>11 .
use Test::Deep;
use Data::Dumper;
use lib qw(../../modules/);


########################################
#Test of the use of pairing modules
########################################
use_ok('localConfig') or exit;
use_ok('pairing') or exit;

can_ok('pairing','pairRecognition');
can_ok('pairing','createDirPerCouple');
can_ok('pairing','repairing');
can_ok('pairing','extractName');

use localConfig;
use pairing;

my $expectedData="$toggle/data/expectedData/";


#########################################
#Remove files and directory created by previous test
#########################################
my $testingDir="$toggle/dataTest/pairingTestDir";
my $creatingDirCom="rm -Rf $testingDir ; mkdir -p $testingDir";                                    #Allows to have a working directory for the tests
system($creatingDirCom) and die ("ERROR: $0 : Cannot execute the command $creatingDirCom\n$!\n");

chdir $testingDir or die ("ERROR: $0 : Cannot go into the new directory with the command \"chdir $testingDir\"\n$!\n");

my $makeDirCmd = "mkdir pairingDir";
system ($makeDirCmd) and die ("ERROR: $0 : Cannot create the new directory with the command $makeDirCmd\n$!\n");


#######################################
#Creating the IndividuSoft.txt file
#######################################
my $creatingCommand="echo \"pairing\nTEST\" > individuSoft.txt";
system($creatingCommand) and die ("ERROR: $0: Cannot create the individuSoft.txt file with the command $creatingCommand \n$!\n");


#######################################
#Cleaning the logs for the test
#######################################
my $cleaningCommand="rm -Rf pairing_TEST_log.*";
system($cleaningCommand) and die ("ERROR: $0: Cannot clean the previous log files for this test with the command $cleaningCommand \n$!\n");

########################################
##### pairing::extractName
########################################

my $expectName1=("RC3_1");
my $expectRG1=("RC3");

my ($obsName1, $obsRG1)=pairing::extractName('RC3_1.fastq');
is_deeply($obsName1,$expectName1,'pairing::extractName - individu RC3_1');
is_deeply($obsRG1,$expectRG1,'pairing::extractName - RG RC3');

my $expectName2=("RC3_2");
my $expectRG2=("RC3");

my ($obsName2, $obsRG2)=pairing::extractName('RC3_2.fastq');
is_deeply($obsName2,$expectName2,'pairing::extractName - individu RC3_2');
is_deeply($obsRG2,$expectRG2,'pairing::extractName - RG RC3');


########################################
##### pairing::pairRecognition
########################################

# input file
my $checkFastq = 1;
my $fastqFile1 = $expectedData."RC*_1.fastq";     # fastq file 
my $copyCmd= "cp $fastqFile1 pairingDir";           # command to copy the original fastq file into the test directory
system ($copyCmd) and die ("ERROR: $0 : Cannot copy the file $fastqFile1 in the test directory with the command $copyCmd\n$!\n");    # RUN the copy command

my $fastqFile2 = $expectedData."RC*_2.fastq";     # fastq file 
$copyCmd= "cp $fastqFile2 pairingDir";           # command to copy the original fastq file into the test directory
system ($copyCmd) and die ("ERROR: $0 : Cannot copy the file $fastqFile2 in the test directory with the command $copyCmd\n$!\n");    # RUN the copy command


my $expectedOutput={
          '@H3:C39R6ACXX:3:1101:1215:1877' => {
                                                         'ReadGroup' => 'RC1',
                                                         'forward' => 'pairingDir/RC1_1.fastq',
                                                         'reverse' => 'pairingDir/RC1_2.fastq'
                                                       },
          '@H2:C381HACXX:5:1101:1359:1908' => {
                                                 'ReadGroup' => 'RC3',
                                                 'forward' => 'pairingDir/RC3_1.fastq',
                                                 'reverse' => 'pairingDir/RC3_2.fastq'
                                               },
          '@H3:C39R6ACXX:3:1101:1192:1848' => {
                                                   'ReadGroup' => 'RC2',
                                                   'forward' => 'pairingDir/RC2_1.fastq',
                                                } 
        };

my $observedOutput=pairing::pairRecognition("pairingDir",$checkFastq);
##DEBUG print "pairRecognition Expected :\n"; print Dumper ($expectedOutput);print "pairRecognition Observed:\n"; print Dumper ($observedoutput);
is_deeply($observedOutput,$expectedOutput,'pairing::pairRecognition - output list');


########################################
##### pairing::createDirPerCouple
########################################

my $checkValue3=pairing::createDirPerCouple($observedOutput,"pairingDir");
is ($checkValue3,1,'pairing::createDirPerCouple - running');

# Filetree expected
my $expectedFileTree = "pairingDir/:
RC1
RC2
RC3

pairingDir/RC1:
RC1_1.fastq
RC1_2.fastq

pairingDir/RC2:
RC2_1.fastq

pairingDir/RC3:
RC3_1.fastq
RC3_2.fastq
";

my $observedFileTree = `ls -R pairingDir/`;

##DEBUG print "Expected: \n"; print Dumper ($expectedFileTree);print "Observed: \n"; print Dumper ($observedFileTree);
is_deeply($observedFileTree,$expectedFileTree,'pairing::pairRecognition - Filetree created');

########################################
##### pairing::repairing
########################################

# input file
my $rmDirCmd= "rm -r pairingDir";           # command to copy the original fastq file into the test directory
system ($rmDirCmd) and die ("ERROR: $0 : Cannot removed pairingDir in the test directory with the command $rmDirCmd\n$!\n");    # RUN the rm command

$fastqFile1 = $expectedData."RC3_1.CUTADAPT.fastq";     # fastq file 
$copyCmd= "ln -s $fastqFile1 ./";           # command to copy the original fastq file into the test directory
system ($copyCmd) and die ("ERROR: $0 : Cannot link the file $fastqFile1 in the test directory with the command $copyCmd\n$!\n");    # RUN the copy command

$fastqFile2 = $expectedData."RC3_2.CUTADAPT.fastq";     # fastq file 
$copyCmd= "ln -s $fastqFile2 ./";           # command to copy the original fastq file into the test directory
system ($copyCmd) and die ("ERROR: $0 : Cannot link the file $fastqFile2 in the test directory with the command $copyCmd\n$!\n");    # RUN the copy command


#Check if running
my $checkValue=pairing::repairing('RC3_1.CUTADAPT.fastq','RC3_2.CUTADAPT.fastq',".",$checkFastq);
is ($checkValue,'1','pairing::repairing - running');

#Check if files created
$expectedFileTree = ".:
individuSoft.txt
pairing_TEST_log.e
pairing_TEST_log.o
RC3_1.CUTADAPT.fastq
RC3_1.REPAIRING.fastq
RC3_2.CUTADAPT.fastq
RC3_2.REPAIRING.fastq
RC3_Single

./RC3_Single:
RC3Single.fastq
";

$observedFileTree = `ls -R`;

##DEBUG print "Expected: \n"; print Dumper ($expectedFileTree);print "Observed: \n"; print Dumper ($observedFileTree);
is_deeply($observedFileTree,$expectedFileTree,'pairing::pairRecognition - Filetree created');

#Check if working
my $numberOfLinesObserved=`wc -l RC3_Single/RC3Single.fastq`;
chomp $numberOfLinesObserved;
is ($numberOfLinesObserved,'8 '.'RC3_Single/RC3Single.fastq','pairing::repairing - single file');

#Check if the files created are the same as planned
my $diffForward=`diff -q RC3_1.REPAIRING.fastq $toggle/data/expectedData/RC3_1.REPAIRING.fastq`;
is ($diffForward,'','pairing::repairing - forward file');

#Check if the files created are the same as planned
my $diffReverse=`diff -q RC3_2.REPAIRING.fastq $toggle/data/expectedData/RC3_2.REPAIRING.fastq`;
is ($diffReverse,'','pairing::repairing - reverse file');

