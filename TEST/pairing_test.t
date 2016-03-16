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

#Will test if pairing.pm works correctly

use strict;
use warnings;
use Test::More  'no_plan';
use Data::Dumper;

use lib qw(../Modules/);


########################################
#use of pairing module ok
########################################

use_ok('pairing');

use pairing;

#########################################
#Remove files and directory created by previous test
#########################################
my $testingDir="../DATA-TEST/pairingTestDir";
my $cleaningCmd="rm -Rf $testingDir"; 
system ($cleaningCmd) and die ("ERROR: $0 : Cannot remove the previous test directory with the command $cleaningCmd \n$!\n");

my $expectedData="../../DATA/expectedData/";

########################################
#Creation of test directory
########################################
my $makeDirCmd = "mkdir $testingDir";
system ($makeDirCmd) and die ("ERROR: $0 : Cannot create the new directory with the command $makeDirCmd\n$!\n");
chdir $testingDir or die ("ERROR: $0 : Cannot go into the new directory with the command \"chdir $testingDir\"\n$!\n");

$makeDirCmd = "mkdir pairingDir";
system ($makeDirCmd) and die ("ERROR: $0 : Cannot create the new directory with the command $makeDirCmd\n$!\n");

#######################################
#Creating the IndividuSoft.txt file
#######################################
my $creatingCmd="echo \"pairing\nTEST\" > individuSoft.txt";
system($creatingCmd) and die ("ERROR: $0 : Cannot create the individuSoft.txt file with the command $creatingCmd\n$!\n");

#######################################
#Cleaning the logs for the test
#######################################
$cleaningCmd="rm -Rf pairing_TEST_log.*";
system($cleaningCmd) and die ("ERROR: $0 : Cannot remove the previous log files with the command $cleaningCmd \n$!\n");


########################################
#Input files
########################################
my $originalFastqcFile = $expectedData."RC*_1.fastq";     # fastq file 
my $copyCmd= "cp $originalFastqcFile pairingDir";           # command to copy the original fastq file into the test directory
system ($copyCmd) and die ("ERROR: $0 : Cannot link the file $originalFastqcFile in the test directory with the command $copyCmd\n$!\n");    # RUN the copy command

$originalFastqcFile = $expectedData."RC*_2.fastq";     # fastq file 
$copyCmd= "cp $originalFastqcFile pairingDir";           # command to copy the original fastq file into the test directory
system ($copyCmd) and die ("ERROR: $0 : Cannot link the file $originalFastqcFile in the test directory with the command $copyCmd\n$!\n");    # RUN the copy command

$originalFastqcFile = $expectedData."/*.CUTADAPT.fastq";     # fastq file 
my $lnCmd= "ln -s $originalFastqcFile .";           # command to copy the original fastq file into the test directory
system ($lnCmd) and die ("ERROR: $0 : Cannot copy the file $originalFastqcFile in the test directory with the command $lnCmd\n$!\n");    # RUN the copy command


########################################
#extractName
########################################
my $expectName1=("RC3_1");
my $expectRG1=("RC3");

my ($obsName1, $obsRG1)=pairing::extractName('RC3_1.fastq');
is_deeply($obsName1,$expectName1,'pairing::extractName... individu RC3_1');
is_deeply($obsRG1,$expectRG1,'pairing::extractName... RG RC3');

my $expectName2=("RC3_2");
my $expectRG2=("RC3");

my ($obsName2, $obsRG2)=pairing::extractName('RC3_2.fastq');
is_deeply($obsName2,$expectName2,'pairing::extractName... individu RC3_2');
is_deeply($obsRG2,$expectRG2,'pairing::extractName... RG RC3');


########################################
#pairRecognition
########################################
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

my $observedoutput=pairing::pairRecognition("pairingDir");
##DEBUG print "pairRecognition Expected :\n"; print Dumper ($expectedOutput);print "pairRecognition Observed:\n"; print Dumper ($observedoutput);
is_deeply($observedoutput,$expectedOutput,'pairing::pairRecognition');



#########################################
##createDirPerCouple
#########################################
my $checkValue3=pairing::createDirPerCouple(pairing::pairRecognition("pairingDir"),"pairingDir");
is ($checkValue3,1,'pairing::createDirPerCouple... running');

# Filetree expected
my $expectedFileTree = 
        [
            'pairingDir/RC1:',
            'RC1_1.fastq',
            'RC1_2.fastq',
            '',
            'pairingDir/RC2:',
            'RC2_1.fastq',
            '',
            'pairingDir/RC3:',
            'RC3_1.fastq',
            'RC3_2.fastq'
        ];

my $observedFileTree=toolbox::readDir("pairingDir");
##DEBUG print "Expected: \n"; print Dumper ($expectedFileTree);print "Observed: \n"; print Dumper ($observedFileTree);
is_deeply($observedFileTree,$expectedFileTree,'pairing::pairRecognition... Filetree created');


########################################
#repairing 
########################################

#Check if running
my $checkValue=pairing::repairing('RC3_1.CUTADAPT.fastq','RC3_2.CUTADAPT.fastq',".");
is ($checkValue,'1','pairing::repairing... running');

#Check if working
my $numberOfLinesObserved=`wc -l RC3_Single/RC3Single.fastq`;
chomp $numberOfLinesObserved;
is ($numberOfLinesObserved,'4 '.'RC3_Single/RC3Single.fastq','pairing::repairing... single file');

#Check if the files created are the same as planned
my $diffForward=`diff -q RC3_1.REPAIRING.fastq ../../DATA/expectedData/RC3_1.REPAIRING.fastq`;
is ($diffForward,'','pairing::repairing... forward file');

#Check if the files created are the same as planned
my $diffReverse=`diff -q RC3_2.REPAIRING.fastq ../../DATA/expectedData/RC3_2.REPAIRING.fastq`;
is ($diffReverse,'','pairing::repairing... reverse file');

