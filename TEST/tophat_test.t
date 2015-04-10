#!/usr/bin/perl -w

###################################################################################################################################
#
# Copyright 2014 IRD-CIRAD
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
# Intellectual property belongs to IRD, CIRAD and South Green developpement plateform
# Written by Cecile Monat, Christine Tranchant, Ayite Kougbeadjo, Cedric Farcy, Mawusse Agbessi, Marilyne Summo, and Francois Sabot
#
###################################################################################################################################


use strict;

#Will test if tophat works correctly
use warnings;
use Test::More 'no_plan'; #Number of tests, to modify if new tests implemented. Can be changed as 'no_plan' instead of tests=>11 .
use Test::Deep;
use lib qw(../Modules/);
use Data::Dumper;

#######################################
#Creating the IndividuSoft.txt file
#######################################
my $creatingCommand="echo \"tophat\nTEST\" > individuSoft.txt";
system($creatingCommand) and die ("ERROR: $0: Cannot create the individuSoft.txt file with the command $creatingCommand \n$!\n");


#######################################
#Cleaning the logs for the test
#######################################
my $cleaningCommand="rm -Rf tophat_TEST_log.*";
system($cleaningCommand) and die ("ERROR: $0: Cannot clean the previous log files for this test with the command $cleaningCommand \n$!\n");

########################################
#initialisation and setting configs
########################################
my $testingDir="../DATA-TEST/tophatTestDir";
my $creatingDirCom="rm -Rf $testingDir ; mkdir -p $testingDir";                                    #Allows to have a working directory for the tests
system($creatingDirCom) and die ("ERROR: $0 : Cannot execute the command $creatingDirCom\n$!\n");

my $OriginalFastaRef="../DATA/expectedData/Reference.fasta";
my $fastaRef="$testingDir/Reference.fasta";
my $refCopyCom="cp $OriginalFastaRef $fastaRef";
system($refCopyCom) and die ("ERROR: $0 : Cannot copy the Reference $OriginalFastaRef with the command $refCopyCom\n$!\n");     #Now we have a ref to be tested

my $originalFastqFile1="../DATA/expectedData/RC3_1.REPAIRING.fastq";
my $originalFastqFile2="../DATA/expectedData/RC3_2.REPAIRING.fastq";
my $fastqFile1="$testingDir/RC3_1.REPAIRING.fastq";
my $fastqFile2="$testingDir/RC3_2.REPAIRING.fastq";
my $seqCopyCom1="cp $originalFastqFile1 $fastqFile1";
my $seqCopyCom2="cp $originalFastqFile2 $fastqFile2";
system($seqCopyCom1) and die ("ERROR: $0 : Cannot copy the Fastq file $fastqFile1 for test with the command $seqCopyCom1 \n$!\n");    #The sequences are copied for testing
system($seqCopyCom2) and die ("ERROR: $0 : Cannot copy the Fastq file $fastqFile2 for test with the command $seqCopyCom2 \n$!\n");    #The sequences are copied for testing

########################################
#use of module ok
########################################
use_ok('toolbox') or exit;
use_ok('tophat') or exit;
can_ok( 'tophat','bowtieBuild');

use toolbox;
use tophat;

################################################################################################
###tophat::bowtieBuild
################################################################################################
my %optionsHachees = ('');        # Hash containing informations
my $optionHachees = \%optionsHachees;                           # Ref of the hash
my $expectedIndexPrefix=$testingDir."/Reference";
is(tophat::bowtieBuild($fastaRef,$optionHachees),$expectedIndexPrefix, 'OK for bowtieBuild RUNNING');

###Checking the correct structure for the output file using md5sum
my $expectedMD5sum="de1ef57892bd9f508fb466521bd5a5b6";
my $observedMD5sum=`md5sum $expectedIndexPrefix.1.ebwt`;# structure of the test file
my @withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'Ok for the content of the bowtie build 1.ebwt structure');

$expectedMD5sum="5fe542df841de8685b4ee1c694b52f64";
$observedMD5sum=`md5sum $expectedIndexPrefix.2.ebwt`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'Ok for the content of the bowtie build 2.ebwt structure');

$expectedMD5sum="dc12cca8433dfb22df23bc78bc6aeef6";
$observedMD5sum=`md5sum $expectedIndexPrefix.3.ebwt`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'Ok for the content of the bowtie build 3.ebwt structure');

$expectedMD5sum="3d11892beee30c866ee5e2a06bbbc3d8";
$observedMD5sum=`md5sum $expectedIndexPrefix.4.ebwt`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'Ok for the content of the bowtie build 4.ebwt structure');

$expectedMD5sum="cdf0694f4adfc7c5773f59c234081e98";
$observedMD5sum=`md5sum $expectedIndexPrefix.rev.1.ebwt`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'Ok for the content of the bowtie build rev.1.ebwt structure');

$expectedMD5sum="f55fc9bd3bc5298fb0946289db6cff66";
$observedMD5sum=`md5sum $expectedIndexPrefix.rev.2.ebwt`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'Ok for the content of the bowtie build rev.2.ebwt structure');



exit;

################################################################################################
###tophat::bowtie2Build
################################################################################################
is(tophat::bowtie2Build($fastaRef),$expectedIndexPrefix, 'OK for bowtie2Build RUNNING');

###Checking the correct structure for the output file using md5sum
$expectedMD5sum="b7a6d65d4bcefe2332dcdc8e9c0cb9c1";
$observedMD5sum=`md5sum $expectedIndexPrefix.1.bt2`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'Ok for the content of the bowtie build 1.bt2 structure');

$expectedMD5sum="481b0055258e98825bb4a8c52c3e90c0";
$observedMD5sum=`md5sum $expectedIndexPrefix.2.bt2`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'Ok for the content of the bowtie build 2.bt2 structure');

$expectedMD5sum="dc12cca8433dfb22df23bc78bc6aeef6";
$observedMD5sum=`md5sum $expectedIndexPrefix.3.bt2`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'Ok for the content of the bowtie build 3.bt2 structure');

$expectedMD5sum="3d11892beee30c866ee5e2a06bbbc3d8";
$observedMD5sum=`md5sum $expectedIndexPrefix.4.bt2`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'Ok for the content of the bowtie build 4.bt2 structure');

$expectedMD5sum="c5c13aab02f5bf0d2701b8de21df32ec";
$observedMD5sum=`md5sum $expectedIndexPrefix.rev.1.bt2`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'Ok for the content of the bowtie build rev.1.bt2 structure');

$expectedMD5sum="f53340fee1bdbd14d0da74565975c29d";
$observedMD5sum=`md5sum $expectedIndexPrefix.rev.2.bt2`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'Ok for the content of the bowtie build rev.2.bt2 structure');

################################################################################################
###tophat::tophat2
################################################################################################
$optionsHachees=$configInfos->{'tophat'};
is(tophat::tophat2($testingDir, $expectedIndexPrefix, $fastqFile1, $fastqFile2, $gffFile, , 'OK for bowtie2Build RUNNING');
