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
my $cleaningCommand="rm -rf tophat_TEST_log.*";
system($cleaningCommand) and die ("ERROR: $0: Cannot clean the previous log files for this test with the command $cleaningCommand \n$!\n");

########################################
#initialisation and setting configs
########################################
my $testingDir="../DATA-TEST/tophatTestDir";
my $creatingDirCom="rm -rf $testingDir ; mkdir -p $testingDir";                                    #Allows to have a working directory for the tests
system($creatingDirCom) and die ("ERROR: $0 : Cannot execute the command $creatingDirCom\n$!\n");

my $OriginalFastaRef="../DATA/expectedData/referenceRNASeq.fa";
my $fastaRef="$testingDir/referenceRNASeq.fa";
my $refCopyCom="cp $OriginalFastaRef $fastaRef";
system($refCopyCom) and die ("ERROR: $0 : Cannot copy the Reference $OriginalFastaRef with the command $refCopyCom\n$!\n");     #Now we have a ref to be tested

my $OriginalGffRef="../DATA/expectedData/referenceRNASeq.gff3";
my $gffRef="$testingDir/referenceRNASeq.gff3";
$refCopyCom="cp $OriginalGffRef $gffRef";
system($refCopyCom) and die ("ERROR: $0 : Cannot copy the gff Reference $OriginalGffRef with the command $refCopyCom\n$!\n");     #Now we have a ref to be tested

my $originalFastqFile1="../DATA/expectedData/RNASeq_1.fastq";
my $originalFastqFile2="../DATA/expectedData/RNASeq_2.fastq";
my $fastqFile1="$testingDir/RNASeq_1.fastq";
my $fastqFile2="$testingDir/RNASeq_2.fastq";
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
can_ok( 'tophat','bowtie2Build');
can_ok( 'tophat','tophat2');

use toolbox;
use tophat;

################################################################################################
###tophat::bowtieBuild
################################################################################################
my %optionsHachees = ();        # Hash containing informations
my $optionHachees = \%optionsHachees;                           # Ref of the hash
my $expectedIndexPrefix=$testingDir."/referenceRNASeq";
is(tophat::bowtieBuild($fastaRef,$optionHachees),$expectedIndexPrefix, 'OK for bowtieBuild RUNNING');

###Checking the correct structure for the output file using md5sum
my $expectedMD5sum="167cd0622bda91392673aaf207255d2b";
my $observedMD5sum=`md5sum $expectedIndexPrefix.1.ebwt`;# structure of the test file
my @withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'Ok for the content of the bowtie build 1.ebwt structure');

$expectedMD5sum="dd30c97b610f5dc53cf6f02123fcf807";
$observedMD5sum=`md5sum $expectedIndexPrefix.2.ebwt`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'Ok for the content of the bowtie build 2.ebwt structure');

$expectedMD5sum="8aa5c56a0ba0b0ab7e9e7f3fb7ee4a76";
$observedMD5sum=`md5sum $expectedIndexPrefix.3.ebwt`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'Ok for the content of the bowtie build 3.ebwt structure');

$expectedMD5sum="a4ebbf39ff457e410253b571ee79088d";
$observedMD5sum=`md5sum $expectedIndexPrefix.4.ebwt`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'Ok for the content of the bowtie build 4.ebwt structure');

$expectedMD5sum="4bd4f23dc8b98a5dc4b56d7f4d89a9b5";
$observedMD5sum=`md5sum $expectedIndexPrefix.rev.1.ebwt`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'Ok for the content of the bowtie build rev.1.ebwt structure');

$expectedMD5sum="619322d189d42f4eaede8aaaedf9890e";
$observedMD5sum=`md5sum $expectedIndexPrefix.rev.2.ebwt`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'Ok for the content of the bowtie build rev.2.ebwt structure');



################################################################################################
###tophat::bowtie2Build
################################################################################################
%optionsHachees = ("-T" => "RealignerTargetCreator","-nt" => "4");        # Hash containing informations
$optionHachees = \%optionsHachees;                           # Ref of the hash
is(tophat::bowtie2Build($fastaRef),$expectedIndexPrefix, 'OK for bowtie2Build RUNNING');

###Checking the correct structure for the output file using md5sum
$expectedMD5sum="4e0329b55cd2a67490ef96b7b2567e5d";
$observedMD5sum=`md5sum $expectedIndexPrefix.1.bt2`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'Ok for the content of the bowtie build 1.bt2 structure');

$expectedMD5sum="4ad45a523ecaef6310bb6f7f608eb311";
$observedMD5sum=`md5sum $expectedIndexPrefix.2.bt2`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'Ok for the content of the bowtie build 2.bt2 structure');

$expectedMD5sum="8aa5c56a0ba0b0ab7e9e7f3fb7ee4a76";
$observedMD5sum=`md5sum $expectedIndexPrefix.3.bt2`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'Ok for the content of the bowtie build 3.bt2 structure');

$expectedMD5sum="a4ebbf39ff457e410253b571ee79088d";
$observedMD5sum=`md5sum $expectedIndexPrefix.4.bt2`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'Ok for the content of the bowtie build 4.bt2 structure');

$expectedMD5sum="beeb0f44030b631af6f182f4a85f045d";
$observedMD5sum=`md5sum $expectedIndexPrefix.rev.1.bt2`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'Ok for the content of the bowtie build rev.1.bt2 structure');

$expectedMD5sum="df753d8d0522ee01e2479f99a04525dd";
$observedMD5sum=`md5sum $expectedIndexPrefix.rev.2.bt2`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'Ok for the content of the bowtie build rev.2.bt2 structure');

################################################################################################
###tophat::tophat2
################################################################################################
%optionsHachees = ("-i" => "30", "-I" => "20000",
                      "-a" => "8",
                      "-m" => "1",
                      "--no-coverage-search" => '',
                      "-g" => "10",
                      "--bowtie-n" => '',
                      "--library-type" => 'fr-secondstrand',
                      "--microexon-search" => '');
$optionHachees = \%optionsHachees;                           # Ref of the hash

is(tophat::tophat2($testingDir, $expectedIndexPrefix, $fastqFile1, $fastqFile2, $gffRef, $optionHachees), 'OK for tophat2 RUNNING');
