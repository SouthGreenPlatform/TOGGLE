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

#Will test if Cutadapt works correctly
use strict;
use warnings;
use Test::More 'no_plan'; #Number of tests, to modify if new tests implemented. Can be changed as 'no_plan' instead of tests=>11 .
use Test::Deep;
use Data::Dumper;
use lib qw(../Modules/);


########################################
#Test of the use of cutadapt modules
########################################
use_ok('toolbox') or exit;
use_ok('cutadapt') or exit;

can_ok('cutadapt','execution');

use toolbox;
use cutadapt;

my $expectedData="../../DATA/expectedData/";

#########################################
#Remove files and directory created by previous test
#########################################
my $testingDir="../DATA-TEST/cutadaptTestDir";
my $creatingDirCom="rm -Rf $testingDir ; mkdir -p $testingDir";                                    #Allows to have a working directory for the tests
system($creatingDirCom) and die ("ERROR: $0 : Cannot execute the command $creatingDirCom\n$!\n");

chdir $testingDir or die ("ERROR: $0 : Cannot go into the new directory with the command \"chdir $testingDir\"\n$!\n");


#######################################
#Creating the IndividuSoft.txt file
#######################################
my $creatingCommand="echo \"cutadapt\nTEST\" > individuSoft.txt";
system($creatingCommand) and die ("ERROR: $0: Cannot create the individuSoft.txt file with the command $creatingCommand \n$!\n");


#######################################
#Cleaning the logs for the test
#######################################
my $cleaningCommand="rm -Rf cutadapt_TEST_log.*";
system($cleaningCommand) and die ("ERROR: $0: Cannot clean the previous log files for this test with the command $cleaningCommand \n$!\n");


########################################
##### cutadapt::execution Single
########################################

# input file
my $fastqFile = $expectedData."RC3_2.fastq";     		# fastq file

# output file
my $fastqFileOut = "RC3_2.CUTADAPT.fastq";                   # Output file without adaptators sequences

# execution test
my %optionsHachees = (
						"-O" => 10,
						"-m" => 35,
						"-q" => "20",
						"--overlap" => 7,
						"-b" => "GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG -b GTTCGTCTTCTGCCGTATGCTCTAGCACTACACTGACCTCAAGTCTGCACACGAGAAGGCTAG",	
					);        # Hash containing informations
my $optionsHachees = \%optionsHachees;   
is ((cutadapt::execution($fastqFile,$fastqFileOut,undef, undef, $optionsHachees)),1, 'cutadapt::execution Single');    # TEST IF FONCTION WORKS

# expected output test
my $observedOutput = `ls`;
my @observedOutput = split /\n/,$observedOutput;
my @expectedOutput = ('cutadapt_TEST_log.e','cutadapt_TEST_log.o','individuSoft.txt','RC3_2.CUTADAPT.fastq');

is_deeply(\@observedOutput,\@expectedOutput,'cutadapt::execution Single - output list');

# expected content test
my $expectedMD5sum = "5a74f7c91eb01b46c4002a584d48bf6f";                                            # structure of the ref file for checking
my $observedMD5sum = `md5sum $fastqFileOut`;                                                        # structure of the test file for checking
my @withoutName = split (" ", $observedMD5sum);                                                     # to separate the structure and the name of file
$observedMD5sum = $withoutName[0];     										                        # just to have the md5sum result
is($observedMD5sum, $expectedMD5sum, "cutadapt::execution Single - output content");               # TEST IF THE STRUCTURE OF THE FILE OUT IS GOOD
##############################

########################################
##### cutadapt::execution Paired
########################################

# input file
my $fastqFile1 = $expectedData."RC3_1.fastq";     		# fastq file
my $fastqFile2 = $expectedData."RC3_2.fastq";     		# fastq file

# output file
my $fastqFileOut1 = "RC3_1.CUTADAPT.fastq";                   # Output file without adaptators sequences
my $fastqFileOut2 = "RC3_2.CUTADAPT.fastq";                   # Output file without adaptators sequences

# execution test
%optionsHachees = (
						"-O" => 10,
						"-m" => 35,
						"-q" => "20,20",
						"--overlap" => 7,
						"-b" => "GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG -b GTTCGTCTTCTGCCGTATGCTCTAGCACTACACTGACCTCAAGTCTGCACACGAGAAGGCTAG",
						"-B" => "GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG -B GTTCGTCTTCTGCCGTATGCTCTAGCACTACACTGACCTCAAGTCTGCACACGAGAAGGCTAG"	
					);          # Hash containing informations
$optionsHachees = \%optionsHachees;   
is ((cutadapt::execution($fastqFile1,$fastqFileOut1, $fastqFile2, $fastqFileOut2, $optionsHachees)),1, 'cutadapt::execution Paired');    # TEST IF FONCTION WORKS

# expected output test
$observedOutput = `ls`;
@observedOutput = split /\n/,$observedOutput;
@expectedOutput = ('cutadapt_TEST_log.e','cutadapt_TEST_log.o','individuSoft.txt','RC3_1.CUTADAPT.fastq','RC3_2.CUTADAPT.fastq');

is_deeply(\@observedOutput,\@expectedOutput,'cutadapt::execution Paired - output list');

# expected content test
$expectedMD5sum = "9055c369fed4d016212caca9d750f6b6";                                            # structure of the ref file for checking
$observedMD5sum = `md5sum $fastqFileOut1`;                                                        # structure of the test file for checking
@withoutName = split (" ", $observedMD5sum);                                                     # to separate the structure and the name of file
$observedMD5sum = $withoutName[0];     										                        # just to have the md5sum result
is($observedMD5sum, $expectedMD5sum, "cutadapt::execution Paired - output content file 1");               # TEST IF THE STRUCTURE OF THE FILE OUT IS GOOD

# expected content test
$expectedMD5sum = "548103973bb70ca479b08e529d05bdcb";                                            # structure of the ref file for checking
$observedMD5sum = `md5sum $fastqFileOut2`;                                                        # structure of the test file for checking
@withoutName = split (" ", $observedMD5sum);                                                     # to separate the structure and the name of file
$observedMD5sum = $withoutName[0];     										                        # just to have the md5sum result
is($observedMD5sum, $expectedMD5sum, "cutadapt::execution Paired - output content file 2");  
##############################

exit;