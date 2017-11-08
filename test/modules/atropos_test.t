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

#Will test if Atropos works correctly
use strict;
use warnings;

use Test::More 'no_plan'; #Number of tests, to modify if new tests implemented. Can be changed as 'no_plan' instead of tests=>11 .
use Test::Deep;
use Data::Dumper;
use lib qw(../../modules/);

########################################
#Test of the use of atropos modules
########################################
use_ok('localConfig') or exit;
use_ok('atropos') or exit;
can_ok('atropos','execution');

use localConfig;
use atropos;

my $fastqData="$toggle/data/testData/fastq/pairedTwoIndividusIrigin/";

#########################################
#Remove files and directory created by previous test
#########################################
my $testingDir="$toggle/dataTest/atroposTestDir";
my $creatingDirCom="rm -Rf $testingDir ; mkdir -p $testingDir";                                    #Allows to have a working directory for the tests
system($creatingDirCom) and die ("ERROR: $0 : Cannot execute the command $creatingDirCom\n$!\n");

chdir $testingDir or die ("ERROR: $0 : Cannot go into the new directory with the command \"chdir $testingDir\"\n$!\n");


#######################################
#Creating the IndividuSoft.txt file
#######################################
my $creatingCommand="echo \"atropos\nTEST\" > individuSoft.txt";
system($creatingCommand) and die ("ERROR: $0: Cannot create the individuSoft.txt file with the command $creatingCommand \n$!\n");


#######################################
#Cleaning the logs for the test
#######################################
my $cleaningCommand="rm -Rf atropos_TEST_log.*";
system($cleaningCommand) and die ("ERROR: $0: Cannot clean the previous log files for this test with the command $cleaningCommand \n$!\n");


########################################
##### atropos::execution Single
########################################

# input file

my $fastqFile=$fastqData."irigin1_2.fastq";

# output file
my $fastqFileOut = "irigin1_2.ATROPOS.fastq";                   # Output file without adaptators sequences

# execution test
my %optionsHachees = (
						"-O" => 10,
						"-m" => 35,
						"-q" => "20",
						"--overlap" => 7,
						"-b" => "GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG -b GTTCGTCTTCTGCCGTATGCTCTAGCACTACACTGACCTCAAGTCTGCACACGAGAAGGCTAG",
					);        # Hash containing informations
my $optionsHachees = \%optionsHachees;
is ((atropos::execution($fastqFile,$fastqFileOut,undef, undef, $optionsHachees)),1, 'atropos::execution Single');    # TEST IF FONCTION WORKS

# expected output test
my $observedOutput = `ls`;
my @observedOutput = split /\n/,$observedOutput;
my @expectedOutput = ('atropos_TEST_log.e','atropos_TEST_log.o','individuSoft.txt','irigin1_2.ATROPOS.fastq');

is_deeply(\@observedOutput,\@expectedOutput,'atropos::execution Single - output list');

# expected content test
my $expectedMD5sum = "aeb4479faca4ce706f82ef6d752e89d1";                                            # structure of the ref file for checking
my $observedMD5sum = `md5sum $fastqFileOut`;                                                        # structure of the test file for checking
my @withoutName = split (" ", $observedMD5sum);                                                     # to separate the structure and the name of file
$observedMD5sum = $withoutName[0];     										                        # just to have the md5sum result
is($observedMD5sum, $expectedMD5sum, "atropos::execution Single - output content");               # TEST IF THE STRUCTURE OF THE FILE OUT IS GOOD
##############################

########################################
##### atropos::execution Paired
########################################

# input file
my $forwardFastq=$fastqData."irigin1_1.fastq";
my $reverseFastq=$fastqData."irigin1_2.fastq";

# output file
my $fastqFileOut1 = "irigin1_1.ATROPOS.fastq";                   # Output file without adaptators sequences
my $fastqFileOut2 = "irigin1_2.ATROPOS.fastq";                   # Output file without adaptators sequences

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
is ((atropos::execution($forwardFastq,$fastqFileOut1, $reverseFastq, $fastqFileOut2, $optionsHachees)),1, 'atropos::execution Paired');    # TEST IF FONCTION WORKS

# expected output test
$observedOutput = `ls`;
@observedOutput = split /\n/,$observedOutput;
@expectedOutput = ('atropos_TEST_log.e','atropos_TEST_log.o','individuSoft.txt','irigin1_1.ATROPOS.fastq','irigin1_2.ATROPOS.fastq');

is_deeply(\@observedOutput,\@expectedOutput,'atropos::execution Paired - output list');

# expected content test
$expectedMD5sum = "0307339b270fe6d67a940938cb4f3990";                                            # structure of the ref file for checking
$observedMD5sum = `md5sum $fastqFileOut1`;                                                        # structure of the test file for checking
@withoutName = split (" ", $observedMD5sum);                                                     # to separate the structure and the name of file
$observedMD5sum = $withoutName[0];     										                        # just to have the md5sum result
is($observedMD5sum, $expectedMD5sum, "atropos::execution Paired - output content file 1");               # TEST IF THE STRUCTURE OF THE FILE OUT IS GOOD

# expected content test
$expectedMD5sum = "5b1528f6d369af2726326e1eeebf3f3f";                                            # structure of the ref file for checking
$observedMD5sum = `md5sum $fastqFileOut2`;                                                        # structure of the test file for checking
@withoutName = split (" ", $observedMD5sum);                                                     # to separate the structure and the name of file
$observedMD5sum = $withoutName[0];     										                        # just to have the md5sum result
is($observedMD5sum, $expectedMD5sum, "atropos::execution Paired - output content file 2");
##############################

exit;
