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

use strict;
use warnings;

use Test::More tests => 5;
use lib qw(../Modules/);


########################################
#use of cutadapt module ok
########################################
use_ok('toolbox') or exit;                                                                          # Check if toolbox is usable
use_ok('cutadapt') or exit;                                                                         # Check if cutadapt is usable
can_ok('cutadapt','execution');                                                                     # Check if cutadapt::execution is find

use toolbox;
use cutadapt;


#########################################
#Remove files and directory created by previous test
#########################################
my $testingDir="../DATA-TEST/cutadaptTestDir";
my $cleaningCmd="rm -Rf $testingDir"; 
system ($cleaningCmd) and die ("ERROR: $0 : Cannot remove the previous test directory with the command $cleaningCmd \n$!\n");

my $expectedData="../../DATA/expectedData/";

########################################
#Creation of test directory
########################################
my $makeDirCmd = "mkdir $testingDir";
system ($makeDirCmd) and die ("ERROR: $0 : Cannot create the new directory with the command $makeDirCmd\n$!\n");
chdir $testingDir or die ("ERROR: $0 : Cannot go into the new directory with the command \"chdir $testingDir\"\n$!\n");

#######################################
#Creating the IndividuSoft.txt file
#######################################
my $creatingCmd="echo \"cutadapt\nTEST\" > individuSoft.txt";
system($creatingCmd) and die ("ERROR: $0 : Cannot create the individuSoft.txt file with the command $creatingCmd\n$!\n");

#######################################
#Cleaning the logs for the test
#######################################
$cleaningCmd="rm -Rf fastqc_TEST_log.*";
system($cleaningCmd) and die ("ERROR: $0 : Cannot remove the previous log files with the command $cleaningCmd \n$!\n");

########################################
#Input files
########################################
my $fastqFile = "RC3_2.fastq";                      # fastq file for test
my $originalFastqFile = $expectedData."RC3_2.fastq";     # fastq file
my $lnCmd = "ln -s $originalFastqFile .";          # command to copy the original fastq file into the test directory
system ($lnCmd) and die ("ERROR: $0 : Cannot copy the file $originalFastqFile in the test directory with the command $lnCmd\n$!\n");    # RUN the copy command

my $adaptatorFile = "adaptators.txt";                         # adaptator file for test
my $originalAdaptatorFile = $expectedData.$adaptatorFile;     # adaptator file
$lnCmd = "ln -s $originalAdaptatorFile .";             # command to copy the original adaptator file into the test directory
system ($lnCmd) and die ("ERROR: $0 : Cannot copy the file $originalAdaptatorFile in the test directory with the command $lnCmd\n$!\n");    # RUN the copy command

my $fileOut = "RC3_2.CUTADAPT.fastq";                                                  # Output file without adaptators sequences
######################




### Test of cutadapt::exec ###
my %optionsHachees = (
                    "-O" => 10,
                    "-m" => 35,
                    "-q" => 20,
                    "--overlap" => 7
);        # Hash containing informations
my $optionsHachees = \%optionsHachees;   
is ((cutadapt::execution($fastqFile,$fileOut,undef, undef, $optionsHachees)),1, 'cutadapt::execution');                      # TEST IF FONCTION WORKS
my $md5sumOfRefOut = "3275afc598641cd5ed97ab21d371194b";                                            # structure of the ref file for checking
my $md5sumOfFileOut = `md5sum $fileOut`;                                                            # structure of the test file for checking
my @nameless = split (" ", $md5sumOfFileOut);                                                       # to separate the structure and the name of file
$md5sumOfFileOut = $nameless[0];                                                                    # just to have the md5sum result
is_deeply ($md5sumOfFileOut, $md5sumOfRefOut, "cutadapt::execution... Cutadapt out file checkout");                             # TEST IF THE STRUCTURE OF THE FILE OUT IS GOOD
##############################

exit;