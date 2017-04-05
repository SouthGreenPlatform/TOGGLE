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

#Will test if radseq works correctly
use strict;
use warnings;
use Test::More 'no_plan'; #Number of tests, to modify if new tests implemented. Can be changed as 'no_plan' instead of tests=>11 .
use Test::Deep;
use Data::Dumper;
use lib qw(../../modules/);


########################################
#use of radseq module ok
########################################
use_ok('localConfig') or exit;                                                                  # Check if toolbox is usable
use_ok('radseq') or exit;                                                                   # Check if radseq is usable
can_ok('radseq','processRadtags');                                                          # Check if radseq::processRadtags is find
can_ok('radseq','parseKeyFile');                                                            # Check if radseq::parseKeyFile is find
can_ok('radseq','executeDemultiplexing');                                                   # Check if radseq::executeDemultiplexing is find
can_ok('radseq','checkOrder');                                                              # Check if radseq::checkOrder is find
can_ok('toolbox','rmHashOrder');                                                             # Check if radseq::rmHashOrder is find

use localConfig;
use radseq;

my $expectedData="$toggle/data/expectedData/";

#########################################
#Remove files and directory created by previous test
#########################################
my $testingDir="$toggle/dataTest/radseqTestDir";
my $creatingDirCom="rm -Rf $testingDir ; mkdir -p $testingDir";                                    #Allows to have a working directory for the tests
system($creatingDirCom) and die ("ERROR: $0 : Cannot execute the command $creatingDirCom\n$!\n");

chdir $testingDir or die ("ERROR: $0 : Cannot go into the new directory with the command \"chdir $testingDir\"\n$!\n");


#######################################
#Creating the IndividuSoft.txt file
#######################################
my $creatingCommand="echo \"radseq\nTEST\" > individuSoft.txt";
system($creatingCommand) and die ("ERROR: $0 : Cannot create the individuSoft.txt file with the command $creatingCommand\n$!\n");


#######################################
#Cleaning the logs for the test
#######################################
my $cleaningCommand="rm -f radseq_TEST_log.*";
system($cleaningCommand) and die ("ERROR: $0 : Cannot remove the previous log files with the command $cleaningCommand \n$!\n");



##########################################
##### radseq::rmHashOrder
##########################################

my %hashOrder = (
				'1' => 'processRadtags',
				'2' => 'fastqc'
				);
my $hashOrder = \%hashOrder;

# execution test
my $hashOrderRM;
$hashOrderRM = toolbox::rmHashOrder($hashOrder, "processRadtags");

# expected content test
is_deeply($hashOrder, $hashOrderRM, "toolbox::rmHashOrder - cleanning step one OK");               # TEST IF THE STRUCTURE OF THE FILE OUT IS GOOD


##########################################
##### radseq::checkOrder
##########################################


my $outputDir = "./";
my $fileConf = "$toggle/exampleConfigs/radseqSingle.config.txt";
my $initialDir = "$toggle/data/testData/radseq/single/";
my $keyfile = "$toggle/data/testData/radseq/keyfileTestSingle";
my $checkFastq = 1;


# execution test
my $initialDirContent;
#$initialDirContent = radseq::checkOrder($outputDir,%param);
#print Dumper($initialDirContent);

my @expectedOutput = ('.//demultiplexedFiles/33-16.fq', './/demultiplexedFiles/A239.fq', './/demultiplexedFiles/A272.fq', './/demultiplexedFiles/A554.fq', './/demultiplexedFiles/A619.fq', './/demultiplexedFiles/A632.fq', './/demultiplexedFiles/A654.fq', './/demultiplexedFiles/A659.fq', './/demultiplexedFiles/A680.fq', './/demultiplexedFiles/A682.fq', './/demultiplexedFiles/B103.fq', './/demultiplexedFiles/B104.fq', './/demultiplexedFiles/B109.fq', './/demultiplexedFiles/B10.fq', './/demultiplexedFiles/B73Htrhm.fq', './/demultiplexedFiles/B76.fq', './/demultiplexedFiles/B77.fq', './/demultiplexedFiles/B84.fq', './/demultiplexedFiles/B97.fq', './/demultiplexedFiles/C103.fq', './/demultiplexedFiles/CH701-30.fq', './/demultiplexedFiles/CI28A.fq', './/demultiplexedFiles/CI3A.fq', './/demultiplexedFiles/CI66.fq', './/demultiplexedFiles/CI-7.fq', './/demultiplexedFiles/CM105.fq', './/demultiplexedFiles/CM174.fq', './/demultiplexedFiles/CML14.fq', './/demultiplexedFiles/CML157Q.fq', './/demultiplexedFiles/CML220.fq', './/demultiplexedFiles/CML228.fq', './/demultiplexedFiles/CML238.fq', './/demultiplexedFiles/CML258.fq', './/demultiplexedFiles/CML277.fq', './/demultiplexedFiles/CML281.fq', './/demultiplexedFiles/CML311.fq', './/demultiplexedFiles/CML323.fq', './/demultiplexedFiles/CML332.fq', './/demultiplexedFiles/CML333.fq', './/demultiplexedFiles/CML45.fq', './/demultiplexedFiles/CML52.fq', './/demultiplexedFiles/CML91.fq', './/demultiplexedFiles/CML92.fq', './/demultiplexedFiles/CO106.fq', './/demultiplexedFiles/CO125.fq', './/demultiplexedFiles/DE-2.fq', './/demultiplexedFiles/DE-3.fq', './/demultiplexedFiles/EMPTY.fq', './/demultiplexedFiles/EP1.fq', './/demultiplexedFiles/F6.fq', './/demultiplexedFiles/F7.fq', './/demultiplexedFiles/H91.fq', './/demultiplexedFiles/H95.fq', './/demultiplexedFiles/I137TN.fq', './/demultiplexedFiles/I29.fq', './/demultiplexedFiles/IDS28.fq', './/demultiplexedFiles/K55.fq', './/demultiplexedFiles/Ki11.fq', './/demultiplexedFiles/Ki14.fq', './/demultiplexedFiles/L317.fq', './/demultiplexedFiles/M14.fq', './/demultiplexedFiles/NC230.fq', './/demultiplexedFiles/NC250.fq', './/demultiplexedFiles/NC262.fq', './/demultiplexedFiles/NC290A.fq', './/demultiplexedFiles/NC298.fq', './/demultiplexedFiles/NC300.fq', './/demultiplexedFiles/NC318.fq', './/demultiplexedFiles/NC320.fq', './/demultiplexedFiles/NC328.fq', './/demultiplexedFiles/NC338.fq', './/demultiplexedFiles/NC356.fq', './/demultiplexedFiles/NC364.fq', './/demultiplexedFiles/NC368.fq', './/demultiplexedFiles/Oh40B.fq', './/demultiplexedFiles/Oh43E.fq', './/demultiplexedFiles/Oh603.fq', './/demultiplexedFiles/Os420.fq', './/demultiplexedFiles/Pa875.fq', './/demultiplexedFiles/R109B.fq', './/demultiplexedFiles/R168.fq', './/demultiplexedFiles/SC213R.fq', './/demultiplexedFiles/Sg18.fq', './/demultiplexedFiles/tripsacum.fq', './/demultiplexedFiles/Tzi25.fq', './/demultiplexedFiles/Va22.fq', './/demultiplexedFiles/Va35.fq', './/demultiplexedFiles/Va59.fq', './/demultiplexedFiles/Va99.fq', './/demultiplexedFiles/W117Ht.fq', './/demultiplexedFiles/W153R.fq', './/demultiplexedFiles/W182B.fq', './/demultiplexedFiles/W22.fq', './/demultiplexedFiles/WD.fq' );
is_deeply(radseq::checkOrder($outputDir,$fileConf,$initialDir,$checkFastq,$keyfile),\@expectedOutput,'radseq::checkOrder - running');

# expected output test
#Check if files created
@expectedOutput = ('33-16.fq', 'A239.fq', 'A272.fq', 'A554.fq', 'A619.fq', 'A632.fq', 'A654.fq', 'A659.fq', 'A680.fq', 'A682.fq', 'B103.fq', 'B104.fq', 'B109.fq', 'B10.fq', 'B73Htrhm.fq', 'B76.fq', 'B77.fq', 'B84.fq', 'B97.fq', 'C103.fq', 'CH701-30.fq', 'CI28A.fq', 'CI3A.fq', 'CI66.fq', 'CI-7.fq', 'CM105.fq', 'CM174.fq', 'CML14.fq', 'CML157Q.fq', 'CML220.fq', 'CML228.fq', 'CML238.fq', 'CML258.fq', 'CML277.fq', 'CML281.fq', 'CML311.fq', 'CML323.fq', 'CML332.fq', 'CML333.fq', 'CML45.fq', 'CML52.fq', 'CML91.fq', 'CML92.fq', 'CO106.fq', 'CO125.fq', 'DE-2.fq', 'DE-3.fq', 'EMPTY.fq', 'EP1.fq', 'F6.fq', 'F7.fq', 'H91.fq', 'H95.fq', 'I137TN.fq', 'I29.fq', 'IDS28.fq', 'K55.fq', 'Ki11.fq', 'Ki14.fq', 'L317.fq', 'M14.fq', 'NC230.fq', 'NC250.fq', 'NC262.fq', 'NC290A.fq', 'NC298.fq', 'NC300.fq', 'NC318.fq', 'NC320.fq', 'NC328.fq', 'NC338.fq', 'NC356.fq', 'NC364.fq', 'NC368.fq', 'Oh40B.fq', 'Oh43E.fq', 'Oh603.fq', 'Os420.fq', 'Pa875.fq', 'R109B.fq', 'R168.fq', 'SC213R.fq', 'Sg18.fq', 'tripsacum.fq', 'Tzi25.fq', 'Va22.fq', 'Va35.fq', 'Va59.fq', 'Va99.fq', 'W117Ht.fq', 'W153R.fq', 'W182B.fq', 'W22.fq', 'WD.fq' );
my $observedOutput = `ls demultiplexedFiles`;
my @observedOutput = split /\n/,$observedOutput;

is_deeply(\@observedOutput,\@expectedOutput,'radseq::checkOrder - Demultiplexed files');



