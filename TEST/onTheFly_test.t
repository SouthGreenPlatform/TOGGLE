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

#Will test if onTheFly works correctly

use strict;
use warnings;
use Test::More 'no_plan'; #tests => 19; #Number of tests, to modify if new tests implemented. Can be changed as 'no_plan' instead of tests=>11 .
use Data::Dumper;
use Test::Deep;
use lib qw(../Modules/);

use localConfig;
my $configFile='software.config.txt';

use_ok('onTheFly');
can_ok('onTheFly','checkOrder');
can_ok('onTheFly','generateScript');
can_ok('onTheFly','indexCreator');

use onTheFly;


#######################################
#Creating the IndividuSoft.txt file
#######################################
#system("mkdir LOGS");
my $creatingCommand="echo \"onTheFly\nTEST\" > individuSoft.txt";
system($creatingCommand) and die ("ERROR: $0: Cannot create the individuSoft.txt file with the command $creatingCommand \n$!\n");

#######################################
#Cleaning the logs for the test
#######################################
my $cleaningCommand="rm -Rf onTheFly_TEST_log.*";
system($cleaningCommand) and die ("ERROR: $0: Cannot clean the previous log files for this test with the command $cleaningCommand \n$!\n");



#########################################
#Remove the files and directory created by the previous test
#########################################
my $testingDir="../DATA-TEST/onTheFlyTestDir";
$cleaningCommand="rm -Rf ../DATA-TEST/$testingDir"; 
system ("rm -Rf $testingDir") and die ("ERROR: $0 : Cannot remove the previous test directory with the command $cleaningCommand \n$!\n");


########################################
#Creation of test directory
########################################
my $makeDirCom = "mkdir $testingDir";
system ($makeDirCom) and die ("ERROR: $0 : Cannot create the new directory with the command $makeDirCom\n$!\n");

########################################
#Creation of test files
########################################
my $originalFastaRef="../DATA/expectedData/Reference.fasta";
my $fastaRef="$testingDir/Reference.fasta";
my $refCopyCom="cp $originalFastaRef $fastaRef";
system($refCopyCom) and die ("ERROR: $0 : Cannot copy the Reference for test with the command $refCopyCom \n$!\n");
  #Now we have a ref to be tested

########################################
#checkOrder test
#######################################

#testing the correct rendering
#Adding a configHash
my $hashOrderOk =   {
                        "order"=>   {
                                        "1" => "bwaSampe",
                                        "2" => "samtoolsView",
                                        "3" => "samtoolsSort"
                                    }
                    };

my @output = onTheFly::checkOrder($hashOrderOk);
my @expected=('1','3');
is_deeply(\@output,\@expected,'Test for correct pipeline onTheFly::checkOrder');

#testi;ng the uncorrect rendering TEST Ok in dev, but die so cannot be tested...
#Adding a configHash
#my $hashOrderNotOk =   {
#                        "order"=>   {
#                                        "2" => "bwaSampe",
#                                        "1" => "samtools view",
#                                        "3" => "samtools sort"
#                                    }
#                        };
#is (onTheFly::checkOrder($hashOrderNotOk),'0','Test for uncorrect pipeline onTheFly::checkOrder');

#testing for dead-end program beginning
my $hashOrderNAOk =   {
                        "order"=>   {
                                        "1" => "fastqc",
                                        "2" => "bwaSampe",
                                        "3" => "samtoolsView",
                                        "4" => "samtoolsSort"
                                    }
                    };

@output = onTheFly::checkOrder($hashOrderNAOk);
@expected=('1','4');
is_deeply(\@output,\@expected,'Test for correct pipeline with \'dead-end\' software beginning onTheFly::checkOrder');

#Testing for dead-end program in the middle
my $hashOrderNAOkBis =   {
                        "order"=>   {
                                        "3" => "fastqc",
                                        "1" => "bwaSampe",
                                        "2" => "samtoolsView",
                                        "4" => "samtoolsSort"
                                    }
                    };
@output = onTheFly::checkOrder($hashOrderNAOkBis);
@expected=('1','4');
is_deeply(\@output,\@expected,'Test for correct pipeline with \'dead-end\' software in the middle onTheFly::checkOrder');

#testing the single program
my $hashOrderSingle =   {
                        "order"=>   {
                                        "3" => "bwaSampe"
                                    }
                    };
@output = onTheFly::checkOrder($hashOrderSingle);
@expected=('3','3');
is_deeply(\@output,\@expected,'Test for correct pipeline with a single software onTheFly::checkOrder');


#testing multiple call of the same program
my $hashOrderMultiple =   {
                        "order"=>   {
                                        "1" => "bwaSampe",
                                        "2" => "samtoolsView",
                                        "3" => "samtoolsSort",
                                        "4" => "samtoolsView"
                                    }
                    };
@output = onTheFly::checkOrder($hashOrderMultiple);
@expected=('1','4');
is_deeply(\@output,\@expected,'Test for correct pipeline with multiple call of the same software onTheFly::checkOrder');

########################################
#generateScript test
#######################################

#testing the correct rendering
#Adding a configHash
my $hashConf =   {
                        "order"=>   {
                                        "1" => "bwaAln",
                                        "2" => "bwaSampe",
                                        "3" => "samtoolsView"
                                    }
                    };
is (onTheFly::generateScript($hashConf,"$testingDir/ToggleBzzz.pl"),'1','Test for correct pipeline onTheFly::generateScript');

########################################
#indexCreator test
#######################################

#testing the correct rendering
#Adding a configHash
$hashConf =   {
                        "order"=>   {
                                        "1" => "bwaSampe",
                                        "2" => "gatkHaplotypeCaller",
                                        "3" => "samtoolsSort"
                                    }
                    };
is (onTheFly::indexCreator($hashConf,$fastaRef),'1','Test for correct onTheFly::indexCreator running');

#Testing if creating in case of existing refs.
is (onTheFly::indexCreator($hashConf,$fastaRef),'1','Test for not recreating the ref index for onTheFly::indexCreator running');

#Testing fro re-creating forced of index
$hashConf =   {
                        "order"=>   {
                                        "1" => "bwaSampe",
                                        "2" => "gatkHaplotypeCaller",
                                        "3" => "bwaIndex"
                                    }
                    };
is (onTheFly::indexCreator($hashConf,$fastaRef),'1','Test for forced recreating the ref index for onTheFly::indexCreator running');