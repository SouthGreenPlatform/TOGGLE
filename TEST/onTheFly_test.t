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
can_ok('onTheFly','generateGraphviz');

use onTheFly;

#########################################
#Remove the files and directory created by the previous test
#########################################
my $testingDir="../DATA-TEST/onTheFlyTestDir";
my $rmDirCommand="rm -Rf $testingDir; mkdir $testingDir"; 
system ("$rmDirCommand") and die ("ERROR: $0 : Cannot remove the previous test directory with the command $rmDirCommand \n$!\n");
chdir $testingDir or die ("ERROR: $0 : Cannot create $testingDir\n$!\n");

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


########################################
#Creation of test files
########################################
my $originalFastaRef="../../DATA/expectedData/Reference.fasta";
my $fastaRef="Reference.fasta";
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
                                        "2" => "samToolsView",
                                        "3" => "samToolsSort"
                                    }
                    };

my $refFastaFile="../../DATA/Bank/referenceIrigin.fasta";
my $gffFile="None";
my $keyfile="None";

my @output = onTheFly::checkOrder($hashOrderOk,$refFastaFile,$gffFile,$keyfile);
my @expected=('1','3');
is_deeply(\@output,\@expected,'onTheFly::checkOrder - normal order');

#testing the uncorrect rendering TEST Ok in dev, but die so cannot be tested...
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
                                        "3" => "samToolsView",
                                        "4" => "samToolsSort"
                                    }
                    };

@output = onTheFly::checkOrder($hashOrderNAOk,$refFastaFile,$gffFile,$keyfile);
@expected=('1','4');
is_deeply(\@output,\@expected,'onTheFly::checkOrder - dead-end beginning');

#Testing for dead-end program in the middle
my $hashOrderNAOkBis =   {
                        "order"=>   {
                                        "3" => "fastqc",
                                        "1" => "bwaSampe",
                                        "2" => "samToolsView",
                                        "4" => "samToolsSort"
                                    }
                    };
@output = onTheFly::checkOrder($hashOrderNAOkBis,$refFastaFile,$gffFile,$keyfile);
@expected=('1','4');
is_deeply(\@output,\@expected,'onTheFly::checkOrder - dead-end intermediate');

#testing the single program
my $hashOrderSingle =   {
                        "order"=>   {
                                        "3" => "bwaSampe"
                                    }
                    };
@output = onTheFly::checkOrder($hashOrderSingle,$refFastaFile,$gffFile,$keyfile);
@expected=('3','3');
is_deeply(\@output,\@expected,'onTheFly::checkOrder - Single software');


#testing multiple call of the same program
my $hashOrderMultiple =   {
                        "order"=>   {
                                        "1" => "bwaSampe",
                                        "2" => "samToolsView",
                                        "3" => "samToolsSort",
                                        "4" => "samToolsView"
                                    }
                    };
@output = onTheFly::checkOrder($hashOrderMultiple,$refFastaFile,$gffFile,$keyfile);
@expected=('1','4');
is_deeply(\@output,\@expected,'onTheFly::checkOrder - multiple calls of the same software');

########################################
#generateScript test
#######################################

#testing the correct rendering
#Adding a configHash
my $hashConf =  {
                "1" => "bwaAln",
                "2" => "bwaSampe",
                "3" => "samToolsView"
                };
#output file
my $outputScript = "toggleBzzz.pl";

#execution test
is (onTheFly::generateScript($hashConf,$outputScript),'1','onTheFly::generateScript');

# expected output test
my $observedOutput = `ls`;
my @observedOutput = split /\n/,$observedOutput;
my @expectedOutput = ('individuSoft.txt','onTheFly_TEST_log.e','onTheFly_TEST_log.o','Reference.fasta','toggleBzzz.pl');
#
is_deeply(\@observedOutput,\@expectedOutput,'onTheFly::generateScript - output list');

# expected content test

my $expectedMD5sum="ef3cc44ab074e4a66b23558de171f0c9";
my $observedMD5sum=`md5sum $outputScript`;# structure of the test file
my @withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'onTheFly::generateScript - output structure');

########################################
#indexCreator test
#######################################

#testing the correct rendering
#Adding a configHash
$hashConf =   {
                        "order"=>   {
                                        "1" => "bwaSampe",
                                        "2" => "gatkHaplotypeCaller",
                                        "3" => "samToolsSort"
                                    }
                    };
is (onTheFly::indexCreator($hashConf,$fastaRef),'1','onTheFly::indexCreator - creation');

# expected output test
$observedOutput = `ls`;
@observedOutput = split /\n/,$observedOutput;
@expectedOutput = ('individuSoft.txt','onTheFly_TEST_log.e','onTheFly_TEST_log.o','Reference.dict','Reference.fasta','Reference.fasta.amb','Reference.fasta.ann','Reference.fasta.bwt','Reference.fasta.fai','Reference.fasta.pac','Reference.fasta.sa','toggleBzzz.pl');
#
is_deeply(\@observedOutput,\@expectedOutput,'onTheFly::indexCreator - output list');

# expected content test

$expectedMD5sum="4b9a4431e72c9db7e5c1f2153eba9fe7";
$observedMD5sum=`md5sum $fastaRef.fai`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
#is($observedMD5sum,$expectedMD5sum,'onTheFly::generateScript - output structure');

#Testing if creating in case of existing refs.
#Input file
my $initialModifDate = `stat -c %y Reference.fasta.pac`; #type such as 2016-05-11 15:49:19.696893241 +0200
chomp $initialModifDate;

#execution test
is (onTheFly::indexCreator($hashConf,$fastaRef),'1','onTheFly::indexCreator no creation');

# expected output test
$observedOutput = `ls`;
@observedOutput = split /\n/,$observedOutput;
@expectedOutput = ('individuSoft.txt','onTheFly_TEST_log.e','onTheFly_TEST_log.o','Reference.dict','Reference.fasta','Reference.fasta.amb','Reference.fasta.ann','Reference.fasta.bwt','Reference.fasta.fai','Reference.fasta.pac','Reference.fasta.sa','toggleBzzz.pl');
#
is_deeply(\@observedOutput,\@expectedOutput,'onTheFly::indexCreator no creation - output list');

# expected content test
$observedOutput = `stat -c %y Reference.fasta.pac`;
chomp $observedOutput;
is($observedOutput,$initialModifDate,'onTheFly::indexCreator no creation - output structure');


#Testing fro re-creating forced of index
$hashConf =   {
                        "order"=>   {
                                        "1" => "bwaSampe",
                                        "2" => "gatkHaplotypeCaller",
                                        "3" => "bwaIndex"
                                    }
                      };
#execution test
is (onTheFly::indexCreator($hashConf,$fastaRef),'1','onTheFly::indexCreator forcing creation');

# expected output test
$observedOutput = `ls`;
@observedOutput = split /\n/,$observedOutput;
@expectedOutput = ('individuSoft.txt','onTheFly_TEST_log.e','onTheFly_TEST_log.o','Reference.dict','Reference.fasta','Reference.fasta.amb','Reference.fasta.ann','Reference.fasta.bwt','Reference.fasta.fai','Reference.fasta.pac','Reference.fasta.sa','toggleBzzz.pl');
#
is_deeply(\@observedOutput,\@expectedOutput,'onTheFly::indexCreator forcing creation - output list');

# expected content test
$observedOutput = `stat -c %y Reference.fasta.pac`;
chomp $observedOutput;
isnt($observedOutput,$initialModifDate,'onTheFly::indexCreator forcing creation - output structure');


########################################
#generateGraphviz test
#######################################

#Generating normal hash for order
$hashConf = {
            "1" => "bwaSampe",
            "2" => "gatkHaplotypeCaller",
            "3" => "bwaIndex"
            };

#execution test
is (onTheFly::generateGraphviz($hashConf,$toggle."/DATA-TEST/onTheFlyTestDir"),'1','onTheFly::generateGraphviz');

# expected output test
$observedOutput = `ls`;
@observedOutput = split /\n/,$observedOutput;
@expectedOutput = ('individuSoft.txt','onTheFly_TEST_log.e','onTheFly_TEST_log.o','Reference.dict','Reference.fasta','Reference.fasta.amb','Reference.fasta.ann','Reference.fasta.bwt','Reference.fasta.fai','Reference.fasta.pac','Reference.fasta.sa','toggleBzzz.pl','togglePipeline.dot','togglePipeline.png');
#
is_deeply(\@observedOutput,\@expectedOutput,'onTheFly::generateGraphviz - output list');

# expected content test

$expectedMD5sum="14";
$observedMD5sum=`wc -l togglePipeline.dot`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'onTheFly::generateGraphviz - output structure');