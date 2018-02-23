#!/usr/bin/env perl

###################################################################################################################################
#
# Copyright 2014-2018 IRD-CIRAD-INRA-ADNid
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

use strict;
use warnings;
use Test::More 'no_plan';
use Test::Deep;
use fileConfigurator;
use localConfig;

#####################
## PATH for datas test
#####################
my $vcfData="$toggle/data/testData/vcf/duplicVCF/";

#########################################
#Remove files and directory created by previous test
#########################################
my $testingDir="$toggle/dataTest/duplicationDetectorBlocks";
my $creatingDirCom="rm -Rf $testingDir ; mkdir -p $testingDir";                                    #Allows to have a working directory for the tests
system($creatingDirCom) and die ("ERROR: $0 : Cannot execute the command $creatingDirCom\n$!\n");


print "\n\n#################################################\n";
print "#### TEST DUPLICATIONDETECTOR \n";
print "#################################################\n";

#Creating config file for this test
my @listSoft = ("duplicationDetector");
fileConfigurator::createFileConf(\@listSoft,"blockTestConfig.txt");

my $runCmd = "toggleGenerator.pl -c blockTestConfig.txt -d ".$vcfData." -o ".$testingDir;
print "\n### Toggle running : $runCmd\n";
system("$runCmd") and die "#### ERROR : Can't run TOGGLE for duplicationDetector";

# check final results
my $bedFileOut = $testingDir."/finalResults/smallDuplic.DUPLICATIONDETECTOR.bed";

#Check if files created
my @expectedOutput = ("smallDuplic.DUPLICATIONDETECTOR.bed");
my $observedOutput = `ls $testingDir/finalResults/`;
my @observedOutput = split /\n/,$observedOutput;
is_deeply(\@observedOutput,\@expectedOutput,'duplicationDetector::execute - Files created');

# expected content test FileOut
my $expectedLineNumber = "3 $bedFileOut";                                            # structure of the ref file for checking
my $observedLineNumber = `wc -l $bedFileOut`;                                                        # structure of the test file for checking
chomp $observedLineNumber;                                                     # to separate the structure and the name of file
is($observedLineNumber, $expectedLineNumber, "duplicationDetector::execute - output content file Bed");               # TEST IF THE STRUCTURE OF THE FILE OUT IS GOOD

###Test for correct file value of duplicationDetector
#GREP command result
my $grepResult=`grep "359593" $bedFileOut`;
chomp $grepResult;
is($grepResult,"Reference	359593	359695	102	11	1.00",'duplicationDetector::execute - output grep in file Bed');
