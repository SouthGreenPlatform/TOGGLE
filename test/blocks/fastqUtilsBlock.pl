#!/usr/bin/env perl

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

use strict;
use warnings;
use Test::More 'no_plan';
use Test::Deep;
use lib qw(../../modules/);
use fileConfigurator;
use localConfig;

#####################
## PATH for datas test
#####################

# references files
my $dataFastq33 = "$toggle/data/testData/fastq/pairedTwoIndividusIrigin/";
my $dataFastq64 = "$toggle/data/testData/fastq/singleIndividualPHRED64/";

print "\n\n#################################################\n";
print "#### TEST checkEncodeByASCIIcontrol - PHRED+33\n";
print "#################################################\n";

# Remove files and directory created by previous test
my $testingDir="$toggle/dataTest/checkEncodeByASCIIcontrol33";
my $cleaningCmd="rm -Rf $testingDir";
system ($cleaningCmd) and die ("ERROR: $0 : Cannot remove the previous test directory with the command $cleaningCmd \n$!\n");

#Creating config file for this test
my @listSoft = ("checkEncodeByASCIIcontrol");
fileConfigurator::createFileConf(\@listSoft,"blockTestConfig.txt");

my $runCmd = "toggleGenerator.pl -c blockTestConfig.txt -d ".$dataFastq33." -o ".$testingDir;
print "\n### Toggle running : $runCmd\n";
system("$runCmd") and warn "#### ERROR : Can't run TOGGLE for checkEncodeByASCIIcontrol33";

# check final results
print "\n### TEST Ouput list & content : $runCmd\n";
my $observedOutput = `ls $testingDir/finalResults`;
my @observedOutput = split /\n/,$observedOutput;
my @expectedOutput = ('irigin1.checkEncodeByASCIIcontrol.log','irigin3.checkEncodeByASCIIcontrol.log');

# expected output test
is_deeply(\@observedOutput,\@expectedOutput,'toggleGenerator - fastqUtils checkEncodeByASCIIcontrol33 list ');

# expected output content
$observedOutput=`grep -c \"is encoded as a PHRED+33 file\" $testingDir/finalResults/irigin1.checkEncodeByASCIIcontrol.log`; 
chomp $observedOutput;
my $expectedOutput = "2";
is($observedOutput,$expectedOutput, 'toggleGenerator - fastqUtils checkEncodeByASCIIcontrol33 content ');

print "\n\n#################################################\n";
print "#### TEST checkEncodeByASCIIcontrol - PHRED+64+\n";
print "#################################################\n";

# Remove files and directory created by previous test
$testingDir="$toggle/dataTest/checkEncodeByASCIIcontrol64";
$cleaningCmd="rm -Rf $testingDir";
system ($cleaningCmd) and die ("ERROR: $0 : Cannot remove the previous test directory with the command $cleaningCmd \n$!\n");

#Creating config file for this test
@listSoft = ("checkEncodeByASCIIcontrol");
fileConfigurator::createFileConf(\@listSoft,"blockTestConfig.txt");

$runCmd = "toggleGenerator.pl -c blockTestConfig.txt -d ".$dataFastq64." -o ".$testingDir;
print "\n### Toggle running : $runCmd\n";
system("$runCmd") and warn "#### ERROR : Can't run TOGGLE for checkEncodeByASCIIcontrol64";

# check final results
print "\n### TEST Ouput list & content : $runCmd\n";
$observedOutput = `ls $testingDir/finalResults`;
@observedOutput = split /\n/,$observedOutput;
@expectedOutput = ('RC1.checkEncodeByASCIIcontrol.log');

# expected output test
is_deeply(\@observedOutput,\@expectedOutput,'toggleGenerator - fastqUtils checkEncodeByASCIIcontrol64 list ');

# expected output content
$observedOutput=`grep -c \"is not encoded as a PHRED+33 file\" $testingDir/finalResults/RC1.checkEncodeByASCIIcontrol.log`; 
chomp $observedOutput;
$expectedOutput = "1";
is($observedOutput,$expectedOutput, 'toggleGenerator - fastqUtils checkEncodeByASCIIcontrol64 content ');
