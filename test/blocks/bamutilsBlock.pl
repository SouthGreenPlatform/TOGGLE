
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
use fileConfigurator;
use localConfig;


## PATH for datas test


# input file
my $dataOneBam = "$toggle/data/testData/samBam/oneBam/";


#####################
## TOGGLE bamutilsFilter
#####################

print "\n\n#################################################\n";
print "#### TEST bamutilsFilter \n";
print "#################################################\n";

# Remove files and directory created by previous test
my $testingDir="$toggle/dataTest/bamutilsFilter-noSGE-Blocks";
my $cleaningCmd="rm -Rf $testingDir";
system ($cleaningCmd) and die ("ERROR: $0 : Cannot remove the previous test directory with the command $cleaningCmd \n$!\n");

#Creating config file for this test
my @listSoft = ("bamutilsFilter", "bamutilsRemoveclipping");
fileConfigurator::createFileConf(\@listSoft,"blockTestConfig.txt");

my $runCmd = "toggleGenerator.pl -c blockTestConfig.txt -d ".$dataOneBam." -o ".$testingDir;
print "\n### Toggle running : $runCmd\n";
system("$runCmd") and die "#### ERROR : Can't run TOGGLE for bamutilsFilter and bamutilsRemoveclipping";


## check final results
print "\n### TEST Ouput list & content : $runCmd\n";
my $observedOutput = `ls $testingDir/finalResults/`;
my @observedOutput = split /\n/,$observedOutput;
my @expectedOutput = ('RC3-SAMTOOLSVIEW.bamutilsremoveclipping.bam');

## expected output test
is_deeply(\@observedOutput,\@expectedOutput,'toggleGenerator - One Bam (no SGE) bamutilsFilter and bamutilsRemoveclipping file list');

# expected output structure
my $expectedMD5sum = "40d83df65cee1917dfb847c3c4ab0bd3";
my $observedMD5sum=`md5sum $testingDir/finalResults/RC3-SAMTOOLSVIEW.bamutilsremoveclipping.bam`;;# structure of the test file
my @withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'toggleGenerator - One Bam (no SGE) bamutilsFilter and bamutilsRemoveclipping content');

