
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


## PATH for datas test


# input file
my $bamSorted = "$toggle/data/testData/rnaseq/bamSorted/";
my $gtf="$toggle/data/testData/rnaseq/gtf/";

#####################
## TOGGLE stringtie
#####################

print "\n\n#################################################\n";
print "#### TEST stringtie \n";
print "#################################################\n";

# Remove files and directory created by previous test
my $testingDir="$toggle/dataTest/stringtie-noSGE-Blocks";
my $cleaningCmd="rm -Rf $testingDir";
system ($cleaningCmd) and die ("ERROR: $0 : Cannot remove the previous test directory with the command $cleaningCmd \n$!\n");

#Creating config file for this test
my @listSoft = ("stringtie");
fileConfigurator::createFileConf(\@listSoft,"blockTestConfig.txt");
my $runCmd = "toggleGenerator.pl -c blockTestConfig.txt -d ".$bamSorted." -o ".$testingDir;
print "\n### Toggle running : $runCmd\n";
system("$runCmd") and die "#### ERROR : Can't run TOGGLE for stringtie";

## check final results
print "\n### TEST Ouput list & content : $runCmd\n";
my $observedOutput = `ls $testingDir/finalResults/`;
my @observedOutput = split /\n/,$observedOutput;
my @expectedOutput = ('RNASeq2.STRINGTIE.gtf','RNASeq.STRINGTIE.gtf');

## expected output test
is_deeply(\@observedOutput,\@expectedOutput,'toggleGenerator - Two Sorted Bam (no SGE) stringtie file list');

# expected output structure
my $expectedMD5sum = "c3e6f50be438ed87fba88911aec66a3e";
my $observedMD5sum=`md5sum $testingDir/finalResults/RNASeq.STRINGTIE.gtf`;;# structure of the test file
my @withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'toggleGenerator - Two Sorted Bam (no SGE) stringtie file content');


print "\n\n#################################################\n";
print "#### TEST stringtie with --merge option  TODO \n";
print "#################################################\n";

# Remove files and directory created by previous test
$testingDir="$toggle/dataTest/stringtie-noSGE-Blocks";
$cleaningCmd="rm -Rf $testingDir";
system ($cleaningCmd) and die ("ERROR: $0 : Cannot remove the previous test directory with the command $cleaningCmd \n$!\n");

#Creating config file for this test
@listSoft = ("1000","stringtie");
fileConfigurator::createFileConf(\@listSoft,"blockTestConfig.txt");
my $sed="sed -i -e 's|\$stringtie|\$stringtie\\n--merge\\n|'  blockTestConfig.txt"; #adding --merge option in blocktestconfig
## DEBUG print $sed;
system($sed) and die ("#### ERROR  SED COMMAND: $sed\n");

$runCmd = "toggleGenerator.pl -c blockTestConfig.txt -d ".$gtf." -o ".$testingDir;
print "\n### Toggle running : $runCmd\n";
system("$runCmd") and die "#### ERROR : Can't run TOGGLE for stringtie";

## check final results
print "\n### TEST Ouput list & content : $runCmd\n";
$observedOutput = `ls $testingDir/finalResults/`;
@observedOutput = split /\n/,$observedOutput;
@expectedOutput = ('globalAnalysis.STRINGTIEMERGE.gtf');

## expected output test
is_deeply(\@observedOutput,\@expectedOutput,'toggleGenerator - Two GTF (no SGE) stringtie file list');

# expected output structure
$expectedMD5sum = "981d3db4929a0bfbe75ccf575d1c7169";
$observedMD5sum=`md5sum $testingDir/finalResults/globalAnalysis.STRINGTIEMERGE.gtf`;;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'toggleGenerator - Two GTF (no SGE) stringtie file content');

