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

# references files
my $dataRefIrigin = "$toggle/data/Bank/referenceIrigin.fasta";
# input file
my $dataFastq="$toggle/data/testData/fastq/pairedTwoIndividusIrigin";


print "\n\n#################################################\n";
print "#### TEST BOWTIE bowtie\n";
print "#################################################\n";

# Remove files and directory created by previous test
my $testingDir="$toggle/dataTest/bowtie-noSGE-Blocks";
my $cleaningCmd="rm -Rf $testingDir";
system ($cleaningCmd) and die ("ERROR: $0 : Cannot remove the previous test directory with the command $cleaningCmd \n$!\n");

#Creating config file for this test
my @listSoft = ("bowtie");
fileConfigurator::createFileConf(\@listSoft,"blockTestConfig.txt");

my $runCmd = "toggleGenerator.pl -c blockTestConfig.txt -d ".$dataFastq." -r ".$dataRefIrigin." -o ".$testingDir;
print "\n### Toggle running : $runCmd\n";
system("$runCmd") and die "#### ERROR : Can't run TOGGLE for bowtie";

# check final results

# expected output content
my $observedOutput = `ls $testingDir/finalResults`;
my @observedOutput = split /\n/,$observedOutput;
my @expectedOutput = ('irigin1.BOWTIE.sam','irigin3.BOWTIE.sam');

is_deeply(\@observedOutput,\@expectedOutput,'toggleGenerator - Two fastq (no SGE) BOWTIE file list ');

# expected output value
my $grepResult= `grep -c 'H3:C39R6ACXX' $testingDir/finalResults/irigin1.BOWTIE.sam`;
chomp $grepResult;
is($grepResult,8,'toggleGenerator - Two fastq (no SGE) BOWTIE result of bowtie irigin1.BOWTIE.sam');
# expected output value
$grepResult= `grep -c 'H2:C381HACXX' $testingDir/finalResults/irigin3.BOWTIE.sam`;
chomp $grepResult;
is($grepResult,456,'toggleGenerator - Two fastq (no SGE) BOWTIE result of bowtie irigin3.BOWTIE.sam');


print "\n\n#################################################\n";
print "#### TEST BOWTIE2 bowtie2\n";
print "#################################################\n";

# Remove files and directory created by previous test
$testingDir="$toggle/dataTest/bowtie2-noSGE-Blocks";
$cleaningCmd="rm -Rf $testingDir";
system ($cleaningCmd) and die ("ERROR: $0 : Cannot remove the previous test directory with the command $cleaningCmd \n$!\n");

#Creating config file for this test
@listSoft = ("bowtie2");
fileConfigurator::createFileConf(\@listSoft,"blockTestConfig.txt");

$runCmd = "toggleGenerator.pl -c blockTestConfig.txt -d ".$dataFastq." -r ".$dataRefIrigin." -o ".$testingDir;
print "\n### Toggle running : $runCmd\n";
system("$runCmd") and die "#### ERROR : Can't run TOGGLE for bowtie2";

# check final results

# expected output content
$observedOutput = `ls $testingDir/finalResults`;
@observedOutput = split /\n/,$observedOutput;
@expectedOutput = ('irigin1.BOWTIE2.sam','irigin3.BOWTIE2.sam');

is_deeply(\@observedOutput,\@expectedOutput,'toggleGenerator - Two fastq (no SGE) BOWTIE2 file list ');

# expected output value
$grepResult= `grep -c 'H3:C39R6ACXX' $testingDir/finalResults/irigin1.BOWTIE2.sam`;
chomp $grepResult;
is($grepResult,2000,'toggleGenerator - Two fastq (no SGE) BOWTIE result of bowtie2 irigin1.BOWTIE2.sam');
# expected output value
$grepResult= `grep -c 'H2:C381HACXX' $testingDir/finalResults/irigin3.BOWTIE2.sam`;
chomp $grepResult;
is($grepResult,2000,'toggleGenerator - Two fastq (no SGE) BOWTIE result of bowtie2 irigin3.BOWTIE2.sam');