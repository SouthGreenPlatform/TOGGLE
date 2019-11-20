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
#use Test::More 'no_plan';
#use Test::Deep;
use fileConfigurator;
use localConfig;

#####################
## PATH for datas test
#####################

# references files
my $dataRefIrigin = "$toggle/data/Bank/referenceIrigin.fasta";
my $dataFasta = "$toggle/data/testData/fasta/TGICL/";
my $dataOneBam = "$toggle/data/testData/samBam/oneBam/";
my $dataFastq = "$toggle/data/testData/fastq/pairedOneIndividuArcad/";
my $dataVcf = "$toggle/data/testData/vcf/singleVCF";

print "\n\n#################################################\n";
print "#### TEST checkFormatFasta\n";
print "#################################################\n";

# Remove files and directory created by previous test
my $testingDir="$toggle/dataTest/checkFormatFasta";
my $cleaningCmd="rm -Rf $testingDir";
system ($cleaningCmd) and die ("ERROR: $0 : Cannot remove the previous test directory with the command $cleaningCmd \n$!\n");

#Creating config file for this test
my @listSoft = ("checkFormatFasta");
fileConfigurator::createFileConf(\@listSoft,"blockTestConfig.txt");

my $runCmd = "toggleGenerator.pl -c blockTestConfig.txt -d ".$dataFasta." -r ".$dataRefIrigin." -o ".$testingDir;
print "\n### Toggle running : $runCmd\n";
system("$runCmd") and warn "#### ERROR : Can't run TOGGLE for checkFormatFasta";

print "\n\n#################################################\n";
print "#### TEST checkFormatFastq\n";
print "#################################################\n";

# Remove files and directory created by previous test
$testingDir="$toggle/dataTest/checkFormatFastq";
$cleaningCmd="rm -Rf $testingDir";
system ($cleaningCmd) and die ("ERROR: $0 : Cannot remove the previous test directory with the command $cleaningCmd \n$!\n");

#Creating config file for this test
@listSoft = ("checkFormatFastq");
fileConfigurator::createFileConf(\@listSoft,"blockTestConfig.txt");

$runCmd = "toggleGenerator.pl -c blockTestConfig.txt -d ".$dataFastq." -r ".$dataRefIrigin." -o ".$testingDir;
print "\n### Toggle running : $runCmd\n";
system("$runCmd") and warn "#### ERROR : Can't run TOGGLE for checkFormatFastq";

print "\n\n#################################################\n";
print "#### TEST checkFormatVcf\n";
print "#################################################\n";

# Remove files and directory created by previous test
$testingDir="$toggle/dataTest/checkFormatVcf";
$cleaningCmd="rm -Rf $testingDir";
system ($cleaningCmd) and die ("ERROR: $0 : Cannot remove the previous test directory with the command $cleaningCmd \n$!\n");

#Creating config file for this test
@listSoft = ("checkFormatVcf");
fileConfigurator::createFileConf(\@listSoft,"blockTestConfig.txt");

$runCmd = "toggleGenerator.pl -c blockTestConfig.txt -d ".$dataVcf." -r ".$dataRefIrigin." -o ".$testingDir;
print "\n### Toggle running : $runCmd\n";
system("$runCmd") and warn "#### ERROR : Can't run TOGGLE for checkFormatVcf";

print "\n\n#################################################\n";
print "#### TEST checkFormatSamOrBam\n";
print "#################################################\n";

# Remove files and directory created by previous test
$testingDir="$toggle/dataTest/checkFormatSamOrBam";
$cleaningCmd="rm -Rf $testingDir";
system ($cleaningCmd) and die ("ERROR: $0 : Cannot remove the previous test directory with the command $cleaningCmd \n$!\n");

#Creating config file for this test
@listSoft = ("checkFormatSamOrBam");
fileConfigurator::createFileConf(\@listSoft,"blockTestConfig.txt");

$runCmd = "toggleGenerator.pl -c blockTestConfig.txt -d ".$dataOneBam." -r ".$dataRefIrigin." -o ".$testingDir;
print "\n### Toggle running : $runCmd\n";
system("$runCmd") and warn "#### ERROR : Can't run TOGGLE for checkFormatSamOrBam";