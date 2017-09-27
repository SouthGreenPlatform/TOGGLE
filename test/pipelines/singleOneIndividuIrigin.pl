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
use lib qw(./);
use sedToggle;

#####################
## PIPELINE TEST
#####################

## PATH for datas test
# references files
my $dataRefIrigin = "$toggle/data/Bank/referenceIrigin.fasta";

#####################
## TOGGLE fastq singleOneIndividuIrigin
#####################


print "\n\n#################################################\n";
print "#### TEST SNPdiscoverySingle Irigin (one individu) / no SGE mode\n";
print "#################################################\n";

# Copy file config
my $fileSNPSingleIni="$toggle/exampleConfigs/SNPdiscoverySingle.config.txt";          # Path of the SNPdiscoveryPaired.config.txt
my $fileSNPSingleNoSGE="$toggle/test/pipelines/SNPdiscoverySingleTest.config.txt";

my $cmd="cp $fileSNPSingleIni $fileSNPSingleNoSGE";
## DEBUG print "\n### COPY conf file SNPdiscoveryPaired : $cmd\n";
system($cmd) and die ("#### ERROR COPY CONFIG FILE: $cmd\n");     # Copy into TEST

# Change the TOGGLE addaptator configuration file
sedToggle::sedFunction($fileSNPSingleNoSGE);

# data
my $dataFastqsingleOneIndividuIrigin = "$toggle/data/testData/fastq/singleOneIndividuIrigin";

# Remove files and directory created by previous test
my $testingDir="$toggle/dataTest/singleOneIndividuIrigin-noSGE";
my $cleaningCmd="rm -Rf $testingDir";
system ($cleaningCmd) and die ("ERROR: $0 : Cannot remove the previous test directory with the command $cleaningCmd \n$!\n");

my $runCmd = "toggleGenerator.pl -c ".$fileSNPSingleNoSGE." -d ".$dataFastqsingleOneIndividuIrigin." -r ".$dataRefIrigin." -o ".$testingDir;
print "\n### Toggle running : $runCmd\n";
system("$runCmd") and die "#### ERROR : Can't run TOGGLE for singleOneIndividusIrigin no SGE mode";


# check final results
print "\n### TEST Ouput list & content : $runCmd\n";
my $observedOutput = `ls $testingDir/finalResults`;
my @observedOutput = split /\n/,$observedOutput;
my @expectedOutput = ('multipleAnalysis.GATKSELECTVARIANT.vcf','multipleAnalysis.GATKSELECTVARIANT.vcf.idx');

# expected output test
is_deeply(\@observedOutput,\@expectedOutput,'toggleGenerator - singleOneIndividu (no SGE) list ');

# expected output content
$observedOutput=`tail -n 1 $testingDir/finalResults/multipleAnalysis.GATKSELECTVARIANT.vcf`;
chomp $observedOutput;
my $expectedOutput="#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	irigin2";
is($observedOutput,$expectedOutput, 'toggleGenerator - singleOneIndividu (no SGE) content ');

