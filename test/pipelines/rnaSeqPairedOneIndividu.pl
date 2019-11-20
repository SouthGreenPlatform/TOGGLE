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
use lib qw(./);
use sedToggle;

#####################
## PIPELINE TEST
#####################

## PATH for datas test
# references files
my $dataRefRnaseq = "$toggle/data/Bank/referenceRnaseq.fa";
my $dataRefRnaseqGFF = "$toggle/data/Bank/referenceRnaseqGFF.gff3";

#####################
## TOGGLE RNASeq pairedOneIndividu
#####################


print "\n\n#################################################\n";
print "#### TEST RNASEQPaired  (one individu) / no SGE mode\n";
print "#################################################\n";

# Copy file config
my $fileRNAPairedIni="$toggle/exampleConfigs/RNASeq.config.txt";
my $fileRNAPairedNoSGE="$toggle/test/pipelines/RNASeqNoSGE.config.txt";

my $cmd="cp $fileRNAPairedIni $fileRNAPairedNoSGE";
## DEBUG print "\n### COPY conf file SNPdiscoveryPaired : $cmd\n";
system($cmd) and die ("#### ERROR COPY CONFIG FILE: $cmd\n");     # Copy into TEST

# Change the TOGGLE addaptator configuration file
sedToggle::sedFunction($fileRNAPairedNoSGE);

# data
my $dataRNAseqPairedOneIndividu = "$toggle/data/testData/rnaseq/pairedOneIndividu";

# Remove files and directory created by previous test
my $testingDir="$toggle/dataTest/RNAseq-pairedOneIndividu-noSGE";
my $cleaningCmd="rm -Rf $testingDir";
system ($cleaningCmd) and die ("ERROR: $0 : Cannot remove the previous test directory with the command $cleaningCmd \n$!\n");

my $runCmd = "toggleGenerator.pl -c ".$fileRNAPairedNoSGE." -d ".$dataRNAseqPairedOneIndividu." -r ".$dataRefRnaseq." -g ".$dataRefRnaseqGFF." -o ".$testingDir;
print "\n### Toggle running : $runCmd\n";
system("$runCmd") and die "#### ERROR : Can't run TOGGLE for pairedOneIndividuRNASEQ no SGE mode";

# check final results
print "\n### TEST Ouput list & content : $runCmd\n";
my $observedOutput = `ls $testingDir/finalResults`;
my @observedOutput = split /\n/,$observedOutput;
my @expectedOutput = ('RNASeq.accepted_hits.HTSEQCOUNT.txt');

# expected output test
is_deeply(\@observedOutput,\@expectedOutput,'toggleGenerator - pairedOneIndividuRNASEQ (no SGE) list ');

# expected output content
$observedOutput=`wc -l $testingDir/finalResults/RNASeq.accepted_hits.HTSEQCOUNT.txt`;
chomp $observedOutput;
my $expectedOutput="2388 $testingDir/finalResults/RNASeq.accepted_hits.HTSEQCOUNT.txt";
is($observedOutput,$expectedOutput, 'toggleGenerator - pairedOneIndividuRNASEQ (no SGE) content ');
