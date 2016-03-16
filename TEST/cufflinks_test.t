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

use strict;

#Will test if bwa works correctly
use warnings;
use lib qw(../Modules/);
use Test::More 'no_plan'; #Number of tests, to modify if new tests implemented. Can be changed as 'no_plan' instead of tests=>11 .
use Test::Deep;
use Data::Dumper;


use_ok('toolbox') or exit;
use_ok('cufflinks') or exit;
can_ok( 'cufflinks','cufflinks');
can_ok( 'cufflinks','cuffmerge');
can_ok( 'cufflinks','cuffdiff');
use toolbox;
use cufflinks;
toolbox::readFileConf("software.config.txt");
#
#
#
####Test for cufflinks:cufflinks
my $optionsHachees=$configInfos->{'cufflinks cufflinks'};
my $refFasta="/home/al/Documents/test_module_cufflinks/ALL_MSU7.fa";
my $annotation="/home/al/Documents/test_module_cufflinks/MSU7.gff3";
print toolbox::readFile($annotation);
my $alignementFile="/home/al/Documents/test_module_cufflinks/AYR_HOSW_4_1/accepted_hits.bam";
#is(cufflinks::cufflinks($refFasta, $annotation,$alignementFile, $optionsHachees),1,'Ok for cufflinks::cufflinks');
####
###
###
#####Test for cufflinks::cuffmerge
$optionsHachees=$configInfos->{'cufflinks cuffmerge'};
my $inputFile="";
my $inputDir="/home/al/Documents/test_module_cufflinks/";
my $outdir="/home/al/Documents/test_module_cufflinks/cuffmerge";
is(cufflinks::cuffmerge($refFasta, $annotation, $inputDir, $inputFile, $outdir),1,'Ok for cufflinks::cuffmerge');
#####
####
####
#######Test for cufflinks::cuffdiff
#$optionsHachees=$configInfos->{'cufflinks cuffdiff'};
#my $mergedFile="";
#my $bamFileIn="";
#my $outdir="";
#is(cufflinks::cuffdiff($refFasta, $mergedFile, $bamFileIn, $annotation, $outdir, $optionsHachees),1,'Ok for gatk Variant Filtratrion');
