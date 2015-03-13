#!/usr/bin/perl -w

###################################################################################################################################
#
# Licencied under CeCill-C (http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html) and GPLv3
#
# Intellectual property belongs to IRD, CIRAD and South Green developpement plateform 
# Written by Cecile Monat, Ayite Kougbeadjo, Mawusse Agbessi, Christine Tranchant, Marilyne Summo, Cedric Farcy, Francois Sabot
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
