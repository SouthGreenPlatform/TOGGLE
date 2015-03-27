#!/opt/perl-5.16.2/bin/perl
# script to execute tophhat module
###################################################################################################################################
#
# Copyright 2014 IRD-CIRAD
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
# Intellectual property belongs to IRD, CIRAD and South Green developpement plateform
# Written by Cecile Monat, Christine Tranchant, Ayite Kougbeadjo, Cedric Farcy, Mawusse Agbessi, Marilyne Summo, and Francois Sabot
#
###################################################################################################################################

use strict;
use warnings;

use lib qw(/data/projects/Floripalm/STAGE-SOUHILA/TOGGLE/Modules);
use localConfig;
use toolbox;
use tophat;
use Data::Dumper;
my $fileConf = $ARGV[0];
my $refFile = $ARGV[1];
my $forwardFile = $ARGV[2];
my $reverseFile = $ARGV[3];
my $gffFile = $ARGV[4];

$forwardFile =~ /^(.+)\/.+\./; #$forwardFile =~ /^.+\/(.+)\./; 
my $dirOut = $1."/tophat";
#######################################
#Creating the IndividuSoft.txt file
#######################################
my $creatingCommand="echo \"tophat\nTEST\" > individuSoft.txt";
system($creatingCommand) and die ("ERROR: $0: Cannot create the individuSoft.txt file with the command $creatingCommand \n$!\n");

#######################################
#Cleaning the logs for the test
#######################################
my $cleaningCommand="rm -Rf tophat_TEST_log.*";
system($cleaningCommand) and die ("ERROR: $0: Cannot clean the previous log files for this test with the command $cleaningCommand \n$!\n");

#######################################
#Gestion des options 
#######################################
#Chargement des options de chaque logiciel d<8e>finies dans le fichier softwareConfig dans une table de hachage
my $optionref = toolbox::readFileConf($fileConf);

#Recuperer les options propres <88> bowtie-build
my $softParameters = toolbox::extractHashSoft($optionref,"bowtie-build");   # get the bowtie-index option defined in software file and stored it into the variable $softParameters
tophat::bowtieBuild($refFile);

#Recuperer les options propres <88> bowtie2-build
$softParameters = toolbox::extractHashSoft($optionref,"bowtie2-build");   # get the bowtie2-index option defined in software file and stored it into the variable $softParameters
my $refIndex=tophat::bowtie2Build($refFile);

$softParameters = toolbox::extractHashSoft($optionref,"tophat2");
#tophat::tophat2($dirOut, $refFile, $forwardFile, $reverseFile, $gffFile, $softParameters);    # get the tophat option defined in software file and stored it into the variable $softParameters
tophat::tophat2($dirOut, $refIndex, $forwardFile, $reverseFile, $gffFile,$softParameters);# get the tophat option defined in software file and stored it into the variable $softParameters