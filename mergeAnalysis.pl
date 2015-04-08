#!/opt/perl-5.16.2/bin/perl



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
use lib qw(./Modules);
use localConfig;
use Data::Dumper;

use gatk;
use samTools;
use toolbox;




##########################################
# recovery of parameters/arguments given when the program is executed
##########################################
my $cmd_line=$0." @ARGV";
unless ($#ARGV>=0)                                                                                          # if no argument given
{
  my ($nomprog)=$0=~/([^\/]+)$/;
  print <<"Mesg";

  perldoc $nomprog display the help

Mesg

  exit;
}

my %param = @ARGV;                                                                                          # get the parameters 


##########################################
# recovery of initial informations/files
##########################################
my $initialDir = $param{'-d'};                                                                              # recovery of the name of the directory to analyse
my $fileConf = $param{'-c'};                                                                                # recovery of the name of the software.configuration.txt file
my $refFastaFile = $param{'-r'};                                                                            # recovery of the reference file
toolbox::existsDir($initialDir);                                                                            # check if this directory exists



##########################################
# Creation of IndividuSoft.txt for creation of logs files later
##########################################
my @pathIndividu = toolbox::extractPath($initialDir);
my $individu = $pathIndividu[0];
chdir "$initialDir";
my $infosFile = "individuSoft.txt";
#my $infosFile = "$pathIndividu[1]/individuSoft.txt";
open (F1, ">",$infosFile) or die ("ERROR: $0 : Cannot open the file $infosFile\n$!\n");
print F1 "$individu\n";
print F1 "Initialization\n";


my $indivName = `head -n 1 individuSoft.txt`;
chomp $indivName;

my $logFile=$indivName."_Global"."_log";
open (LOG, ">",$logFile) or die ("ERROR: $0 : Cannot open the file $logFile\n$!\n");
print LOG "#########################################\nINFOS: Merge analysis started\n#########################################\n\n";

toolbox::checkFile($fileConf);                                                                              # check if this file exists
toolbox::checkFile($refFastaFile);                                                                          # check if the reference file exists

##########################################
# gatk HAPLOTYPE CALLER, VARIANT FILTRATION and SELECT VARIANTS
##########################################
print LOG "INFOS: $0 : Start GATK Haplotype Caller, Select Variants and Variant Filtration\n";
print F1 "gatk\n";

my $listOfBam = toolbox::readDir($initialDir,"bam");                                                        # read the BAM directory to recover files in it
my $optionref = toolbox::readFileConf($fileConf);                                                           # read the configuration file
my $vcfCalled = "$initialDir"."/GATKHAPLOTYPECALLER.vcf";                                                   # name of the first VCF file
my $softParameters = toolbox::extractHashSoft($optionref, "GATK gatkHaplotypeCaller");                         # recovery of specific parameters of gatk haplotype caller
##DEBUG print "###DEBUG : $softParameters";
gatk::gatkHaplotypeCaller($refFastaFile, $vcfCalled, $listOfBam, $softParameters);

my $vcfFiltered= "$initialDir"."/GATKVARIANTFILTRATION.vcf";
$softParameters = toolbox::extractHashSoft($optionref, "GATK gatkVariantFiltration");
gatk::gatkVariantFiltration($refFastaFile, $vcfFiltered, $vcfCalled, $softParameters);

my $vcfVariantsSelected= "$initialDir"."/GATKSELECTVARIANTS.vcf";
$softParameters = toolbox::extractHashSoft($optionref, "GATK gatkSelectVariants");
gatk::gatkSelectVariants($refFastaFile, $vcfFiltered, $vcfVariantsSelected, $softParameters);


print LOG "#########################################\nINFOS: Merge analysis done correctly\n#########################################\n";


close LOG;
close F1;

exit;