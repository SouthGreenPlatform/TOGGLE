#!/usr/bin/perl -w


###################################################################################################################################
#
# Licencied under CeCill-C (http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html) and GPLv3
#
# Intellectual property belongs to IRD, CIRAD and SouthGreen developpement plateform 
# Written by CŽcile Monat, Ayite Kougbeadjo, Mawusse Agbessi, Christine Tranchant, Marilyne Summo, CŽdric Farcy, Franois Sabot
#
###################################################################################################################################


=head1 Name

NGS.pl

=head1 Usage

NGS.pl -d directory with fastq -c configuration file -a adaptator file (cutadapt step)

=head1 Synopsis

Juste step fastqc + cutadapt

=head2 Inputs

directory with one or several fastq sequence

=head2 Outputs

fastqc outut
cutadapt output

=head1 Notes
version 0.0

=head1  Author
GGR TEAM, UMR DIADE - IRD

=cut

# perl NGS.pl -d /scratch/tranchant/TEST/ -c /scratch/tranchant/TEST/software.config.txt -a /scratch/tranchant/TEST/adaptator.txt

use strict;
use diagnostics;
use Data::Dumper;
use Getopt::Std;

use lib qw (../Modules);
use fastqc;
use cutadapt;
use local_config;
use toolbox;
use pairing;


#####################################################
my $cmd_line=$0." @ARGV  \n"; # recovery of program name with options

# Start of Log file
toolbox::exportLog("\n###### commande: $cmd_line",1);

unless ($#ARGV>=0)  # Check if there is one arguments at least
{
  my ($nomprog)=$0=~/([^\/]+)$/;
  print <<"Mesg";

   perldoc $nomprog display the help

Mesg

  exit;
}
#####################################################





#####################################################
############ VARIABLES DECLARATION ##################

#### Command-line arguments/option
#### VOIR SI METTRE OPTION PLUS LONGUE
my %param = @ARGV;  # recovery of options
@param{ map { lc $_ } keys %param } = values %param;  # make a list of your options

## Data directory with one ou several fastq
my $fastqDir=defined($param{'-d'})?$param{'-d'}:'';   # define the value of fastqDir according to options previously defined, if not defined, the value of fastqDir is empty

## File with option software
my $softConfFile=defined ($param{'-c'})? $param{'-c'} : $fastqDir."/software.config"; #define the value of softConfFile according to option previously defined, if not defined, search the softwre.config file in the data directory (default value)
toolbox::readFileConf("$softConfFile"); # read the confFile and create the variable $configInfos
print Dumper(\$configInfos);  # to informe user, and check if options are correct 
#### TODO: MAKE TEST IF THE VARIABLE $configInfos ISN'T EMPTY #####

### File with adaptator sequences / CUTADAPT STEP
my $adaptatorFile=defined ($param{'-a'})?$param{'-a'}:$fastqDir."/adaptator.txt"; ##### TODO: TO MOVE IN CONFIG FILE ######

## Other variable
my %stat;                       ## read info after fastqc step

my @fqList;                     ## fastq filename list stored in the data directory

#my $fqFile;

my $erase=1;                    ## Manage if output directory must be remove


######################################################
##############  PARAMETERS TEST ######################

toolbox::exportLog("\n###### INITIALIZATION STEP\n",1);
toolbox::exportLog("\n## REPERTOIRE DE DONNEES: $fastqDir\n",1);

## Test data directory existence
toolbox::existsDir($fastqDir); # directory exists?


##### TEST droit acces repertoire?

## Select fastq file in the input directory
my $fqList_ref = toolbox::readDir($fastqDir,"fastq");
@fqList=@{$fqList_ref};  # fastq filename list stored in the data directory

print STDERR Data::Dumper::Dumper(\@fqList);  ########## A COMMENTER

my $pairFile_ref= pairing::pairRecognition($fastqDir);
my %pairFile= %{$pairFile_ref};
print STDERR Data::Dumper::Dumper(\%pairFile);  ########## A COMMENTER

## Test adaptator file existence
toolbox::exportLog("\n## FICHIER ADAPTATEUR: $adaptatorFile\n",1);
toolbox::existsFile($adaptatorFile) and toolbox::readFile($adaptatorFile,"check");



######################################################
##############  OUTPUT DIRECTORY CREATION ############
toolbox::exportLog("\n###### OUTPUT DIRECTORY CREATION\n",1);


# create the differents directory necessary to the analysis
  my $fastqcDir=$fastqDir."1_FASTQC";
  toolbox::makeDir($fastqcDir, $erase);
  
  my $cutadaptDir=$fastqDir."2_CUTADAPT";
  toolbox::makeDir($cutadaptDir, $erase);


############### STEP CUTADAPT
## To create cutadapt configuration file
  my $cutadaptConf=$cutadaptDir."/cutadapt.conf";
  toolbox::exportLog("## Creation du fichier configuration $cutadaptConf\n",1);
  cutadapt::createConfFile($adaptatorFile,$cutadaptConf,$configInfos->{'cutadapt'});
  
  


foreach my $fqFile (@fqList)
{

  ## TODO : INTEGRER LE CHECK FORMAT FASTQC

  ######################################################
  ############### STEP FASTQC ##########################
  toolbox::exportLog("\n###### FASTQC STEP\n",1);
  toolbox::exportLog("## fastqc $fqFile\n",1);  # print the current file
  
  fastqc::exec($fqFile,$fastqcDir); # execute fastqc

  
  my ($fileName) = toolbox::extractPath($fqFile); # recovery of fastqFile name

  my $fastqcDirSeq = ($fileName =~ /^(.*)\.fastq$/) ? $fastqcDir."/".$1."_fastqc" : $fastqcDir."/".$fileName."_fastqc"; # define the name of new directory created by fastqc for the current fastq file
  my ($statRef)=fastqc::parse($fastqcDirSeq);      # recovery of fastqc stats informations
  $stat{$fqFile}={%$statRef};  # associate the fastqc stats informations to the current file


  #####################################################
  ############## STEP CUTADAPT ########################
  toolbox::exportLog("\n###### CUTADAPT STEP\n",1);

  toolbox::exportLog("## fastqc $fqFile\n",1);
  my $fastqCutadapt=$fqFile."_cutadapt.fq";
  cutadapt::exec($fqFile,$cutadaptConf,$fastqCutadapt); # execute cutadapt  
}

print STDERR Data::Dumper::Dumper(\%stat); ### TO COMMENT AFTER 

pairing::repairing("/teams/ggr/pipelineNGS/DATA-TEST/RC1_1.CUTADAPT.fastq","/teams/ggr/pipelineNGS/DATA-TEST/RC1_2.CUTADAPT.fastq");

#### STEP PAIRING ####
