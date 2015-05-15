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


use bwa;
use cutadapt;
use fastqc;
use fastxToolkit;
use fastqUtils;
use toolbox;
use tophat;
use samTools;
use HTSeq;


##########################################
# recovery of initial informations/files
##########################################
my $initialDir = $ARGV[0];                                                                                  # recovery of the name of the directory to analyse
my $fileConf = $ARGV[1];                                                                                    # recovery of the name of the software.configuration.txt file
my $refFastaFile = $ARGV[2];                                                                                # recovery of the reference file
my $gffFile = $ARGV[3];

toolbox::existsDir($initialDir);                                                                            # check if this directory exists




##########################################
# Creation of IndividuSoft.txt for creation of logs files later
##########################################
my @pathIndividu = toolbox::extractPath($initialDir);
my $individu = $pathIndividu[0];
chdir "$initialDir";
my $infosFile = "individuSoft.txt";
#my $infosFile = "$pathIndividu[1]/individuSoft.txt";
open (F1, ">$infosFile") or die ("ERROR: $0 : Cannot open the file $infosFile\n$!\n");
print F1 "$individu\n";
print F1 "Initialization\n";



my $indivName = `head -n 1 individuSoft.txt`;
chomp $indivName;

my $logFile=$indivName."_Global"."_log";
open (LOG, ">",$logFile) or die ("ERROR: $0 : Cannot open the file $logFile\n$!\n");
print LOG "#########################################\nINFOS: Single sequence analysis started\n#########################################\n\n";



toolbox::checkFile($fileConf);                                                                              # check if this file exists
toolbox::checkFile($refFastaFile);                                                                          # check if the reference file exists

my $optionref = toolbox::readFileConf($fileConf);                                                           # recovery of option for all softwares written in the $fileConf


### Create the Arborescence
toolbox::makeDir("$initialDir/0_PAIRING_FILES/");
toolbox::makeDir("$initialDir/1_FASTQC/");
toolbox::makeDir("$initialDir/11_FASTXTRIMMER/");
toolbox::makeDir("$initialDir/2_CUTADAPT/");
toolbox::makeDir("$initialDir/4_MAPPING/");



### Copy fastq into 0_PAIRING_FILES/    
my $copyCom = "cp $initialDir/*.fastq $initialDir/0_PAIRING_FILES/.";                               # command to move the initial fastq files into the directory appropriate for the pipeline
toolbox::run($copyCom);                                                                             # move the files
  
my $removeCom = "rm $initialDir/*.fastq";                                                            # command to remove the files in the main directory
toolbox::run($removeCom);


### Check the number of fastq
my $listOfFiles = toolbox::readDir($initialDir."/0_PAIRING_FILES/");                               # read it to recover files in it
##DEBUG print LOG "INFOS toolbox::ReadDir : @$listOfFiles\n";
my @listOfFiles = @$listOfFiles;

if (@$listOfFiles != 1)                                                                            # check if the number of files in the direcory is two as excpected for pair analysis
{
    toolbox::exportLog("ERROR: $0 : The directory ".$initialDir."/0_PAIRING_FILES/ don't contain the right number of files\n",1);
}




##########################################
# CHECK ENCODE OF FILES TO ANALYSE
##########################################
print LOG "----------------------------------------\n";
print LOG "INFOS: $0 : Checking sequences number on file $listOfFiles[0]\n";
print F1 "checkNumberLines\n";
toolbox::checkFile($listOfFiles[0]);                                                                        # check that the file exists, is readble, writable and has something in
my $numberOfReads = (toolbox::checkNumberLines($listOfFiles[0]))/4;                                           # check the number of sequences
print LOG "INFOS: $0 : Number of reads: $numberOfReads\n";
my $phred33Control = fastqUtils::checkEncodeByASCIIcontrol($listOfFiles[0]);                                # check the encode format (PHRED 33 or 64)
print LOG "INFOS: $0 : Return 1 if PHRED33, 0 if PHRED64: $phred33Control\n";
#### a tester avec fichier PHRED33 
    if ( $phred33Control == 1)                                                                              # if encode format is PHRED 33 then convert it in PHRED 64
    {
        fastqUtils::convertLinePHRED33ToPHRED64($listOfFiles[0]);                                           # change encode from PHRED 33 to PHRED 64
    }
####

#########################################
# FASTQC
#########################################
print LOG "----------------------------------------\n";
print LOG "INFOS: $0 : Start FASTQC\n";
print F1 "FastQC\n";
my @fileAndPath = toolbox::extractPath($listOfFiles[0]);                                                    # recovery of file name and path to have it
print LOG "INFOS: $0 : File: $fileAndPath[0]\n";
##DEBUG print LOG "INFOS extract path: $fileAndPath[1]\n";
my $newDir = toolbox::changeDirectoryArbo($initialDir,1);                                                   # change for the FASTQC directory
##DEBUG print LOG "CHANGE DIRECTORY TO $newDir\n";
fastqc::execution($listOfFiles[0],$newDir);                                                                 # run fastQC program on current file
my $fastqcStats = fastqc::parse($newDir);                                                                   # parse fastQC file to get statistics of fastqc
print LOG "INFOS: $0 : Statistics of fastqc:\n";
print LOG Dumper ($fastqcStats);

#########################################
# fastxToolkit::fastxTrimmer                                                                        
#########################################
print LOG "INFOS: $0 : start fastx_trimmer\n";
print F1 "fastxTrimmer\n";                                             
print LOG "INFOS: $0 : File: $listOfFiles[0]\n";
##DEBUG print LOG "INFOS extract path: $listOfFiles[0]\n";

                         
my $softParameters = toolbox::extractHashSoft($optionref,"fastx_trimmer");                            # get options for fastqx_trimmer
                               
$newDir = toolbox::changeDirectoryArbo($initialDir,11);                                                # change for the trimmer directory
##DEBUG print LOG "CHANGE DIRECTORY TO $newDir\n";
my $fileBasename=toolbox::extractName($listOfFiles[0]);                                            # get the fastq filename without the complete path and the extension
my $fileTrimmed= $newDir."/".$fileBasename.".FASTXTRIMMER.fastq";                                     # define the filename generated by trimmer 
fastxToolkit::fastxTrimmer($listOfFiles[0],$fileTrimmed,$softParameters);

#########################################
# CUTADAPT CREATE CONF FILE and EXECUTION
#########################################
print LOG "----------------------------------------\n";
print LOG "INFOS: $0 : Start cutadapt create configuration file\n";
print F1 "cutadapt\n";
$newDir = toolbox::changeDirectoryArbo($initialDir,2);                                                   # change for the cutadapt directory
##DEBUG print LOG "CHANGE DIRECTORY TO $newDir\n";
my $fileAdaptator = "$toggle/adaptator.txt";     # /!\ ARGV[3] et si non reseigné ce fichier là, mais on le place où ?
toolbox::checkFile($fileAdaptator);
my $cutadaptSpecificFileConf = "$newDir"."/cutadapt.conf";                                                  # name for the cutadapt specific configuration file
$softParameters = toolbox::extractHashSoft($optionref,"cutadapt");
##DEBUG print LOG "DEBUG: optionref\n";
##DEBUG print LOG Dumper ($optionref);

cutadapt::createConfFile($fileAdaptator, $cutadaptSpecificFileConf, $softParameters);                            # create the configuration file specific to cutadapt software
my $trimmedFiles=toolbox::readDir($initialDir."/11_FASTXTRIMMER/");
my @trimmedFiles=@$trimmedFiles;
#DEBUG print LOG Dumper(@trimmedFiles);

my $fileWithoutExtention = toolbox::extractName($trimmedFiles[0]);                                           # extract name of file without the extention
my $fileCutadaptOut = "$newDir"."/"."$fileWithoutExtention".".CUTADAPT.fastq";                                      # name for the output file of cutadapt execution
print LOG "INFOS: $0 : Start cutadapt execution on file $trimmedFiles[0]\n";
cutadapt::execution($trimmedFiles[0],$cutadaptSpecificFileConf,$fileCutadaptOut);                                    # run cutadapt program on current file


##########################################
# tophat::bowtieBuild
##########################################
print LOG "----------------------------------------\n";
print LOG "INFOS: $0 : Start BOWTIE inde\n";
print F1 "BOWTIEBUILD\n";
$newDir = toolbox::changeDirectoryArbo($initialDir,4);                                                      # change for the tophat direcotry
my $tophatDir = $newDir;

##DEBUG
print LOG "CHANGE DIRECTORY TO $newDir\n";
$softParameters = toolbox::extractHashSoft($optionref, "bowtieBuild");                              # recovery of specific parameters of bowtiebuild index
tophat::bowtieBuild($refFastaFile,$softParameters);                                           # indexation of Reference sequences file


##########################################
# tophat::bowtie2Build
##########################################
print LOG "----------------------------------------\n";
print LOG "INFOS: $0 : Start BOWTIE2-BUILD\n";
print F1 "BOWTIE2BUILD\n";
                          
##DEBUG
print LOG "CHANGE DIRECTORY TO $newDir\n";                          
$softParameters = toolbox::extractHashSoft($optionref, "bowtie2-build");                                      # recovery of specific parameters of tophat index
my $refIndex=tophat::bowtie2Build($refFastaFile,$softParameters);                                            # indexation of Reference sequences file               


##########################################
# tophat::tophat2
##########################################
print LOG "INFOS: $0 : start tophat2\n";
print F1 "tophat2\n";

my $tophatdirOut = $newDir;   #créer le répertoire des résultats de tophat
$softParameters = toolbox::extractHashSoft($optionref,"tophat2");  

print LOG "INFOS tophats argument: $tophatdirOut,$refIndex,$fileCutadaptOut,$gffFile";
my $fileReverse;
tophat::tophat2($tophatdirOut,$refIndex,$fileCutadaptOut,undef $fileReverse, $gffFile,$softParameters);            # generate alignement in SAM format


##########################################
# samTools::sort
##########################################
print LOG "INFOS: $0 : start samtools sort\n";
print F1 "samtools\n";
$softParameters = toolbox::extractHashSoft($optionref,"samtools sort");  
my $tophatBam=$tophatdirOut."/accepted_hits.bam";
samTools::samToolsSort($tophatBam,$softParameters);


##########################################
# HTSeq::htseqCount
##########################################
print LOG "INFOS: $0 : start htseq-count\n";
print F1 "htseqcount\n";
$softParameters = toolbox::extractHashSoft($optionref,"htseqcount");

my $htseqcountBam=$tophatdirOut."/accepted_hits.SAMTOOLSSORT.bam";
my $htseqcountOut=$tophatdirOut."/accepted_hits.HTSEQCOUNT.txt";
HTSeq::htseqCount($htseqcountBam,$htseqcountOut,$gffFile,$softParameters);


print LOG "#########################################\nINFOS: Single sequence analysis done correctly\n#########################################\n";
close F1;
close LOG;
exit;
