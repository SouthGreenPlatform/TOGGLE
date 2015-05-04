#!/usr/bin/env perl



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
use fastqUtils;
use gatk;
use pairing;
use picardTools;
use samTools;
use toolbox;


##########################################
# recovery of initial informations/files
##########################################
my $initialDir = $ARGV[0];                                                                                  # recovery of the name of the directory to analyse
my $fileConf = $ARGV[1];                                                                                    # recovery of the name of the software.configuration.txt file
my $refFastaFile = $ARGV[2];                                                                                # recovery of the reference file
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



### Create the Arborescence
toolbox::makeDir("$initialDir/0_PAIRING_FILES/");
toolbox::makeDir("$initialDir/1_FASTQC/");
toolbox::makeDir("$initialDir/2_CUTADAPT/");
toolbox::makeDir("$initialDir/4_BWA/");
toolbox::makeDir("$initialDir/5_PICARDTOOLS/");
toolbox::makeDir("$initialDir/6_SAMTOOLS/");
toolbox::makeDir("$initialDir/7_GATK/");


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
my $numberOfReads = toolbox::checkNumberLines($listOfFiles[0])/4;                                           # check the number of sequences
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
my $optionref = toolbox::readFileConf($fileConf);                                                           # recovery of option for cutadapt
my $softParameters = toolbox::extractHashSoft($optionref,"cutadapt");
##DEBUG print LOG "DEBUG: optionref\n";
##DEBUG print LOG Dumper ($optionref);
cutadapt::createConfFile($fileAdaptator, $cutadaptSpecificFileConf, $softParameters);                            # create the configuration file specific to cutadapt software
my $fileWithoutExtention = toolbox::extractName($listOfFiles[0]);                                           # extract name of file without the extention
my $fileOut = "$newDir"."/"."$fileWithoutExtention".".CUTADAPT.fastq";                                      # name for the output file of cutadapt execution
print LOG "INFOS: $0 : Start cutadapt execution on file $listOfFiles[0]\n";
cutadapt::execution($listOfFiles[0],$cutadaptSpecificFileConf,$fileOut);                                    # run cutadapt program on current file

#########################################
# BWA INDEX, ALN and SAMSE
#########################################
print LOG "----------------------------------------\n";
print LOG "INFOS: $0 : Start BWA index and aln\n";
print F1 "BWA\n";
$newDir = toolbox::changeDirectoryArbo($initialDir,4);                                                      # change for the bwa direcotry
##DEBUG print LOG "CHANGE DIRECTORY TO $newDir\n";
$softParameters = toolbox::extractHashSoft($optionref,"BWA index");                                      # recovery of specific parameters of bwa index
bwa::bwaIndex($refFastaFile,$softParameters);                                                               # indexation of Reference sequences file
$softParameters = toolbox::extractHashSoft($optionref,"BWA aln");                                           # recovery of specific parameters of bwa aln
$fileWithoutExtention = toolbox::extractName($listOfFiles[0]);                                              # extract name of file without the extention
my $saiFileOut = "$newDir"."/"."$fileWithoutExtention".".BWAALN.sai";                                       # name for the output file of bwa aln
bwa::bwaAln($refFastaFile,$fileOut,$saiFileOut,$softParameters);                                     # find the SA coordinates of the current file
my $samFileOut = "$newDir"."/"."$fileWithoutExtention".".BWASAMSE.sam";                                     # name for the output file of bwa samse
$softParameters = toolbox::extractHashSoft($optionref,"BWA samse");                                         # recovery of specific parameters of bwa samse

@fileAndPath = toolbox::extractPath($listOfFiles[0]); 
#my $infoForRG = pairing::pairRecognition($fileAndPath[1]);                                                      # recovery of ReadGroup information
bwa::bwaSamse($samFileOut,$refFastaFile,$saiFileOut,$listOfFiles[0],$fileWithoutExtention,$softParameters);            # generate alignement in SAM format

#########################################
# PICARD CREATE SEQUENCE DICTIONARY and SORT SAM
#########################################
print LOG "----------------------------------------\n";
print LOG "INFOS: $0 : Start picard create sequence dictionary and sort sam\n";
print F1 "picardTools\n";
$newDir = toolbox::changeDirectoryArbo($initialDir,5);                                                      # change for the picardtools directory
##DEBUG print LOG "CHANGE DIRECTORY TO $newDir\n";
my $dictFileOut=$refFastaFile;                                                                              # name of dictionary file    
$dictFileOut =~ s/\.[^\.]*$/.dict/;
##DEBUG print LOG "INFOS: $0 : dict filename: $dictFileOut\n";
$softParameters = toolbox::extractHashSoft($optionref,"picardTools createSequenceDictionary");              # recovery of specific parameters of picard create sequence dictionary
picardTools::picardToolsCreateSequenceDictionary($refFastaFile,$dictFileOut,$softParameters);               # create ".dict" of the reference file
$fileWithoutExtention = pairing::extractName($samFileOut);                                                  # extract name of file without the extention
my $bamFileOut = "$newDir"."/"."$fileWithoutExtention".".PICARDTOOLSSORT.bam";                              # name for the output file of picardtools sort sam
$softParameters = toolbox::extractHashSoft($optionref,"picardTools sortsam single");                        # recovery of specific parameters of picard sort sam
picardTools::picardToolsSortSam($samFileOut,$bamFileOut,$softParameters);                                   # convert from SAM format to BAM format

#########################################
# SAMtools INDEX, VIEW, INDEX and FAIDX
#########################################
print LOG "----------------------------------------\n";
print LOG "INFOS: $0 : Start SAMtools index1, view and index2\n";
print F1 "SAMtools\n";
$newDir = toolbox::changeDirectoryArbo($initialDir,6);                                                      # change for the samtools directory
##DEBUG print LOG "CHANGE DIRECTORY TO $newDir\n";
samTools::samToolsIndex($bamFileOut);                                                                       # indexation of BAM file
$softParameters = toolbox::extractHashSoft($optionref,"samtools view single");                              # recovery of specific parameters of samtools view pair
my $bamFileIn = $bamFileOut;                                                                                # passing input for this software from the output of the previous one
$bamFileOut = "$newDir"."/"."$fileWithoutExtention".".SAMTOOLSVIEW.bam";                                    # name for the output file of samtools view
samTools::samToolsView($bamFileIn,$bamFileOut,$softParameters);                                             # extraction of not clean mapped alignement
samTools::samToolsIndex($bamFileOut);                                                                       # indexation of BAM file after cleaning
samTools::samToolsFaidx($refFastaFile);                                                                     # create ".fai" of the reference file

#########################################
# GATK REALIGNER TARGET CREATOR and INDEL REALIGNER
#########################################
print LOG "----------------------------------------\n";
print LOG "INFOS: start gatk realigner target creator and indel realigner\n";
print F1 "gatk\n";
$newDir = toolbox::changeDirectoryArbo($initialDir,7);                                                      # change for the gatk directory
##DEBUG print LOG "CHANGE DIRECTORY TO $newDir\n";
$softParameters = toolbox::extractHashSoft($optionref,"GATK gatkRealignerTargetCreator");                   # recovery of specific parameters of GATK Realigner Target Creator
$fileWithoutExtention = pairing::extractName($bamFileOut);                                                  # extract name of file without extention
my $intervalsFile = "$newDir"."/"."$fileWithoutExtention".".GATKREALIGNERTARGETCREATOR.intervals";          # name for the output file of gatk realigner target creator
gatk::gatkRealignerTargetCreator($refFastaFile, $bamFileOut, $intervalsFile, $softParameters);              # determine (small) suspicious intervals which are likely in need of realignment
$softParameters = toolbox::extractHashSoft($optionref,"GATK gatkIndelRealigner");                           # recovery of specific parameters of GATK Indel Realigner
my $bamRealigned = "$newDir"."/"."$fileWithoutExtention".".GATKINDELREALIGNER.bam";                         # name for the output file of gatk indel realigner
gatk::gatkIndelRealigner($refFastaFile, $bamFileOut, $intervalsFile, $bamRealigned, $softParameters);       # run the realigner over the intervals producted by gatk::gatkRealignerTargetCreator (see above)

#########################################
# SAMtools INDEX
#########################################
print LOG "----------------------------------------\n";
print LOG "INFOS: $0 : Start SAMtools index3\n";
print F1 "SAMtools\n";
$newDir = toolbox::changeDirectoryArbo($initialDir,6);                                                      # change for the samtools directory
##DEBUG print LOG "CHANGE DIRECTORY TO $newDir\n";
samTools::samToolsIndex($bamRealigned);                                                                     # indexation of BAM file after local realignment

#########################################
# PICARD MARK DUPLICATES
#########################################
print LOG "----------------------------------------\n";
print LOG "INFOS: $0 : Start picard mark duplicates\n";
print F1 "picardTools\n";
$newDir = toolbox::changeDirectoryArbo($initialDir,5);                                                      # change for the picardtools directory
##DEBUG print LOG "CHANGE DIRECTORY TO $newDir\n";
$softParameters = toolbox::extractHashSoft($optionref,"picardTools markDuplicates");                        # recovery of specific parameters of picard mark duplicates
$fileWithoutExtention = pairing::extractName($bamRealigned);                                                # extract name of file without extention
my $bamAnalyzed = "$newDir"."/"."$fileWithoutExtention".".PICARDTOOLSMARKDUPLICATES.bam";                   # name for the output file of picard mark duplicates
my $bamDuplicates = "$newDir"."/"."$fileWithoutExtention".".PICARDTOOLSMARKDUPLICATES.bamDuplicates";       # name for the output file of duplicates finded by picard mark duplicates
picardTools::picardToolsMarkDuplicates($bamRealigned, $bamAnalyzed, $bamDuplicates, $softParameters);       # examines aligned records in the supplied BAM file to locate duplicate molecules. All records are then written to the output file with the duplicate records flagged

#########################################
# SAMtools INDEX
#########################################
print LOG "----------------------------------------\n";
print LOG "INFOS: $0 : Start SAMtools index4\n";
print F1 "SAMtools\n";
$newDir = toolbox::changeDirectoryArbo($initialDir,6);                                                      # change for the samtools directory
##DEBUG print LOG "CHANGE DIRECTORY TO $newDir\n";
samTools::samToolsIndex($bamAnalyzed);                                                                      # indexation of BAM file after duplicate molecule flagged


print LOG "#########################################\nINFOS: Single sequence analysis done correctly\n#########################################\n";
close F1;
close LOG;
exit;
