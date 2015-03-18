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
open (F1, ">",$infosFile) or die ("ERROR: $0 : Cannot open the file $infosFile\n$!\n");
print F1 "$individu\n";
print F1 "Initialization\n";


my $indivName = `head -n 1 individuSoft.txt`;
chomp $indivName;

my $logFile=$indivName."_Global"."_log";
open (LOG, ">",$logFile) or die ("ERROR: $0 : Cannot open the file $logFile\n$!\n");
print LOG "#########################################\nINFOS: Paired sequences analysis started\n#########################################\n\n";

toolbox::checkFile($fileConf);                                                                              # check if this file exists
toolbox::checkFile($refFastaFile);                                                                          # check if the reference file exists




### Create the Arborescence
toolbox::makeDir("$initialDir/0_PAIRING_FILES/");
toolbox::makeDir("$initialDir/1_FASTQC/");
toolbox::makeDir("$initialDir/2_CUTADAPT/");
toolbox::makeDir("$initialDir/3_PAIRING_SEQUENCES/");
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
my $initialsFiles = toolbox::readDir($initialDir."/0_PAIRING_FILES/");                               # read it to recover files in it
##DEBUG print LOG "INFOS: toolbox::ReadDir : @$initialsFiles\n";
##DEBUG print LOG "----------------------------------------\n";
my @initialsFiles = @$initialsFiles;

if (@$initialsFiles != 2)                                                                            # check if the number of files in the direcory is two as excpected for pair analysis
{
    toolbox::exportLog("ERROR: $0 : The directory ".$initialDir."/0_PAIRING_FILES/ don't contain the right number of files\n",1);
}





# to all files on
for (my $i=0; $i<=$#initialsFiles; $i++)
{
    ##########################################
    # fastqUtils::checkNumberByWC
    ##########################################
    print LOG "INFOS: $0 : Checking sequences number on file $initialsFiles[$i]\n";

    #print F1 "checkNumberByWC\n";
    my $numberOfReads = fastqUtils::checkNumberByWC($initialsFiles[$i]);                                      # check the number of sequences
    print LOG "INFOS: $0 : Number of reads: $numberOfReads\n";
    ##########################################
    # fastqUtils::checkEncodeByASCIIcontrol
    ##########################################
    my $phred33Control = fastqUtils::checkEncodeByASCIIcontrol($initialsFiles[$i]);                           # check the encode format (PHRED 33 or 64)
    print LOG "INFOS: $0 : Return 1 if PHRED33, 0 if PHRED64: $phred33Control\n";
    print LOG "----------------------------------------\n";
    ##########################################
    # fastqUtils::changeEncode
    ##########################################
    if ( not $phred33Control )                                                                          # if encode format is PHRED 64 then convert it in PHRED 33
    {
        print LOG "INFOS: $0 : Change encoding quality from 64 to 33\n";
        print LOG "----------------------------------------\n";
        my $tmpFile=$initialsFiles[$i]."sanger.tmp";
        fastqUtils::changeEncode($initialsFiles[$i],$tmpFile,64,33);                                    # change encode from PHRED 64 to PHRED 33
        my $moveCom="mv $tmpFile $initialsFiles[$i]";
        system($moveCom) and die ("ERROR: $0 : Cannot create the file with the good quality with the command $moveCom \n$!\n");     
    }
    ##########################################
    # fastqc::execution
    ##########################################
    print LOG "INFOS: $0 : Start FASTQC\n";
    print F1 "FastQC\n";
    my @fileAndPath = toolbox::extractPath($initialsFiles[$i]);                                               # recovery of file name and path to have it
    print LOG "INFOS: $0 : File: $fileAndPath[0]\n";
    ##DEBUG print LOG "INFOS extract path: $fileAndPath[1]\n";
    my $newDir = toolbox::changeDirectoryArbo($initialDir,1);                                                  # change for the FASTQC directory
    ##DEBUG print LOG "CHANGE DIRECTORY TO $newDir\n";
    fastqc::execution($initialsFiles[$i],$newDir);                                                            # run fastQC program on current file
    ##########################################
    # fastqc::parse
    ##########################################
    my $fastqcStats = fastqc::parse($newDir);                                                               # parse fastQC file to get statistics of fastqc
    print LOG "INFOS statistics of fastqc: \n";
    print LOG Dumper ($fastqcStats);
    print LOG "----------------------------------------\n";
}

##########################################
# cutadapt::createConfFile
##########################################
print LOG "INFOS: $0 : Start cutadapt create configuration file\n";
print F1 "cutadapt\n";
my $newDir = toolbox::changeDirectoryArbo($initialDir,2);                                                   # change for the cutadapt directory
my $cutadaptDir = $newDir;                                                                                  # to keep the information of cutadapt folder for the following analysis
##DEBUG print LOG "CHANGE DIRECTORY TO $newDir\n";
my $fileAdaptator = "$toggle/adaptator.txt";     # /!\ ARGV[3] et si non reseignÃ© ce fichier lÃ , mais on le place oÃ¹ ?
toolbox::checkFile($fileAdaptator);
my $cutadaptSpecificFileConf = "$newDir"."/cutadapt.conf";                                                  # name for the cutadapt specific configuration file
my $optionref = toolbox::readFileConf($fileConf);                                                           # recovery of option for cutadapt
my $softParameters = toolbox::extractHashSoft($optionref, "cutadapt");                                         # recovery of specific informations for cutadapt
##DEBUG print LOG "DEBUG: optionref\n";
##DEBUG print LOG Dumper ($optionref);
cutadapt::createConfFile($fileAdaptator, $cutadaptSpecificFileConf, $softParameters);                            # create the configuration file specific to cutadapt software



# A FAIRE SUR TOUT LES FICHIERS
for (my $i=0; $i<=$#initialsFiles; $i++)
{
    ##########################################
    # cutadapt::execution
    ##########################################
    print LOG "INFOS: $0 : Start cutadapt execution on file $initialsFiles[$i]\n";
    $newDir = toolbox::changeDirectoryArbo($initialDir,2);                                                  # change for the cutadapt directory
    ##DEBUG print LOG "CHANGE DIRECTORY TO $newDir\n";
    my $fileWithoutExtention = toolbox::extractName($initialsFiles[$i]);                                      # extract name of file without the extention
    my $fileOut = "$newDir"."/"."$fileWithoutExtention".".CUTADAPT.fastq";                                  # name for the output file of cutadapt execution
    cutadapt::execution($initialsFiles[$i],$cutadaptSpecificFileConf,$fileOut);                               # run cutadapt program on current file
}

##########################################
# pairing::repairing
##########################################
print LOG "----------------------------------------\n";
print LOG "INFOS: $0 : Start pairing sequences\n";
print F1 "PairingSequences\n";
$newDir = toolbox::changeDirectoryArbo($initialDir,3);                                                      # change for the pairing_sequences directory
my $repairingDir = $newDir;                                                                                 # to keep the information of repairing folder for the following analysis
##DEBUG print LOG "CHANGE DIRECTORY TO $newDir\n";
my $cutadaptList = toolbox::readDir($cutadaptDir);                                                          # read it to recover files in it
##DEBUG print LOG "INFOS toolbox ReadDir: @$cutadaptList\n";
my @cutadaptList = @$cutadaptList;
my @finalCutadaptList;


for (my $i=0; $i<=$#cutadaptList; $i++)
{
    if ($cutadaptList[$i]=~m/cutadapt.conf/)
    {
        next;
    }
    else
    {
        push (@finalCutadaptList,"$cutadaptList[$i]");
    }
}

pairing::repairing($finalCutadaptList[0],$finalCutadaptList[1],$newDir);                                                    # from two de-paired files, forward + reverse, will create three files, forward, reverse and singl

##########################################
# bwa::bwaIndex
##########################################
print LOG "----------------------------------------\n";
print LOG "INFOS: $0 : Start BWA index and aln\n";
print F1 "BWA\n";
$newDir = toolbox::changeDirectoryArbo($initialDir,4);                                                  # change for the bwa direcotry
my $bwaDir = $newDir;                                                                                   # to keep the information of bwa folder for the following analysis
##DEBUG print LOG "CHANGE DIRECTORY TO $newDir\n";
$softParameters = toolbox::extractHashSoft($optionref,"BWA index");                                  # recovery of specific parameters of bwa index
bwa::bwaIndex($refFastaFile,$softParameters);                                                           # indexation of Reference sequences file

# A FAIRE SUR TOUT LES FICHIERS
my $repairingList = toolbox::readDir($repairingDir);
##DEBUG print LOG "INFOS toolbox ReadDir: @$repairingList\n";
my @repairingList = @$repairingList;
for (my $i=0; $i<=$#repairingList; $i++)
{
    ##########################################
    # bwa::bwaAln
    ##########################################
    $softParameters = toolbox::extractHashSoft($optionref,"BWA aln");                                       # recovery of specific parameters of bwa aln
    my $fileWithoutExtention = toolbox::extractName($repairingList[$i]);                                         # extract name of file without the extention
    my $saiFileOut = "$newDir"."/"."$fileWithoutExtention".".BWAALN.sai";                                   # name for the output file of bwa aln
    bwa::bwaAln($refFastaFile,$repairingList[$i],$saiFileOut,$softParameters);                                # find the SA coordinates of the current file
}


##########################################
# bwa::bwaSampe
##########################################
print LOG "INFOS: $0 : Start BWA sampe\n";
$newDir = toolbox::changeDirectoryArbo($initialDir,4);                                                      # change for the bwa direcotry
##DEBUG print LOG "CHANGE DIRECTORY TO $newDir\n";
my $listOfSai = toolbox::readDir($bwaDir);                                                                  # read it to recover files in it
my @listOfSai = @$listOfSai;
my $fileWithoutExtention = pairing::extractName($listOfSai[0]);                                             # extract name of file without the extention
my $samFileOut = "$newDir"."/"."$fileWithoutExtention".".BWASAMPE.sam";                                     # name for the output file of bwa sampe
$softParameters = toolbox::extractHashSoft($optionref,"BWA sampe");                                      # recovery of specific parameters of bwa sampe
bwa::bwaSampe($samFileOut,$refFastaFile,$listOfSai[0],$listOfSai[1],$repairingList[0],$repairingList[1],$fileWithoutExtention,$softParameters); # generate alignement in SAM format

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
$softParameters = toolbox::extractHashSoft($optionref,"picardTools sortsam pair");                           # recovery of specific parameters of picard sort sam
picardTools::picardToolsSortSam($samFileOut,$bamFileOut,$softParameters);                                   # convert from SAM format to BAM format

#########################################
# SAMtools INDEX 1, VIEW, INDEX 2 and FAIDX
#########################################
print LOG "----------------------------------------\n";
print LOG "INFOS: $0 : Start SAMtools index1, view and index2\n";
print F1 "SAMtools\n";
$newDir = toolbox::changeDirectoryArbo($initialDir,6);                                                      # change for the samtools directory
##DEBUG print LOG "CHANGE DIRECTORY TO $newDir\n";
samTools::samToolsIndex($bamFileOut);                                                                       # indexation of BAM file
$softParameters = toolbox::extractHashSoft($optionref,"samtools view pair");                                # recovery of specific parameters of samtools view pair
my $bamFileIn = $bamFileOut;                                                                                # passing input for this software from the output of the previous one
$bamFileOut = "$newDir"."/"."$fileWithoutExtention".".SAMTOOLSVIEW.bam";                                    # name for the output file of samtools view
samTools::samToolsView($bamFileIn,$bamFileOut,$softParameters);                                             # extraction of not clean mapped alignement
samTools::samToolsIndex($bamFileOut);                                                                       # indexation of BAM file after cleaning
samTools::samToolsFaidx($refFastaFile);                                                                     # create ".fai" of the reference file


#########################################
# GATK REALIGNER TARGET CREATOR and INDEL REALIGNER
#########################################
print LOG "----------------------------------------\n";
print LOG "INFOS: $0 : Start gatk realigner target creator and indel realigner\n";
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
# SAMtools INDEX 3 
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
# SAMtools INDEX 4 
#########################################
print LOG "----------------------------------------\n";
print LOG "INFOS: $0 : Start SAMtools index4\n";
print F1 "SAMtools\n";
$newDir = toolbox::changeDirectoryArbo($initialDir,6);                                                      # change for the samtools directory
##DEBUG print LOG "CHANGE DIRECTORY TO $newDir\n";
samTools::samToolsIndex($bamAnalyzed);                                                                      # indexation of BAM file after duplicate molecule flagged




print LOG "#########################################\nINFOS: Paired sequences analysis done correctly\n#########################################\n";

close F1;
close LOG;

exit;
