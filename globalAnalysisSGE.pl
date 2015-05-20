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

use pairing;
use toolbox;

my $initialDir = $ARGV[0];                                                                                  # recovery of the name of the directory to analyse
my $fileConf = $ARGV[1];                                                                                    # recovery of the name of the software.configuration.txt file
my $refFastaFile = $ARGV[2];                                                                                # recovery of the reference file
my $cmd_line=$0." @ARGV";

my $optionref = toolbox::readFileConf($fileConf);
my $softParameters = toolbox::extractHashSoft($optionref, "SGE");   
my $options="";
$options=toolbox::extractOptions($softParameters); ##Get given options



my $jobList;

my $infosFile = "individuSoft.txt";

my ($sec, $min, $h, $mois_jour, $mois, $an, $sem_jour, $cal_jour, $heure_ete) = localtime(time);
$mois+=1;
$mois = $mois < 10 ? $mois = "0".$mois : $mois;
$mois_jour = $mois_jour < 10 ? $mois_jour = "0".$mois_jour : $mois_jour;
$h = $h < 10 ? $h = "0".$h : $h;
$min = $min < 10 ? $min = "0".$min : $min;
$an+=1900;
my $date="$mois_jour-$mois-$an-$h"."_"."$min";

#my $infosFile = "$pathIndividu[1]/individuSoft.txt";
open (F1, ">",$infosFile) or die ("$0 : open error of $infosFile .... $!\n");
print F1 "GLOBAL\n";
print F1 "ANALYSIS_$date\n";


toolbox::exportLog("#########################################\nINFOS: Global analysis \n#########################################\n",1);
toolbox::exportLog("----------------------------------------",1);
toolbox::exportLog("INFOS: $0 : Command line : $cmd_line\n",1);
toolbox::exportLog("----------------------------------------",1);
toolbox::checkFile($fileConf);                                                                              # check if this file exists
toolbox::existsDir($initialDir);                                                                            # check if this directory exists
toolbox::checkFile($refFastaFile);                                                                          # check if the reference file exists
my $loop = 0;                                                                                               # for the second loop



my $listOfFiles = toolbox::readDir($initialDir);                                                            # read it to recover files in it
##DEBUG toolbox::exportLog("INFOS toolbox ReadDir: @$listOfFiles\n",1);
my @listOfFiles = @$listOfFiles;



#########################################
# check if initial directory contain folder or not
#########################################
toolbox::exportLog("----------------------------------------",1);
my $folder = toolbox::checkInitialDirContent($initialDir);
if ($folder == 0)                                                                                           # if folder = 0, it's mean that there is only files in initial directory
{
    #########################################
    # recognition of pairs of files and create a folder for each pair
    #########################################
    my $pairsInfos = pairing::pairRecognition($initialDir);                                                 # from files fasta recognition of paired files
    pairing::createDirPerCouple($pairsInfos,$initialDir);                                                   # from infos of pairs, construction of the pair folder

    $listOfFiles = toolbox::readDir($initialDir);                                                           # read it to recover files in it
    toolbox::exportLog("INFOS: $0 : toolbox::readDir : $initialDir after create dir per couple: @$listOfFiles\n",1);
    @listOfFiles = @$listOfFiles;
}
else
{
    my @listOfFolder = @$folder;                                                                            # if contain folder, recovery if the list of them in this table
}

LOOP:
#########################################
# Creation of the numerous arborescence and lunch analysis single or pair as appropriate
#########################################
toolbox::exportLog("----------------------------------------",1);
for (my $i=0; $i<=$#listOfFiles; $i++)                                                                      # for each folder, create directories for different step of analysis
{
    if ($listOfFiles[$i]=~m/.+\..+/)                                                                        # if it's a file and not a folder
    {
        ##DEBUG toolbox::exportLog("FILES: $listOfFiles[$i]\n",1);
        next;
    }
    elsif ($listOfFiles[$i]=~m/.+\/.+:$/)
    {
        
        
      
        ##DEBUG toolbox::exportLog("FOLDER: $listOfFiles[$i]\n",1);
        my @fileAndPath = toolbox::extractPath($listOfFiles[$i]);                                           # recovery of file name and path to have it
        ##DEBUG toolbox::exportLog("INFOS extract file: $fileAndPath[0]\n",1);
        ##DEBUG toolbox::exportLog("INFOS extract path: $fileAndPath[1]\n",1);
        my @splitName = split (":", $fileAndPath[0]);
     
        my $firstDir = "$fileAndPath[1]$splitName[0]/";
        my $listOfFastq = toolbox::readDir($firstDir);                                                      # recovery of fastq file(s)
        ##DEBUG toolbox::exportLog("INFOS toolbox ReadDir: @$listOfFastq\n",1);
        my @listOfFastq = @$listOfFastq;
        if ($#listOfFastq == 0)                                                                             # if 1 file --> single analysis to do
        {
            toolbox::exportLog("INFOS: $0 : Run singleAnalysis.pl on $firstDir\n",1);
            my $singleCom = 'qsub -N singleAnalysis '.$options.' "singleAnalysis.pl '.$firstDir.' '.$fileConf.' '.$refFastaFile.'"';
            ##DEBUG
            toolbox::exportLog("DEBUG: $0 : qsub singleAnalysis command : $singleCom\n",1);
            my $job_id = `$singleCom`;
            
            toolbox::run("sleep 50");
            
            if ($job_id =~ /^[^\d]+(\d+)\s/)
            {
                $jobList.=$1."|";
                ##DEBUG
                toolbox::exportLog("DEBUG: $0 : $jobList\n",1);
            }
            else { toolbox::exportLog("ERROR: $0 : Problem with list of job id : $job_id.\n",0); }
        }
        elsif ($#listOfFastq == 1)                                                                          # if 2 files --> pair analysis to do
        {
            toolbox::exportLog("INFOS: $0 : Run pairAnalysis.pl on $firstDir\n",1);
            my $pairCom = 'qsub -N pairAnalysis '.$options.' "pairAnalysis.pl '.$firstDir.' '.$fileConf.' '.$refFastaFile.'"';
            ##DEBUG
            toolbox::exportLog("DEBUG: $0 : qsub pairAnalysis command : $pairCom\n",1);
            my $job_id = `$pairCom`;
            
            toolbox::run("sleep 50");
            
            if ($job_id =~ /^[^\d]+(\d+)\s/)
            {
                $jobList.=$1."|";
                ##DEBUG
                toolbox::exportLog("DEBUG: $0 : $jobList\n",1);
            }
            else { toolbox::exportLog("ERROR: $0 : Problem with list of job id : $job_id.\n",0); }
        
        }
        else                                                                                                # if more than 2 files, there is a problem
        {
            my $nbFiles = ($#listOfFastq) +1;
            toolbox::exportLog("ERROR: $0 : There is $nbFiles files in your directory $firstDir, this is not appropriate\n",0);
        }
    }
    else
    {
        next;
    }
}

chop $jobList if ($jobList =~/\|$/);
my $jobRunning=`qstat | egrep -c '$jobList'`; #count row with right job id
chomp $jobRunning;
toolbox::exportLog("INFOS: $0 : There are $jobRunning sub jobs running with job id $jobList\n",1);

while ( $jobRunning > 0 ) 	# when no more matching row, job is over
{
	toolbox::exportLog("INFOS: $0 : There are $jobRunning sub jobs running with job id $jobList\n",1);
        toolbox::run("sleep 3");			#go to bed for 3 seconds
	$jobRunning=`qstat | egrep -c '$jobList'`;
        chomp $jobRunning;
}

if ($loop == 1)
{
    goto END;
}

my $listOfFilesV2 = toolbox::readDir($initialDir);                                                            # read it to recover files in it
##DEBUG toolbox::exportLog("DEBUG: $0 : toolbox::readDir: @$listOfFilesV2\n",1);
my @listOfFilesV2 = @$listOfFilesV2;
my @toRedoLoop = ();

for (my $j=0; $j<=$#listOfFilesV2; $j++)
{
    if ($listOfFilesV2[$j]=~m/.+\/\w+_Single:$/)
    {
        push (@toRedoLoop,"$listOfFilesV2[$j]");
    }
}

@listOfFiles = ();
@listOfFiles = @toRedoLoop;
unless (@listOfFiles == 0)
{
    toolbox::exportLog("----------------------------------------",1);
    toolbox::exportLog("INFOS: $0 : Some single folder have been created by repairing step, analysis of them start here\n",1);
    $loop = 1;
    goto LOOP;
}

END:
toolbox::exportLog("#########################################\nCONGRATS: Mapping analysis done correctly !\n#########################################\n",1);

#########################################
# copy the final "BAM" file into the BAM directory
#########################################
my @fileAndPath = toolbox::extractPath($initialDir);
my $bamDirPath = "$fileAndPath[1]"."BamDirectory/";                                                         # name of the BAM directory
toolbox::makeDir("$bamDirPath");                                                                            # to create the BAM directory

$listOfFiles = toolbox::readDir2($initialDir);                                                              # read the initial dir
@listOfFiles = @$listOfFiles;
my @listOfGATK;
my @listOfBAM;
my $gatkPath;
my $okFinal;

for (my $i=0; $i<=$#listOfFiles; $i++)                                                                      # for each folder...
{
    ##DEBUG toolbox::exportLog("DEBUG: $0 : File directory created @ listOfFiles : $listOfFiles[$i]\n",1);
    if ($listOfFiles[$i]=~m/^(.+):$/)
    {
        $gatkPath=$1;                                                                                       # recovery of the current GATK's folder path
    }
    elsif ($listOfFiles[$i]=~m/(7_GATK)$/)
    {
        $gatkPath.="/7_GATK";
        ##DEBUG toolbox::exportLog("DEBUG: $0 : if File directory eq 7_GATK : $listOfFiles[$i]\n",1);
        my $copyCom = "cp $gatkPath/*.GATKINDELREALIGNER.ba* $bamDirPath.";                                 # command to move bam files just created into the directory appropriate for the pipeline
        ##DEBUG toolbox::exportLog("DEBUG: $0 : Check bam cp command: $copyCom\n",1);
        $okFinal = toolbox::run($copyCom);                                                                  # move the files
        undef($gatkPath);
    }
}

if ($okFinal == 1)
{
    toolbox::exportLog("#########################################\nYour final BAM files are into the $bamDirPath directory\n#########################################\n",1);
}



toolbox::exportLog("INFOS: $0 : Run mergeAnalysis.pl on $bamDirPath\n",1);
my $mergeCom = 'qsub -N mergeAnalysis '.$options.'  "mergeAnalysis.pl '.$bamDirPath.' '.$fileConf.' '.$refFastaFile.'"';
##DEBUG
toolbox::exportLog("DEBUG: $0 : qsub pairAnalysis command : $mergeCom\n",1);
toolbox::run($mergeCom);

toolbox::exportLog("#########################################\nCONGRATS: SNP calling done correctly !\n#########################################\n",1);


close F1;
exit;
