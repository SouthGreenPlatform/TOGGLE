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

##import of libraries
use strict;
use warnings;
use lib qw(../Modules/);
use toolbox;
use fastqc;
use cutadapt;
use pairing;
use bwa;
use samTools;
use gatk;
use Data::Dumper;
use picardTools;
##initialization of parameters

###fastq files in arg

toolbox::readFileConf("software.config.txt");
my $optionsIndex=$configInfos->{'BWA index'};
my $optionsAln=$configInfos->{'BWA aln'};
my $optionsCutadapt=$configInfos->{'cutadapt'};
my $optionsGatkUnifiedgenotyper=$configInfos->{'GATK UnifiedGenotyper'};
my $optionsGatkSelectVariant=$configInfos->{'GATK selectVariants'};
my $optionsGatkReadBackedPhasing=$configInfos->{'GATK ReadBackedPhasing'};
my $optionsSamToolsViewSingle=$configInfos->{'samtools view single'};
my $optionsSamToolsViewPair=$configInfos->{'samtools view pair'};
my $optionsGatkRealignerTargetCreator=$configInfos->{'GATK gatkRealignerTargetCreator'};
my $optionsgatkIndelRealigner=$configInfos->{'GATK gatkIndelRealigner'};
my $optionsPicarToolSortSamPair=$configInfos->{'picardTools sortsam pair'};
my $optionsPicarToolSortSamSingle=$configInfos->{'picardTools sortsam single'};
###################################
#print toolbox::checkFormatFastq($outputDir."reverse.fastq");
###Etape 1
###Fastqc
my $outputDir=""; ##Variable qui contient le chemin du dossier principale
if (@ARGV==0) { ## on vérifie si le dossier est bien passé en paramêtre
    print "I don't have any file to work on!!! Please give me some :-) \n";
    toolbox::exportLog("I don't have any file to work on!!! Please give me some  :-) \n",0);
}else{
    $outputDir=shift @ARGV;
    if (toolbox::existsDir($outputDir)==1) {
        print "Demarrage du pipeLine";
        
    }else{
        die "je ne trouve pas le dossier";
    }
}

###BWA INDEX
my $refFastaFileIn=shift toolbox::readDir($outputDir,".fasta");
bwa::bwaIndex($refFastaFileIn,$optionsIndex);
samTools::samToolsFaidx($refFastaFileIn);
picardTools::picardToolsCreateSequenceDictionary($refFastaFileIn);

##Pairing and repairer
my $dictionnairePair=pairing::pairRecognition($outputDir);
print Dumper($dictionnairePair);
pairing::createDirPerCouple($dictionnairePair);

#
my @pairFolders=split "\n", `ls -d $outputDir/*/ `;
print Dumper(@pairFolders);
#
foreach my $pairFolder (@pairFolders){
   
    ##FastqC
    my $fastqFiles=toolbox::readDir($pairFolder,".fastq");
    foreach my $fastqfile (@{$fastqFiles}){
        if(toolbox::sizeFile($fastqfile)==1){
            if (toolbox::checkFormatFastq($fastqfile)==1) {
                fastqc::exec($fastqfile,$pairFolder);
            }else{
                toolbox::exportLog("$fastqfile  is not a fastq file",0)
            }
        }else{
            toolbox::exportLog("$fastqfile doesn't exist or is empty ",0)
        }
    }   



    ###Cut Adapt
    my $fileAdaptator= shift toolbox::readDir($outputDir,".conf");
    print Dumper($fileAdaptator);
    my $confFile=toolbox::extractName($fileAdaptator)."Conf.conf";
    if (cutadapt::createConfFile($fileAdaptator,$confFile , $optionsCutadapt)==1){
    
        my $fastqFiles=toolbox::readDir($pairFolder,".fastq");
        foreach my $fastqfile (@{$fastqFiles}){
            if(toolbox::sizeFile($fastqfile)==1){
                if (toolbox::checkFormatFastq($fastqfile)==1) {
                    my $fileOut=toolbox::extractName($fastqfile)."_cutadapt.fastq";
                    cutadapt::exec( $fastqfile,$confFile,$fileOut);
                }else{
                    toolbox::exportLog("$fastqfile  is not a fastq file",0)
                }
            }else{
                toolbox::exportLog("$fastqfile doesn't exist or is empty ",0)
            }
        }
    }
    
    ###Repairing
    $fastqFiles=toolbox::readDir($pairFolder,"_cutadapt.fastq");
    my $forward=shift @{$fastqFiles};
    my $reverse=shift @{$fastqFiles};
    pairing::repairing($forward,$reverse);
    
    
    ##Pair::Forward and reverse
    ###########################
    ##############################################################
    ###Align
    ###Aln des deux fichiers pairer
    $fastqFiles=toolbox::readDir($pairFolder,"_repaired.fastq");
    foreach my $fastqfile (@{$fastqFiles}){
        if(toolbox::sizeFile($fastqfile)==1){
                if (toolbox::checkFormatFastq($fastqfile)==1) {
                    bwa::bwaAln($refFastaFileIn,$fastqfile,$optionsAln);
                }else{
                    toolbox::exportLog("$fastqfile  is not a fastq file",0)
                }
            }else{
                toolbox::exportLog("$fastqfile doesn't exist or is empty ",0)
            }
    }
    ###SAMPE
    $forward=shift @{$fastqFiles};
    $reverse=shift @{$fastqFiles};
    my $forwardSaiFileIn=toolbox::extractName($forward).".sai";
    my $reverseSaiFileIn=toolbox::extractName($reverse).".sai";
    my $samFileOut=toolbox::extractName($forward).".sam";
    my $readGroup=$pairFolder;
    $readGroup=~ s/\///g;
    $readGroup=~ s/\.\.//g;
    bwa::bwaSampe($samFileOut,$refFastaFileIn,$forwardSaiFileIn,$reverseSaiFileIn,$forward,$reverse,$pairFolder);
    
    ###SortSam
    my $samFileIn=$samFileOut;
    my $bamFileOut=toolbox::extractName($samFileIn)."_sorted.bam";
    picardTools::picardToolsSortSam($samFileIn,$bamFileOut,$optionsPicarToolSortSamPair);
     
    ###Samtools View
    my $bamFileIn=$bamFileOut;
    samTools::samToolsView($bamFileIn,$optionsSamToolsViewPair);
       
    
    ##Single
    ###########################
    ##############################################################
    ###Aln du fichier single
    ### A supprimer si le fichier single s appelle nom_single_repaired
    $fastqFiles=toolbox::readDir($pairFolder,"_single.fastq");
    foreach my $fastqfile (@{$fastqFiles}){
        if(toolbox::sizeFile($fastqfile)==1){
                if (toolbox::checkFormatFastq($fastqfile)==1) {
                    bwa::bwaAln($refFastaFileIn,$fastqfile,$optionsAln);
                }else{
                    toolbox::exportLog("$fastqfile  is not a fastq file",0)
                }
            }else{
                toolbox::exportLog("$fastqfile doesn't exist or is empty ",0)
            }
    }
    
     ###SAMSE
    my $single=shift @{$fastqFiles};
    my $singleSaiFileIn=toolbox::extractName($single).".sai";
    $samFileOut=toolbox::extractName($single).".sam";
    $readGroup=$pairFolder;
    $readGroup=~ s/\///g;
    $readGroup=~ s/\.\.//g;
    bwa::bwaSamse($samFileOut,$refFastaFileIn,$singleSaiFileIn,$single,$readGroup);
    
    ######Sort sam
    $samFileIn=$samFileOut;
    $bamFileOut=toolbox::extractName($samFileIn)."_sorted.bam";
    picardTools::picardToolsSortSam($samFileIn,$bamFileOut,$optionsPicarToolSortSamSingle);
 
    ###Samtools View
    $bamFileIn=$bamFileOut;
    samTools::samToolsView($bamFileIn,$optionsSamToolsViewSingle);
       
    

    ###Samtools Index
    my $bamFiles=toolbox::readDir($pairFolder,"_view.bam");
    foreach my $bamFile (@{$bamFiles}){
        if(toolbox::sizeFile($bamFile)==1){
            samTools::samToolsIndex($bamFile);
        }else{
            toolbox::exportLog("$bamFile doesn't exist or is empty ",0);
        }
    }
    
    
    foreach my $bamFile (@{$bamFiles}){
        my $intervalsFile=toolbox::extractName($bamFile).".list";
        if(toolbox::sizeFile($bamFile)==1){
            gatk::gatkRealignerTargetCreator($refFastaFileIn, $bamFile, $intervalsFile, $optionsGatkRealignerTargetCreator);
        }else{
            toolbox::exportLog("$bamFile doesn't exist or is empty ",0);
        }
    }
    
    
    
    
    foreach my $bamFile (@{$bamFiles}){
        my $intervalsFile=toolbox::extractName($bamFile).".list";
        my $bamRealigned=toolbox::extractName($bamFile)."_realigned.bam";
        if(toolbox::sizeFile($bamFile)==1){
            gatk::gatkIndelRealigner($refFastaFileIn, $bamFile, $intervalsFile, $bamRealigned, $optionsgatkIndelRealigner);
        }else{
            toolbox::exportLog("$bamFile doesn't exist or is empty ",0);
        }
    }
    
    ###Samtools Index
    $bamFiles=toolbox::readDir($pairFolder,"_realigned.bam");
    foreach my $bamFile (@{$bamFiles}){
        if(toolbox::sizeFile($bamFile)==1){
            samTools::samToolsIndex($bamFile);
        }else{
            toolbox::exportLog("$bamFile doesn't exist or is empty ",0);
        }
    }
}


