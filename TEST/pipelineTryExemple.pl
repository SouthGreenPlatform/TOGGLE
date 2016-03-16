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
use Getopt::Long;

my %Options;
my $ok=GetOptions(\%Options,"c=s","f=s");

##initialization of parameters
my $configFile=$Options{"c"};
toolbox::readFileConf($configFile);
my $optionsIndex=$configInfos->{'BWA index'};
my $optionsAln=$configInfos->{'BWA aln'};
my $optionsCutadapt=$configInfos->{'cutadapt'};
my $optionsGatkUnifiedgenotyper=$configInfos->{'GATK UnifiedGenotyper'};
my $optionsGatkSelectVariant=$configInfos->{'GATK selectVariants'};
my $optionsGatkReadBackedPhasing=$configInfos->{'GATK ReadBackedPhasing'};
my $optionsSamToolsViewSingle=$configInfos->{'samtools view single'};
my $optionsSamToolsViewPair=$configInfos->{'samtools view pair'};
my $optionsGatkRealignerTargetCreator=$configInfos->{'GATK gatkRealignerTargetCreator'};
my $optionsGatkHaplotypeCaller=$configInfos->{'GATK HaplotypeCaller'};
my $optionsGatkVariantFiltration=$configInfos->{'GATK VariantFiltration'};
my $optionsGatkSelectVariants=$configInfos->{'GATK SelectVariants'};
my $optionsgatkIndelRealigner=$configInfos->{'GATK gatkIndelRealigner'};
my $optionsgatkReadBackedPhasing=$configInfos->{'GATK ReadBackedPhasing'};
my $optionsPicarToolSortSamPair=$configInfos->{'picardTools sortsampair'};
my $optionsPicarToolSortSamSingle=$configInfos->{'picardTools sortsamsingle'};
my $optionsPicarToolMarkDuplicate=$configInfos->{'picardTools markduplicate'};


###################################
#print toolbox::checkFormatFastq($outputDir."reverse.fastq");
###Etape 1
###Fastqc
my $outputDir=""; ##Variable qui contient le chemin du dossier principale
$outputDir=$Options{"f"};
if ($outputDir eq "") { ## on vérifie si le dossier est bien passé en paramêtre
    print "I don't have any file to work on!!! Please give me some :-) \n";
    toolbox::exportLog("I don't have any file to work on!!! Please give me some  :-) \n",0);
}else{
    if (toolbox::existsDir($outputDir)==1) {
        print "Demarrage du pipeLine";
        
    }else{
        die "je ne trouve pas le dossier";
    }
}



   
#########################Copie des fichiers fastq dans un nouveau dossier
##création du nouveau dossier
my $fastqFolder="$outputDir/fastqFolder";
toolbox::makeDir($fastqFolder);
my $fastqFiles=toolbox::readDir("$outputDir","fastq");
foreach my $readfile (@{$fastqFiles}){
    toolbox::run("cp $readfile $fastqFolder");
    toolbox::run("rm $readfile"); 
}

####BWA INDEX on reference
my $refFastaFileIn=shift toolbox::readDir($outputDir,".fasta");
bwa::bwaIndex($refFastaFileIn,$optionsIndex);

###samTools Faidx (index) on reference
samTools::samToolsFaidx($refFastaFileIn);


###picardTools Create Sequence Dictionnary from reference
my $refdict=$outputDir."/".toolbox::extractName($refFastaFileIn).".dict";
picardTools::picardToolsCreateSequenceDictionary($refFastaFileIn,$refdict);



############################################################################
############################################################################
##Pairing and repairer
my $dictionnairePair=pairing::pairRecognition($fastqFolder);
pairing::createDirPerCouple($dictionnairePair,$outputDir);


while ((my $cle,my $valeur)=each %{$dictionnairePair}) {
    my %valeur=%{$valeur};
    my $pairFolder="$outputDir/$valeur{'ReadGroup'}";

#

my $fastqCDir="$pairFolder/1_FastqC";
my $CutAdaptDir="$pairFolder/2_CutAdapt";
my $repairingDir="$pairFolder/3_repairing";
my $bwaAlnDir="$pairFolder/4_bwaAln";
my $bwaSampeDir="$pairFolder/5_bwaSampe";

#    ##FastqC
   
    toolbox::makeDir($fastqCDir);
    my $fastqFiles=toolbox::readDir($pairFolder,".fastq");
    foreach my $fastqfile (@{$fastqFiles}){
        if(toolbox::sizeFile($fastqfile)==1){
            if (toolbox::checkFormatFastq($fastqfile)==1) {
                fastqc::execution($fastqfile,$fastqCDir);
            }else{
                toolbox::exportLog("$fastqfile  is not a fastq file",0)
            }
        }else{
            toolbox::exportLog("$fastqfile doesn't exist or is empty ",0)
        }
    }   


######################################################################
    ###Cut Adapt
   
    toolbox::makeDir($CutAdaptDir);
    my $fileAdaptator= shift toolbox::readDir($outputDir,".conf");
    print Dumper($fileAdaptator);
    my $confFile=toolbox::extractName($fileAdaptator)."Conf.conf";
    if (cutadapt::createConfFile($fileAdaptator,$confFile , $optionsCutadapt)==1){
    
        my $fastqFiles=toolbox::readDir($pairFolder,".fastq");
        foreach my $fastqfile (@{$fastqFiles}){
            if(toolbox::sizeFile($fastqfile)==1){
                if (toolbox::checkFormatFastq($fastqfile)==1) {
                    my $fileOut="$CutAdaptDir/".toolbox::extractName($fastqfile)."_cutadapt.fastq";
                    cutadapt::execution( $fastqfile,$confFile,$fileOut);
                }else{
                    toolbox::exportLog("$fastqfile  is not a fastq file",0)
                }
            }else{
                toolbox::exportLog("$fastqfile doesn't exist or is empty ",0)
            }
        }
    }
    
    ###Repairing
    $fastqFiles=toolbox::readDir($CutAdaptDir,"_cutadapt.fastq");
    my $forward=shift @{$fastqFiles};
    my $reverse=shift @{$fastqFiles};
    toolbox::makeDir($repairingDir);
    pairing::repairing($forward,$reverse,$repairingDir);
    
    
    ##
    ###########################
    ##############################################################
    ###Align
    ###Aln des deux fichiers pairer
    $fastqFiles=toolbox::readDir($repairingDir,".REPAIRING.fastq");
  
    toolbox::makeDir($bwaAlnDir);
    foreach my $fastqfile (@{$fastqFiles}){
        if(toolbox::sizeFile($fastqfile)==1){
                if (toolbox::checkFormatFastq($fastqfile)==1) {
                    my $saifileOut="$bwaAlnDir/".toolbox::extractName($fastqfile)."_aln.sai";
                    bwa::bwaAln($refFastaFileIn,$fastqfile,$saifileOut,$optionsAln);
                }else{
                    toolbox::exportLog("$fastqfile  is not a fastq file",0)
                }
            }else{
                toolbox::exportLog("$fastqfile doesn't exist or is empty ",0)
            }
    }
    ###SAMPE
    #
    toolbox::makeDir($bwaSampeDir);
    $forward=shift @{$fastqFiles};
    $reverse=shift @{$fastqFiles};
    my $forwardSaiFileIn=toolbox::extractName($forward).".sai";
    my $reverseSaiFileIn=toolbox::extractName($reverse).".sai";
    my $samFileOut="$bwaSampeDir/".toolbox::extractName($forward).".sam";
    my $readGroup=$pairFolder;
    bwa::bwaSampe($samFileOut,$refFastaFileIn,$forwardSaiFileIn,$reverseSaiFileIn,$forward,$reverse,$readGroup);
    
    ####SortSam
    #my $samFileIn=$samFileOut;
    #my $bamFileOut=toolbox::extractName($samFileIn)."_sorted.bam";
    #picardTools::picardToolsSortSam($samFileIn,$bamFileOut,$optionsPicarToolSortSamPair);
    # 
    ####Samtools View
    #my $bamFileIn=$bamFileOut;
    #samTools::samToolsView($bamFileIn,$optionsSamToolsViewPair);
       
#    
  
#    
 
#    
#    ##Single
#    ###########################
#    ##############################################################
#    ###Aln du fichier single
#    ### A supprimer si le fichier single s appelle nom_single_repaired
#    $fastqFiles=toolbox::readDir($pairFolder,"_single.fastq");
#    foreach my $fastqfile (@{$fastqFiles}){
#        if(toolbox::sizeFile($fastqfile)==1){
#                if (toolbox::checkFormatFastq($fastqfile)==1) {
#                    bwa::bwaAln($refFastaFileIn,$fastqfile,$optionsAln);
#                }else{
#                    toolbox::exportLog("$fastqfile  is not a fastq file",0)
#                }
#            }else{
#                toolbox::exportLog("$fastqfile doesn't exist or is empty ",0)
#            }
#    }
#    
#     ###SAMSE
#    my $single=shift @{$fastqFiles};
#    my $singleSaiFileIn=toolbox::extractName($single).".sai";
#    $samFileOut=toolbox::extractName($single).".sam";
#    $readGroup=$pairFolder;
#    $readGroup=~ s/\///g;
#    $readGroup=~ s/\.\.//g;
#    bwa::bwaSamse($samFileOut,$refFastaFileIn,$singleSaiFileIn,$single,$readGroup);
#    
#    ######Sort sam
#    $samFileIn=$samFileOut;
#    $bamFileOut=toolbox::extractName($samFileIn)."_sorted.bam";
#    picardTools::picardToolsSortSam($samFileIn,$bamFileOut,$optionsPicarToolSortSamSingle);
# 
#    ###Samtools View
#    $bamFileIn=$bamFileOut;
#    samTools::samToolsView($bamFileIn,$optionsSamToolsViewSingle);
#    ##################################################################################################################
#    
#    
#    #####Suite du mapping
#    ###Samtools Index
#    my $bamFiles=toolbox::readDir($pairFolder,"_view.bam");
#    foreach my $bamFile (@{$bamFiles}){
#        if(toolbox::sizeFile($bamFile)==1){
#            samTools::samToolsIndex($bamFile);
#        }else{
#            toolbox::exportLog("$bamFile doesn't exist or is empty ",0);
#        }
#    }
#    
#    ########################################################################
#    ##Target CReator
#    foreach my $bamFile (@{$bamFiles}){
#        my $intervalsFile=toolbox::extractName($bamFile).".list";
#        if(toolbox::sizeFile($bamFile)==1){
#            gatk::gatkRealignerTargetCreator($refFastaFileIn, $bamFile, $intervalsFile, $optionsGatkRealignerTargetCreator);
#        }else{
#            toolbox::exportLog("$bamFile doesn't exist or is empty ",0);
#        }
#    }
#    
#    
#    ########################################################################
#    ##Target Indel REaligner
#    foreach my $bamFile (@{$bamFiles}){
#        my $intervalsFile=toolbox::extractName($bamFile).".list";
#        my $bamRealigned=toolbox::extractName($bamFile)."_realigned.bam";
#        if(toolbox::sizeFile($bamFile)==1){
#            gatk::gatkIndelRealigner($refFastaFileIn, $bamFile, $intervalsFile, $bamRealigned, $optionsgatkIndelRealigner);
#        }else{
#            toolbox::exportLog("$bamFile doesn't exist or is empty ",0);
#        }
#    }
#    
#    ###Samtools Index
#    $bamFiles=toolbox::readDir($pairFolder,"_realigned.bam");
#    foreach my $bamFile (@{$bamFiles}){
#        if(toolbox::sizeFile($bamFile)==1){
#            samTools::samToolsIndex($bamFile);
#        }else{
#            toolbox::exportLog("$bamFile doesn't exist or is empty ",0);
#        }
#    }
#    
#    ####Remove duplicates
#    foreach my $bamFile (@{$bamFiles}){
#        my $bamOut=toolbox::extractName($bamFile)."_duplicateRemoved.bam";
#        my $metrics=toolbox::extractName($bamFile).".metrics";
#        if(toolbox::sizeFile($bamFile)==1){
#            picardTools::picardToolsMarkDuplicates($bamFile, $bamOut, $metrics, $optionsPicarToolMarkDuplicate);
#        }else{
#            toolbox::exportLog("$bamFile doesn't exist or is empty ",0);
#        }
#    }
#    
#    ###Samtools Index
#    $bamFiles=toolbox::readDir($pairFolder,"_duplicateRemoved.bam");
#    foreach my $bamFile (@{$bamFiles}){
#        if(toolbox::sizeFile($bamFile)==1){
#            samTools::samToolsIndex($bamFile);
#        }else{
#            toolbox::exportLog("$bamFile doesn't exist or is empty ",0);
#        }
#    }  
}
#
####################################################################################################
#    ###Merge bam
#    my $bamOutFile=$outputDir."/merged.bam";
#    my @bamFiles;
#    foreach my $pairFolder (@pairFolders){
#        foreach my $bamfile (toolbox::readDir($pairFolder,"_duplicateRemoved.bam")){
#            push @bamFiles, $bamfile;
#        }
#    }
#   
#   
#   my $bamFiles= shift \@bamFiles;
####   
#    samTools::samToolsMerge($bamFiles,$bamOutFile);
#    
#
#########
#
#
#    ##Samtools Index
#    $bamFiles=toolbox::readDir($outputDir,"merged.bam");
#    foreach my $bamFile (@{$bamFiles}){
#        if(toolbox::sizeFile($bamFile)==1){
#            samTools::samToolsIndex($bamFile);
#        }else{
#            toolbox::exportLog("$bamFile doesn't exist or is empty ",0);
#        }
#    }
#    
#    ###################################################################################################
#    ####################################################################################################
#    ################################### Calling ########################################################
#    
#    ######Happlotype Caller
#    my $vcfCalled=$outputDir."/haplotype.vcf";
#    my $bamsToCall=toolbox::readDir($outputDir,"merged.bam");
#    
#    gatk::gatkHaplotypeCaller($refFastaFileIn, $vcfCalled, $bamsToCall,$optionsGatkHaplotypeCaller);
#    
#    
#    
#    ######Variant filtrator
#    my $vcfToFilter=$vcfCalled;
#    my $vcfFiltered=toolbox::extractName($vcfToFilter)."_filtred.vcf";
#    gatk::gatkVariantFiltration($refFastaFileIn, $vcfFiltered, $vcfToFilter, $optionsGatkVariantFiltration);
#    
#    
#    
#    ###############Select Variant
#    my $vcfSnpKnownFile=$vcfFiltered;
#    my $vcfVariantsSelected=$vcfFiltered=toolbox::extractName($vcfSnpKnownFile)."_SNP.vcf";
#    gatk::gatkSelectVariants($refFastaFileIn, $vcfSnpKnownFile, $vcfVariantsSelected, $optionsGatkSelectVariants);
#    
#    my $vcfVariant=$vcfVariantsSelected;
#    my $vcfFileOut=toolbox::extractName($vcfVariant)."_phased.vcf";
#    my $bamFileIn=shift toolbox::readDir($outputDir,"merged.bam");
#    #gatk::gatkReadBackedPhasing($refFastaFileIn, $bamFileIn,$vcfVariant, $vcfFileOut, $optionsgatkReadBackedPhasing);


