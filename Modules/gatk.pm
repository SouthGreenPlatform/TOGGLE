package gatk;



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
use localConfig;
use toolbox;
use Data::Dumper;

# GATK Base Recalibrator: recalibrate the quality score of bases from informations stemming from SNP VCF file
sub gatkBaseRecalibrator
{
    my ($refFastaFileIn, $bamToRecalibrate, $vcfSnpKnownFile, $tableReport, $optionsHachees) = @_;      # recovery of information
    if ((toolbox::checkSamOrBamFormat($bamToRecalibrate)==2) and (toolbox::sizeFile($refFastaFileIn)==1) and (toolbox::sizeFile($vcfSnpKnownFile)==1) and (toolbox::sizeFile($bamToRecalibrate)==1))     # check if files exists and arn't empty and stop else
    {
        my $options=toolbox::extractOptions($optionsHachees);       # extraction of options parameters
        my $comGatkBaseRecalibrator = "$GATK"."$options"." -I $bamToRecalibrate -R $refFastaFileIn -knownSites $vcfSnpKnownFile -o $tableReport";       # command line
        toolbox::run($comGatkBaseRecalibrator);
        if(toolbox::run($comGatkBaseRecalibrator)==1)
        {
            toolbox::exportLog("INFOS: gatk::gatkBaseRecalibrator : Correctly done\n",1);
            return 1;
        }
        else        # if one or some previous files doesn't exist or is/are empty or if gatkBaseRecalibrator failed
        {
            toolbox::exportLog("ERROR: gatk::gatkBaseRecalibrator : Uncorrectly done\n", 0);        # returns the error message
            return 0;
        }
    }
    else        # if one or some previous files doesn't exist or is/are empty or if gatkBaseRecalibrator failed
    {
        toolbox::exportLog("ERROR: gatk::gatkBaseRecalibrator : The files $refFastaFileIn, $vcfSnpKnownFile or/and $bamToRecalibrate are incorrects\n", 0);     # returns the error message
        return 0;
    }
}
# GATK Realigner Target Creator: determine (small) suspicious intervals which are likely in need of realignment
sub gatkRealignerTargetCreator
{
    my ($refFastaFileIn, $bamToRealigne, $intervalsFile, $optionsHachees) = @_;     # recovery of information
    if ((toolbox::checkSamOrBamFormat($bamToRealigne)==2) and (toolbox::sizeFile($refFastaFileIn)==1) and (toolbox::sizeFile($bamToRealigne)==1))     # check if files exists and arn't empty and stop else
    {
        my $options=toolbox::extractOptions($optionsHachees);       # extraction of options parameters
        my $comGatkRealignerTargetCreator = "$GATK"."$options"." -R $refFastaFileIn -I $bamToRealigne -o $intervalsFile ";#--fix_misencoded_quality_scores -fixMisencodedQuals";        # command line
        if(toolbox::run($comGatkRealignerTargetCreator)==1)     # command line execution
        {
            toolbox::exportLog("INFOS: gatk::gatkRealignerTargetCreator : Correctly done\n",1);
            return 1;
        }
    }
    else        # if one or some previous files doesn't exist or is/are empty or if gatkRealignerTargetCreator failed
    {
        toolbox::exportLog("ERROR: gatk::gatkRealignerTargetCreator : The files $refFastaFileIn or/and $bamToRealigne are incorrects\n", 0);        # returns the error message
        return 0;
    }
}
# GATK Indel Realigner: run the realigner over the intervals producted by gatk::gatkRealignerTargetCreator (see above)
sub gatkIndelRealigner
{
    my ($refFastaFileIn, $bamToRealigne, $intervalsFile, $bamRealigned, $optionsHachees) = @_;      # recovery of information
    if ((toolbox::checkSamOrBamFormat($bamToRealigne)==2) and (toolbox::sizeFile($refFastaFileIn)==1) and (toolbox::sizeFile($bamToRealigne)==1) and (toolbox::readFile($intervalsFile)==1))      # check if files exists and arn't empty and stop else
    {
        my $options=toolbox::extractOptions($optionsHachees);       # extraction of options parameters
        my $comGatkIndelRealigner = "$GATK"."$options"." -R $refFastaFileIn -I $bamToRealigne -targetIntervals $intervalsFile -o $bamRealigned";# --fix_misencoded_quality_scores -fixMisencodedQuals";     # command line
        if(toolbox::run($comGatkIndelRealigner)==1)
        {                                                                                                                                                                               # command line execution
            toolbox::exportLog("INFOS: gatk::gatkIndelRealigner : Correctly done\n",1);
            return 1;
        }
    }
    else        # if one or some previous files doesn't exist or is/are empty or if gatkIndelRealigner failed
    {
        toolbox::exportLog("ERROR: gatk:gatkIndelRealigner : The files $refFastaFileIn, $bamToRealigne or/and $intervalsFile are incorrects\n", 0);     # returns the error message
        return 0;
    }
}
# GATK Haplotype Caller: Haplotypes are evaluated using an affine gap penalty Pair HMM.
sub gatkHaplotypeCaller
{
    my ($refFastaFileIn, $vcfCalled, $listOfBam, $optionsHachees, $vcfSnpKnownFile, $intervalsFile) = @_;      # recovery of information
 
    #TODO adding VCF control
    if (toolbox::sizeFile($refFastaFileIn)==1)     # check if files exist and isn't empty and stop else
    {
        ## informations about SNP know file
        my $dbsnp="";
        if (($vcfSnpKnownFile) and (toolbox::sizeFile($vcfSnpKnownFile)==0))        # If specified, check if the file of known SNP is correct
        {
            toolbox::exportLog("ERROR: gatk::gatkHaplotypeCaller : The file $vcfSnpKnownFile is uncorrect\n", 0);       # returns the error message
            return 0;
        }
        elsif (($vcfSnpKnownFile) and (toolbox::sizeFile($vcfSnpKnownFile)==1))     # if specified and file correct ...
        {
            $dbsnp=" --dbsnp $vcfSnpKnownFile";     # ... recovery of informations for command line used later
        }
        
        ## informations about intervals file
        my $intervals="";
        if (($intervalsFile) and (toolbox::sizeFile($intervalsFile)==0))        # if specified, check if the file of intervals is correct
        {
            toolbox::exportLog("ERROR: gatk::gatkHaplotypeCaller : The file $intervalsFile is uncorrect\n", 0);     # returns the error message
            return 0;
        }
        elsif (($intervalsFile) and (toolbox::sizeFile($intervalsFile)==1))     # if specified and file of intervals is correct ...
        {
            $intervalsFile="-L $intervalsFile";     # ... recovery of informations for command line used later
        }
        
        my $bamFiles_names="";
        foreach my $file (@{$listOfBam})       # for each BAM file(s)
        {
            if (toolbox::checkSamOrBamFormat($file)==2 and toolbox::sizeFile($file)==1)        # if current file is not empty
            {
                $bamFiles_names.="-I ".$file." ";       # recovery of informations fo command line used later
            }
            else        # if current file is empty
            {
                toolbox::exportLog("ERROR: gatk::gatkHaplotypeCaller : The file $file is uncorrect\n", 0);      # returns the error message
                return 0;
            }
        }
        
        my $options="";
        if ($optionsHachees)
        {
            $options=toolbox::extractOptions($optionsHachees);      ##Get given options
        }
        my $comGatkHaplotypeCaller = "$GATK"."$options"." -R $refFastaFileIn $bamFiles_names $dbsnp $intervals -o $vcfCalled";      # command line
        if(toolbox::run($comGatkHaplotypeCaller)==1)        # command line execution
        {
            toolbox::exportLog("INFOS: gatk::gatkHaplotypeCaller : Correctly done\n",1);
            return 1;
        }
        else
        {
            toolbox::exportLog("ERROR: gatk::gatkHaplotypeCaller : Uncorrectly done\n", 0);     # returns the error message
            return 0;
        }
    }
    else        # if one or some previous files doesn't exist or is/are empty or if gatkHaplotypeCaller failed
    {
        toolbox::exportLog("ERROR: gatk::gatkHaplotypeCaller : The file $refFastaFileIn is uncorrect\n", 0);        # returns the error message
        return 0;
    }
}
# GATK Select Variants: Selects variants from a VCF source.
sub gatkSelectVariants
{
    my ($refFastaFileIn, $vcfSnpKnownFile, $vcfVariantsSelected, $optionsHachees) = @_;     # recovery of information    
    #TODO adding the VCF type control
    if ((toolbox::sizeFile($refFastaFileIn)==1)  and  (toolbox::sizeFile($vcfSnpKnownFile)==1))     # check if ref file exist and isn't empty and stop else
    {
        my $options=toolbox::extractOptions($optionsHachees);       # extraction of options parameters
        my $comGatkSelectVariants = "$GATK"."$options"." -R $refFastaFileIn --variant $vcfSnpKnownFile -o $vcfVariantsSelected";        # command line
        if(toolbox::run($comGatkSelectVariants)==1)     # command line execution
        {
            toolbox::exportLog("INFOS: gatk::gatkSelectVariants : Correctly done\n",1);
            return 1;
        }
        else        # if one or some previous files doesn't exist or is/are empty or if gatkSelectVariants failed
        {
            toolbox::exportLog("ERROR: gatk::gatkSelectVariants : Uncorrectly done\n", 0);      # returns the error message
            return 0;
        }
    }
    else        # if one or some previous files doesn't exist or is/are empty or if gatkSelectVariants failed
    {
        toolbox::exportLog("ERROR: gatk::gatkSelectVariants : The files $refFastaFileIn or/and $vcfSnpKnownFile are incorrects\n", 0);      # returns the error message
        return 0;
    }
}
# GATK Variant Filtration: filter variant calls using a number of user-selectable, parameterizable criteria.
sub gatkVariantFiltration
{
    my ($refFastaFileIn, $vcfFiltered, $vcfToFilter, $optionsHachees) = @_;     # recovery of information
    #TODO adding the VCF type control
    if ((toolbox::sizeFile($refFastaFileIn)==1) and (toolbox::sizeFile($vcfToFilter)==1))       # check if ref file exist and isn't empty and stop else
    {
        my $options="";
        if ($optionsHachees)
        {
            $options=toolbox::extractOptions($optionsHachees);      ##Get given options
        }
        print Dumper($options);     # extraction of options parameters
        my $comGatkVariantFiltration = "$GATK"."$options"." -R $refFastaFileIn -o $vcfFiltered --variant $vcfToFilter";     # command line
        if(toolbox::run($comGatkVariantFiltration)==1)      # command line execution
        {
            toolbox::exportLog("INFOS: gatk::gatkVariantFiltration : Correctly done\n",1);
            return 1;
        }
        else        # if one or some previous files doesn't exist or is/are empty or if gatkVariantFiltration failed
        {
            toolbox::exportLog("ERROR: gatk::gatkVariantFiltration : Uncorrectly done\n", 0);       # returns the error message
            return 0;
        }
    }
    else        # if one or some previous files doesn't exist or is/are empty or if gatkVariantFiltration failed
    {
        toolbox::exportLog("ERROR: gatk::gatkVariantFiltration :  The files $refFastaFileIn or/and $vcfToFilter are incorrects\n", 0);      # returns the error message
        return 0;
    }
}
# GATK Unified Genotyper: A variant caller which unifies the approaches of several disparate callers -- Works for single-sample and multi-sample data.
sub gatkUnifiedGenotyper
{
    my ($refFastaFileIn, $bamFileIn, $vcfFileOut, $optionsHachees) = @_;        # recovery of information  
    #TODO adding VCF type control
    if ((toolbox::checkSamOrBamFormat($bamFileIn)==2) and (toolbox::sizeFile($refFastaFileIn)==1) and (toolbox::sizeFile($bamFileIn)==1))     # check if ref file exist and isn't empty and stop else
    {
        my $options=toolbox::extractOptions($optionsHachees);
        my $comGatkUnifiedGenotyper = "$GATK"."$options"." -R $refFastaFileIn -I $bamFileIn -o $vcfFileOut";        # command line
        if(toolbox::run($comGatkUnifiedGenotyper)==1)       # command line execution
        {
            toolbox::exportLog("INFOS: gatk::gatkUnifiedGenotyper : Correctly done\n",1);
            return 1;
        }
        else        # if one or some previous files doesn't exist or is/are empty or if gatkVariantFiltration failed
        {
            toolbox::exportLog("ERROR: gatk::gatkUnifiedGenotyper : Uncorrectly done\n", 0);        # returns the error message
            return 0;
        }
    }
    else        # if one or some previous files doesn't exist or is/are empty or if gatkVariantFiltration failed
    {
        toolbox::exportLog("ERROR: gatk::gatkUnifiedGenotyper : The files $refFastaFileIn or/and $bamFileIn are incorrects\n", 0);      # returns the error message
        return 0;
    }
}
# GATK Read backedPhasing: Walks along all variant ROD loci, caching a user-defined window of VariantContext sites, and then finishes phasing them when they go out of range (using upstream and downstream reads).
sub gatkReadBackedPhasing
{
    my ($refFastaFileIn, $bamFileIn,$vcfVariant, $vcfFileOut, $optionsHachees) = @_;         # recovery of information
    #TODO adding VCF type control
    if ((toolbox::checkSamOrBamFormat($bamFileIn)==2) and (toolbox::sizeFile($refFastaFileIn)==1) and (toolbox::sizeFile($bamFileIn)==1) and (toolbox::sizeFile($vcfVariant)==1))     # check if ref file exist and isn't empty and stop else
    {
        my $options=toolbox::extractOptions($optionsHachees);
        my $comGatkReadBackedPhasing = "$GATK"."$options"." -R $refFastaFileIn -I $bamFileIn --variant $vcfVariant -o $vcfFileOut";     # command line
        if(toolbox::run($comGatkReadBackedPhasing)==1)      # command line execution
        {
            toolbox::exportLog("INFOS: gatk::gatkReadBackedPhasing : Correctly done\n",1);
            return 1;
        }
    }
    else        # if one or some previous files doesn't exist or is/are empty or if gatkVariantFiltration failed
    {
        toolbox::exportLog("ERROR: gatk::gatkReadBackedPhasing : The files $refFastaFileIn, $bamFileIn or/and $vcfVariant are incorrects\n", 0);        # returns the error message
        return 0;
    }
}
1;

=head1 NAME

    Package I<gatk> 

=head1 SYNOPSIS

        use gatk;
    
        gatk::gatkBaseRecalibrator ($refFastaFileIn, $bamToRecalibrate, $vcfSnpKnownFile, $tableReport, $option_prog{'GATK gatkBaseRecalibrator'});
    
        gatk::gatkRealignerTargetCreator ($refFastaFileIn, $bamToRealigne, $intervalsFile, $option_prog{'GATK gatkRealignerTargetCreator'});
    
        gatk::gatkIndelRealigner ($refFastaFileIn, $bamToRealigne, $intervalsFile, $bamRealigned, $option_prog{'GATK gatkIndelRealigner'});
    
        gatk::gatkHaplotypeCaller ($refFastaFileIn, $vcfCalled, $bamsToCall,$option_prog{'GATK gatkHaplotypeCaller'}, $vcfSnpKnownFile, $intervalsFile);
    
        gatk::gatkSelectVariants ($refFastaFileIn, $vcfSnpKnownFile, $vcfVariantsSelected, $option_prog{'GATK selectVariants'});
    
        gatk::gatkVariantFiltration ($refFastaFileIn, $vcfFiltered, $vcfToFilter, $option_prog{'GATK VariantFiltration'});
    
        gatk::gatkUnifiedGenotyper ($refFastaFileIn, $bamFileIn, $vcfFileOut, $option_prog{'GATK UnifiedGenotyper'});
    
        gatk::gatkReadBackedPhasing ($refFastaFileIn, $bamFileIn,$vcfVariant, $vcfFileOut, $option_prog{'GATK ReadBackedPhasing'});

=head1 DESCRIPTION

    Package GATK (McKenna et al, 2010) The Genome Analysis Toolkit or GATK is a software package developed at the Broad Institute to analyze high-throughput sequencing data. 

=head2 FUNCTIONS


=head3 gatk::gatkBaseRecalibrator

This module recalibrate the quality score of bases from informations stemming from SNP VCF file
It takes at least three arguments: the database indexed, the file ".bam" to recalibrate, the name of table of report which will be created by this module
The file of already known SNP is not mandatory
The last argument is the options of gatk baseRecalibrator, for more informations see https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_bqsr_BaseRecalibrator.php



=head3 gatk::gatkRealignerTargetCreator

This module determine (small) suspicious intervals which are likely in need of realignment
It takes at least three arguments: the database indexed, the file ".bam" to realigne, the name of the intervals file created by this module
The last argument is the options of gatk realignerTargetCreator, for more informations see https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_indels_RealignerTargetCreator.php


=head3 gatk::gatkIndelRealigner

/!\ IMPORTANT TO NOTE: This module requires the use of the previous one /!\ The input BAM(s), reference, and known indel file(s) should be the same in both modules.
This module run the realigner over the intervals producted by gatk::gatkRealignerTargetCreator (see above)
It takes four arguments: the database indexed, the file ".bam" to realigne, the intervals file generated previously, the name of the file realigned which will be created bu this module
The last arguement is the options of gatk indelRealigner, for more informations see https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_indels_IndelRealigner.php



=head3 gatk::gatkHaplotypeCaller

This module evaluate haplotypes using an affine gap penalty Pair HMM
It takes at least three arguments: the database indexed, the file(s) ".bam" to evaluate, the name of the ouput file ".vcf"
The fourth argument is the options of gatk haplotypeCaller, for more informations see https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php
The fifth argument is the file of already known SNPs at VCF format
The last argument is the intervals file created by gatkRealignerTargetCreator (see above)



=head3 gatk::gatkSelectVariants

This module selects variants from a VCF source
It takes at least three arguments: the database indexed, the file of already known SNPs, the name of the output file
The last argument is the options of gatk selectVariants, for more informations see https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_variantutils_SelectVariants.php



=head3 gatk::gatkVariantFiltration

This module filter variant calls using a number of user-selectable, parameterizable criteria
It takes at least three arguments: the database indexed, the file, the name of the output file, the file to filter in ".vcf" format
The last argument is the options of gatk variantFiltration, for more informations see https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_filters_VariantFiltration.php



=head3 gatk::gatkUnifiedGenotyper

This module is a variant caller which unifies the approaches of several disparate callers -- Works for single-sample and multi-sample data
It takes at least three arguments: the database indexed, the ".bam" file in, the name of the output file in ".vcf" format
The last argument is the options of gatk unifiedGenotyper, for more informations see https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_genotyper_UnifiedGenotyper.php



=head3 gatk::gatkReadBackedPhasing

This module walks along all variant ROD loci, caching a user-defined window of VariantContext sites, and then finishes phasing them when they go out of range (using upstream and downstream reads)
It takes at least four arguments: the database indexed, the ".bam" file in, the already known variant in ".vcf" format, the name of the output file in ".vcf" format
The last argument is the options of gatk readBackPhasing, for more informations see https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_phasing_ReadBackedPhasing.php

=head1 AUTHORS

Intellectual property belongs to IRD, CIRAD and South Green developpement plateform 
Written by Cecile Monat, Ayite Kougbeadjo, Marilyne Summo, Cedric Farcy, Mawusse Agbessi, Christine Tranchant and Francois Sabot

=head1 SEE ALSO

L<http://www.southgreen.fr/>

=cut
