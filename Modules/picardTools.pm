package picardTools;



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
use lib qw(.);
use localConfig;
use toolbox;

######################################
#PicardTools MarkDuplicates
######################################
# This module examine aligned records in the supplied BAM file to locate duplicate molecules. All records are then written to the output file with the duplicate records flagged
sub picardToolsMarkDuplicates
{
   
    my ($bamToAnalyze, $bamAnalyzed, $bamDuplicates, $optionsHachees) = @_;     # recovery of information
    if (toolbox::checkSamOrBamFormat($bamToAnalyze)==2 && toolbox::sizeFile($bamToAnalyze)==1)      # check if file are really bam file; if file exists and isn't empty
    {
        my $options="";
        if ($optionsHachees) 
	{
            $options=toolbox::extractOptions($optionsHachees,"=");      # recovery of options if they are provided
        }
        my $comPicardToolsMarkDuplicates = "$picard/picard.jar MarkDuplicates $options INPUT=$bamToAnalyze OUTPUT=$bamAnalyzed METRICS_FILE=$bamDuplicates ";      #command line 
        toolbox::run($comPicardToolsMarkDuplicates);        # command line execution                                                                                                                                                                                                                    
    }
    else        # if something wrong (size, format) in the file to examine, don't run the module ...                                                                                                                                                                                                                                                                 # if previous files doesn't exists or are empty or if picardToolsMarkDuplicates failed
    {
        toolbox::exportLog("ERROR: picardTools::picardToolsMarkDuplicates : The file $bamToAnalyze is incorrect\n", 0);     # ... and report an error                                                                                                                                                                              # returns error message
    }
}
######################################
#picardToolsCreateSequenceDictionary
######################################
# This module read fasta or fasta.gz containing reference sequences, and write as a SAM or BAM file with only sequence dictionary.
sub picardToolsCreateSequenceDictionary
{
    my($refFastaFile,$dictFileOut,$optionsHachees)= @_;     # recovery of informations
    if (toolbox::sizeFile($refFastaFile)==1)        # check if the reference fasta file is not empty 
    {
	if (toolbox::existsFile($dictFileOut,0)==0)
	{
	    my $options="";
	    if ($optionsHachees)
	    {
	        $options=toolbox::extractOptions($optionsHachees);      # recovery of options if they are provided
	    }
	    my $command="$picard/picard.jar CreateSequenceDictionary $options REFERENCE=$refFastaFile OUTPUT=$dictFileOut";        #creation of the command line
	    if(toolbox::run($command)==1)       #Execution of the command line
	    {
	        toolbox::exportLog("INFOS: picardTools::picardToolsCreateSequenceDictionary : Correctly run\n",1);
	        return 1;
	    }
	else 
	{
	    toolbox::exportLog("INFOS: picardTools::picardToolsCreateSequenceDictionary : The file $dictFileOut already exists\n",1);
	    return 1;
	}
    }
    else        # if the reference fasta file is empty, don't run the module ...
    {
        toolbox::exportLog("ERROR: picardTools::picardToolsCreateSequenceDictionary : The file $refFastaFile is incorrect\n",0);        # ... and return and error message
        return 0;
    }
}    
######################################
#PicardToolsSortSam
######################################
# This module sorts the input SAM or BAM.
sub picardToolsSortSam
{
    my($bamOrSamFileIn,$bamOrSamFileOut,$optionsHachees)= @_;       # recovery of informations
    if ((toolbox::checkSamOrBamFormat($bamOrSamFileIn)==1) && (toolbox::sizeFile($bamOrSamFileIn)==1))        # check if the file to sort is a ".sam" one and is not empty
    {
        my $options="";
        if ($optionsHachees)
        {
            $options=toolbox::extractOptions($optionsHachees,"=");      # recovery of options if they are provided
        }
        my $command="$picard/picard.jar SortSam $options INPUT=$bamOrSamFileIn OUTPUT=$bamOrSamFileOut";       #creation of the command line
        if(toolbox::run($command)==1)       #Execute command
        {
            toolbox::exportLog("INFOS: picardTools::picardToolsSortSam : Correctly run\n",1);
            return 1;
        }
    }
    else        # if the file is not a ".bam" one or is empty don't run the module ...
    {
        toolbox::exportLog("ERROR: picardTools::picardToolsSortSam : The file $bamOrSamFileIn is incorrect\n",0);       # ... and return an error message
        return 0;
    }
}
1;

=head1 NAME

    Package I<picardTools> 

=head1 SYNOPSIS

        use picardTools;
    
        picardTools::picardToolsMarkDuplicates ($bamToAnalyze, $bamAnalyzed, $bamDuplicates, $option_prog{'picardTools markDuplicates'});
    
        picardTools::picardToolsCreateSequenceDictionary ($refFastaFile,$dictFileOut,$option_prog{'picardTools createSequenceDictionary'});
    
        picardTools::picardToolsSortSam ($bamOrSamFileIn,$bamOrSamFileOut,$option_prog{'picardTools sortsam single/pair'});

=head1 DESCRIPTION

    Package picardTools is a set of Java command line tools for manipulating high-throughput sequencing data (HTS) data and formats.

=head2 FUNCTIONS


=head3 picardTools::picardToolsMarkDuplicates

This module examines aligned records in the supplied SAM or BAM file to locate duplicate molecules. All records are then written to the output file with the duplicate records flagged
It takes at least three arguments: the ".bam" file user wants to examine, the name of the output file,the name of file to write duplication metrics
The last argument is the options of picardTools markDuplicates, for more informations see http://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates



=head3 picardTools::picardToolsCreateSequenceDictionary

This module read fasta or fasta.gz containing reference sequences, and write as a SAM or BAM file with only sequence dictionary
It takes at least two arguments: the database indexed, the name of the output file
The last argument is the options of picardTools createSequenceDictionnary, for more informations see http://broadinstitute.github.io/picard/command-line-overview.html#CreateSequenceDictionary



=head3 picardTools::picardToolsSortSam

This module sorts the input SAM or BAM
It takes at least two arguments: the file to sort, the name of the output file
The last argument is the options of picardTools sortSam, for more informations see http://broadinstitute.github.io/picard/command-line-overview.html#SortSam


=head1 AUTHORS

Intellectual property belongs to IRD, CIRAD and South Green developpement plateform 
Written by Cecile Monat, Ayite Kougbeadjo, Marilyne Summo, Cedric Farcy, Mawusse Agbessi, Christine Tranchant and Francois Sabot

=head1 SEE ALSO

L<http://www.southgreen.fr/>

=cut
