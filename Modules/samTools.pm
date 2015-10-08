package samTools;



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

##############################################
##SamTools
##Module containing SamTools functions
##############################################
##
##
##SAMTOOLS Faidx
#Index reference sequence in the FASTA format or extract subsequence from indexed reference sequence.
sub samToolsFaidx{
     my($refFastaFileIn)=@_;
     if (toolbox::sizeFile($refFastaFileIn)==1){ ##Check if entry file exist and is not empty
                    
          #Execute command
          my $command=$samtools." faidx ".$refFastaFileIn;
          if(toolbox::run($command)==1){
               toolbox::exportLog("INFOS: samTools::samToolsFaidx : Correctly done\n",1);
               return 1;#Command Ok
          }else{
               toolbox::exportLog("ERROR: samTools::samToolsFaidx : Uncorrecly done\n",0);
               return 0;#Error in command
          }
     }else{
        toolbox::exportLog("ERROR: samTools::samToolsFaidx : The file $refFastaFileIn is incorrect;failed\n",0);
        return 0;#Error in the file itself
     }
}
##
##
##SAMTOOLS INDEX
#Index sorted alignment for fast random access.
sub samToolsIndex{
     my($bamFileIn)=@_;
     if (toolbox::sizeFile($bamFileIn)==1){ ##Check if entry file exist and is not empty
          
          #Check if the format is correct
          if (toolbox::checkSamOrBamFormat($bamFileIn)==0) {#The file is not a BAM/SAM file
               toolbox::exportLog("ERROR: samTools::samToolsIndex : The file $bamFileIn is not a SAM/BAM file\n",0);
               return 0;
          }
          
          #Execute the command
          my $command=$samtools." index ".$bamFileIn;
          #Execute command
          if(toolbox::run($command)==1){
               toolbox::exportLog("INFOS: samTools::samToolsIndex : Correctly done\n",1);
               return 1;#Command Ok
          }else{
               toolbox::exportLog("ERROR: samTools::samToolsIndex : Uncorrectly done\n",0);
               return 0;#Command not Ok
               }
     }else{
        toolbox::exportLog("ERROR: samTools::samToolsIndex : The file $bamFileIn is incorrect\n",0);
        return 0;#The file is incorrect
    }
}

##
##
##SAMTOOLS VIEW
#Extract/print all or sub alignments in SAM or BAM format.
sub samToolsView{
     my($bamFileIn,$bamFileOut,$optionsHachees)=@_;
     if (toolbox::sizeFile($bamFileIn)==1){ ##Check if entry file exist and is not empty
          
          #Check if the format is correct
          if (toolbox::checkSamOrBamFormat($bamFileIn)==0) {#The file is not a BAM/SAM file
               toolbox::exportLog("ERROR: samTools::samToolsView : The file $bamFileIn is not a SAM/BAM file\n",0);
               return 0;
          }

          my $options="";
          if ($optionsHachees) {
               $options=toolbox::extractOptions($optionsHachees); ##Get given options
          }
          my $command=$samtools." view ".$options." -o ".$bamFileOut." ".$bamFileIn;
          #toolbox::exportLog($command."\n",1);
          #Execute command
          if(toolbox::run($command)==1){
               toolbox::exportLog("INFOS: samTools::samToolsView : Correctly done\n",1);
               return 1;#Command Ok
          }else{
               toolbox::exportLog("ERROR: samTools::samToolsView : Uncorrectly done\n",0);
               return 0;#Command not Ok
               }
     }else{
        toolbox::exportLog("ERROR: samTools::samToolsView : The file $bamFileIn is uncorrect\n",0);
        return 0;#File not Ok
    }
}

##
##
##SAMTOOLS SORT
#Sort alignments by leftmost coordinates.
sub samToolsSort{
     my($bamFileIn,$optionsHachees)=@_;
     if (toolbox::sizeFile($bamFileIn)==1){ ##Check if entry file exist and is not empty
          
          #Check if the format is correct
          if (toolbox::checkSamOrBamFormat($bamFileIn)==0) {#The file is not a BAM/SAM file
               toolbox::exportLog("ERROR: samTools::samToolsSort : The file $bamFileIn is not a SAM/BAM file\n",0);
               return 0;
          }
          
          my $bamPath=toolbox::extractPath($bamFileIn);
          my $bamFileName=toolbox::extractName($bamFileIn).".SAMTOOLSSORT";
          my $bamFileOut = $bamPath.$bamFileName;
          my $options="";
          if ($optionsHachees) {
               $options=toolbox::extractOptions($optionsHachees);
          }
          my $command=$samtools." sort ".$options." ".$bamFileIn." ".$bamFileOut;
          #Execute command
          if(toolbox::run($command)==1){
               toolbox::exportLog("INFOS: samTools::samToolsSort : Correctly done\n",1);
               return 1;#Command Ok
          }else{
                    toolbox::exportLog("ERROR: samTools::samToolsSort : Uncorrectly done\n",0);
                    return 0;#Command not Ok
               }
     }else{
        toolbox::exportLog("ERROR: samTools::samToolsSort : The file $bamFileIn is uncorrect\n",0);
        return 0;#File not Ok
     }
}
##MergeHeader
##Create Header containing all bam  headers (for merging), by adding the new RG in a currently existing header
sub mergeHeader{
     my($bamFiles,$outHeader)=@_;
     open(FILE, ">> $outHeader") or toolbox::exportLog("ERROR: samTools::mergeHeader : Cannot open $outHeader\n$!\n",0); #Create the outfile (or re-open it)
     my @listOfBam=@{$bamFiles};#De-reference the list of BAMs
     foreach my $bamfile (@listOfBam){
            toolbox::checkSamOrBamFormat($bamfile);#Check if the file is a correct one (SAM or BAM)
            my $command=$samtools." view -H ".$bamfile." | grep '\@RG'"; #Just pick up the ReadGroup line
            my $header_char=`$command`;#Compute the command
            print FILE $header_char;#Add the RG info in the new header file
        }
     return 1;
}


##SAMTOOLS MERGE
#Merge multiple sorted alignments.
sub samToolsMerge{
     my($bamFiles,$bamOutFile,$header,$optionsHachees)=@_;
    
     my $bamToMerge; # set of correct bam, ready for command launch
     foreach my $putativeBam (@{$bamFiles})
          {
          if (toolbox::sizeFile($putativeBam)==1)
               {#The file exists and is not empty
               #check if it is a BAM
               if (toolbox::checkSamOrBamFormat($putativeBam)==2)
                    {
                    #The file is a BAM
                    $bamToMerge.=" ".$putativeBam;
                    }
               else
                    {#The file is not a BAM
                    toolbox::exportLog("ERROR: samTools::samToolsDepth : The file $putativeBam is not a BAM file\n",0);
                    return 0;  
                    }
               }
          else #The file does not exist
               {
               toolbox::exportLog("ERROR: samTools::samToolsDepth : The file $putativeBam does not exist\n",0);
               return 0;
               }
          }
     #PickUp options
     my $options="";
     if ($optionsHachees) {
          $options=toolbox::extractOptions($optionsHachees);
     }
         
     my $command=$samtools." merge -f ".$options."-h ".$header." ".$bamOutFile." ".$bamToMerge;#Using the SAMtools merge command
     #Execute command
     if(toolbox::run($command)==1){
          toolbox::exportLog("INFOS: samTools::samToolsMerge : Correctly done\n",1);
          return 1;#Command Ok
          }
     else{
          toolbox::exportLog("ERROR: samTools::samToolsMerge : Uncorrectly done\n",0);
          return 0;#Command not Ok
          }
}


##SAMTOOLS IDXSTATS
# samtools idxstats <aln.bam>
#Retrieve and print stats in the index file. The output is TAB delimited with each line consisting of reference sequence name, sequence length, # mapped reads and # unmapped reads. 
    sub samToolsIdxstats{
	my($bamFileIn,$idxstatsOutput)=@_;
	if (toolbox::sizeFile($bamFileIn)==1){ ##Check if entry file exists and is not empty
          #Check if the format is correct
          if (toolbox::checkSamOrBamFormat($bamFileIn)==0) {#The file is not a BAM/SAM file
               toolbox::exportLog("ERROR: samTools::samToolsIdxstats : The file $bamFileIn is not a SAM/BAM file\n",0);
               return 0;
          }
          my $command="$samtools idxstats $bamFileIn > $idxstatsOutput";#Command to be executed
          #Execute command
          if(toolbox::run($command)==1){
               toolbox::exportLog("INFOS: samTools::samToolsIdxstats : Correctly done\n",1);
               return 1;#Command Ok
          }else{
                    toolbox::exportLog("ERROR: samTools::samToolsIdxstats : Uncorrectly done\n",0);
                    return 0;#Command not ok
               }
     }else{
        toolbox::exportLog("ERROR: samTools::samToolsIdxstats : The file $bamFileIn is uncorrect\n",0);
        return 0;#Not a good file
     }
}

##
##
##SAMTOOLS DEPTH
#Compute the depth on a BAM file (or a set of bam) and provide a tabular file with the number of reads per positions.
sub samToolsDepth {
     my ($bamFilesList, $depthOutputFile,$optionsHachees)=@_;
     #Testing if sent bam are true bam and recover them
     my $bamFilesToCompute; # set of correct bam, ready for command launch
     foreach my $putativeBam (@{$bamFilesList})
          {
          if (toolbox::sizeFile($putativeBam)==1)
               {#The file exists and is not empty
               #check if it is a BAM
               if (toolbox::checkSamOrBamFormat($putativeBam)==2)
                    {
                    #The file is a BAM
                    $bamFilesToCompute.=" ".$putativeBam;
                    }
               else
                    {#The file is not a BAM
                    toolbox::exportLog("ERROR: samTools::samToolsDepth : The file $putativeBam is not a BAM file\n",0);
                    return 0;  
                    }
               }
          else #The file does not exist
               {
               toolbox::exportLog("ERROR: samTools::samToolsDepth : The file $putativeBam does not exist\n",0);
               return 0;
               }
          }
     #We have the complete BAM list validated and arranged, let's launch the command!!
     my $options="";
     if ($optionsHachees)#Picking up the options
          {
          $options=toolbox::extractOptions($optionsHachees);
          }
     
     my $command = "$samtools depth $options $bamFilesToCompute > $depthOutputFile";#Command to be launched
     if(toolbox::run($command)==1)
          {
          toolbox::exportLog("INFOS: samTools::samToolsDepth : Correctly done\n",1);
          return 1;#Command Ok
          }
     else
          {
          toolbox::exportLog("ERROR: samTools::samToolsDepth : Uncorrectly done\n",0);
          return 0;#Command not Ok
          } 
}

##
##
##SAMTOOLS FLAGSTATS
#Provide simple stats on a BAM/SAM file.
sub samToolsFlagstat
{
	my ($bamFile,$flagstatsOutputFile)=@_;
	if (toolbox::sizeFile($bamFile)==1)##Check if entry file exist and is not empty
	{ 
		#Check if the format is correct
		if (toolbox::checkSamOrBamFormat($bamFile)!=2)#The file is not a BAM file
		{
			if (toolbox::checkSamOrBamFormat($bamFile)==1)#The file is a SAM file
			{
				toolbox::exportLog("ERROR: samTools::samToolsFlagstat : The file $bamFile is a SAM file and samTools::samToolsFlagstat accept only BAM file\n",0);
				return 0;
			}
			elsif (toolbox::checkSamOrBamFormat($bamFile)==0)#The file is not a SAM/BAM file
			{
				toolbox::exportLog("ERROR: samTools::samToolsFlagstat : The file $bamFile is not either a SAM or BAM file and samTools::samToolsFlagstat accept only BAM file\n",0);
				return 0;
			}
		}
		else ##We have the BAM validated
		{
			toolbox::exportLog("INFOS: samTools::samToolsFlagstat : The file $bamFile is a BAM file\n",1);
			my $command = "$samtools flagstat $bamFile > $flagstatsOutputFile";#Command to launch
			if(toolbox::run($command)==1)
			{
				toolbox::exportLog("INFOS: samTools::samToolsFlagstat : Correctly done\n",1);
				return 1;
			}
			else
			{
				toolbox::exportLog("ERROR: samTools::samToolsFlagstat : Uncorrectly done\n",0);
				return 0;
			} 
		}
	}
	else ##The file does not exist or is empty
	{
		toolbox::exportLog("ERROR: samTools::samToolsFlagstat : The file $bamFile does not exist or is empty\n",0);
		return 0;
	}
}

1;

=head1 NAME

    Package I<samtools> 

=head1 SYNOPSIS

    This package contains the whole modules for the SAMtools software

=head1 DESCRIPTION

    Package SAMtools (Li et al, 2009, http://http://www.htslib.org/ ) is a software package for working on SAM and BAM files (sorting, selection, merging...)

=head2 Functions

=over 4

=item mergeHeader (Merges the headers from various sources in a single one)

=item samToolsDepth (Computes the depth coverage using a BAM file)

=item samToolsFaidx (Generates an index for a given Fasta file)

=item samToolsFlagstat (Generates simple statistics for a BAM file. Provides the number of correctly mapped, etc...)

=item samToolsIdxstats (Generates more complex BAM statistics)

=item samToolsIndex (Generates an index for a BAM file)

=item samToolsMerge (Merges two or more BAMs in a single resulting BAM. The header will be the same as of the first BAM.)

=item samToolsSort (Sorts in different ways a BAM file: coordinates, random, reads)

=item samToolsView (Allows to see a BAM in SAM, or to extract a subregion or with specific qualifiers - correctly mapped...)

=back

=head1 AUTHORS

Intellectual property belongs to IRD, CIRAD and South Green developpement plateform 
Written by Cecile Monat, Ayite Kougbeadjo, Marilyne Summo, Cedric Farcy, Mawusse Agbessi, Christine Tranchant and Francois Sabot

=head1 SEE ALSO

L<http://www.southgreen.fr/>		# SOUTH GREEN

=cut
