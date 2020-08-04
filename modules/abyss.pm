package abyss;

###################################################################################################################################
#
# Copyright 2014-2018 IRD-CIRAD-INRA-ADNid
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
# Version 3 written by Cecile Monat, Christine Tranchant, Laura Helou, Abdoulaye Diallo, Julie Orjuela-Bouniol, Sebastien Ravel, Gautier Sarah, and Francois Sabot
#
###################################################################################################################################

use strict;
use warnings;
use localConfig;
use toolbox;
use checkFormat;
use Data::Dumper;
use Switch;

#This function will validate that the given file is at least in one of the accepted format (fastq/fasta/sam/bam)
sub localFormatCheck
{
	my ($file)=@_;
	my $validation =0;
	switch (1)
	{
		case ($file =~ m/fastq$/i){$validation = 1 if (checkFormat::checkFormatFastq($file) == 1)}
		case ($file =~ m/fasta$/i){$validation = 1 if (checkFormat::checkFormatFasta($file) == 1)}
		case ($file =~ m/am$/i){$validation = 1 if (checkFormat::checkFormatSamOrBam($file) > 0)}
		
		else {toolbox::exportLog("ERROR: abyss : The file $file is not a FASTQ/FASTA/BAM/SAM file\n",0);}
	}
	
	return $validation;
	
}

#sub transAbyss
#{
#    my($outputDir,$readGroup,$firstListOfFile,$secondListOfFile,$optionsHachees)=@_;
#    ## if paired reads give $single as --left and $paired as --right; else give only $single
#    my $orderedList=""; 
#    my $command="";
#    my $options="";
#
#    if ($optionsHachees)
#    {
#         $options=toolbox::extractOptions($optionsHachees); ##Get given options
#    }
#
#    if ($options !~ m/--outdir\s\w+|--name\s\w+/) # The type of walker is not informed in the options
#    {
#        $options .= "--outdir ".$outputDir." --name ".$readGroup;
#    }
#    else
#    {
#        toolbox::exportLog("WARNING: abyss::transAbyss : You do not need to specify name or outdir options\n",2);
#    }
#    
#    if (scalar (@{$secondListOfFile}) > 0) # if is not empty,we have paired reads
#    {
#        if ((scalar (@{$firstListOfFile}) > 0 ) and (scalar (@{$firstListOfFile}) == scalar (@{$secondListOfFile})) )  # if is not empty and same length
#        {
#            my $i = 0;
#            foreach my $localFile (@{$firstListOfFile})       # for each pair of Fastq file(s)
#            {
#                if ($localFile ne "NA" and &localFormatCheck($localFile) == 1)
#                { ##Check if entry file exist and is not empty
#                    if (@{$secondListOfFile}[$i] ne "NA" and &localFormatCheck(@{$secondListOfFile}[$i])==1 )
#                    { ##Check if entry file exist and is not empty
#                        $orderedList .= " --pe ".$localFile." ".@{$secondListOfFile}[$i];    
#                        
#                    }
#                    else
#                    {
#                        toolbox::exportLog("ERROR: abyss::transAbyss : The file ".@{$secondListOfFile}[$i]." is uncorrect or is not a FASTQ/FASTA/SAM/BAM\n",0);
#                        return 0;#File not Ok
#                    }
#                }
#                else
#                {
#                    toolbox::exportLog("ERROR: abyss::transAbyss : The file ".$localFile." is uncorrect or is not a FASTQ/FASTA/SAM/BAM\n",0);
#                    return 0;#File not Ok
#                }
#                $i++;
#            }   
#            #print @{$firstListOfFastq},"\n";
#        }
#        else
#        {
#            toolbox::exportLog("ERROR: abyss::transAbyss : Incorrect list of files ".@{$firstListOfFile}." \n",0);
#            return 0;
#        }
#    }
#    else  # Single mode
#    {
#        if (scalar (@{$firstListOfFile}) > 0)  # if is not empty
#        {
#            $orderedList .= "--se ";
#            foreach my $localFile (@{$firstListOfFile})       # for each Fastq file(s)
#            {                 
#                if ($localFile ne "NA" and &localFormatCheck($localFile) == 1 )
#                { ##Check if entry file exist and is not empty
#                    $orderedList .= $localFile." ";
#                }
#                else
#                {
#                    toolbox::exportLog("ERROR: abyss::transAbyss : The file $localFile is uncorrect or is not a SAM/BAM/FASTA/FASTQ\n",0);
#                    return 0;#File not Ok
#                }
#            }
#        }
#        else
#        {
#            toolbox::exportLog("ERROR: abyss::transAbyss : Incorrect list of files ".@{$firstListOfFile}." \n",0);
#            return 0;
#        }
#    }
#           
#    $command=$transAbyss." ".$orderedList." ".$options;
#
#    #DEBUG toolbox::exportLog("RUN: transabyss::transabyssRun : CMD:: $command\n",1);
#    
#    #Execute command
#    if(toolbox::run($command)==1)
#    {
#        
#         return 1;#Command Ok
#    }
#    else
#    {
#         toolbox::exportLog("ERROR: abyss::transAbyss: Uncorrectly done\n",0);
#         return 0;#Command not Ok
#    }
#}

sub abyssSimple
{
    #Will use FASTQ/FASTA/SAM/BAM data for assembly, single-end and pair-end FASTQ/FASTA data or single library BAM/SAM
    my($outputDir,$readGroup,$forwardFile,$reverseFile,$optionsHachees)=@_;
    my $orderedList=""; 
    my $options="";
    
    if (ref $reverseFile)
    {
        #No reverse file, probably a BAM/SAM file
        $optionsHachees = $reverseFile;
        $reverseFile = "NA";
    }

    if ($optionsHachees)
    {
         $options=toolbox::extractOptions($optionsHachees,"="); ##Get given options
    }
    
    if (&localFormatCheck($forwardFile) != 1 )
    { #Not the good format
        toolbox::exportLog("ERROR: abyss::abyssSimple : The file $forwardFile is uncorrect or is not a SAM/BAM/FASTA/FASTQ\n",0);
        return 0;#File not Ok
    }
    if ($reverseFile ne "NA" && &localFormatCheck($reverseFile) != 1 )
    { #Not the good format
        toolbox::exportLog("ERROR: abyss::abyssSimple : The file $reverseFile is uncorrect or is not a SAM/BAM/FASTA/FASTQ\n",0);
        return 0;#File not Ok
    }
    
	#De-gzip files
	if ($forwardFile =~ m/\.gz$/)
	{
		my $plainFile = toolbox::degzip($forwardFile);
		toolbox::exportLog("ERROR: abyss::abyssSimple: Cannot decompress the gzip file $forwardFile", 0) if $plainFile == 0;
		$forwardFile = $plainFile;
	}
	if ($reverseFile ne "NA" && $reverseFile =~ m/\.gz$/)
	{
		my $plainFile = toolbox::degzip($reverseFile);
		toolbox::exportLog("ERROR: abyss::abyssSimple: Cannot decompress the gzip file $reverseFile", 0) if $plainFile == 0;
		$reverseFile = $plainFile;
	}
	
	#Checking Mandatory options
	
    if ($options !~ m/ k=\d+/)
	{
		#No k-mer value is affected, will fix it at 65
		toolbox::exportLog("WARNING: abyss::abyssSimple : No k-mer value (k option) provided, fixing it at 65\n", 2);
		$options .= " k=65 ";
		$options =~ s/\s+/ /g; #Removing extraspaces
	}
	if ($options =~ m/name=/)
	{
		#The name is fixed by the user, but no need
		toolbox::exportLog("WARNING: abyss::abyssSimple : The name option (name=) will be fixed by TOGGLe automatically\n", 2);
		$options =~ s/name=\w+//;
	}
	
    my $command = "$abyss -C $outputDir $options name=".$readGroup.".ABYSS in='".$forwardFile."";
    
    if ($reverseFile ne "NA")
    {
        $command .=" ".$reverseFile."'";
    }
    else
    {
        $command .="'";
    }
    
    #Execute command
    if(toolbox::run($command)==1)
    {
        return 1;#Command Ok
    }
    else
    {
         toolbox::exportLog("ERROR: abyss::abyssSimple : Uncorrectly done\n",0);
         return 0;#Command not Ok
    }
}

1;

=head1 NAME

Package I<abyss> 

=head1 SYNOPSIS

    use abyss;
	
	abyss::transAbyss($outputDir,$readGroup,$firstListOfFile,$secondListOfFile,$optionsHachees);
	
	abyss::abyssSimple ($outputDir,$readGroup,$forwardFile,$reverseFile,$optionsHachees);
    
    localFormatCheck($file);
    
        
=head1 DESCRIPTION

    this module allows to run Abyss and Trans-Abyss
    
=head2 FUNCTIONS

=head3 abyss::transAbyss($folderIn,tmpRoot)

	Transcriptome assembly
	
=head3  abyss::abyssSimple ($localFolder,$folderIn)

	Genome assembly from single-end and pair-end data, single library

=head3  localFormatCheck($file)

    Local format check tester
    
=head1 AUTHORS

Intellectual property belongs to IRD, CIRAD and South Green developpement plateform for all versions also for ADNid for v2 and v3 and INRA for v3
Written by Christine Tranchant, Cecile Monat, Laura Helou, Abdoulaye Diallo, Julie Orjuela-Bouniol, Sebastien Ravel, Gautier Sarah, and Francois Sabot

=head1 SEE ALSO

L<http://toggle.southgreen.fr/>

=cut