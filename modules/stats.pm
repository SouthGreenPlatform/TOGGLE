package stats;

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
use Data::Dumper;
use checkFormat;

use picardTools;
use samTools;
use pairing;
use eaUtils;

##############################################
##stats
##Module containing basic stats functions
##used to write reports
##############################################
##
##

## creatingMappingStatFileRaw
# Execute samtools flagstat / idxstats to generate raw files stat related to sam/bam files
sub creatingMappingStatFileRaw
{
     
     # Getting arguments
     my($bamFileIn,$bamFileOut)=@_;
     
     my $softParameters;
     $softParameters->{"CREATE_INDEX"}="TRUE";
     $softParameters->{"SORT_ORDER"}="coordinate";
     $softParameters->{"VALIDATION_STRINGENCY"}="SILENT";     
        
     # Check if the sam or bam file exists and is not empty
     if (toolbox::sizeFile($bamFileIn)==1)
     {     
          #Check if the format is not correct neither sam or ban
          if (checkFormat::checkFormatSamOrBam($bamFileIn)==0)
          {
               #The file is not a BAM/SAM file
               toolbox::exportLog("ERROR: stats::creatingMappingStatFileRaw : The file $bamFileIn is not a SAM/BAM file\n",0);
               return 0;
          }

          # Automatic sorting sam/bam file, indexing it and generating a bam if sam file          
          picardTools::picardToolsSortSam($bamFileIn,$bamFileOut,$softParameters);
     
          # flagstat
          my $flagstatsOutputFile=$bamFileOut;
          $flagstatsOutputFile =~ s/\.bam$/\.flagstat\.mapping\.stat/;
          samTools::samToolsFlagstat($bamFileOut,$flagstatsOutputFile);
          
          # idxstats
          my $idxstatsOutputFile=$bamFileOut;
          $idxstatsOutputFile =~ s/\.bam$/\.idxstats\.mapping\.stat/;
          samTools::samToolsIdxstats($bamFileOut,$idxstatsOutputFile);

          # Delete bam generated
          my $rmCommand = "rm $bamFileOut*";
          toolbox::run($rmCommand,"noprint");
          
          # Delete bai generated
          $bamFileOut =~ s/\.bam$/\.bai/;
          $rmCommand = "rm $bamFileOut*";
          toolbox::run($rmCommand,"noprint");
     }
     else
     {
        toolbox::exportLog("ERROR: stats::mappingStat : The file $bamFileIn is uncorrect\n",0);
        return 0;#File not Ok
     }
     
    # creatingMappingStatFileTex
}



# Execute grep to count polymorphism on a vcf
sub creatingCallingStatFileRaw
{  
     # Getting arguments
     my ($vcfFileIn)=@_;
         		
     # Check if the vcf file exists and is not empty
     if (toolbox::sizeFile($vcfFileIn)==1)
     {      
          my $vcfOutputFile=$vcfFileIn;
          $vcfOutputFile =~ s/\.vcf$/\.calling\.stat/;
         
          my $grepcmd='grep "#" '.$vcfFileIn ." -v -c > $vcfOutputFile " ;
          system($grepcmd) and die "ERROR: stats::creatingCallingStatFileRaw : The command $grepcmd is failed\n";  #toolbox::run($grepcmd); #,"noprint");        
     }
     else
     {
        toolbox::exportLog("ERROR: stats::creatingCallingStatFileRaw : The file $vcfFileIn is uncorrect\n",0);
        return 0;#File not Ok
     }
     
}

# Execute eaUtils::fastqStats on fastq files into directory
sub creatingFastqStatFileRaw
{  
     # Getting arguments
     my ($fastqDir )=@_;
    
     my $softParameters;
     $softParameters->{"-D"}="";
     
     my $fastqFiles=toolbox::readDir($fastqDir);
     foreach my $fastqFile (@{$fastqFiles})
     {
          my $statOutputFile = $fastqFile;
          $statOutputFile =~ s/\.fastq$|\.fastq.gz$/\.fastq\.stat/;
          eaUtils::fastqStats($fastqFile, $statOutputFile, $softParameters); 
     }

}




# Generate tex files from stats files in a statDir
sub creatingStatFileTex
{

     my ($statDir)=@_;		# get directory with stats file
	my $fileList = toolbox::readDir($statDir);		# get stat files list
	
     # variable initialisation
     my $texMapping='NA';
     my $texCalling='NA'; # store tex file generated on the fly by every part
     my $texFastq='NA';
     
     # open stat file tex  
     my $statTexFile = "stats.tex";
     open(my $texFh, ">>", $statTexFile) or toolbox::exportLog("$0 : open error of $statTexFile .... $!\n",0);
	     
	# Parsing raw stat files
     foreach my $file (@{$fileList}) #Copying the final data in the final directory
	{
          
          # open raw stat files generated by parallel and global analysis (finalresults) ### TODO : INTEGRATE INPUT ALSO
          open (my $fh, "<", $file) or toolbox::exportLog("$0 : open error of $file .... $!\n",0);
          
          # extract sample name (global or sample)    
          my $sample = pairing::extractName($file);
          
          ########### MAPPING PART - FLAGSTAT
		if ($file =~ /\.flagstat.mapping.stat$/)
		{
               # tex header (table header)
               $texMapping = "\\subsection{Mapping}
	\\begin{table}[ht]
         \\resizebox{\\textwidth}{!}{%
		\\centering
		\\begin{tabularx}{18cm}{X|c|c|c}
			Samples & Raw sequences & Mapped sequences & Properly mapped  \\\\\\hline  \n" if  ($texMapping eq 'NA');
          
			# variable initialisation (number of sequences)
			my $raw = 0;
               my $mapped = 0;
               my $properly = 0;
               
               # reading raw stat files generated by flagstat
			while (my $line = <$fh>)
			{    
                    # Extract #mapped / properly / raw read number from flasgtat file 
                    my (@split) = split /\s\+/ , $line;
                    my $val=$split[0];
                    
				if ($line =~ /in\stotal/) { $raw = $val; }
				elsif ($line =~ /\sproperly\spaired/) { $properly = $val;  }
				elsif ($line =~ /\smapped\s\(/) { $mapped = $val; }
			}

			$texMapping .= " $sample & $raw & $mapped (" . sprintf("%.2f",($mapped/$raw*100)) . " \\%) & $properly (". sprintf("%.2f",($properly/$raw*100)) . " \\%) \\\\ \n";
		
		}
          
          ########### CALLING PART - NB POLYMORPHISMs
          if ($file =~ /\.calling.stat$/)
		{
               # tex header (table header)
               $texCalling = "\\subsection{Calling}
	\\begin{table}[ht]
         \\resizebox{\\textwidth}{!}{%
		\\centering
		\\begin{tabularx}{15cm}{X|c}
			Samples & Polymorphisms detected   \\\\\\hline \n  " if  ($texCalling eq 'NA');
               
               # variable initialisation (number of polymorphism)
			my $polymorphism = 0;

			while (my $line = <$fh>)
			{    
                    chomp $line;
                    $texCalling .= " $sample & $line  \\\\ \n";
			}
			
		}
          
          ########### FASTQ PART - fastqd
		if ($file =~ /\.fastq.stat$/)
		{
               
               # tex header (table header)
               $texFastq = "\\subsection{Fastq stat}
	\\begin{table}[ht]
     \\resizebox{\\textwidth}{!}{%
		\\centering
		\\begin{tabularx}{18cm}{X|c|c|c|c|c|c|c|c|c|c|c|c|X}
Samples & seq & \\multicolumn{3}{c|}{Length} & \\multicolumn{3}{c|}{Qual} & \\%A & \\%C & \\%G & \\%T & \\%N & Total Bases    \\\\ \\hline
& & Len. & Mean & Min & Min & Max & Mean & & & & & & \\\\" if  ($texFastq eq 'NA');

			# variable initialisation (number of sequences)
			my $read = 0;
               my $len = 0;
               my $lenMean = 0;   
               my $lenMin = 0;
               my $qualMin = 0;
               my $qualMax = 0;
               my $qualMean = 0;
               my $A=0;
               my $C=0;
               my $G=0;
               my $T=0;
               my $N=0;
               my $totBase=0;
               
               # reading raw stat files generated by fastq stats
			while (my $line = <$fh>)
			{
                    chomp $line;
                    
                    # Extract value from fastqstat file 
                    my (@split) = split /\s+/ , $line;
                    my $val=$split[$#split];
                    
				if ($line =~ /^reads/) { $read = $val; }
				elsif ($line =~ /^len\t/) { $len = $val;  }
				elsif ($line =~ /^len\smean/) { $lenMean = sprintf("%.2f",($val)); }
                    elsif ($line =~ /^len\smin/) { $lenMin = $val; }
                    elsif ($line =~ /^qual\smin/) { $qualMin = $val; }
                    elsif ($line =~ /^qual\smax/) { $qualMax = $val; }
                    elsif ($line =~ /^qual\smean/) { $qualMean = sprintf("%.2f",($val)); }
                    elsif ($line =~ /^%A/) { $A = sprintf("%.2f",($val)); }
                    elsif ($line =~ /^%C/) { $C = sprintf("%.2f",($val)); }
                    elsif ($line =~ /^%G/) { $G = sprintf("%.2f",($val)); }
                    elsif ($line =~ /^%T/) { $T = sprintf("%.2f",($val)); }
                    elsif ($line =~ /^%N/) { $N = sprintf("%.2f",($val)); }
                    elsif ($line =~ /^total\sbases/) { $totBase = $val; }
			}

			$texFastq .= " $sample & $read & $len & $lenMean & $lenMin & $qualMin & $qualMax & $qualMean & $A & $C  & $G & $T & $N & $totBase \\\\ \n";
		
		}
          close $fh;
	}
     
	my $texFooter =  " \\end{tabularx}}
\\end{table}";
	
     if ($texFastq ne 'NA') { $texFastq .= $texFooter } else { $texFastq = ''; }
     if ($texCalling ne 'NA') { $texCalling .= $texFooter } else { $texCalling = ''  ; }
     if ($texMapping ne 'NA') { $texMapping .= $texFooter } else { $texMapping = ''  ; }
     
     print $texFh $texFastq."\n" . $texMapping ."\n" . $texCalling;
	close $texFh;
}





1;

=head1 NAME

    Package I<stats> 

=head1 SYNOPSIS

    This package contains the whole modules for the ...

=head1 DESCRIPTION

   ....

=head2 Functions

=over 4

=item mergeHeader (Merges the headers from various sources in a single one)


=back

=head1 AUTHORS

Intellectual property belongs to IRD, CIRAD and South Green developpement plateform 
Written by Cecile Monat, Ayite Kougbeadjo, Marilyne Summo, Cedric Farcy, Mawusse Agbessi, Christine Tranchant and Francois Sabot

=head1 SEE ALSO

L<http://www.southgreen.fr/>		# SOUTH GREEN

=cut
