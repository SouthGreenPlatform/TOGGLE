package bedTools;

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
use Switch;

#This function will validate that the given file is at least in one of the accepted format (vcf/bam/gff/bed)
sub localFormatCheck{
	my ($file)=@_;
	my $validation =0;
	switch (1)
	{
		case ($file =~ m/vcf$/i){$validation = 1 if (checkFormat::checkFormatVcf($file) == 1)}
		case ($file =~ m/bam$/i){$validation = 1 if (checkFormat::checkFormatSamOrBam($file) == 2)}
		case ($file =~ m/gff$/i){$validation = 1 if (checkFormat::checkFormatGff($file) == 1)}
		case ($file =~ m/bed$/i){$validation = 1 if (checkFormat::checkFormatBed($file) == 1)}
		else {toolbox::exportLog("ERROR: bedtools : The file $file is not a BAM/GFF/VCF/BED file\n",0);}
	}
	
	return $validation;
	
}

#BEDtools intersectBed/intersect
sub intersectBed {
	my($fileIn,$fileOut,$optionsHachees)=@_;
	if (toolbox::sizeFile($fileIn)==1)
	{ ##Check if entry file exist and is not empty
		
		#Check if the format is correct
		if (&localFormatCheck($fileIn) == 0)
		{#The file is not a BAM/VCF/GFF/BED file
			toolbox::exportLog("ERROR: bedtools::intersectBed : The file $fileIn is not a BAM/GFF/VCF/BED file\n",0);
			return 0;
		}
		
		my $options="";
		
		if ($optionsHachees)
		{
			$options=toolbox::extractOptions($optionsHachees);
		}
		else
		{
				toolbox::exportLog("ERROR: bedtools::intersectBed : No options provided, we cannot execute any command\n",0);
				return 0;
		}
		
		#Picking up the second file
		my $crossFile = $options;
		$crossFile =~ s/.*-[a|b] ([A-Za-z0-9\/-_\.~]+).*/$1/;
		$crossFile = `readlink -f $crossFile`;
		chomp $crossFile;
		##DEBUG toolbox::exportLog("file crossed: $crossFile",2);
		if (&localFormatCheck($crossFile) == 0)
		{#The file is not a BAM/VCF/GFF/BED file
			toolbox::exportLog("ERROR: bedtools::intersectBed : The file $crossFile is not a BAM/GFF/VCF/BED file\n",0);
			return 0;
		}
		
		#The second file and other options are within the global file options
		
		my $variableQualifier = "-a";
				
		if ($options =~ m/ -a /) #The user has specified another file as the -a file
		{
			$variableQualifier = "-b"
		}
		
		my $command=$bedtools." intersect ".$options." ".$variableQualifier." ".$fileIn." > ".$fileOut;
		
		#Execute command
		if(toolbox::run($command)==1)
		{
			return 1;#Command Ok
		}
		else
		{
			toolbox::exportLog("ERROR: bedtools::intersectBed : Uncorrectly done\n",0);
			return 0;#Command not Ok
		}
	}
	else
	{
	   toolbox::exportLog("ERROR: bedtools::intersectBed : The file $fileIn is uncorrect\n",0);
	   return 0;#File not Ok
	}
}

#BEDtools intersectBed/intersect
sub windowBed {
	my($fileIn,$fileOut,$optionsHachees)=@_;
	if (toolbox::sizeFile($fileIn)==1)
	{ ##Check if entry file exist and is not empty
		
		#Check if the format is correct
		if (&localFormatCheck($fileIn) == 0)
		{#The file is not a BAM/VCF/GFF/BED file
			toolbox::exportLog("ERROR: bedtools::windowBed : The file $fileIn is not a BAM/GFF/VCF/BED file\n",0);
			return 0;
		}
		
		my $options="";
		
		if ($optionsHachees)
		{
			$options=toolbox::extractOptions($optionsHachees);
		}
		else
		{
				toolbox::exportLog("ERROR: bedtools::windowBed : No options provided, we cannot execute any command\n",0);
				return 0;
		}
		
		#Picking up the second file
		my $crossFile = $options;
		$crossFile =~ s/.*-[a|b] ([A-Za-z0-9\/-_\.~]+).*/$1/;
		$crossFile = `readlink -f $crossFile`;
		chomp $crossFile;
		if (&localFormatCheck($crossFile) == 0)
		{#The file is not a BAM/VCF/GFF/BED file
			toolbox::exportLog("ERROR: bedtools::windowBed : The file $crossFile is not a BAM/GFF/VCF/BED file\n",0);
			return 0;
		}
		
		#The second file and other options are within the global file options
		
		my $variableQualifier = "-a";
		
		if ($options =~ m/ -a /) #The user has specified another file as the -a file
		{
			$variableQualifier = "-b"
		}
		
		my $command=$bedtools." window ".$options." ".$variableQualifier." ".$fileIn." > ".$fileOut;
		
		#Execute command
		if(toolbox::run($command)==1)
		{
			return 1;#Command Ok
		}
		else
		{
			toolbox::exportLog("ERROR: bedtools::windowBed : Uncorrectly done\n",0);
			return 0;#Command not Ok
		}
	}
	else
	{
	   toolbox::exportLog("ERROR: bedtools::windowBed : The file $fileIn is uncorrect\n",0);
	   return 0;#File not Ok
	}
}

#This command allows user to launch ANY bedtools command using the bedtools subprogram
sub generic {
	my($fileIn,$fileOut,$optionsHachees)=@_;
	if (toolbox::sizeFile($fileIn)==1)
	{ ##Check if entry file exist and is not empty
		
		#Check if the format is correct
		if (&localFormatCheck($fileIn) == 0)
		{#The file is not a BAM/VCF/GFF/BED file
			toolbox::exportLog("ERROR: bedtools::generic : The file $fileIn is not a BAM/GFF/VCF/BED file\n",0);
			return 0;
		}
		
		my $options="";
		
		if ($optionsHachees)
		{
			$options=toolbox::extractOptions($optionsHachees);
		}
		else
		{
				toolbox::exportLog("ERROR: bedtools::generic : No options and commands provided, we cannot execute any command\n",0);
				return 0;
		}
		
		#The generic command system will transform the FILEIN text by the correct FILENAME
		$options =~ s/FILEIN/$fileIn/i;
		
		#The generic command system will transform the FILEOUT text by the correct FILENAME
		$options =~ s/FILEOUT/$fileOut/i;
		my $command;
		if ($options =~ m/$fileOut/)
		{
			$command=$bedtools." ".$options;
		}
		else
		{
			$command=$bedtools." ".$options." > ".$fileOut;
		}
		
		#Execute command
		if(toolbox::run($command)==1)
		{
			return 1;#Command Ok
		}
		else
		{
			toolbox::exportLog("ERROR: bedtools::generic : Uncorrectly done\n",0);
			return 0;#Command not Ok
		}
	}
	else
	{
	   toolbox::exportLog("ERROR: bedtools::generic : The file $fileIn is uncorrect\n",0);
	   return 0;#File not Ok
	}
}

1;
=head1 NAME

    Package I<bedtools>

=head1 SYNOPSIS

	use bedtools;

	

=head1 DESCRIPTION

    Package BEDtools (http:// ) is a software package for 

=head2 FUNCTIONS

=head3 bedtools::generic



=head1 AUTHORS

Intellectual property belongs to IRD, CIRAD and South Green developpement plateform for all versions also for ADNid for v2 and v3 and INRA for v3
Written by Christine Tranchant, Cecile Monat, Laura Helou, Abdoulaye Diallo, Julie Orjuela-Bouniol, Sebastien Ravel, Gautier Sarah, and Francois Sabot

=head1 SEE ALSO

L<http://toggle.southgreen.fr/>

=cut