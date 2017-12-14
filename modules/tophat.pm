package tophat;

###################################################################################################################################
#
# Copyright 2014-2017 IRD-CIRAD-INRA-ADNid
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
use pairing;
use Data::Dumper;

#################################################################################################
## sub bowtieBuild : builds a Bowtie index from a set of DNA sequences.
#################################################################################################
## arguments : fasta file to index and options used for bowtieBuild running
## Returns prefixname of the database created  (1 if the execution is correctly done else 0)
#################################################################################################
#sub bowtieBuild
#{
#	my ($refFastaFileIn,$optionsHachees)=@_;
#	##DEBUG toolbox::exportLog("DEBUG: tophat::bowtieBuild : $prefixRef\n",1);
#
#	if (toolbox::sizeFile($refFastaFileIn)==1)						# check if the reference file exist and is not empty
#	{
#		my $options=toolbox::extractOptions($optionsHachees, " ");			# Get given options
#		my $command=$bowtieBuild.$options." ".$refFastaFileIn." ".$refFastaFileIn;		# command
#		##DEBUG toolbox::exportLog("DEBUG: tophat::bowtieBuild : $command\n",1);
#		# Execute command
#		if(toolbox::run($command)==1)							# The command should be executed correctly (ie return) before exporting the log
#	{
#			toolbox::exportLog("INFOS: tophat::bowtieBuild : correctly done\n",1);	# bowtiebuild have been correctly done
#		}
#		else
#		{
#			toolbox::exportLog("ERROR: tophat::bowtieBuild : ABORTED\n",0);		# bowtiebuild have not been correctly done
#		}
#	}
#	else
#	{
#		toolbox::exportLog("ERROR: tophat::bowtiebBuild : Problem with the file $refFastaFileIn\n",0);		# bowtiebuild can not function because of wrong/missing reference file
#	}
#	return $refFastaFileIn;
#}
#################################################################################################
## END sub bowtieBuild
#################################################################################################
#
#
#
#################################################################################################
## sub bowtie2Build : builds a Bowtie index from a set of DNA sequences.
#################################################################################################
## arguments : fasta file to index and options used for bowtie2Build running
## Returns prefixname of the database created  (1 if the execution is correctly done else 0)
#################################################################################################
#sub bowtie2Build
#{
#	my($refFastaFileIn,$optionsHachees)=@_;
#	##DEBUG toolbox::exportLog("DEBUG: tophat::bowtie2Build : $prefixRef\n",1);
#
#	if (toolbox::sizeFile($refFastaFileIn)==1)						# Check if the reference file exist and is not empty
#	{
#		my $options=toolbox::extractOptions($optionsHachees, " ");			# Get given options
#		my $command=$bowtie2Build.$options." ".$refFastaFileIn." ".$refFastaFileIn;		# command
#		##DEBUG
#		toolbox::exportLog("INFOS: tophat::bowtie2Build : $command\n",1);
#		#Execute command
#		if(toolbox::run($command)==1)							# The command should be executed correctly (ie return) before exporting the log
#	{
#			toolbox::exportLog("INFOS: tophat::bowtie2Build : correctly done\n",1);	# bowtie2build have been correctly done
#		}
#		else
#		{
#			toolbox::exportLog("ERROR: tophat::bowtie2Build : ABORTED\n",0);		# bowtie2build have not been correctly done
#			return 0;
#		}
#	}
#	else
#	{
#		toolbox::exportLog("ERROR: tophat::bowtie2bBuild : Problem with the file $refFastaFileIn\n",0);		# bowtie2build can not function because of wrong/missing reference file
#		return 0;
#	}
#	return $refFastaFileIn;
#}
#################################################################################################
## END sub bowtie2Build
#################################################################################################



################################################################################################
# sub tophat2  : realize mapping with tophat
################################################################################################
sub tophat2
{
	my ($tophatDirOut, $prefixRef, $forwardFastqFileIn,$reverseFastqFileIn,$gffFile,$optionsHachees)=@_;
	my $options="";
	if ($optionsHachees)
	{
		$options=toolbox::extractOptions($optionsHachees);		##Get given options
	}
	my $command = "";
	if ((toolbox::sizeFile($forwardFastqFileIn)==1) and not (defined $reverseFastqFileIn))		##Check if entry files exist and are not empty / single mode
	{
		if ($gffFile ne "None" and toolbox::sizeFile($gffFile)==1)
		{
			$command=$tophat2.$options." -p 8 -G ".$gffFile." -o ".$tophatDirOut." ".$prefixRef." ".$forwardFastqFileIn;		# command line
			toolbox::exportLog("INFOS: tophat::topHat2 : $command\n",1);
		}
		else
		{
			$command=$tophat2.$options." -p 8 -o ".$tophatDirOut." ".$prefixRef." ".$forwardFastqFileIn;		# command line
			toolbox::exportLog("INFOS: tophat::topHat2 : $command\n",1);
		}

		# Command is executed with the run function (package toolbox)
		if (toolbox::run($command)==1)
		{
		# Parse the tophat directory and rename all the files by adding the readgroup and the log subdirectory
		my ($fileName,$readGroup) = pairing::extractName($forwardFastqFileIn);
		my $fileList=toolbox::readDir($tophatDirOut);
		my @fileList=@$fileList;
		## DEBUG toolbox::exportLog("DEBUG: tophat : @fileList\n",1);
		for (my $i=0; $i<=$#fileList;$i ++)
		{
		next if ($fileList[$i] eq '');

		my ($file,$path)=toolbox::extractPath($fileList[$i]);
		$file =~ s/://g;
		##DEBUG toolbox::exportLog("DEBUG: tophat : $fileList[$i] -$file-$readGroup-$fileName- \n",1);
		my $command="mv $tophatDirOut/$file $tophatDirOut/$readGroup".".".$file;   # rename the file by adding read group
		toolbox::run($command);
		last if ($file =~ m/log/);		#rename the lof directory but not the file into log directory
		}

			toolbox::exportLog("INFOS: tophat : correctly done\n",1);
			return 1;
		}
		else
		{
			toolbox::exportLog("ERROR: tophat : ABBORTED\n",0);
			return 0;
		}

	}
	elsif ((toolbox::sizeFile($forwardFastqFileIn)==1) and (toolbox::sizeFile($reverseFastqFileIn)==1) )		##Check if entry files exist and are not empty / paired mode
	{
		if ($gffFile ne "None" and toolbox::sizeFile($gffFile)==1)
		{
			$command=$tophat2.$options." -G ".$gffFile." -o ".$tophatDirOut." ".$prefixRef." ".$forwardFastqFileIn." ".$reverseFastqFileIn;		# command line
			toolbox::exportLog("INFOS: tophat::topHat2 : $command\n",1);
		}
		else
		{
			$command=$tophat2.$options." -o ".$tophatDirOut." ".$prefixRef." ".$forwardFastqFileIn." ".$reverseFastqFileIn;		# command line
			toolbox::exportLog("INFOS: tophat::topHat2 : $command\n",1);
		}

		# Command is executed with the run function (package toolbox)
		if (toolbox::run($command)==1)
		{
		# Parse the tophat directory and rename all the files by adding the readgroup and the log subdirectory
		my ($fileName,$readGroup) = pairing::extractName($forwardFastqFileIn);
		my $fileList=toolbox::readDir($tophatDirOut);
		my @fileList=@$fileList;
		toolbox::exportLog("DEBUG: tophat : @fileList\n",1);
		for (my $i=0; $i<=$#fileList;$i ++)
		{
		next if ($fileList[$i] eq '');

		my ($file,$path)=toolbox::extractPath($fileList[$i]);
		$file =~ s/://g;
		##DEBUG toolbox::exportLog("DEBUG: tophat : $fileList[$i] -$file-$readGroup-$fileName- \n",1);
		my $command="mv $tophatDirOut/$file $tophatDirOut/$readGroup".".".$file;   # rename the file by adding read group
		toolbox::run($command);
		last if ($file =~ m/log/);		#rename the lof directory but not the file into log directory
		}

			toolbox::exportLog("INFOS: tophat : correctly done\n",1);
			return 1;
		}
		else
		{
			toolbox::exportLog("ERROR: tophat : ABBORTED\n",0);
			return 0;
		}
	}
	else
	{
		toolbox::exportLog("ERROR: tophat::tophat2 : Problem with the files\n",0);
		return 0;
	}

}


1;
