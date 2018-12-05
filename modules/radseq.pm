package radseq;

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
use Data::Dumper;

use lib qw(.);
use localConfig;
use toolbox;
use checkFormat;




################################################################################################
# sub radseq::executeDemultiplexing => to run Demultiplexing with Stacks::process_rad
################################################################################################
# arguments :
#	- keyFile : the keyfile containing barcodes and sample names
#	- initialDir: the directory with fastq files
#	- workingDirRadseq: the output directory with fastq demultiplexed
#	- optionsRadtags: options for radseq::processRadtags
################################################################################################
# return boolean :
#	- 1 if radseq::executeDemultiplexing has runned correctly
################################################################################################

sub executeDemultiplexing
{
	toolbox::exportLog("ERROR: radseq::executeDemultiplexing : should done at least five arguments\n",0) if (@_ < 5);
	my ($keyFile,$initialDir,$workingDirRadseq,$optionsRadtags, $checkFastq)=@_;
	toolbox::exportLog("#########################################\nINFOS: Demultiplexing Step \n#########################################\n",1);

	my $initialDirContent=toolbox::readDir($initialDir);
	
	if ( $checkFastq == 1 )
	{
		foreach my $file (@{$initialDirContent})
		{
			checkFormat::checkFormatFastq($file);
		}
	}

	radseq::processRadtags($keyFile, $initialDir, $workingDirRadseq, $optionsRadtags);	# run demultiplexing

	# Reinitialise the initialDir with output of demultiplexing and add output for toggle
	$initialDirContent=toolbox::readDir($workingDirRadseq);

	#check if data files are not empty because some files can be empty after demultiplex step.
	foreach my $file (@{$initialDirContent})
	{
	  if (toolbox::sizeFile($file) == 0)   #file empty
	  {
		my $rmFileEmpty = "rm $file";
		system ("$rmFileEmpty") and toolbox::exportLog("ERROR: radseq::executeDemultiplexing : cannot remove the file $file\nExiting...\n",0);
		toolbox::exportLog(">>>>>>>>>>>>>>>> WARNINGS: radseq::executeDemultiplexing : Some individuals are empty and were removed: $file\n",2);
	  }
	}
	toolbox::exportLog("#########################################\nINFOS: Demultiplexing Data checking \n#########################################\n",1);

	return 1;
}
#################################################################################################
## END sub radseq::executeDemultiplexing
#################################################################################################



################################################################################################
# sub radseq::parseKeyFile => to parse radseq keyFile in
################################################################################################
# arguments :
# 	- fileIn : the Keyfile file to split
################################################################################################
# return boolean :
#	- 1 if radseq::parseKeyFile has runned correctly
################################################################################################

sub parseKeyFile
{
	toolbox::exportLog("ERROR: radseq::parseKeyFile : should done at least one arguments\n",0) if (@_ < 1);
	my ($keyFile)=@_;

	#open temp keyfile
	my $outputTmpKeyfile ="$keyFile.tmp";
	open(my $TEMPFILE, ">",$outputTmpKeyfile) or toolbox::exportLog("ERROR: radseq::parseKeyFile : Can't open the file $outputTmpKeyfile $!\n",0);

	open(my $KEYFILEOPEN, "<", $keyFile) or toolbox::exportLog ("ERROR: radseq::parseKeyFile : Can't open the file $keyFile $!\n",0); # Reading keyfile input
	while (<$KEYFILEOPEN>)								# for each line of the keyfile
	{
		my $line = $_;								# read line
		chomp $line;
		$line =~ s/\.|_/-/g; 		# modifie "." or "_" by "-"
		print $TEMPFILE $line."\n";
	}
	close $TEMPFILE;
	close $KEYFILEOPEN;

	#mv temp keyfile to replace keyfile
	my $mvTempFile = "mv $outputTmpKeyfile $keyFile";
    system ("$mvTempFile") and die ("\nERROR: radseq::parseKeyFile : cannot rename file $outputTmpKeyfile by $keyFile\nExiting...\n");

	return 1;
}

#################################################################################################
## END sub radseq::parseKeyFile
#################################################################################################



################################################################################################
# sub radseq::processRadtags => to run radseq
################################################################################################
# arguments :
# 	- $keyFile : the keyFile to analyze
#	- $initialDir : the directory that contains all the files fastq
#   - $outDir : output directory with demultiplexed fastq files
#	- $optionsHachees : at least the options -e (enzyme) is necessary; use the option -P if paired-end or others options from process_radtags of stacks
################################################################################################
# return boolean :
#	- 1 if processRadtags runned correctly
#	- else 0
################################################################################################
sub processRadtags
{
	toolbox::exportLog("ERROR: radseq::processRadtags : should get at least four arguments\n",0) if (@_ < 4);
	my ($keyFile,$initialDir,$outDir,$optionsHachees)=@_;

	my $options=toolbox::extractOptions($optionsHachees);		##Get given options
	toolbox::exportLog("ERROR: radseq::processRadtags : you need to specify at least -e process_radtags option of stacks\n",0) if ($options eq "");

	my $pathDir = `dirname $initialDir` or toolbox::exportLog("ERROR: radseq::processRadtags : error in dirname unix command using $initialDir directory\n",0);
	chomp($pathDir);

	#modifie keyfile to change "." or "_" by "-"
	radseq::parseKeyFile($keyFile);

	#my $command = "$stacks/process_radtags --retain_header -p $initialDir -o $outDir -b $keyFile $options "; 		# running radseq
	my $command = "$stacks -p $initialDir -o $outDir -b $keyFile $options "; 		# running radseq


	toolbox::exportLog("INFOS: radseq::processRadtags : $command\n",1);

	if(toolbox::run($command)==1)		## if the command has been excuted correctly, export the log
	{
		my $logBrut = "$outDir/process_radtags.log";
		$command = "mv $logBrut ../";
		toolbox::run($command,"noprint");
	}
	else								## else erreur
	{
	   toolbox::exportLog("ERROR: radseq::processRadtags : ABORTED\n",0);
	}

	return 1;
}
################################################################################################
# END sub radseq::processRadtags
################################################################################################


################################################################################################
# sub radseq::checkOrder => to check if radseq is first step
################################################################################################
# arguments :
# 	- $outputDir : outputDir pass to toggle
#	- $fileConf: toogle file config
#	- $initialDir : initialDir
#	- $checkFastq : arg value of toggle -nocheckFastq option
################################################################################################
# return  :
#	- $initialDirContent: udpate -d argument with path of demultiplexed files generate by stacks
################################################################################################

sub checkOrder
{
	toolbox::exportLog("ERROR: radseq::checkOrder : should get at least five arguments\n",0) if (@_ < 5);
	my ($outputDir,$fileConf,$initialDir,$checkFastq, $keyfile)=@_;

	my $configInfo=toolbox::readFileConf($fileConf);
	my $optionsRadtags=toolbox::extractHashSoft($configInfo,"processRadtags"); 		# Picking up the options for radseq::processRadtags
	my $hashOrder=toolbox::extractHashSoft($configInfo,"order");					#Picking up the options for the order of the pipeline

	my $resultsDirRadseq = "demultiplexedFiles";
	my $workingDirRadseq = $outputDir."/$resultsDirRadseq";
	
	my @firstKeys = sort {$a <=> $b} (keys(%{$hashOrder}));
	if ($$hashOrder{$firstKeys[0]} ne "processRadtags")			#verify if processRadtags is the first step	(if not step 1 error)
	{
		toolbox::exportLog("ERROR  radseq::executeDemultiplexing : demultiplexing must been the first step,\n",0);
	}
	else
	{
		toolbox::makeDir($workingDirRadseq);
		radseq::executeDemultiplexing($keyfile,$initialDir,$workingDirRadseq,$optionsRadtags,$checkFastq);
	}


	# Reinitialise the initialDir with output of demultiplexing and add output for toggle
	my $initialDirContent=toolbox::readDir($workingDirRadseq);

	return $initialDirContent
}
################################################################################################
# END sub radseq::checkOrder
################################################################################################


1;


=head1 NOM

package I<radseq>

=head1 SYNOPSIS

	use radseq;

	radseq::executeDemultiplexing($keyFile,$initialDir,$workingDirRadseq,$optionsRadtags,$checkFastq);

	radseq::parseKeyFile($keyFile);

	radseq::processRadtags($keyFile,$initialDir,$outDir,$optionsHachees);
	
	radseq::checkOrder($outputDir,$fileConf,$initialDir,$checkFastq,%param);


#=head1 DESCRIPTION

This module is a set of functions related to radseq module of Stacks software,  L<http://catchenlab.life.illinois.edu/stacks/>


=head2 Functions


=head3 processRadtags()

This function execute the Stacks::process_radtags software to demultiplex a directory that contains fastq files.
Four arguments are required : a key file, a input directory, a output directory, option of process_radtags.

C<stacks::processRadtags($keyFile,$initialDir,$outDir,$optionsHachees);>

Return 1 if processRadtags runned correctly else 0.

=head3 parseKeyFile()

This function analyze the keyfile to modifie and change "." or "_" by "-" in order to avoid error from Stack or TOGGLE.

One argument is required: the keyfile

C<radseq::parseKeyFile($keyFile);>

=head3 executeDemultiplexing()

This function run stacks::processRadtags and check if data files are not empty because some files can be empty after demultiplex step.

Return 1 if radseq::executeDemultiplexing has ran correctly else 0.

Example :
C<radseq::executeDemultiplexing($keyFile,$initialDir,$workingDirRadseq,$optionsRadtags);>

=head3 checkOrder()

This function check if processRadtags is the first step.

Return initialDirContent: udpate -d argument with path of demultiplexed files generate by stacks.

Example :
C<radseq::checkOrder($outputDir,$fileConf,$initialDir,$checkFastq,%param);>


=head1 AUTHORS

Cecile Monat, Ayite Kougbeadjo, Marilyne Summo, Cedric Farcy, Mawusse Agbessi, Julie Orjuela, Sebastien Ravel Christine Tranchant and Francois Sabot

L<http://www.southgreen.fr/>

=cut
