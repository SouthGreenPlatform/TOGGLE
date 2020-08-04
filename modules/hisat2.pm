package hisat2;

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


################################################################################################
# sub hisat2-build : builds a hisat2  index from a set of DNA sequences.
################################################################################################
# arguments : fasta file to index and options used for hisat2-build running
# Returns prefixname of the database created  (1 if the execution is correctly done else 0)
################################################################################################
sub hisat2Build
{
	my ($refFastaFileIn,$optionsHachees)=@_;
	##DEBUG toolbox::exportLog("DEBUG: hisat2::hisat2-build : $prefixRef\n",1);

	if (toolbox::sizeFile($refFastaFileIn)==1)						# check if the reference file exist and is not empty
	{
		my $options=toolbox::extractOptions($optionsHachees, " ");			# Get given options
		my $command=$hisat2."hisat2-build ".$options." ".$refFastaFileIn." ".$refFastaFileIn;		# command
		##DEBUG toolbox::exportLog("DEBUG: hisat2::hisat2-build : $command\n",1);
		# Execute command
		if(toolbox::run($command)==1)							# The command should be executed correctly (ie return) before exporting the log
	{
			toolbox::exportLog("INFOS: hisat2::hisat2-build : correctly done\n",1);	# hisat2-build have been correctly done
		}
		else
		{
			toolbox::exportLog("ERROR: hisat2::hisat2-build : ABORTED\n",0);		# hisat2-build have not been correctly done
		}
	}
	else
	{
		toolbox::exportLog("ERROR: hisat2::hisat2-build : Problem with the file $refFastaFileIn\n",0);		# hisat2-build can not function because of wrong/missing reference file
	}
	return $refFastaFileIn;
}
################################################################################################
# END sub hisat2-build
################################################################################################



##hisat2 Mapping
sub hisat2
{
    my($samFileOut,$readGroup,$refFastaFileIn,$forwardFastqFile,$reverseFastqFile,$optionsHachees)=@_;

    if (ref $reverseFastqFile) #No mate file, sequences are single ends
    {
        $optionsHachees = $reverseFastqFile;
        $reverseFastqFile = "NA";
    }

    #Options
    my $options = "";
    if ($optionsHachees)
    {
        $options=toolbox::extractOptions($optionsHachees);		##Get given options
    }

    #Command
    if ((toolbox::sizeFile($refFastaFileIn)==1) and (toolbox::sizeFile($forwardFastqFile)==1))		##Check if entry files exist and are not empty
    {

        #Basic command line
        my $command = $hisat2."/hisat2 --rg-id ".$readGroup." --rg SM:".$readGroup." ".$options." -x ".$refFastaFileIn;

        if ($reverseFastqFile ne "NA" && toolbox::sizeFile($reverseFastqFile) == 1) #Mate sequences
        {
            $command .= " -1 ".$forwardFastqFile." -2 ".$reverseFastqFile;
        }
        else
        {
            $command .= " -U ".$forwardFastqFile;
        }
        $command .= " -S ".$samFileOut;

        #Execute command
        if(toolbox::run($command)==1)		## if the command has been executed correctly, export the log
        {
            return 1;
        }
        else
        {
            toolbox::exportLog("ERROR: hisat2::hisat2 : ABORTED\n",0);
            return 0;
        }
    }
    else
    {
        toolbox::exportLog("ERROR: hisat2::hisat2 : Problem with the files $refFastaFileIn or/and $forwardFastqFile\n",0);
        return 0;
    }
}


1;
=head1 NAME

    Package I<hisat2>

=head1 SYNOPSIS

	use hisat2;

    hisat2::hisat2Build

    hisat2::hisat2


=head1 DESCRIPTION

    Package hisat2 ( http:// ) is a software package for mapping low-divergent sequences against a large reference genome, such as the human genome.

=head2 FUNCTIONS

=head3 hisat2::hisat2Build

This module indexes database sequences in the FASTA format.
It takes at least one argument: the name of the database to index

=head3 hisat2::hisat2

This module will map the FASTQ data (single or paired) on the reference

=head1 AUTHORS

Intellectual property belongs to IRD, CIRAD, ADNid and South Green development platform
Written by Christine Tranchant, Julie Orjuela, Sebastien Ravel and Francois Sabot

=head1 SEE ALSO

L<http://www.southgreen.fr/>

=cut
