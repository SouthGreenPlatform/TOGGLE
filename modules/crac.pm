package crac;

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
# sub cracIndex : builds a CRAC index from a set of DNA sequences.
################################################################################################
# arguments : fasta file to index and options used for bowtieBuild running
# Returns prefixname of the database created  (1 if the execution is correctly done else 0)
################################################################################################
sub cracIndex
{
	my ($refFastaFileIn,$optionsHachees)=@_;
	##DEBUG toolbox::exportLog("DEBUG: bowtie::bowtieBuild : $prefixRef\n",1);

	if (toolbox::sizeFile($refFastaFileIn)==1)						# check if the reference file exist and is not empty
	{
		my $options=toolbox::extractOptions($optionsHachees, " ");			# Get given options

		my $command=$cracIndex.$options." index ".$refFastaFileIn.".CRAC.index ".$refFastaFileIn;		# command
		##DEBUG toolbox::exportLog("DEBUG: bowtie::bowtieBuild : $command\n",1);
		# Execute command
		if(toolbox::run($command)!=1)							# The command should be executed correctly (ie return) before exporting the log
		{
			toolbox::exportLog("ERROR: crac::cracIndex : ABORTED\n",0);		# crac-index have not been correctly done
		}
	}
	else
	{
		toolbox::exportLog("ERROR: crac::cracIndex : Problem with the file $refFastaFileIn\n",0);		# crac-index can not function because of wrong/missing reference file
	}
	return 1;
}
################################################################################################
# END sub crac-index
################################################################################################



##crac Mapping
sub crac
{
    my($samFileOut,$refIndex,$forwardFastqFile,$reverseFastqFile,$optionsHachees)=@_;
    
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
    
	unless  ($options =~ m/-k /)
	{
		toolbox::exportLog("WARNING: crac::crac: no kmers length provided, we fix it arbitrary at 22.\n",2);
		$options .= " -k 22"; 
	}
	
    #Command
    if (toolbox::sizeFile($forwardFastqFile)==1)		##Check if entry files exist and are not empty
    {
        
        #Basic command line
        my $command = $crac." ".$options." -i ".$refIndex." -o ".$samFileOut." ";
        
        if (toolbox::sizeFile($reverseFastqFile) == 1) #Mate sequences
        {
            $command .= "-r ".$forwardFastqFile." ".$reverseFastqFile;
        }
        else
        {
            $command .= "-r ".$forwardFastqFile;
        }
        
        #Execute command
        if(toolbox::run($command)==1)		## if the command has been executed correctly, export the log
        {
            return 1;
        }
        else
        {
            toolbox::exportLog("ERROR: crac::crac : ABORTED\n",0);
            return 0;
        }
    }
    else
    {
        toolbox::exportLog("ERROR: crac::crac : Problem with the file $forwardFastqFile\n",0);
        return 0;
    }
}



1;
=head1 NAME

    Package I<bowtie>

=head1 SYNOPSIS

	use crac
	
	crac::cracIndex
	
	crac::crac

=head1 DESCRIPTION

    Package CRAC ( http:// ) is a software package for mapping RNAseq data on large genome such as the human genome.

=head2 FUNCTIONS

=head3 crac::cracIndex

This indexes database sequences in the FASTA format.
It takes at least one argument: the name of the database to index

=head3 crac::crac

This module will map the FASTQ data (single or paired) on the reference

=head1 AUTHORS

Intellectual property belongs to IRD, CIRAD, ADNid and South Green development platform
Written by Christine Tranchant, Julie Orjuela, Sebastien Ravel and Francois Sabot

=head1 SEE ALSO

L<http://www.southgreen.fr/>

=cut
