package bamutils;

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
use Data::Dumper;
use checkFormat;

##############################################
##bamutils
##Module containing bamutils functions
##############################################

##bamutils

sub bamutilsTool
{
     toolbox::exportLog("ERROR: bamutils::bamutilsTool : should done at least four arguments\n",0) if (@_ < 4);
     my($toolname, $bamFileIn,$bamFileOut,$optionsHachees)=@_;
 
     $toolname =~ s/bamutils//g;
     $toolname = lc($toolname);

     if (toolbox::sizeFile($bamFileIn)==1)
     { ##Check if entry file exist and is not empty
          
          #Check if the format is correct
          if (checkFormat::checkFormatSamOrBam($bamFileIn)==0)
          {
               #The file is not a BAM/SAM file
               toolbox::exportLog("ERROR: bamutils::bamutils$toolname : The file $bamFileIn is not a SAM/BAM file\n",0);
               return 0;
          }

          my $options="";
          if ($optionsHachees)
          {
               $options=toolbox::extractOptions($optionsHachees); ##Get given options
          }
          
          #check l'extention
          my $command;
          if ( $bamFileOut =~ m/bam$|sam$/)
          {
               $command=$bamutils." ".$toolname." ".$bamFileIn." ".$bamFileOut." ".$options;
          }
          elsif ( $bamFileOut =~ m/bed$|fastq$|fasta$/)
          {
               $command=$bamutils." ".$toolname." ".$bamFileIn." ".$options." > ".$bamFileOut;
          }

          #toolbox::exportLog($command."\n",1);
          #Execute command
          if(toolbox::run($command)==1)
          {
               return 1;#Command Ok
          }
          else
          {
               toolbox::exportLog("ERROR: bamutils::bamutils$toolname : Uncorrectly done\n",0);
               return 0;#Command not Ok
          }
     }
     else
     {
        toolbox::exportLog("ERROR: bamutils::bamutils$toolname : The file $bamFileIn is uncorrect\n",0);
        return 0;#File not Ok
     }
}

1;

=head1 NAME

    Package I<bamutils> 

=head1 SYNOPSIS

    This package contains the whole modules for the bamutils software

=head1 DESCRIPTION

    Package bamutils (Marcus R. Breese and Yunlong Liu, 2013, http://ngsutils.org/modules/bamutils/ ) is a software package for working on SAM and BAM files (sorting, selection, merging...)

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
