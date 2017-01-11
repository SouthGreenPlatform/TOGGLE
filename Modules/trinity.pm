package trinity;

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


sub trinityRun
{
     my($outputDir,$readGroup,$firstListOfFastq,$secondListOfFastq,$optionsHachees)=@_;  ## if paired reads give $single as --left and $paired as --right; else give only $single
     #print "taille: ", scalar (@{$firstListOfFastq}),"\n";
     if (scalar (@{$firstListOfFastq}) > 0)  # if is not empty
     {
          my $singleFiles=""; my $pairedFiles="";
          my $command="";
          my $options="";
          my $i=0;
          foreach my $file (@{$firstListOfFastq})       # for each Fastq file(s)
          {
               if ($file ne "NA" and toolbox::checkFormatFastq($file)==1 )
               { ##Check if entry file exist and is not empty
                    if (++$i < scalar (@{$firstListOfFastq}))
                    {
                         $singleFiles .= $file."," ;       # recovery of informations fo command line used later
                    }
                    else 
                    {
                         $singleFiles .= "".$file;
                    }
               }
               else
               {
                   toolbox::exportLog("ERROR: trinity::trinityRun : The file $file is uncorrect or is not a Fastq\n",0);
                   return 0;#File not Ok
               }
          }
          #print @{$firstListOfFastq},"\n";
          if ($optionsHachees)
          {
               $options=toolbox::extractOptions($optionsHachees); ##Get given options
          }
          #print $options,"\n";
          if ($options !~ m/--max_memory\s\w+/) # The type of walker is not informed in the options
          {
              $options .= "--max_memory 20G";
          }      
          
          if (scalar (@{$secondListOfFastq}) > 0)  # if is not empty,we have paired reads
          { ##Check if entry file exist and is not empty
              $i=0;
              foreach my $file (@{$secondListOfFastq})       # for each Fastq file(s)
               {
                    if ($file ne "NA" and toolbox::checkFormatFastq($file)==1 )              
                    {
                         if (++$i < scalar (@{$firstListOfFastq}))
                         {
                              $pairedFiles .=$file."," ;       # recovery of informations fo command line used later
                         }
                         else 
                         {
                              $pairedFiles .=$file;
                         }
                    }
                    else
                    {
                        toolbox::exportLog("ERROR: trinity::trinityRun : The file $file is not a Fastq\n",0);
                        return 0;#File not Ok
                    }
               }
               if ($options !~ m/--seqType fq/) # The type of walker is not informed in the options
               {
                   $options .= " --seqType fq";
               }
               $command=$trinity." ".$options." --left ".$singleFiles." --right ".$pairedFiles.' --output '. $outputDir;
               
          }
          else  ## $paired is empty or doesn't exist
          {
              if ($options !~ m/--seqType fq/) # The type of walker is not informed in the options
              {
                  $options .= " --seqType fq";
              }
              $command=$trinity." ".$options." --single ".$singleFiles.' --output '.$outputDir;
              
          }
          #Execute command
          if(toolbox::run($command)==1)
          {
               #my $moveCmd='mv '.$outputDir.'/*Trinity.fasta '.$outputDir.$readGroup; #rename the output file
               chdir "$outputDir";
               #rename the output files
               my $moveCmd="ls | sed -rn \"s/(.*)\$/mv '&' '".$readGroup."_\\1'/ p\" | sh"; 
               toolbox::run($moveCmd,"noprint");
               #remove the sub-repositories to avoid errors during copy in the finalResults repository.     
               my $rmCmd='find . -maxdepth 1 -mindepth 1 -type d -exec rm -r {} \;';
               toolbox::run($rmCmd,"noprint");
               # come back to the working directory.
               chdir "../";
               
               return 1;#Command Ok
          }
          else
          {
               toolbox::exportLog("ERROR: trinity::trinityRun : Uncorrectly done\n",0);
               return 0;#Command not Ok
          }
     }
     else
     {
        toolbox::exportLog("ERROR: trinity::trinityRun : The list of fastq files is empty\n",0);
        return 0;#File not Ok
     }
}
1;