package cufflinks;



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
use lib qw(.);
use localConfig;
use toolbox;
use Data::Dumper;

##CUFFLINKS
##Module containing CUFFLINKS functions
##
###
##Adding the tag XS 
##XS tag is founded in bam file  coming from TopHat or Bowtie and not in other bam file such as Bwa
sub addingXStag
{

    my ($file,$out) = @_;
    
    if ($file !~ m/\.bam$/)	##Check if the extension of the file in iput is ".bam"
    {
        warn ("\n$file is not a BAM file in adding XS tag\n");
        return 0; 		###### 30-09 -CD : REVUE CODE pas gestion exportLog?
    }
    else
    {
        ##Check if the bam file is correct
	if (toolbox::checkSamOrBamFormat($file)==1)
	{
	    toolbox::exportLog("The bam file is a real binary file \n",1);
            return 1;
	}
    
        ##catch the name of the mapping tools in the bam file
        my $test=system("samtools view -H $file | grep \@PG");
        if ($test !~m/TopHat$/)
	{
            my $temp = "temp.sam";
            my $xscom = "$samtools view -h $file | awk '{if(\$0 ~ /XS:A:/ || \$1 ~ /^@/) print \$0; else {if(and(\$2,0x10)) print \$0\"\tXS:A:-\"; else print \$0\"\tXS:A:+\";}}' > $temp ";
            system("$xscom") and return 0; 			###### 30-09 -CD : REVUE CODE pas gestion par run?

            my $bamcom = "$samtools view -bS -o $out $temp";
            system("$bamcom") and return 0; 			###### 30-09 -CD : REVUE CODE pas gestion par run?
            return 1;
        }
        
    }
    
}



##Assembly with Cufflinks
sub cufflinks
{
    my($refFasta, $annotationGff, $alignementFile, $optionsHachees)=@_;
    
    #if (toolbox::sizeFile($refFastaFileIn)==1){     ##Check if the reference file exist and is not empty
	my $options=toolbox::extractOptions($optionsHachees, " "); 
        my $command= "$cufflinks/cufflinks"." $options"." -b $refFasta -g $annotationGff $alignementFile";
        
        toolbox::run($command);
        if(toolbox::run($command)==1)
        {
            toolbox::exportLog("Assembly with Cufflinks done \n",1);
            return 1;
        }
        else                                                                                                                                                                                                                    # if one or some previous files doesn't exist or is/are empty or if gatkBaseRecalibrator failed
        {
            toolbox::exportLog("Error when running cufflinks", 0);                                                                                                                                      # returns the error message
        }

}
    


##
##
##Merging with Cuffmerge
sub cuffmerge
{

    my ($refFasta, $annotationGff, $inputDir, $inputFile, $outdir)=@_;

    toolbox::makeDir($outdir);	    #create the outputdir
    
    system  ("ls $inputDir | grep '*.gtf' > list.txt");	#creation the assemblies.txt file ###### 30-09 -CD : REVUE CODE pas gestion par readDir?
    
    open (IN, "list.txt") or toolbox::exportLog("ERROR: cufflinks::cuffmerge : Cannot open the file list.txt $!\n",0); ###### 30-09 -CD : REVUE CODE fichier en dur?
    open (OUT, ">assemblies.txt") or toolbox::exportLog("ERROR: cufflinks::cuffmerge : Cannot open the file assemblies.txt $!\n",0);   ###### 30-09 -CD : REVUE CODE fichier en dur?
    while (<IN>)
    {
	if ($_ =~/^(.*)$/)
	{
	    my $file=$1;
	    print OUT "./$file\n" ;
	}
    }
    
    my $cuffmergecommand =  $cufflinks."/cuffmerge". " -o $outdir  -s $refFasta -g $annotationGff  assemblies.txt";
    #
    #execution of the command
    toolbox::run($cuffmergecommand);
    
    if(toolbox::run($cuffmergecommand)==1)
    {
        toolbox::exportLog("Mering with Cuffmerge done \n",1);
        return 1;
    }
    else                                                                                                                                                                                                                    # if one or some previous files doesn't exist or is/are empty or if gatkBaseRecalibrator failed
    {
        toolbox::exportLog("Error when running Merging files", 0);                                                                                                                                      # returns the error message
    }

}


##
##
##Cuffdiff
#Index database sequences in the FASTA format.
sub cuffdiff
{
    my($refFasta, $mergedFile, $bamFileIn, $annotation, $outdir, $optionsHachees)=@_;
    
    #if (toolbox::sizeFile($refFastaFileIn)==1){     ##Check if the reference file exist and is not empty
        my $options=toolbox::extractOptions($optionsHachees, " "); 
        # my $command=$bwa." index ".$options." ".$refFastaFileIn; ##command
        my $command=$cufflinks."/cuffdiff". $options." -o $outdir -b $refFasta -u $mergedFile" ; 

        toolbox::run($command);
        if(toolbox::run($command)==1)
        {
            toolbox::exportLog("Assembly with Cufflinks done \n",1);
            return 1;
        }
        else                                                                                                                                                                                                                    # if one or some previous files doesn't exist or is/are empty or if gatkBaseRecalibrator failed
        {
            toolbox::exportLog("Error when running cufflinks", 0);                                                                                                                                      # returns the error message
            return 0;
        }
}

    
1;