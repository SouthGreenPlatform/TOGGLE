package onTheFly;

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
use List::Compare;

use localConfig;
use toolbox;
use bwa;
use samTools;
use picardTools;
use tophat;
use bowtie;
use crac;
use stats;

################################################################################################
# sub checkOrder =>  Will verify the order of the softwares in the pipeline
#                                           are Ok (ie outfile 1 in Ok as infile 2)
################################################################################################
# arguments :
# 	- hash of complete configuration
################################################################################################
sub checkOrder
{
	toolbox::exportLog("ERROR: onTheFly::checkOrder : should done at least four arguments\n",0) if (@_ < 4);
    my ($hashConf, $refFastaFile,$gffFile,$keyfile)=@_;
    ##DEBUG print Dumper $hashConf;

    my $hashOrder=toolbox::extractHashSoft($hashConf,"order"); #Picking up the options for the order of the pipeline

    #Picking up input output for each software
    my $hashInOut=toolbox::readFileConf("$toggle/softwareFormats.txt");

	# checking MANDATORY file for each software ( defined in softwareFormats.txt)
	checkMandatory($hashOrder,$hashInOut,$refFastaFile,$gffFile,$keyfile);
    ##DEBUG print Dumper $hashInOut;

    #Verifying the coherence of input/output
    my ($previousSoft,$previousFormat,$currentFormat,$initialStep,$lastStep);
    foreach my $step (sort {$a<=> $b} keys %{$hashOrder})
    {
	my $currentSoft=$$hashOrder{$step};
        $currentSoft =~ s/ \d+$//; # Removing number if a software is used more than once with different options
		$currentSoft =~ s/bamutils.*/bamutilsTool/g; #Rename special for bamutils tools
	##DEBUG print $previousSoft,"->",$currentSoft,"\n";
        ##DEBUG print "**",$$hashInOut{$currentSoft}{"OUT"},"\n";
	#if first round
	if (!defined $previousFormat && $$hashInOut{$currentSoft}{"OUT"} ne "NA")
	{
	    $previousFormat=$$hashInOut{$currentSoft}{"OUT"};
	    $previousSoft=$currentSoft;
            $initialStep=$step unless $initialStep;
            $lastStep = $step;
	    #print "prout\n";
	    next;
	}
	elsif (!defined $previousFormat && $$hashInOut{$currentSoft}{"OUT"} eq "NA")
	{
            $initialStep = $step;
	    ##DEBUG print "pas prout\n";
	    next;
	}

	#For other rounds
	$currentFormat=$$hashInOut{$currentSoft}{"IN"};

	#Comparison IN/OUT
	my @listCurrentFormat = split (",", $currentFormat);
	my @listPreviousFormat = split (",", $previousFormat);

	## DEBUG print "++",@listCurrentFormat,"++",@listPreviousFormat,"\n";

	my $compareList = List::Compare->new(\@listCurrentFormat,\@listPreviousFormat);
	my @intersection = $compareList->get_intersection;
	##DEBUG print "**",@intersection,"\n";
	unless (scalar (@intersection)) #No element are Ok...
	{
	    #Print Error
	    toolbox::exportLog("ERROR: onTheFly::checkOrder : The $previousSoft outputs ($previousFormat) are not compatible with the $currentSoft inputs ($currentFormat).\nPlease check your pipeline order.\n",0);
	}

	#Preparing for the next round

	next if ($$hashInOut{$currentSoft}{"OUT"} eq "NA"); #for a brick such as FastQC which is a 'dead-end'

	$previousSoft=$currentSoft;
	$previousFormat=$$hashInOut{$currentSoft}{"OUT"};
        $lastStep = $step;
        ##DEBUG print $lastStep,"\n";

    }
    ##DEBUG print $initialStep,"--",$lastStep,"\n";
    return ($initialStep,$lastStep); #Will return the last step number
}

################################################################################################
# sub checkMandatory =>  Will verify the argument of the softwares in the pipeline
#                                           are Ok (ie processRadtags needs -k argument)
################################################################################################
# arguments :
# 	- hash of IN, OUT and MANDATORY
# ->- hash order
#$refFastaFile
#$gffFile,
#$keyfile
################################################################################################
sub checkMandatory
{
	toolbox::exportLog("ERROR: onTheFly::checkMandatory : should done at least five arguments\n",0) if (@_ < 5);
	my ($hashOrder,$hashInOut,$refFastaFile,$gffFile,$keyfile)=@_;

	foreach my $step (sort {$a<=> $b} keys %{$hashOrder})
	{
		my $currentSoft=$$hashOrder{$step};
		$currentSoft =~ s/ \d+$//; # Removing number if a software is used more than once with different options
		if (defined($$hashInOut{$currentSoft}{"MANDATORY"}))
		{
			my $paramMandatory = $$hashInOut{$currentSoft}{"MANDATORY"};

			if ($paramMandatory =~ m/reference/)
			{
				if ($refFastaFile eq 'None')
				{
					toolbox::exportLog("ERROR: onTheFly::checkMandatory : $currentSoft needs a reference file. Use -r option in Toggle command line\nExiting...\n",0);
				}
			}
			if ($paramMandatory =~ m/gff/)
			{
				if ($gffFile eq 'None')
				{
					toolbox::exportLog("ERROR: onTheFly::checkMandatory : $currentSoft needs a gff file. Use -g option in Toggle command line\nExiting...\n",0);
				}
			}
			if ($paramMandatory =~ m/keyfile/)
			{
				if ($keyfile eq 'None')
				{
					toolbox::exportLog("ERROR: onTheFly::checkMandatory : $currentSoft needs a keyfile file . Use -k option in Toggle command line\nExiting...\n",0);
				}
			}
		}
	}
}

################################################################################################
# sub generateScript =>  will generate scripts on the fly
################################################################################################
# arguments :
# 	- hash of complete configuration
#   - script name
################################################################################################
sub generateScript
{
    my ($hashOrder,$script,$hashCleaner,$hashCompressor,$hashmerge)=@_;

    #Picking up input output for each software
    my $hashSoftware=toolbox::readFileConf("$toggle/softwareFormats.txt");

    my $catCommand = "cat $toggle/onTheFly/startBlock.txt"; #Will create the preambule for the pipeline code, including paths, use modules, etc...
    my @stepsList = sort{$a <=> $b} keys %{$hashOrder};
    my $cleanerCounter=1; #
    my $compressorCounter=1;#for compressing previous folder
    foreach my $step (@stepsList)
    {
        my $currentSoft=$$hashOrder{$step}; #Picking up the step name
        $currentSoft =~ s/ \d+$//; #Removing numbers if a soft is called more than once
        $currentSoft =~ s/ /_/g; #Removing extraspace
        $currentSoft =~ s/bamutils.*/bamutilsTool/g; #Rename special for bamutils tools


        $catCommand .= " ".$toggle."/onTheFly/previousBlock.txt"; # adding the infos of previous block
        $catCommand .= " ".$toggle."/onTheFly/".$currentSoft."Block.txt"; #Adding the code block for the corresponding step in the cat command, as all txt files with these blocks are name as softBlock.txt
		if ($$hashSoftware{$currentSoft}{'OUT'} eq "NA")
		{# will not add the switcher of previous directory for 'dead end' protgrams such as fastqc, samtools flagstats....
			$cleanerCounter++;
			$compressorCounter++;
			if (defined $$hashmerge{$step} and $step != $stepsList[-1])
			{
				$catCommand .= " ".$toggle."/onTheFly/mergeBlock.txt";
			}
			$catCommand .= " ".$toggle."/onTheFly/afterBlockNa.txt"; # adding infos for next block

			next;
		}
        if (defined $$hashCompressor{$step-$compressorCounter})
		{# The previous step has to be compressed BUT NOT CLEANED
			unless (defined $$hashCleaner{$step-$cleanerCounter})
			{
			$catCommand .= " ".$toggle."/onTheFly/compressorBlock.txt";
			}
		}

		if (defined $$hashCleaner{$step-$cleanerCounter})
		{# The previous step has to be cleaned
			$catCommand .= " ".$toggle."/onTheFly/cleanerBlock.txt";
		}
		if (defined $$hashmerge{$step} and $step != $stepsList[-1])
		{# The previous step has to be cleaned
			$catCommand .= " ".$toggle."/onTheFly/mergeBlock.txt";
		}
		#Re-initializing the counters for compressing and cleaning data
		$compressorCounter=1;
		$cleanerCounter=1;

		$catCommand .= " ".$toggle."/onTheFly/afterBlock.txt"; # adding infos for next block
    }

	
	
	
	###########################################
	###########################################
	############# EN TEST REPORT STAT MAPPING
	$catCommand .= " ".$toggle."/onTheFly/afterBlockNa.txt";
	$catCommand .= " ".$toggle."/onTheFly/statsMappingBlock.txt";
	###########################################
	###########################################
	
	
	
	
	
    $catCommand .= " $toggle/onTheFly/endBlock.txt > $script && chmod 775 $script"; #Adding the end of the script and rending it executable

    ##DEBUG print $catCommand,"\n";

    if(toolbox::run($catCommand,"noprint")==1)       #Execute command
    {
        toolbox::exportLog("INFOS: onTheFly::generateScript : The script $script has been generated\n",1);
        return 1;
    }

    else
    {
        toolbox::exportLog("ERROR: onTheFly::generateScript : The script $script cannot be created using the following command:\n $catCommand\n",0);       # ... and return an error message
        return 0;
    }
}

################################################################################################
# sub indexCreator =>  will create the different index on the reference needed for the analysis
################################################################################################
# arguments :
# 	- hash of complete configuration
#       - reference file
################################################################################################
sub indexCreator
{
    my ($hashConf,$reference)=@_;
    my $hashOrder=toolbox::extractHashSoft($hashConf,"order"); #Picking up the options for the order of the pipeline

    my @listConfig = keys %{$hashConf}; #Picking up all softwares with any option declared

    foreach my $step (sort {$a <=> $b}  keys %{$hashOrder})
    {
        my $currentSoft = $$hashOrder{$step};

        #INDEXING for BWA
        if ($currentSoft =~ m/bwa/i) #Any step involving BWA
        {
            if ($currentSoft eq "bwaIndex") # If the index is expressely asked
            {
                my $softParameters = toolbox::extractHashSoft($hashConf,"bwaIndex");                                  # recovery of specific parameters of bwa index
                bwa::bwaIndex($reference,$softParameters);
            }
            else #We check if the index is present or not
            {
                my $refIndexedFile = $reference.".ann";
                if (-e $refIndexedFile)
                {# The index is already created
                    toolbox::exportLog("INFOS: onTheFly::indexCreator : The reference index for bwa index already exists, skipped...\n",1);
                    next;
                }
                my $softParameters = toolbox::extractHashSoft($hashConf,"bwaIndex");                                  # recovery of specific parameters of bwa index
                bwa::bwaIndex($reference,$softParameters);
            }
        }
        #INDEXING for PICARDTOOLS
        if ($currentSoft eq "picardToolsCreateSequenceDictionary" or $currentSoft =~ m/gatk/i) #Any step involving GATK
        {
            my $dictFileOut=$reference;   # name of dictionary file
            $dictFileOut =~ s/\.[^\.]*$/.dict/;
            if ($currentSoft eq "picardToolsCreateSequenceDictionary")
            {
                my $softParameters = toolbox::extractHashSoft($hashConf,"picardToolsCreateSequenceDictionary"); # recovery of specific parameters of picardToolsCreateSequenceDictionary
                my $rmCommand = "rm -f $dictFileOut";
                toolbox::run($rmCommand); #Removing of any existing previous dictionary
                picardTools::picardToolsCreateSequenceDictionary($reference,$dictFileOut,$softParameters);
            }
            else #We check if the dict is present or not
            {
                if (-e $dictFileOut)
                {# The index is already created
                    toolbox::exportLog("INFOS: onTheFly::indexCreator : The reference index for picardTools CreateSequenceDictionary already exists, skipped...\n",1);

                }
                else
                {
                    my $softParameters = toolbox::extractHashSoft($hashConf,"picardToolsCreateSequenceDictionary");# recovery of specific parameters of picardToolsCreateSequenceDictionary
                    picardTools::picardToolsCreateSequenceDictionary($reference,$dictFileOut,$softParameters);
                }
            }
        }

        #INDEXING for SAMTOOLS
        if ($currentSoft eq "samToolsFaidx" or $currentSoft =~ m/gatk/i) #Any step involving GATK
        {
            if ($currentSoft eq "samToolsFaidx")
            {
                my $softParameters = toolbox::extractHashSoft($hashConf,"samToolsFaidx"); # recovery of specific parameters of samToolsFaidx
                samTools::samToolsFaidx($reference);
            }
            else #We check if the dict is present or not
            {
                my $indexFileOut=$reference.".fai";
                ##DEBUG print $reference,"--",$indexFileOut,"\n";
                if (-e $indexFileOut)
                {# The index is already created
                    toolbox::exportLog("INFOS: onTheFly::indexCreator : The reference index for samtools faidx already exists, skipped...\n",1);
                    ##DEBUG print "skipped faidx\n";
                }
                else
                {
                    ##DEBUG print "samtools faidx\n";
					my $softParameters = toolbox::extractHashSoft($hashConf,"samToolsFaidx");# recovery of specific parameters of samToolsFaidx
                    samTools::samToolsFaidx($reference);
                }
            }
        }

		#INDEXING for topHat
        if ($currentSoft =~ m/bowtie/i or $currentSoft =~ m/tophat/i) #Any step involving bowtie or tophat
        {
            if ($currentSoft eq "bowtieBuild") # If the index is expressely asked
            {
                my $softParameters = toolbox::extractHashSoft($hashConf,"bowtieBuild");                                  # recovery of specific parameters of bowtieBuild index
                bowtie::bowtieBuild($reference,$softParameters);
            }
            else #We check if the index is present or not
            {
                my $refIndexedFile = $reference.".1.ebwt";
                if (-e $refIndexedFile)
                {# The index is already created
                    toolbox::exportLog("INFOS: onTheFly::indexCreator : The reference index for bowtie-built already exists, skipped...\n",1);
                    next;
                }
                my $softParameters = toolbox::extractHashSoft($hashConf,"bowtieBuild");                                  # recovery of specific parameters of bowtieBuild index
                bowtie::bowtieBuild($reference,$softParameters);
            }

			if ($currentSoft eq "bowtie2-Build" or $currentSoft =~ m/tophat/i) # If the index is expressely asked
            {
                my $softParameters = toolbox::extractHashSoft($hashConf,"bowtie2-Build");                                  # recovery of specific parameters of bowtie2Build index
                bowtie::bowtie2Build($reference,$softParameters);
            }
            else #We check if the index is present or not
            {
                my $refIndexedFile = $reference.".1.bt2";
                if (-e $refIndexedFile)
                {# The index is already created
                    toolbox::exportLog("INFOS: onTheFly::indexCreator : The reference index for bowtie2-built already exists, skipped...\n",1);
                    next;
                }
                my $softParameters = toolbox::extractHashSoft($hashConf,"bowtie2-Build");                                  # recovery of specific parameters of bowtie2Build index
                bowtie::bowtie2Build($reference,$softParameters);
            }


        }
		
		 #INDEXING for CRAC
        if ($currentSoft eq "crac" or $currentSoft =~ m/crac/i)  #Any step involving crac
        {
			my $softParameters = toolbox::extractHashSoft($hashConf,"cracIndex"); # recovery of specific parameters of crac
			crac::cracIndex($reference, $softParameters);
        }
    }


    return 1;
}



################################################################################################
# sub generateGraphviz =>  will generate graphical view of the pipeline on the fly
################################################################################################
# arguments :
# 	- hash of complete configuration
#       - output directory
################################################################################################
sub generateGraphviz
{
    my ($hashOrder,$outDir)=@_;

	# filenames generated with extension . dot and .png
    my $dotFileOut=$outDir."/togglePipeline.dot"; #Creation of the dot file
    my $graphicFileOut=$outDir."/togglePipeline.png"; #Creation of the figure file in png

	# get the input and outpu of the sofwares
    my $hashInOut=toolbox::readFileConf("$toggle/softwareFormats.txt"); #We need the format IN/OUT

	# Get date of analysis
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
	$year += 1900;
	my @months = qw( Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec );
	my @days = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
    my $date = "$days[$wday], the $mday $months[$mon], $year";
	
	# Print Log
	toolbox::exportLog("INFOS : $0 : onTheFly::generateGraphviz is creating the graphical view of the current pipeline.\n",1);     

	# Open dot file to generate
    open(OUT,">", $dotFileOut) or die "Cannot create $dotFileOut: $0\n";
 
    # Dot header
    print OUT "digraph G \n{
    \tgraph [fontsize=24,fontname=Arial, colorscheme=X11,  nodesep=1, compound=true;]
	node [label=\"\\N\",  style=\"rounded,filled\", width=4, fontcolor=white, fontsize=24]
    \tlabel=\"\nGenerated by TOGGLe $date\"\n";

	# Generating progressively workflow by reading  hash of configuration file
    my ($previousSoft,$soft, $softLabel, $input,$output);
	my $wkf;
	
    foreach my $step (sort {$a <=> $b} keys %{$hashOrder})
    {
		# soft name formatting
		$soft=$$hashOrder{$step};	
		$soft =~ s/bamutils.*/bamutilsTool/g; #Rename special for bamutils tools
		$soft =~ s/ .+$//; #Removing anything after a space. E.g a samtoolsview 1 will become samtoolsView
			
		$input=$$hashInOut{$soft}{"IN"};
		$output=$$hashInOut{$soft}{"OUT"};
		
		# Get input and output format from sofware Formats. txt file
		$softLabel="$soft ($step)";
		$soft=$soft."_".$step;
		my $color = "darkorange" ;
		$color = "deeppink1" if ( $step>=1000);
		
		if (not defined ($previousSoft)) #first line, input file 
		{
			print OUT "\nFile_input [shape=record, style=rounded, color=dodgerblue2, fontcolor=dodgerblue2, width=2, label=\" $input \" ]";
			$wkf="\n File_input ";
		}
		$wkf.=" -> $soft";
		
		print OUT "$soft [shape=box, color=$color, label=\" $softLabel \"]\n";
		
		$previousSoft = $soft;
    }

    print OUT "File_output [shape=record, style=rounded, color=dodgerblue2, fontcolor=dodgerblue2, width=2, label= \" $output \" ]; \n";    
	$wkf.=" -> File_output [color=dodgerblue3, arrowhead=none];\n\t}\n";
    print OUT $wkf;
    close OUT;

	
	# Dot command to generate png picture
    my $dotCom="dot -Tpng -o$graphicFileOut $dotFileOut"; #To generate png file

    # Check if dot can work on this cluster
    my $dotHelpCommand = `dot -? 2>/dev/null`;
    if ($dotHelpCommand !~ m/Usage: dot/)
    {
		#The dot soft is not installed on this machine
		toolbox::exportLog("WARNING : $0 - onTheFly::generateGraphviz: Cannot generate graphical view, Graphviz is not installed.\n Only the dot file has been created. \n The command line to create png file is : $dotCom\n ",1);
		return 1;
    }

    toolbox::run("$dotCom");
}



###############################################################################################
# sub generateReport =>  will generate lateX reports of workflow and analysis 
################################################################################################
# arguments :
# - in: outdir
################################################################################################
sub generateReports
{
    my ($outDir, $configInfo)=@_;
	my $reportDirWF="$outDir/texReport/workflow/";
	
	# copying texReport files from 
	my $cpCmd="cp $toggle/texReport $outDir -rf";
	toolbox::run($cpCmd, "noprint");
	
	# Report file name generated by TOGGLe to describe workflow architecture
	my $texWorkflowFile=$reportDirWF."/TOGGLe_Workflow_Report.tex";
	my $pdfWorkflowFile=$reportDirWF."/TOGGLe_Workflow_Report.pdf";
	
	
	#copying differents tex files to input repertory
	
	#Pipeline picture generated by graphivz
	my $mvCmd="cp $outDir/togglePipeline.png $reportDirWF/";
	toolbox::run($mvCmd,"noprint");
	
	#parallel sample created by scheduler::schedulerWait
	if ( -e $outDir."/sample_parallel_table.tex")
	{
		$mvCmd="mv $outDir/sample_parallel_table.tex $reportDirWF/input";
		toolbox::run($mvCmd,"noprint");
	}
	# The file has not been created, no parallel analysis (step number > 1000 uniquely) 
	else
	{
		my $echoCmd="echo 'No single sample analyzed, step number >1000 only' >  $reportDirWF/input/sample_parallel_table.tex";
		toolbox::run($echoCmd,"noprint");
	}
	
	#parallel sample
	if ( -e  $outDir."/sample_global_table.tex")
	{
		$mvCmd="mv $outDir/sample_global_table.tex $reportDirWF/input";
		toolbox::run($mvCmd,"noprint");
	}
	# The file has not been created, no global analysis (step number < 1000 uniquely) 
	else
	{
		my $echoCmd="echo 'No global analysis, step number <1000 only' >  $reportDirWF/input/sample_global_table.tex";
		toolbox::run($echoCmd,"noprint");
	}
	
	# stats
	my $statDir = $outDir."/statsReport";

	my @fileStatList = `ls $statDir` or die "ERROR: onTheFly::generateReports : ls $statDir";
    if (scalar (@fileStatList)>0)
	{		
		stats::creatingStatFileTex($statDir);
		$mvCmd="mv $outDir/stats.tex $reportDirWF/input";
		toolbox::run($mvCmd,"noprint");
	}
	
	#Command line
	$mvCmd="mv $outDir/commandLine.tex $reportDirWF/input";
	toolbox::run($mvCmd);
	
	#Software Version
	$mvCmd="mv $outDir/software.txt $reportDirWF/input";
	toolbox::run($mvCmd,"noprint");
	
	#config file
	$cpCmd="cp $configInfo $reportDirWF/input/configuration.txt";
	toolbox::run($cpCmd,"noprint");
	
	#generating pdf report in $texWorkflowFile 	
	my $texCmd="pdflatex $texWorkflowFile";
	
	#generating biblio
	my $auxFile =$texWorkflowFile;
	my @path = split /\//, $auxFile;
	$auxFile = $path[$#path];
	$auxFile =~ s/\.tex$//;

	my $bibtexCmd="bibtex $auxFile";
	## DEBUG toolbox::exportLog($texCmd,1);
    ## DEBUG toolbox::exportLog($bibtexCmd,1);

	# Moving into directory with .tex files. Requested for pdf lateX to avoid loop.
	chdir "$reportDirWF" or die "\nTEX REPORT ERROR: $0 : cannot cd into the output folder $reportDirWF \nExiting...\n";
	
	#toolbox::run($texCmd);
    `$texCmd` or die "\nTEX REPORT ERROR: $0 : cannot generate pdf from tex document $! \nExiting...\n"; 
	`$bibtexCmd` or die "\nTEX REPORT ERROR: $0 : cannot generate bbl from bib document $! \nExiting...\n";
	`$texCmd` or die "\nTEX REPORT ERROR: $0 : cannot generate pdf from tex document $! \nExiting...\n"; 
	`$texCmd` or die "\nTEX REPORT ERROR: $0 : cannot generate pdf from tex document $! \nExiting...\n";
	
	my $lnCmd="ln -s $pdfWorkflowFile $outDir/.";
	system($lnCmd) and die "\nTEX REPORT ERROR: $0 : cannot generate pdf link $! \nExiting...\n"; 
}




1;

=head1 NAME

package I<onTheFly>

=head1 SYNOPSIS

	use onTheFly;



=head1 DESCRIPTION

	This package is used mainly by toggleGenerator.pl that will allow to generate scripts based on the user definition.

=head2 Functions

=over 4

=item checkOrder (Verifies that the order provided by user is correct in term of input/output)

=item generateScript (Assembles the code block into a script based on the user request)

=item indexCreator (Creates all references index/dictionary requested by the given pipeline)

=back

=head1 AUTHORS

Intellectual property belongs to IRD, CIRAD and South Green developpement plateform for all versions also for ADNid for v2 and v3 and INRA for v3
Version 3 written by Cecile Monat, Christine Tranchant, Laura Helou, Abdoulaye Diallo, Julie Orjuela-Bouniol, Sebastien Ravel, Gautier Sarah, and Francois Sabot

=head1 SEE ALSO

L<http://www.southgreen.fr/>		# SOUTH GREEN

=cut
