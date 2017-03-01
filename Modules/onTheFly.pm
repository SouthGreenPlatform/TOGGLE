package onTheFly;

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
use Data::Dumper;
use List::Compare;

use localConfig;
use toolbox;
use bwa;
use samTools;
use picardTools;
use tophat;

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
    my ($hashOrder,$script,$hashCleaner,$hashCompressor)=@_;
    
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
        $catCommand .= " ".$toggle."/onTheFly/previousBlock.txt"; # adding the infos of previous block
        $catCommand .= " ".$toggle."/onTheFly/".$currentSoft."Block.txt"; #Adding the code block for the corresponding step in the cat command, as all txt files with these blocks are name as softBlock.txt
	if ($$hashSoftware{$currentSoft}{'OUT'} eq "NA")
	{# will not add the switcher of previous directory for 'dead end' protgrams such as fastqc, samtools flagstats....
	    $cleanerCounter++;
	    $compressorCounter++;
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
	#Re-initializing the counters for compressing and cleaning data
	$compressorCounter=1;
	$cleanerCounter=1;


	$catCommand .= " ".$toggle."/onTheFly/afterBlock.txt"; # adding infos for next block
    }
    
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
        if ($currentSoft =~ m/bowtie/i or $currentSoft =~ m/tophat/i) #Any step involving BWA
        {
            if ($currentSoft eq "bowtieBuild") # If the index is expressely asked
            {
                my $softParameters = toolbox::extractHashSoft($hashConf,"bowtieBuild");                                  # recovery of specific parameters of bwa index
                tophat::bowtieBuild($reference,$softParameters);  
            }
            else #We check if the index is present or not
            {
                my $refIndexedFile = $reference.".1.ebwt";
                if (-e $refIndexedFile)
                {# The index is already created
                    toolbox::exportLog("INFOS: onTheFly::indexCreator : The reference index for bowtie-built already exists, skipped...\n",1);
                    next;
                }
                my $softParameters = toolbox::extractHashSoft($hashConf,"bowtieBuild");                                  # recovery of specific parameters of bwa index
                tophat::bowtieBuild($reference,$softParameters);  
            }
	    
	    if ($currentSoft eq "bowtie2-Build" or $currentSoft =~ m/tophat/i) # If the index is expressely asked
            {
                my $softParameters = toolbox::extractHashSoft($hashConf,"bowtie2-Build");                                  # recovery of specific parameters of bwa index
                tophat::bowtie2Build($reference,$softParameters);  
            }
            else #We check if the index is present or not
            {
                my $refIndexedFile = $reference.".1.bt2";
                if (-e $refIndexedFile)
                {# The index is already created
                    toolbox::exportLog("INFOS: onTheFly::indexCreator : The reference index for bowtie2-built already exists, skipped...\n",1);
                    next;
                }
                my $softParameters = toolbox::extractHashSoft($hashConf,"bowtie2-Build");                                  # recovery of specific parameters of bwa index
                tophat::bowtie2Build($reference,$softParameters);  
            }
	    
	    
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
    
    my $dotFileOut=$outDir."/togglePipeline.dot"; #Creation of the dot file
    my $graphicFileOut=$outDir."/togglePipeline.png"; #Creation of the figure file in png
    
    my $hashInOut=toolbox::readFileConf("$toggle/softwareFormats.txt"); #We need the format IN/OUT

    #Log info
    toolbox::exportLog("INFOS : $0 : onTheFly::generateGraphviz is creating the graphical view of the current pipeline.\n",1);
    
    open(OUT,">", $dotFileOut) or die "Cannot create $dotFileOut: $0\n";
    my $date = `date`;
    chomp $date;
    print OUT "digraph G {
    \tgraph [fontsize=18,fontname=\"Arial\"]
    \tlabel=\"TOGGLE pipeline generated on the $date\"
    \tnode [shape=box,style=\"rounded,filled\",color=lightblue,width=3,fontname=\"Arial\",fontsize=12]\n";
    
    my ($previousSoft,$input,$output);
    foreach my $step (sort {$a <=> $b} keys %{$hashOrder})
    {
	my $soft=$$hashOrder{$step};
	##DEBUG print $soft,"-->";
	##DEBUG toolbox::exportLog("DEBUG : $0 : onTheFly::generateGraphviz, step = $step, soft = $soft.\n",2);

	#Removing anything after a space. E.g a samtoolsview 1 will become samtoolsView
	$soft =~ s/ .+$//;
	
	$input=$$hashInOut{$soft}{"IN"};
	##DEBUG print $input,"\n";
	$output=$$hashInOut{$soft}{"OUT"};
	$soft=$soft."_".$step;
	unless ($previousSoft) #first line, initiation
	{
	    $previousSoft=$soft;
	    my $inputLine = "\t\"".$input."\" [shape=record,style=\"rounded\",color=gray,width=2];\n";
	    print OUT $inputLine;
	    my $outline = "\t\"".$input."\"->".$soft." [color=gray];\n";
	    print OUT $outline;
	    next;
	}
	if ($output eq "NA")
	{
	    #The soft is a 'dead-end'
	    print OUT "\tnode [shape=ellipse,style=\"rounded,filled\",color=\".7 .3 1.0\",width=2,fontsize=8];\n";
	    my $outline = "\t\"".$previousSoft."\"->".$soft." [style=dotted,weight=1];\n";
	    $outline .= "\t\"".$soft."\"->".$previousSoft." [style=dotted,weight=1];\n";
	    $outline .= "\tnode [shape=box,style=\"rounded,filled\",color=lightblue,width=3,fontsize=12];\n";
	    print OUT $outline;
	    next;
	}
	
	 
	my $outline = "\t\"".$previousSoft."\"->\"".$soft."\" [weight=500];\n";
	print OUT $outline;
	$previousSoft = $soft;
    }
	
    my $trueName=$previousSoft;
    $trueName =~ s/_\d{1,}$//;
    print OUT "\t\"".$$hashInOut{$trueName}{"OUT"}."\" [shape=record,style=\"rounded\",color=gray,width=2];\n";
    my $lastLine="\t\"".$previousSoft."\"->\"".$$hashInOut{$trueName}{"OUT"}."\"[color=gray];\n\t}\n";
    print OUT $lastLine;
    close OUT;
    
    
    my $dotCom="dot -Tpng -o$graphicFileOut $dotFileOut"; #To generate png file
    
    #Verification if dot can work on this installation
    my $dotHelpCommand = `dot -? 2>/dev/null`;
    if ($dotHelpCommand !~ m/Usage: dot/)
    {
	#The dot soft is not installed on this machine
	toolbox::exportLog(">>>>>>>>>>>>>> WARNING : onTheFly::generateGraphviz: Cannot generate graphical view, Graphviz is not installed. Only the dot file has been created. \n The command line to create png file is : $dotCom\n ",1);
	return 1;
    }
    
    toolbox::run("$dotCom");
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
