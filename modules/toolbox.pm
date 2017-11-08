package toolbox;

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
use Exporter;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

use localConfig;
use namingConvention;
use checkFormat;

#Global infos
#our @ISA=qw(Exporter);
#our @EXPORT=qw($configInfos);
#my $configInfos; #Config informations, ref of an hash



##############################################
#Global functions
##############################################



################################################################################################
# sub toolbox::exportLog => log printing
################################################################################################
# arguments :
# 	- logLines :  lines to print in log file
#	- Boolean :  control the type of log. If 0, print logs into error log file and die.
# if 1, print logs into log file. If 2, print logs into error log file and the program continue.
################################################################################################
# no value returned
################################################################################################
sub exportLog
{

    my ($logLines,$controlValue)=@_;
	my $indivName = "TOGGLE";
	my $currentSoft = "toggleGenerator";
	# test if directory contain individuSoft.txt
	if ((-e "individuSoft.txt"))
	{
		# Get the name of indivudal analysed when this function is called
		$indivName = `head -n 1 individuSoft.txt`;
		chomp $indivName;

		# Get the name of software executed when this function is called
		$currentSoft = `tail -n 1 individuSoft.txt`;
		chomp $currentSoft;
	}

    # Initialization of the log name from the name of the individual and the name of software
    my $logBasicName=$indivName."_".$currentSoft."_log";
	$logBasicName = relativeToAbsolutePath($logBasicName,0);
    # Opening of the log file and erro file
    open (OUT, '>>', $logBasicName.".o") or toolbox::exportLog("Cannot create $logBasicName.o file: $!",0);
    open (ERR, '>>', $logBasicName.".e") or toolbox::exportLog("Cannot create $logBasicName.e file: $!",0);

    # According to the value of the $controlValue, the log message is written either in log file or in log file and error file
    if ($controlValue eq "0")		#Something wrong and die
    {
	print OUT "ERROR: toolbox::exportLog : Look at $logBasicName.e for more infos\n";
        print ERR $logLines;
	print OUT "ANALYSIS ABORTED\n";
	die("ANALYSIS ABORTED... Look at $logBasicName.e\n");   #toolbox::exportLog("ANALYSIS ABORTED... Look at $logBasicName.e\n",0);
    }
    elsif ($controlValue eq "1")	#Everything is Ok
    {
	print OUT $logLines."\n";
    }
    elsif ($controlValue eq "2")	#Something wrong but continue with warning
    {
	print OUT "WARNING: toolbox::exportLog : Look at $logBasicName.e for more infos\n";
	print ERR $logLines."\n";
    }

    close OUT;
    close ERR;
}
################################################################################################
# END sub toolbox::exportLog
################################################################################################




##############################################
#Functions related to file, generic
##############################################


################################################################################################
# sub toolbox::checkFile => Check if a file exists, is readable/writable, and has something in
################################################################################################
# arguments :
# 	- file :  name of the file to check
################################################################################################
# A boolean value is returned. 1 if the file exist, if it isn't empty and if it can be read and
# written. else 0
################################################################################################
sub checkFile
{
    my ($file)=@_;

    #Check existence
    my $existence = existsFile($file);
    if ($existence == 0)		#File does not exist
    {
        return 0;
    }

    #Check size
    my $size = sizeFile($file);
    if ($size == 0)		#Empty file
    {
        my $infoSize ="ERROR: toolbox::checkFile : $file empty\n";
        return $infoSize;
    }

    #Check read and write right
    my $readingRights = readFile($file);

    my $logOut;


    if($readingRights == 0)
    {
        $logOut="ERROR: toolbox::checkFile : The file $file is not readable\n";
        exportLog($logOut,0);
    }

    return 1;
}
################################################################################################
# END sub toolbox::checkFile
################################################################################################



################################################################################################
# sub toolbox::readFile => Check if a file is readable
################################################################################################
# arguments :
# 	- file :  name of the file to check
################################################################################################
# A boolean value is returned. 1 if the file is readable. else 0.
################################################################################################
sub readFile
{
    my ($file)=@_;
    existsFile($file); 		#Check if exists

    if (-r $file){return 1;}  	#Check if readable
    else {return 0;}
}
################################################################################################
# END sub toolbox::readFile
################################################################################################




################################################################################################
# sub toolbox::writeFile => Check if a file is writable
################################################################################################
# arguments :
# 	- file :  name of the file to check
################################################################################################
# A boolean value is returned. 1 if the file is writable. else 0.
################################################################################################
sub writeFile
{
    my ($file)=@_;
    existsFile($file); 		#Check if exists
    if (-w $file){return 1;}	#Check if writable
    else {return 0};
}
################################################################################################
# END sub toolbox::writeFile
################################################################################################




################################################################################################
# sub toolbox::sizeFile => Check if a file is empty
################################################################################################
# arguments :
# 	- file :  name of the file to check
################################################################################################
# A boolean value is returned. 1 if the file isn't empty. else 0.
################################################################################################
sub sizeFile
{
    my ($file) = @_;
    existsFile($file); 		#Check if exists
    if (-s $file) {return 1;}	#Check if the file is more than 0 bytes
    else {return 0;}    	#File does not exists or has a size of 0
}
################################################################################################
# END sub toolbox::sizeFile
################################################################################################




################################################################################################
# sub toolbox::existsFile => Check if a file exists
################################################################################################
# arguments :
# 	- file :  name of the file to check
#	- boolean : if the file doesn't exist and the boolean parameter is equal to 0, print just a log
# if the file doesn't exist and the boolean parameter is equal to 1, print a log/error and die (via toolbox::exportLog)
################################################################################################
# A boolean value is returned. 1 if the file exists. else 0.
################################################################################################
sub existsFile
{
    my ($file, $boolean)=@_;
    $boolean = 1 if (not defined $boolean);

    if ((-e $file))	#file exists
    {return 1;}
    elsif ($boolean)	#file does not exists and boolean==1 means error + die
    {

        my $logOut ="ERROR: toolbox::existsFile : The file $file does not exist or is not a file\n$!\n";
        exportLog($logOut,0); 	#report an error to exportLog and die
    }
    else		#file does not exists and boolean==0 means log
    {
        ##DEBUG my $logOut ="INFOS: toolbox::existsFile : The file $file does not exist or is not a file\n$!\n";
        ##DEBUG exportLog($logOut,1);
	return 0;
    }
}
################################################################################################
# END sub toolbox::existsFile
################################################################################################




################################################################################################
# sub toolbox::existsDir => Check if a directory exists
################################################################################################
# arguments :
# 	- name of the directory to check
#	- boolean : if the directory doesn't exist and the boolean parameter is equal to 0, print just a log
# if the file doesn't exist and the boolean parameter is equal to 1, print a log/error and die (via toolbox::exportLog)
################################################################################################
# A boolean value is returned. 1 if the directory exists. else 0.
################################################################################################
sub existsDir
{
    my ($dir,$boolean)=@_;
    $boolean = 1 if (not defined $boolean);

    if (-e $dir and -d $dir) { return 1; } 	#directory exists
    elsif ($boolean) 				#directory does not exists and boolean==1 means error + die
    {
        toolbox::exportLog("ERROR: toolbox::existsDir : The directory $dir does not exist\n$!\n",0);
    }
    else					#directory does not exists and boolean==0 means log
    {
        toolbox::exportLog("INFOS: toolbox::existsDir : The directory $dir does not exist or is not a file\n$!\n",1);
	return 0;
    }
}
################################################################################################
# END sub toolbox::existsDir
################################################################################################




################################################################################################
# sub toolbox::makeDir => to create a directory
################################################################################################
# arguments :
# 	- name of the directory to check
#	- boolean (optional, 0 by default): if the boolean parameter is equal to 1 and
# the directory already exists, the existing directory is removed before creatong again
#
# A boolean value is returned. 1 if the directory has been created. else 0.
################################################################################################
sub makeDir
{
    toolbox::exportLog("ERROR: toolbox::makeDir : should get one argument at least\n",0) if (@_ < 1 );

    my ($dir, $erase)=@_;

    $erase = 0 if (not defined $erase);

    # if $erase equal to 1 and the directory exists, the directory is removed
    run("rm -rf $dir") if ($erase and existsDir($dir));

    # the directory is created unless there is a problem. Then a error is generated and managed
    # with the function exportLog
    unless (mkdir $dir)
    {
	my $logOut = "ERROR: toolbox::makeDir : Unable to create the directory $dir \n$!\n";
        toolbox::exportLog($logOut,0);
    }

    return 1;

}
################################################################################################
# END sub toolbox::makeDir
################################################################################################




################################################################################################
# sub toolbox::readDir => to read the content of a directory
################################################################################################
# arguments :
# 	- name of the directory to parse
#	- part of the filename (optional): list only the file that end by the word given (part of the filename)
#
# The list of files (table) is returned.
################################################################################################
sub readDir
{
    toolbox::exportLog("ERROR: toolbox::readDir : should get at one argument at least\n",0) if (@_ < 1 );

    my ($dir)= shift @_;

    # if the part of name is given in the second argument, search for files directory/*part_of_name
    # else search for all files directory/*
    my $path = defined ($_[0]) ? $dir.'/*'.$_[0] : $dir."/*";

    # ls run
    my $file=`ls $path` or toolbox::exportLog("ERROR: toolbox::readDir : Can't open the directory $path\n$!\n",0);
    chomp $file;

    # split the list into a table returned after
    my @fileList = split /\n/, $file; #print Dumper(\@fileList);


    #Validation of the complete path provided, to avoid the 'ls *' bug with a single subfolder
    foreach my $files (@fileList)
	{
	##DEBUG print $files," --> ";
	#Check if there is at least a '/' in the name
	if ($files =~ m/\//)
	    { #The file has the full path
	    ##DEBUG print "No changes to $files\n";
	    last; # all must be with full path, so let's leave the correction
	    }
	#The full path is not implemented
	my $subfolder = `ls $dir`;
	chomp $subfolder;
	$files = $dir."/".$subfolder."/".$files; #Adding the complete path
	##DEBUG print $files,"\n";
	}

    return(\@fileList);

}
################################################################################################
# sub toolbox::readDir
################################################################################################




################################################################################################
# sub toolbox::readDir2 => to read the content of a directory
################################################################################################
# arguments :
# 	- name of the directory to parse
#	- part of the filename (optional): list only the file that begin by the word given (part of the filename)
#
# The list of files (table) is returned.
################################################################################################
#sub readDir2
#{
#    toolbox::exportLog("ERROR: toolbox::readDir2 : should get at one argument at least\n",0) if (@_ < 1 );
#
#    my ($dir)= shift @_;
#
#    # if the part of name is given in the second argument, search for files directory/part_of_name*
#    # else search for all files directory/*
#    my $path = defined ($_[0]) ? $dir.'/'.$_[0].'*' : $dir."/*";
#
#    # ls command
#    my $file=`ls $path` or toolbox::exportLog("ERROR: toolbox::readDir : Can't open the directory $path\n$!\n",0);
#    chomp $file;
#
#    # split the list into a table returned after
#    my @fileList = split /\n/, $file;
#    return(\@fileList);
#
#}
#########################################
# END sub toolbox::readDir2
#########################################




################################################################################################
# sub toolbox::extractPath => to extract separately the directory name from the filename
################################################################################################
# arguments :
# 	- pathfile
#
# return separately the filename and the path
################################################################################################
sub extractPath
{

  my ($path) = shift @_;

  # split the filepath into a table
  # the last element only contains the name of the file
  my @pathTab = split /\//, $path;
  my $file = $pathTab[$#pathTab];

  # remove the filename in the path given by argument
  $path =~ s/$file//;

  # return 2 parameters : the filename, the path
  return ($file,$path);
}
################################################################################################
# END sub toolbox::extractPath
################################################################################################






##############################################
#Function related to a conf file
###############################################



################################################################################################
# sub toolbox::readFileConf => Read the Software Configuration File and export the value in a hash of hashes
# ex:
# configuration file
#       ----------------------
#	| $samtools view pair
#	| -F=0x02
#	----------------------
# Hash returned. Ex: $configInfos->{"samtools view pair"}{-F} -> '0x02'
################################################################################################
# arguments :
# 	- name of the software configuration file directory
#
# Returns a hash of hashes
# ex : $configInfos->{$currentProgram}{$optionName}=$optionValue;
################################################################################################
sub readFileConf
{
    my($fileConf) = @_;
    my $configInfos;
    #chech if the configuration file is readable
    readFile($fileConf) or exportLog("ERROR: toolbox::readFileConf : Cannot read the configuration file $fileConf\n$!\n",0);

    #store the file content into the table @lines
    open (FILE,"<",$fileConf) or toolbox::exportLog("ERROR: toolbox::readFileConf  : Cannot open the file $$fileConf\n$!\n",0);;
    my @lines = <FILE>;
    close FILE;

    #Generating the hash with all infos
    my $currentProgram;		#Need to be declared outside the loop

    while (@lines)		#Reading all lines
    {
        my $currentLine=shift @lines;
        chomp $currentLine;

        #Avoided lines
        next if $currentLine =~ m/^$/;#Empty line
        next if $currentLine =~ m/^#/;#Commented line

        if ($currentLine =~ m/^\$/)		#New program to be configured, line starting by $
        {
            $currentProgram=$currentLine;
            $currentProgram =~ s/\$//;#remove the "$" in front of the name

            if ($currentProgram =~ m/sge/i) #even if no options are provided, we must see the SGE key
            {
                $configInfos->{$currentProgram}{' '}=' ';
            }

        }
        else		#Config lines
        {
            if($currentLine =~/^no_option$/)
            {
            	$configInfos->{$currentProgram}{' '}=' '; ## add a key for programme without options to be considered
            	next;
            }
            my($optionName,$optionValue)=split /=/,$currentLine, 2;
            $optionValue = "NA" unless $optionValue; # For option without value, eg -v for verbose
            #Populating the hash
            $configInfos->{$currentProgram}{$optionName}=$optionValue;
        }
    }

    $configInfos=namingConvention::softwareNomenclature($configInfos); # will correct the configInfo names

    return $configInfos;

}
################################################################################################
# END sub toolbox::readFileConf
################################################################################################





################################################################################################
# sub toolbox::extractOptions => for extracting options (parameters that'll be given to the tool
# commands)
################################################################################################
# arguments :
# 	- a hash of hash (reference) with the software name as key, then the option as key and
#         the value as element
#       - separator used to separate an option from its value (\t, =)
#	- concatenator used to concatenate each option/value couples
#
# ex : my $optionLine=toolbox::extractOptions($configInfos->{"BWA aln"}," "," ");
#
# Returns the list of options for a software given.
################################################################################################
sub extractOptions
{
    ##Getting the two parameters, the options hash and the option-value separator
    my($optionsHashees,$separateur,$concatenator)=@_;

    # The concatenation of option is classically a space (eg -b 1 -c2), but sometimes we may need to put them one per line (eg
    # export TOTO=$TOTO:/my/new/path
    # export PATH=$PATH:/my/second/path)
    $concatenator = " " unless $concatenator;

    #The separator between option and its value is generally a space but can also be an equal (as in INPUT=/my/file for picardTools) sign.
    $separateur=" " unless $separateur;

    if ($optionsHashees)			# the option hash isn't empty/is defined => options are extracted
    {
	my %options=%{$optionsHashees};
	my $option=" ";
	 	## if no separator is given, set it as a single space
	try
	{
	    foreach my $cle (keys %options)
	    {
		if ($options{$cle} eq 'NA') 	## The option has no value  => no print $options{$cle}
		{
		    $option=$option.$cle.$separateur.$concatenator;
		}
		else				## The option has value  => print $options{$cle}
		{
		    $option=$option.$cle.$separateur.$options{$cle}.$concatenator;

		    #if option is of type -l mem_free=10G
		    $option=~ s/ =/=/g;
		}
	    }
	    return $option;
	}
	catch
	{
	    exportLog("ERROR: toolbox::extractOptions : $_",0);
	}
    }
    else
    {
	return ""; # the option hash is empty/isn't defined => Empty string is returned # REVIEW CODE quel soft???
    }
}
################################################################################################
# END sub toolbox::extractOptions
################################################################################################






################################################################################################
# sub toolbox::extractName => for extracting Name (remove the extention and return the name)
################################################################################################
# arguments : complete path of a file (path+name of the file)
# Returns just the filename without the path
################################################################################################
sub extractName
{
    my $bruteName=shift @_;	# get the complete path of the file

    $bruteName=`basename $bruteName` or toolbox::exportLog("ERROR : $0 : toolbox::extractName cannot work on $bruteName name\n",0); # get just the name of the file after the last /
    chomp $bruteName;
					    ##### REVUE CODE CD 18-09 : et s'il y a un autre point dans le nom de fichier??? ou pas de point
    ##DEBUG print Dumper($bruteName);
    return $bruteName;
}
################################################################################################
# END sub toolbox::extractName
################################################################################################







################################################################################################
# sub toolbox::run => executes the command you gave it and exports the log via toolbox::exportLog() function
################################################################################################
# arguments : command to execute
# Returns boolean (1 if the command has correctly worked else 0)
################################################################################################
sub run
{
    use Capture::Tiny qw(capture);

    my($command,$print)=@_;
    exportLog("INFOS: toolbox::run : $command\n",1) if (not defined $print);

    ##Execute the command
    my ($result,$stderr)=capture {system( $command)};

    ##Log export according to the error
    if ($?==0)
    {
		##DEBUG exportLog("INFOS: toolbox::run : ".$result."\n--".$stderr."\n",1);
		return 1;
    }
    else
    {
		exportLog("ERROR: toolbox::run : ".$command."\n--STDERR: ".$stderr."\n--STDOUT: ".$result."\n",0);
		return 0;
    }

		##########################################################################
		## CD 18-09: ADD LOG
		##########################################################################

    # restore STDOUT and STDERR
    #open (STDOUT, '>&', $STDOUT_);
    #open (STDERR, '>&', $STDERR_);

}

################################################################################################
# END sub toolbox::run
################################################################################################


#########################################
#checkNumberLines
#check the sequences number using the wc -l command from bash and return the number of lines
#########################################

sub checkNumberLines
{
    my ($fileName)=@_;		# recovery of informations
    #my $nbLineCommand="wc -l ".$fileName; #command to count the line number
    my $nbLine=4000;
    #my $nbLine = `$nbLineCommand` or exportLog("ERROR: toolbox::checkNumberLines : Cannot run $nbLineCommand\n$!\n",0);	# execution of the command or if not possible, return an error message
    chomp $nbLine;

    #Add a split to only keep the number of line without the file name
    my @splitLine = split (" ", $nbLine);
    $nbLine = $splitLine[0];

    return $nbLine;  # return the number of line of file
}









#########################################
# Experimental feature: adding tag infos in a BAM file
#########################################
#Highly experimental for now!!
################################################################################################
# sub addInfoHeader => adding tag infos in a BAM file (@CO field)
################################################################################################
# arguments : filename to analyze, text to add
# Returns boolean (1 if this function has been correctly running else 0)
################################################################################################
sub addInfoHeader
{
    my ($bamFile,$textToAdd)=@_;

    #Is the file sam of bam ?
    my $formatCheck=checkFormat::checkFormatSamOrBam($bamFile);

    if ($formatCheck == 0)	#The file is not a BAM nor a SAM and cannot be treated here
    {
	toolbox::exportLog("ERROR: toolbox::addInfoHeader : The file $bamFile is not a SAM/BAM file\n",0);
	return 0;
    }
    #Not sure it is requested, as checkSamOrBamFormat will already kill the job...
    #TODO: using a warning in the checkSamOrBamFormat and a die here, to make the differences in the logs ?

    #Picking up current header
    my $headerCommand="$samtools view -H $bamFile > headerFile.txt"; #extract the header and put it in the headerFile.txt
    run($headerCommand);

    #Adding the @CO field to the header
    my $addingCoLineCommand = "echo \"\@CO\t$textToAdd\" | cat - >> headerFile.txt";#Adding the text at the end of the file.
    run($addingCoLineCommand);

    #reheading the sam/bam file
    my $reheadedFile = $bamFile.".reheaded.bam";
    my $reheaderCommand = "$samtools reheader headerFile.txt $bamFile > $reheadedFile";
    run($reheaderCommand);

    #copy reheaded file to original bam file
    my $copyCommand = "cp $reheadedFile $bamFile";
    run($copyCommand);

    #returning if OK
		##########################################################################
    return 1    ### CD 18-09 : Return toujours 1
}		##########################################################################
################################################################################################
# END sub addInfoHeader
################################################################################################











################################################################################################
# sub changeDirectoryArbo => change directory into the arborescence that we've defined
################################################################################################
# arguments :
#       - Directory that will contain all the directory created by the pipeline
#       - NumberOfDirectory to define what is the directory name that will be returned
# Returns the directory name that will be created according to the number given as argument
################################################################################################
sub changeDirectoryArbo
{
    my ($inputDir,$numberOfDirectory) = @_;
    # numberOfDirectory could be:
    # 0  for /0_PAIRING_FILES/
    # 1  for /1_FASTQC/
    # 11 for /11_FASTXTRIMMER/
    # 2  for /2_CUTADAPT/
    # 3  for /3_PAIRING_SEQUENCES/
    # 4  for /4_MAPPING/
    # 5  for /5_PICARDTOOLS/
    # 6  for /6_SAMTOOLS/
    # 7  for /7_GATK/

    #TODO can be factorized using a HASH I think

    my $finalDirectory;	# output directory returned by the fonction

    if ($numberOfDirectory == 0)		# if you want to move to 0_PAIRING_FILES/
	{
	    $finalDirectory = "$inputDir"."0_PAIRING_FILES";
	}
	elsif ($numberOfDirectory == 1)		# if you want to move to 1_FASTQC/
	{
	    $finalDirectory = "$inputDir"."1_FASTQC";
	}
	elsif ($numberOfDirectory == 11)		# if you want to move to 1_FASTQC/
	{
	    $finalDirectory = "$inputDir"."11_FASTXTRIMMER";
	}
	elsif ($numberOfDirectory == 2)		# if you want to move to 2_CUTADAPT/
	{
	    $finalDirectory = "$inputDir"."2_CUTADAPT";
	}
	elsif ($numberOfDirectory == 3)		# if you want to move to 3_PAIRING_SEQUENCES/
	{
	    $finalDirectory = "$inputDir"."3_PAIRING_SEQUENCES";
	}
	elsif ($numberOfDirectory == 4)		# if you want to move to 4_BWA/
	{
	    $finalDirectory = "$inputDir"."4_MAPPING";
	}
	elsif ($numberOfDirectory == 5)		# if you want to move to 5_PICARDTOOLS/
	{
	    $finalDirectory = "$inputDir"."5_PICARDTOOLS";
	}
	elsif ($numberOfDirectory == 6)		# if you want to move to 6_SAMTOOLS/
	{
	    $finalDirectory = "$inputDir"."6_SAMTOOLS";
	}
	elsif ($numberOfDirectory == 7)		# if you want to move to 7_GATK/
	{
	    $finalDirectory = "$inputDir"."7_GATK";
	}

    return $finalDirectory;
}
################################################################################################
# END sub changeDirectoryArbo
################################################################################################




################################################################################################
# sub extractHashSoft => from output of readFileConf, give the sub-hash for the soft you want to use
################################################################################################
# arguments :
#       - reference oh the hash that contains all softawre name and their parameters
#       - software name
# Returns the sub hash for the software name given as parameter
################################################################################################
sub extractHashSoft
{
    my ($optionref,$soft) = @_;
    my %optionsRef = %$optionref; #WHY ??? because $optionref is a ref but not a hash.

    my $softInfos;
    foreach my $option (keys %optionsRef)	# for each options in the options hash
    {
	if ($option eq $soft)			# if the current option is the soft the user wants
	{
	    foreach my $parameter (keys %{($optionsRef{$option})})			# recovery of parameter for this soft
	    {
		##DEBUG toolbox::exportLog("DEBUG : toolbox::extractHashSoft: $parameter $optionsRef{$option}{$parameter}\n",1); # print in the configuration file the parameter and the options cooresponding to
		$softInfos->{$parameter}=$optionsRef{$option}{$parameter};		# recovery into a reference for the soft
	    }
	}
    }
    ##DEBUG   toolbox::exportLog("DEBUG : toolbox::extractHashSoft: @$softInfos\n",1);
    return $softInfos;
}


################################################################################################
# END sub extractHashSoft
################################################################################################






################################################################################################
# sub checkInitialDirContent => from initial directory determine if contains onlyfiles or folder with files
################################################################################################
# arguments : directory name to analyze
# Returns 0 if the directory contains only files or the list of folders
################################################################################################
sub checkInitialDirContent
{
    my ($initialDir) = @_;			# recovery of initial directory

    toolbox::existsDir($initialDir);		# check the directory exists
    my $list = toolbox::readDir($initialDir);	# read directory
    my @list = @$list;				# recovery directory contents
    my $arbo = 0;				# count for folder
    my @listOfFolder;
    for (my $i=0; $i<=$#list; $i++)		# for each file/folder contain in the directory
    {
	if ($list[$i]=~m/(.+):$/)		# there is folder(s)
	{
	    $arbo++;
	    my $folder = $1;
	    push (@listOfFolder, $folder);
	}
    }
    if ($arbo == 0)		# return 0 if only files in directory
    {
	toolbox::exportLog("INFOS: toolbox::checkInitialDirContent : The directory $initialDir contain only files\n",1);
	return 0;
    }
    else			# return the list of folders in the directory
    {
	toolbox::exportLog("INFOS: toolbox::checkInitialDirContent : The directory $initialDir contain $arbo folder(s)\n",1);
	return (\@listOfFolder);
    }
}
################################################################################################
# END sub checkInitialDirContent
################################################################################################







################################################################################################
# sub transferDirectoryFromMasterToNode =>  to copy data directory from a master server
#                                           to one local scratch space of the node
################################################################################################
# arguments :
# 	- the directory name that contains the data to transfer
#	- the name of the server that contains the initial data
# Returns the path of the directory created on the scratch space (local) and the node name that contains the data transfered
################################################################################################
#sub transferDirectoryFromMasterToNode
#{
#    my ($localDir,$master) = @_;
#								    ###########################################################################
#    $master = 'nas2' if (not defined $master or $master eq '');     #######SUPPOSE QUE LES DONNEES SONT TOUJOURS DANS NAS2 (donc pas /teams) - AJOUTER DANS fichier configuration?
#								    ###########################################################################
#    # get the SGE user name
#    my $SGE_User = `echo \$USER` or toolbox::exportLog("ERROR: toolbox::transferDirectoryFromMasterToNode : The SGE USER isn't defined\n",0);
#    chomp($SGE_User);
#    ##DEBUG exportLog("INFOS: toolbox::transferDirectoryFromMasterToNode : SGE USER $SGE_User\n",1);
#
#    # get the SGE node on what the script is running
#    my $SGE_Node = `echo \$HOSTNAME` or toolbox::exportLog("ERROR: toolbox::transferDirectoryFromMasterToNode : The SGE node isn't defined\n",0);
#    chomp($SGE_Node);
#    ##DEBUG exportLog("INFOS: toolbox::transferDirectoryFromMasterToNode : SGE NODE $SGE_Node\n",1);
#
#    # get the SGE job id
#    my $SGE_JobId = `echo \$JOB_ID` or toolbox::exportLog("ERROR: toolbox::transferDirectoryFromMasterToNode : The SGE id isn't defined\n",0);
#    chomp($SGE_JobId);
#    ##DEBUG exportLog("DEBUG INFOS: toolbox::transferDirectoryFromMasterToNode : SGE JOB ID $SGE_JobId\n",1);
#
#    # Creating local tmp folder according this convention : /scratch/$user-$jobid-$sge_task_id
#    my $tmpDir = "/scratch/".$SGE_User."-".$SGE_JobId."/";
#    toolbox::makeDir($tmpDir) if (not existsDir($tmpDir,0));      # Creation of the tempory directory
#    ##DEBUG exportLog("DEBUG INFOS: toolbox::transferDirectoryFromMasterToNode : SGE JOB ID $SGE_JobId\n",1);
#
#    # Copy the data from the master server to the local scrtatch space
#    my $cmd="scp -r $SGE_User\@$master:$localDir $tmpDir";
#    toolbox::run($cmd);
#
#    # Extract the name of the directory and return the complete local directory name and the node name
#    my ($file,$path)=toolbox::extractPath($localDir);
#    return($tmpDir.$file,$SGE_Node);
#
#}
################################################################################################
# END sub transferDirectoryFromMasterToNode
################################################################################################




################################################################################################
# sub transferDirectoryFromNodeToMaster =>  to transfer data from one local scratch space
#                                           of the node to a master server
################################################################################################
# arguments :
# 	- the path of the directory on the scratch space (local)
#	- the node name that contains the data to transfer
#	- the name of the server that contains the initial data
#	- a boolean parameter to define if the data wil be removed from the local space after
#         transferring from the local space to the master server. By default equal to 1 (Data removing)
# No parameter returned
################################################################################################
#sub transferDirectoryFromNodeToMaster
#{
#
#    my ($localDir,$distantDir,$erase,$master) = @_;
#								###########################################################################
#    $master = 'nas2' if (not defined $master or $master eq ''); #######SUPPOSE QUE LES DONNEES SONT TOUJOURS DANS NAS2 (donc pas /teams) - AJOUTER DANS fichier configuration?
#								###########################################################################
#
#    $erase=1 if (not defined $erase);
#
#     # get the SGE user name
#    my $SGE_User = `echo \$USER` or toolbox::exportLog("ERROR: toolbox::transferDirectoryFromMasterToNode : The SGE USER isn't defined\n",0);
#    chomp($SGE_User);
#    ##DEBUG exportLog("INFOS: toolbox::transferDirectoryFromMasterToNode : SGE USER $SGE_User $erase $master\n",1);
#
#    # Data transferring from the local directory to the distant directory of the master
#    my $cmd="scp -r $localDir $SGE_User\@$master:$distantDir";
#    toolbox::run($cmd);
#
#    # Removing of the local data directory if $erase equal to 1
#    $cmd="rm -rf $localDir";
#    toolbox::run($cmd) if ($erase);
#
#}
################################################################################################
# END sub transferDirectoryFromNodeToMaster
################################################################################################



################################################################################################
# sub relativeToAbsolutePath => modify a relative path in absolute one
################################################################################################
# arguments :
# 	- relative path
# returns :
#	- absolute path
################################################################################################
sub relativeToAbsolutePath
{
    my ($relative, $print)=@_;
    $print = 1 if (not defined $print);
    my ($absolutePath,$log);

    if ($relative =~ m/None$/) {
        $log = "INFOS : $0 : toolbox::relativeToAbsolutePath : the path $relative is not a path.";
        $absolutePath=$relative;
    }
    elsif ($relative !~ m/^\//) # The relative path is a relative path, ie do not starts with /
    {
        my $com = "readlink -m $relative";
        $absolutePath = `$com`;
        chomp $absolutePath;
        $log = "INFOS : $0 : toolbox::relativeToAbsolutePath : the relative path $relative has been modified as an absolute path in $absolutePath \n";
    }
    else #relative is in fact an absolute path, send a warning
    {
    	$log = "INFOS : $0 : toolbox::relativeToAbsolutePath : the path $relative is not a relative but an absolute. TOGGLE will not modify it \n";
    	$absolutePath = $relative;
    }

    # return only absolutePath if boolean = 0  or both (absolutePath and log) if boolean = 1
    if ($print==0)
    {
        return $absolutePath;
    }
    else
    {
        return ($absolutePath,$log);
    }


}
################################################################################################
# END sub relativeToAbsolutePath
################################################################################################


################################################################################################
# sub rmHashOrder => remove step in hashOrder if stepName if the first key
# NOTE/WARNING: this function was developed to remove only the first step (for processRadstacks)
################################################################################################
# arguments :
# 	- $hashOrder
#	- $stepName
# returns :
#	- $hashOrder with step removed
################################################################################################
sub rmHashOrder
{
	toolbox::exportLog("ERROR: toolbox::rmHashOrder : should get at least one arguments\n",0) if (@_ < 1);
	my ($hashOrder,$stepName)=@_;

	my @firstKeys = sort {$a <=> $b} (keys(%{$hashOrder}));
	if ($hashOrder->{$firstKeys[0]} eq $stepName)
	{
		delete $$hashOrder{$firstKeys[0]};
	}
	return $hashOrder;
}
################################################################################################
# END sub rmHashOrder
################################################################################################



1;

=head1 NAME

package I<toolbox>

=head1 SYNOPSIS

	use toolbox;

	toolbox::exportLog("ERROR: $0 : The directory don't contain the right number of files\n",1);


=head1 DESCRIPTION

This module aims to provide a set of functions related to check file format, to check/list the content of a directory etc.


=head2 Functions

=head3 toolbox::exportLog()

This function aims to manage the  logs and the errors generated by others functions or scripts. Two arguments are required :
- lines to print in a log or error file
- boolean value to control the type of the log and of the error. 0 : print logs into error log file and die.
1 : print logs into log file. 2 : print logs into error log file and the program continue.

No parameters returned

Example :
C<toolbox::exportLog("ERROR: $0 : The directory don't contain the right number of files\n",1);>


=head3 toolbox::checkFile()

This function checks if a file exists, is readable/writable, and isn't empty. One argument is required : the name of the file to analyze.

Returns 1 if the file exists, is readable and writable AND isn't empty.

Example :
C<toolbox::checkFile($fileToCheck);>


=head3 toolbox::readFile()

This function checks if a file is readable. One argument is required : the name of the file to analyze.

Returns 1 if the file is readable else return 0.

Example :
C<toolbox::readFile($fileToCheck);>


=head3 toolbox::writeFile()

This function checks if a file is writable. One argument is required : the name of the file to analyze.

Returns 1 if the file is writable else return 0.

Example :
C<toolbox::writeFile($fileToCheck);>


=head3 toolbox::sizeFile()

This function checks if a file is empty or not. One argument is required : the name of the file to analyze.

Returns 1 if the file isn't empty else return 0.

Example :
C<toolbox::sizeFile($fileToCheck);>



=head3 toolbox::existFile()

This function checks if a file exists. One argument is required : the name of the file to analyze. An other optional argument can be used,
a boolean parameter to determine when the file exists, if it's a error or not.

Returns 1 if the file exists else return 0.

Example :
C<toolbox::existFile($fileToCheck);>



=head3 toolbox::existsDir()

This function checks if a directory exists. One argument is required : the name of the directory to analyze. An other optional argument can be used,
a boolean parameter to determine when the directory exists, if it's a error or not.

Returns 1 if the directory exists else return 0.

Example :
C<toolbox::existsDir($dirToCheck);>


=head3 toolbox::makeDir()

This function allows to create a directory. One argument is required : the name of the directory to create. An other optional argument can be used,
a boolean parameter to control when the directory exists, if it will be removed before creating the directory.
If the boolean parameter is equal to 1 and, the existing directory is removed before creating it again

Returns 1 if the directory has been correctly created else return 0.

Example :
C<toolbox::makeDir($dirToCreate);>


=head3 toolbox::readDir()

This function lists the content of a directory and returns the list of files. One argument is required : the name of the directory.
An other optional argument can be used that allows to filter on the end of a filename.

Returns the list of files.

Example :
C<toolbox::readDir($dirToParse,'fastq');>

=head3 toolbox::readDir2()

This function lists the content of a directory and returns the list of files. One argument is required : the name of the directory.
An other optional argument can be used that allows to filter on the beginning of a filename.

Returns the list of files.

Example :
C<toolbox::readDir2($dirToParse,'RC');>


=head3 toolbox::extractPath()

This function extract the filename and the path from a complete file path. One argument is required : the file path.

Returns separately the filename and the path

Example :
C<toolbox::extractPath('/data/projet/chloro/RC1/RC1_1.fastq');>


=head3 toolbox::readFileConf()

This function reads the Software Configuration File and stores the options of every software into a hash of hashes.
One argument is required : the filename.

Returns a hash of hashes

Example :
C<toolbox::readFileConf('/data/projet/chloro/software.config');>


=head3 toolbox::extractOptions()

This function extracts the options/values for a software given (from a hash of hash) and return them as a string. One argument is required :
a hash of hash (reference) with the software name as key, then the option as key and the value as element. A second optional parameter can be used,as type of separators used to separate an option from its value (space by default)
The third (optional) value is the way to concatenate the different options together (space by default)

Returns a string

Example :
C<my $optionLine=toolbox::extractOptions($configInfos->{"BWA aln"}," ", " ");>


=head3 toolbox::extractName()

This function returns the filename from a complete path given as argument.
Returns a string

Example :
C<my $filename=toolbox::extractPath('/data/projects/oryza/RC1/RC1.fastq');>



=head3 toolbox::run()

This function executes the command you gave it and exports the log via toolbox::exportLog() function.
Returns  boolean (1 if the command has been correctly executed else 0)

Example :
C<toolbox::run("rm /data/projects/oryza/*.fastq"); >






=head3 toolbox::addInfoHeader()

This function adds a information tag in a BAM file (@CO field). Two arguments are required : filename to analyze and text to add.
Returns  boolean (1 if this function has been correctly running else 0)

Example :
C<toolbox::addInfoHeader($bamFile, 'RC1'); >






=head3 toolbox::changeDirectoryArbo()

This function returns a directory name according to the number given as argument. Two arguments are required: directory that will contain all the directory created by the pipeline
and a number to define what is the directory name that will be returned
Returns the directory name that will be created

Example :
C<toolbox::changeDirectoryArbo($initialDir,1); >



=head3 toolbox::extractHashSoft()

This function returns the list of parameters for a software given. Two arguments are required: reference oh the hash that contains all softawre name (from toolbox::readFileConf) and the software name searched for.
Returns the sub hash for the software name given as parameter with the list of parameters.

Example :
C<toolbox::extractHashSoft($optionref,"BWA index"); >



=head3 toolbox::checkInitialDirContent()

This function analyzes a folder given as argument and returns 0 if the directory contains only files or the the list of folders contained.

Example :
C<toolbox::checkInitialDirContent($initialDir); >



=head3 toolbox::transferDirectoryFromMasterToNode()

This function copy data directory from a master server to one local scratch space of the node. One argument is required :
the directory name that contains the data to transfer. The second argument is the name of the server that contains the initial data.
Returns the path of the directory created on the scratch space (local) and the node name that contains the data transfered

Example :
C<( $tmp_refFastaFile,$nodeFinal)=toolbox::transferDirectoryFromMasterToNode($refFastaFile,$nodeInitial); >



=head3 toolbox::transferDirectoryFromNodeToMaster()

This function copy data from one local scratch space of the node into a distant directory of the master server.
Three arguments are required: (1) the path of the directory on the scratch space (local), (2) the node name that contains the data to transfer,
 (3) the name of the master server that contains the initial data.
One argument is optional : a boolean parameter to define if the data wil be removed from the local space after data transferring from the local space to the master server.
No parameter returned.

Example :
C<( toolbox::transferDirectoryFromNodeToMaster($initialDir."/*", $MasterDir, $nodeInitial,1);); >



=head3 toolbox::checkNumberLines

This module check the sequences number of a given file using the wc -l command from bash
It takes only one argument, the file you want to know the number of lines

=head3 toolbox::relativeToAbsolutePath

This module will transform a relative path to the corresponding absolute one using the readlink -m bash command.
Take note that if the relative is wrong, the absolute will be wrong too.
The module will inform for transformation and will not modify already absolute paths (starting by '/')
It takes only one argument, the relative path you want to change

=head3 toolbox::rmHashOrder

This function removes the step given in argument from $hashOrder.
True only if this step is the first one defined in the configuration file

=head1 AUTHORS

Cecile Monat, Ayite Kougbeadjo, Julie Orjuela-Bouniol, Marilyne Summo, Cedric Farcy, Mawusse Agbessi, Christine Tranchant and Francois Sabot

L<http://www.southgreen.fr/>


=cut
