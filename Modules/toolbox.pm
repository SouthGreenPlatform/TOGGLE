package toolbox;


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
use Data::Dumper;
use Exporter;

use lib qw(.);
use localConfig;

#Global infos
our @ISA=qw(Exporter);
our @EXPORT=qw($configInfos);
our $configInfos; #Config informations, ref of an hash



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

    # Get the name of indivudal analysed when this function is called
    my $indivName = `head -n 1 individuSoft.txt`;
    chomp $indivName;
 
    # Get the name of software executed when this function is called 
    my $currentSoft = `tail -n 1 individuSoft.txt`;
    chomp $currentSoft;
    
    # Initialization of the log name from the name of the individual and the name of software
    my $logBasicName=$indivName."_".$currentSoft."_log";
    
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
    my $writingRights = writeFile($file);
    
    my $logOut;
    
    if ($readingRights == 1 and $writingRights == 1)
    {
        $logOut= "INFOS: toolbox::checkFile : The file $file is readable and writable\n";
    }
    elsif($readingRights == 1 and $writingRights == 0)
    {
        $logOut="ERROR: toolbox::checkFile : The file $file is readable but not writable\n";
    }
    elsif($readingRights == 0 and $writingRights ==1)
    {
        $logOut="ERROR: toolbox::checkFile : The file $file is writable but not readable\n";
        exportLog($logOut,0);
        return 1;
    }
    else
    {
        $logOut="ERROR: toolbox::checkFile : The file $file is not readable and writable\n";
        exportLog($logOut,0);
        return 1;
    }
    
    exportLog($logOut,1);   
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
        my $logOut ="INFOS: toolbox::existsFile : The file $file does not exist or is not a file\n$!\n";
        exportLog($logOut,1);
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
sub readDir2
{
    toolbox::exportLog("ERROR: toolbox::readDir2 : should get at one argument at least\n",0) if (@_ < 1 );
    
    my ($dir)= shift @_;
    
    # if the part of name is given in the second argument, search for files directory/part_of_name*
    # else search for all files directory/*
    my $path = defined ($_[0]) ? $dir.'/'.$_[0].'*' : $dir."/*";
    
    # ls command
    my $file=`ls $path` or toolbox::exportLog("ERROR: toolbox::readDir : Can't open the directory $path\n$!\n",0);
    chomp $file;
    
    # split the list into a table returned after
    my @fileList = split /\n/, $file; 
    return(\@fileList);

}
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
#
# ex : my $optionLine=toolbox::extractOptions($configInfos->{"BWA aln"}," ");
#
# Returns the list of options for a software given. 
################################################################################################
sub extractOptions
{
    ##Getting the two parameters, the options hash and the option-value separator
    my($optionsHashees,$separateur)=@_;
    
    if ($optionsHashees)			# the option hash isn't empty/is defined => options are extracted
    {
	my %options=%{$optionsHashees};
	my $option=" ";
	$separateur=" " unless $separateur; 	## if no separator is given, set it as a single space
	try
	{                               
	    foreach my $cle (keys %options)
	    {
		if ($options{$cle} eq 'NA') 	## The option has no value  => no print $options{$cle}
		{
		    $option=$option.$cle.$separateur." ";
		}
		else				## The option has value  => print $options{$cle}
		{
		    $option=$option.$cle.$separateur.$options{$cle}." ";
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
    
    $bruteName=~ s/^.+\/(.+)\..+/$1/; # get just the name of the file after the last /
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
    # copy STDOUT and STDERR to another filehandle
    #open (my $STDOUT_, '>&', STDOUT);
    #open (my $STDERR_, '>&', STDERR);

    # redirect STDOUT and STDERR to log.txt
    #open (STDOUT, '>>', 'log.txt');
    #open (STDERR, '>>', 'log.txt');

    #TODO: Change the adresses for STDOUT and STDERR
    
    my($command)=@_;
    exportLog("INFOS: toolbox::run : $command\n",1);
    
    ##Execute the command
    my ($result,$stderr)=capture {` $command `};
    
    ##Log export according to the error
    if ($?==0)
    {
	##DEBUG
	exportLog("INFOS: toolbox::run : ".$result."\n--".$stderr."\n",1);
	return 1;
    }
    else
    {
	##DEBUG
	exportLog("ERROR: toolbox::run : ".$result."\n--".$stderr."\n",0);
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
    my $nbLineCommand="wc -l ".$fileName; #command to count the line number
    my $nbLine = `$nbLineCommand` or exportLog("ERROR: toolbox::checkNumberLines : Cannot run $nbLineCommand\n$!\n",0);	# execution of the command or if not possible, return an error message
    chomp $nbLine;
    
    #Add a split to only keep the number of line without the file name
    my @splitLine = split (" ", $nbLine);
    $nbLine = $splitLine[0];
    
    return $nbLine;  # return the number of line of file
}


################################################################################################
# sub checkFormatFastq => check if a file is really a FASTQ file
################################################################################################
# arguments : filename to analyze
# Returns boolean (1 if the format is fastq else 0)
################################################################################################

sub checkFormatFastq
{
    
    my $notOk = 0;                                                      # counter of error(s)
    my ($fileToTest) = @_;                                              # recovery of file to test
    my $readOk = readFile($fileToTest);                                 # check if the file to test is readable
    
    my $nbLines = toolbox::checkNumberLines(@_);                    # calculing number lines in file
    my $modulo = ($nbLines % 4);
    my $even   = ($nbLines % 2);
    
    if ( ($nbLines>0) and ($modulo==0) and ($even==0) )                # testing if the number of lines is a multiple of 4
    {
        #print "$nbLines is a multiple of 4\n";
    }
    else {
        toolbox::exportLog("ERROR: toolbox::checkFormatFastq : Number of lines is not a multiple of 4 in file $fileToTest.\n",0);
        return 0;
    }
                                                                        # open and traite the file if the number of lines is a multiple of 4
    open (F1, $fileToTest) or toolbox::exportLog("ERROR: toolbox::checkFormatFastq : Cannot open the file $fileToTest\n$!\n",0); # open the file to test
    
    my  @linesF1=();
    my $comp=0;
    my $countlines=0;
    my $stop=0;
    
    while ((my $line = <F1>))                                           # scanning file and stocking in an array the four lines of a read.
    {
        chomp $line;
        $countlines++;
        
        if ($comp<3)
        {
            $comp++;
            push (@linesF1,$line);
        }
        else                                                            # Completing block, treatment starts here.
        {
            $stop++;
            push (@linesF1,$line);
            
            my $i=0;
            while ( ($i<=$#linesF1) and ($notOk <=1))                 # treatment of a block containing four lines of file and stop if 20 errors found.
            {
                
                my $idLine=$linesF1[$i];
                my $fastaLine=$linesF1[$i+1];
                my $plusLine=$linesF1[$i+2];
                my $qualityLine=$linesF1[$i+3];
                my $nbIDLine=$countlines-3;
                my $nbLineFasta=$countlines-2;
                my $nbPlusLine=$countlines-1;
                my $nbQualityLine=$countlines;
                
                if (($idLine=~m/^$/) and ($plusLine=~m/^$/))            # if the ID line and the "plus" line are not empty ...
                {
                    toolbox::exportLog("ERROR: toolbox::checkFormatFastq : The file $fileToTest is not a FASTQ file => The ID infos line $nbIDLine is not present.\n",0);
                    $notOk++;                                           # one error occured, so count it
                }
                
                elsif ( (($idLine=~m/^\@.*/) or ($idLine=~m/^\>.*/) ) and ($plusLine=~m/^\+$/) )   # if ID ligne is not empty and it starts by "@" or ">" and the
                # plus line has a "+", the block of four lines ($i to $i+3) is traited.
                {
                    if ( length($fastaLine) == length($qualityLine) )   # comparing the fasta line and the quality line lengths.
                    {
                        my @fasta = split //, $fastaLine;
                        foreach my $nucleotide (@fasta)
                        {
                            if ($nucleotide!~m/A|T|G|C|a|g|t|c|N|n/)    # checking nucleotides in the fasta line.
                            {
                                toolbox::exportLog ("ERROR: toolbox::checkFormatFastq : Not basic IUPAC letter, only ATGCNatgcn characters are allowed: unauthorized characters are in the line $nbLineFasta of $fileToTest.\n",0);
                                $notOk++;
                            }
                        }
                    }
                    else 												# error if fasta line length and quality line length are differents.
                    {
                        toolbox::exportLog("ERROR: toolbox::checkFormatFastq : Fasta line $nbLineFasta has not the same lenght than quality line $nbQualityLine in file $fileToTest.\n",0);
                        $notOk++;
                    }
                }
                
                else													#error if the ID line do not start with @ or >.
                {
                    toolbox::exportLog("ERROR: toolbox::checkFormatFastq : ID line has to start with @ or > in line $nbIDLine of file $fileToTest.\n",0);
                    $notOk++;
                }
                $i=$i+4; 												# jumping to next read.
            }
            
            last if ($stop==200000);                                    # stoping treatment if 50000 reads were analysed.
            
            undef @linesF1;
            $comp=0;
        }
        next;
    }
    
    if (($notOk == 0))                    						# if any error occured in the file, the format is right.
    {
        toolbox::exportLog("INFOS: toolbox::checkFormatFastq : The file $fileToTest is a FASTQ file.\n",1);
	return 1;
    }
    else                                						# if one or some error(s) occured on the file, the fastq format is not right.
    {
        toolbox::exportLog("ERROR: toolbox::checkFormatFastq : Invalid FASTQ requirements in file $fileToTest.\n",0);
	return 0;
    }
    
    close F1;
    
}
################################################################################################
# END sub checkFormatFastq
################################################################################################






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
    my $formatCheck=checkSamOrBamFormat($bamFile);
    
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
# sub checkSamOrBamFormat => verifying the SAM/BAM format based on samtools view system
# samtools view gave an error in case of non SAM or BAM format
################################################################################################
# arguments : filename to analyze
# Returns boolean (1 if the fileformat is sam, 2 bam and 0 neither bam or sam)
################################################################################################
sub checkSamOrBamFormat {
    
    my ($samFile)=@_;
    
    existsFile($samFile); # Will check if the submitted file exists
       
    #Is the file sam of bam through the binary mode? Requested for the -S option in samtools view
    my ($inputOption,$binary);
    if (-B $samFile) #The file is a binary BAM file
    {
	$inputOption = ""; #no specific option in samtools view requested
	$binary = 1; # the file is binary
    }
    else #the file is a not binary SAM file
    {
	$inputOption = " -S ";#-S mean input is SAM
	$binary = 0; # the file is not binary
    }
    my $checkFormatCommand="$samtools view $inputOption $samFile -H > /dev/null";
    # The samtools view will ask only for the header to be outputted (-H option), and the STDOUT is redirected to nowher (>/dev/null);
    my $formatValidation=run($checkFormatCommand);
    
    if ($formatValidation == 1)                    # if no error occured in extracting header, ok
    {
        toolbox::exportLog("INFOS: toolbox::checkSamOrBamFormat : The file $samFile is a SAM/BAM file\n",1);
	return 1 if $binary == 0;# Return 1 if the file is a SAM
	return 2 if $binary == 1;# Return 2 if the file is a BAM
    }
    else                                # if one or some error(s) occured in extracting header, not ok
    {
        toolbox::exportLog("ERROR: toolbox::checkSamOrBamFormat : The file $samFile is not a SAM/BAM file\n",0);
	return 0;
    }
}
################################################################################################
# END checkSamOrBamFormat 
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
# sub checkVcfFormat=> checks if the fileformat is VCF
################################################################################################
# arguments : file name to analyze
# Returns 1 if the file format is validated else 0 (error are managed by toolbox::exportLog)
################################################################################################
sub checkVcfFormat
{ 
    #Inspired from The vcftools vcf-validator
    my ($file)=@_;
    
    #Check if we can read the file
    my $rightToRead = readFile($file);
    if ($rightToRead == 0)
    {
	exportLog("ERROR: toolbox::checkVcfFormat : The file $file is not readable\n",0);
	return 0;
    }
    
    #Parsing the file
    my @header;#List to gather the header
    my @listOfFields;
    open(VCF, "<",$file) or toolbox::exportLog("ERROR: toolbox::checkVcfFormat : Cannot open the file $file\n$!\n",0); 
    
    while (my $line=<VCF>)
    {
	chomp $line;
	
	##DEBUG print $line."\n";
	if ($line =~ m/^#/)
	{
	   push (@header, $line);
	}
	else {@listOfFields = split /\t/, $line;}
    }
    
    #Check if the first line of the header is including the version
    my $versionLine=defined $header[0]? $header[0]:undef;
    exportLog("ERROR: toolbox::checkVcfFormat : There is no header of the file $file\n",0) unless (defined $versionLine); #Thrown an error if the $version cannot be obtained (misformatted line)	

    my @version=split /VCFv/,$versionLine;
    exportLog("ERROR: toolbox::checkVcfFormat : Cannot evaluate the VCF version of the file $file file\n",0) if (undef @version); #Thrown an error if the $version cannot be obtained (misformatted line)
    eval ($version[1] == $version[1]) or exportLog("ERROR: toolbox::checkVcfFormat : Cannot obtain the VCF version of $file\n",0); #Verifying if the value obtained is numerical.
    
    # Check the first line format as recommanded
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Individual1
    #Chr1    228069  .       A       .       31.23   LowQual AN=2;DP=1;MQ=37.00;MQ0=0        GT:DP   0/0:1
    if (scalar @listOfFields < 10) #Less than the 10 minimum fields
    {
	exportLog("ERROR: toolbox::checkVcfFormat : The file $file is misformatted (less than 10 colums) ".Dumper(\@listOfFields)."\n",0);
    }
    
    #Check if the second field (the position) is numerical
    eval ($listOfFields[1] == $listOfFields [1]) or exportLog("Cannot confirm that $file is a VCF.\nAborting.\n",0); #Verifying if numerical. Die if not
   
    close VCF;

    return 1; #Return correct if all check are Ok
}
################################################################################################
# END sub checkVcfFormat
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
sub transferDirectoryFromMasterToNode
{
    my ($localDir,$master) = @_;
								    ###########################################################################
    $master = 'nas' if (not defined $master or $master eq '');     #######SUPPOSE QUE LES DONNEES SONT TOUJOURS DANS NAS2 (donc pas /teams) - AJOUTER DANS fichier configuration?
								    ###########################################################################
    # get the SGE user name
    my $SGE_User = `echo \$USER` or toolbox::exportLog("ERROR: toolbox::transferDirectoryFromMasterToNode : The SGE USER isn't defined\n",0);
    chomp($SGE_User);
    ##DEBUG exportLog("INFOS: toolbox::transferDirectoryFromMasterToNode : SGE USER $SGE_User\n",1);
    
    # get the SGE node on what the script is running
    my $SGE_Node = `echo \$HOSTNAME` or toolbox::exportLog("ERROR: toolbox::transferDirectoryFromMasterToNode : The SGE node isn't defined\n",0);
    chomp($SGE_Node);
    ##DEBUG exportLog("INFOS: toolbox::transferDirectoryFromMasterToNode : SGE NODE $SGE_Node\n",1);
    
    # get the SGE job id
    my $SGE_JobId = `echo \$JOB_ID` or toolbox::exportLog("ERROR: toolbox::transferDirectoryFromMasterToNode : The SGE id isn't defined\n",0);
    chomp($SGE_JobId);
    ##DEBUG exportLog("DEBUG INFOS: toolbox::transferDirectoryFromMasterToNode : SGE JOB ID $SGE_JobId\n",1);
    
    # Creating local tmp folder according this convention : /scratch/$user-$jobid-$sge_task_id
    my $tmpDir = "/scratch/".$SGE_User."-".$SGE_JobId."/";
    toolbox::makeDir($tmpDir) if (not existsDir($tmpDir,0));      # Creation of the tempory directory
    ##DEBUG exportLog("DEBUG INFOS: toolbox::transferDirectoryFromMasterToNode : SGE JOB ID $SGE_JobId\n",1);

    # Copy the data from the master server to the local scrtatch space
    my $cmd="scp -r $SGE_User\@$master:$localDir $tmpDir";
    toolbox::run($cmd);
    
    # Extract the name of the directory and return the complete local directory name and the node name
    my ($file,$path)=toolbox::extractPath($localDir);
    return($tmpDir.$file,$SGE_Node);

}
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
sub transferDirectoryFromNodeToMaster
{
    
    my ($localDir,$distantDir,$erase,$master) = @_;
								###########################################################################
    $master = 'nas'; # if (not defined $master or $master eq ''); #######SUPPOSE QUE LES DONNEES SONT TOUJOURS DANS NAS2 (donc pas /teams) - AJOUTER DANS fichier configuration?
								###########################################################################
    
    $erase=1 ; # if (not defined $erase);					
    
     # get the SGE user name
    my $SGE_User = `echo \$USER` or toolbox::exportLog("ERROR: toolbox::transferDirectoryFromMasterToNode : The SGE USER isn't defined\n",0);
    chomp($SGE_User);
    ##DEBUG exportLog("INFOS: toolbox::transferDirectoryFromMasterToNode : SGE USER $SGE_User $erase $master\n",1);
    
    # Data transferring from the local directory to the distant directory of the master
    my $cmd="scp -r $localDir $SGE_User\@$master:$distantDir";
    toolbox::run($cmd);
    
    # Removing of the local data directory if $erase equal to 1
    $cmd="rm -rf $localDir";
    toolbox::run($cmd) if ($erase);   
    
}
################################################################################################
# END sub transferDirectoryFromNodeToMaster 
################################################################################################


################################################################################################
# sub checkFormatFasta => check if a given file is a fasta or not
################################################################################################
# arguments :
# 	- the fasta file
# returns :
#	- a boolean 1 for Ok, 0 for not.
#	- send a warning to the log if not a fasta file
################################################################################################

sub checkFormatFasta{
    my ($file)=@_;

    #Check if we can read the file
    my $rightToRead = toolbox::readFile($file);
    if ($rightToRead == 0)
    {
	toolbox::exportLog("ERROR: toolbox::checkFormatFasta : The file $file is not readable\n",0);
	return 0;
    }
    
    #Opening file
    open(FILE, "<", $file) or toolbox::exportLog("ERROR: toolbox::checkFormatFasta : Cannot open the file $file\n$!\n",0);
   
    #Reading file
    my $initiator = 0; #We need to scan the very first line, that begins with '>', but this is the only one beginning with that we need
    my %errors; # Will score all the error types
    my $lineNumber=0;#Recording line position
    while (my $line = <FILE>)
       {
	chomp $line;
	$lineNumber++;#Counter for scoring lines with errors
    #Checking first line, that must be '>NAME'
	if ($initiator == 0)
	    {
	    #very first line
	    $initiator = 1; #no need for any other lines beginning with ##
	    next if $line =~ m/^>/; #Correct first line, let's go to the following ones
	    #First line is not correct
	    $errors{$lineNumber} = "Not correctly formatted, must be a header name such as '>NAME'.";
            #Will stop immediatly as the file is not correct at all
            last;
	    }
    #Other lines	    
	next if $line =~ m/^$/;# Representing empty lines
	if ($line =~ m/^>/) #A new sequence beginning
	    {
	    next;
	    }
	else # The sequence is read
	    {
	    my $modifiedLine = $line;
	    $modifiedLine =~ s/[A|T|G|C|a|g|t|c|N|n]//g;
	    #DEBUG : print "\n**",$modifiedLine,"**\n";
	    #DEBUG : exit;
	    next if length($modifiedLine eq "");
	    $errors{$lineNumber} = "Not basic IUPAC letter, only ATGCNatgcn characters are allowed: unauthorized characters are $modifiedLine.";
	    }
	#Check if there are too much errors...
        if (scalar (keys %errors) > 19)
          {
          #Will stop after 20 errors
          my $tooMuchErrors = "WARNING : toolbox::checkFormatFasta : Too much errors in the Fasta file $file, only the first 20 errors are shown..."; #Will be printed before all errors
	  toolbox::exportLog($tooMuchErrors,2);
          last;
          }
        next;
    
        }
    close FILE;
    
    #Checking if there are errors
    my $numberOfErrors = scalar(keys %errors);
    
    #Sending the warning message, if any
    if ($numberOfErrors)#There are errors, we will send a warning to the system
	{
	#Formatting errors for the warning
	my $warningLog;
	$warningLog.="WARNING : toolbox::checkFormatFasta : There are $numberOfErrors errors detected for file $file\n";
	foreach my $individualError (sort {$a <=> $b} keys %errors)#Sorting error per line number
	    {
	    $warningLog.="WARNING : toolbox::checkFormatFasta : Line $individualError does not respect FASTA standard : $errors{$individualError}\n";
	    }
	#Sending the warning
	toolbox::exportLog($warningLog,2);
	#returning for failure	
	return 0;
        }
       
    #returning ok
    toolbox::exportLog("INFOS: toolbox::checkFormatFasta : The file $file is a FASTA file\n",1);
    return 1;    
    }
   
   
################################################################################################
# END sub fastaFormatValidator 
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
a hash of hash (reference) with the software name as key, then the option as key and the value as element. A second optional parameter can be used,
it corresponds to the type of separators  used to separate an option from its value (space by default)

Returns a string

Example : 
C<my $optionLine=toolbox::extractOptions($configInfos->{"BWA aln"}," ");>


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



=head3 toolbox::checkFormatFastq()

This function checks if a file is really a FASTQ file. 
Returns  boolean (1 if the format is fastq else 0)

Example : 
C<toolbox::checkFormatFastq("/data/projects/oryza/RC1.fastq"); >


=head3 toolbox::addInfoHeader()

This function adds a information tag in a BAM file (@CO field). Two arguments are required : filename to analyze and text to add.
Returns  boolean (1 if this function has been correctly running else 0)

Example : 
C<toolbox::addInfoHeader($bamFile, 'RC1'); >



=head3 toolbox::checkSamOrBamFormat()

This function verifies the SAM/BAM format based on samtools view system. One argument is required : filename to analyze.
Returns  boolean (1 if the fileformat is sam, 2 bam and 0 neither bam or sam)

Example : 
C<toolbox::checkSamOrBamFormat($samFile); >


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



=head3 toolbox::checkVcfFormat()

This function checks if the format of the given as argument is VCF.
Returns 1  if the file format is validated else 0 (an error is maneged by toolbox::exportLog)

Example : 
C<toolbox::checkVcfForm($file); >



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


=head3 toolbox::checkFormatFasta()

This function checks if a given file is a true DNA fasta file
The only required argument is the filemane. The file can be plain text
Returns a 1 for success, and a 0 for failure. Will send a warning to the log in case of failure.
Will return a maximum of 20 errors.
Will stop immediatly if the first line is misformatted

Example : 
C<( toolbox::checkFormatFasta($fastaFile);); >


=head3 toolbox::checkFormatFastq()
 
This function checks if a given file is a FASTQ file.
The only required argument is the filename.
Returns a 1 for success, and a 0 for failure. Will send a warning to the log in case of failure.
Will return a maximum of 1 errors.
Will stop immediatly if the first line is misformatted
Use toolbox::checkNumberLines to count the number of lines and calculate if this number is a multiple of 4.
Use a 4 lines block to avoid stocking in memory the whole of lines from file.
Check the format of the first 50000 fastq reads.
 
Example :
toolbox::checkFormatFasta($fastaFile);


=head3 toolbox::checkNumberLines

This module check the sequences number of a given file using the wc -l command from bash
It takes only one argument, the file you want to know the number of lines

=head1 AUTHORS
 
Cecile Monat, Ayite Kougbeadjo, Julie Orjuela-Bouniol, Marilyne Summo, Cedric Farcy, Mawusse Agbessi, Christine Tranchant and Francois Sabot

L<http://www.southgreen.fr/>


=cut



