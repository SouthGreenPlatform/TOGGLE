#!/usr/bin/env perl

###################################################################################################################################
#
# Copyright 2014-2015 IRD-CIRAD-INRA-ADNid
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
# Version 3 written by Cecile Monat, Christine Tranchant, Cedric Farcy, Maryline Summo, Julie Orjuela-Bouniol, Sebastien Ravel, Gautier Sarah, and Francois Sabot
#
###################################################################################################################################

use strict;
use warnings;
use localConfig;
use Data::Dumper;

use pairing;
use toolbox;
use onTheFly;
use scheduler;

#For gz files
use IO::Compress::Gzip qw(gzip $GzipError);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);



##########################################
# recovery of parameters/arguments given when the program is executed
##########################################
my $cmd_line=$0." @ARGV";
my ($nomprog)=$0=~/([^\/]+)$/;
unless ($#ARGV>=0)                    # if no argument given
{

  print <<"Mesg";

  perldoc $nomprog display the help

Mesg

  exit;
}

my %param = @ARGV;                     # get the parameters 
if (not defined($param{'-d'}) or not defined($param{'-c'}) or not defined($param{'-r'}) or not defined ($param{'-o'}))
{
  print <<"Mesg";

  ERROR: Parameters -d and -c and -r and -o are required.
  perldoc $nomprog display the help

Mesg
  exit;
}




##########################################
# recovery of initial informations/files
##########################################
##########################################
# transforming relative path in absolute
##########################################
my @logPathInfos;
foreach my $inputParameters (keys %param)
{
  my ($newPath,$log)=toolbox::relativeToAbsolutePath($param{$inputParameters});
  die ("ERROR: $0 : An empty parameter has been given!\n") if ($newPath eq 0);
  $param{$inputParameters}=$newPath;
  push @logPathInfos,$log;
}

my $initialDir = $param{'-d'};        # recovery of the name of the directory to analyse
my $fileConf = $param{'-c'};          # recovery of the name of the software.configuration.txt file
my $refFastaFile = $param{'-r'};      # recovery of the reference file
my $outputDir = $param{'-o'};         # recovery of the output folder

my $gffFile;                          # recovery of the gff file used by topHat and rnaseq analysis
$gffFile = $param{'-g'} if (defined $param{'-g'});





##########################################
# Creation of the output folder
##########################################

if (not -d $outputDir) #if the output folder is not existing yet
{
    #creating the output folder
    my $createOutputDirCommand = "mkdir -p $outputDir";
    system ("$createOutputDirCommand") and die ("\nERROR: $0 : cannot create the output folder $outputDir: $!\nExiting...\n");
}

chdir $outputDir;

#Checking if $outputDir is empty

my $lsOutputDir = `ls`;
chomp $lsOutputDir;
if ($lsOutputDir ne "") # The folder is not empty
{
  die ("\nERROR: $0 : The output directory $outputDir is not empty, TOGGLE will not continue\nPlease provide an empty directory for outputting results.\n\nExiting...\n\n");
}

my $infosFile = "individuSoft.txt";

my ($sec, $min, $h, $mois_jour, $mois, $an, $sem_jour, $cal_jour, $heure_ete) = localtime(time);
$mois+=1;
$mois = $mois < 10 ? $mois = "0".$mois : $mois;
$mois_jour = $mois_jour < 10 ? $mois_jour = "0".$mois_jour : $mois_jour;
$h = $h < 10 ? $h = "0".$h : $h;
$min = $min < 10 ? $min = "0".$min : $min;
$an+=1900;
my $date="$mois_jour-$mois-$an-$h"."_"."$min";


#my $infosFile = "$pathIndividu[1]/individuSoft.txt";
open (F1, ">",$infosFile) or die ("$0 : open error of $infosFile .... $!\n");
print F1 "GLOBAL\n";
print F1 "ANALYSIS_$date\n";

toolbox::exportLog("#########################################\nINFOS: TOGGLE analysis start \n#########################################\n",1);;
toolbox::exportLog("INFOS: $0 : Command line : $cmd_line\n",1);
toolbox::exportLog("INFOS: Your output folder is $outputDir\n",1);

toolbox::exportLog("#########################################\nINFOS: Data checking \n#########################################\n",1);
toolbox::checkFile($fileConf);                              # check if this file exists
toolbox::existsDir($initialDir);                            # check if this directory exists
toolbox::checkFile($refFastaFile);                          #Retriving the configuration

##########################################
# Charging config Infos and copying the software config file
########################################

my $configInfo=toolbox::readFileConf($fileConf);
my $copyCommand = "cp $fileConf $outputDir/.";
toolbox::run($copyCommand,"noprint");

##########################################
# Printing the absolutePath changing logs
#########################################
##DEBUG foreach my $logInfo (@logPathInfos)
##DEBUG {
##DEBUG   toolbox::exportLog($logInfo,1);
##DEBUG }

#Verifying the correct ordering for the experiment, based on input output files and recovering the last step value
my ($firstOrder,$lastOrder) = onTheFly::checkOrder($configInfo);


##########################################
# Transferring data in the output folder and organizing
#########################################

#Checking if inly regular files are present in the initial directory
my $initialDirContent=toolbox::readDir($initialDir);

my $initialDirFolder=toolbox::checkInitialDirContent($initialDir);

if ($initialDirFolder != 0)#The initial dir contains subdirectories, so dying
{
    toolbox::exportLog("ERROR : $0 : The initial data directory $initialDir contains subdirectories and not only regular files.\n",0);
}

#Checking input data homogeneity

my $previousExtension=0;
foreach my $file (@{$initialDirContent})
{
    $file =~ m/\.(\w+)$/;
    my $extension = $1;
    #If the file is a compressed file in gz format
    if ($extension eq "gz")
    {
      if ($file =~ m/fastq\.gz$/ or $file =~ m/fq\.gz$/)#The file is a gz compressed fastq
      {
        $extension = "fastq";
      }
      elsif ($file =~ m/vcf\.gz$/) #The file is a gz compressed vcf
      {
        $extension = "vcf";
      }
      else # The file is neither a fastq.gz nor a vcf.gz file
      {
        toolbox::exportLog("ERROR : $0 : The compressed file $file format in the initial directory is not taken in charge by TOGGLE.\n",0);
      }
    }
    #The file is not a compressed one in gz
    if ($extension eq "fq" or $extension eq "fastq") #homogeneisation for fastq extension
    {
        $extension = "fastq";
    }
    if ($previousExtension eq "0") #first round
    {
        $previousExtension = $extension;
        next;
    }
    if ($previousExtension ne $extension) #not the same extension
    {
        toolbox::exportLog("ERROR : $0 : The file type in the initial directory are not homogeneous : $previousExtension and $extension are not compatible in the same analysis.\n",0);
    }
    next;
}

#Checking if the files are taken in charge by TOGGLE

if ($previousExtension !~ m/fastq|vcf|sam|bam/)
{
    toolbox::exportLog("ERROR : $0 : The filetype $previousExtension is not taken in charge by TOGGLE\n",0);
}


#Linking the original data to the output dir

#Creating a specific name for the working directory depending on the type of analysis

my $resultsDir = "output";

my $workingDir = $outputDir."/$resultsDir";
toolbox::makeDir($workingDir);

foreach my $file (@{$initialDirContent})
{    
    my ($shortName)=toolbox::extractPath($file);
    my $lnCommand = "ln -s $file $workingDir/$shortName";
    toolbox::run($lnCommand,"noprint");      
}

my $listOfFiles = toolbox::readDir($workingDir);                     # read it to recover files in it

if ($previousExtension eq "fastq")               # the data are all in FASTQ format
{
    #########################################
    # recognition of pairs of files and create a folder for each pair
    #########################################
    my $pairsInfos = pairing::pairRecognition($workingDir);            # from files fasta recognition of paired files
    pairing::createDirPerCouple($pairsInfos,$workingDir);              # from infos of pairs, construction of the pair folder
    
    ##DEBUG toolbox::exportLog("INFOS: $0 : toolbox::readDir : $workingDir after create dir per couple: @$listOfFiles\n",1);
    
}

#Other Data are not always treated singlely, but can work together => check if order hash steps higher than 1000 using the $lastStep value
elsif ($firstOrder<1000) #Other types of data requesting a single treatment
{
    #Create a directory per file
    foreach my $file (@$listOfFiles)
    {
        my ($completeName)=toolbox::extractPath($file);
        my ($noNeed,$basicName)=pairing::extractName($completeName);
        my $dirName=$workingDir."/".$basicName;
        toolbox::makeDir($dirName);
        my $mvCommand = "mv $file $dirName/$completeName";
        if (toolbox::run($mvCommand,"noprint") == 1)
        {
            toolbox::exportLog("INFOS : $0 : Transferring $file to $dirName\n",1);
        }
        
    }
}


#Generation of Index required for the analysis to work (on the reference only)
toolbox::exportLog("#########################################\nINFOS: Generating reference index if requested \n#########################################\n",1);

#Linking of the reference file and of already existing index in the output folder to avoid writing rights limitations
##DEBUG print $refFastaFile,"\n";
my $referenceShortName = $refFastaFile;
my $shortRefFileName = toolbox::extractName($refFastaFile); #We have only the reference name, not the full path
$referenceShortName =~ s/\.\w+$//; #Ref can finish with .fa or .fasta, and we need also .dict file
my $refLsCommand = " ls ".$referenceShortName.".*";
my $refLsResults = `$refLsCommand` or die ("ERROR : $0 : Cannot obtain the list of reference associated files with the command $refLsCommand: $!\n");
chomp $refLsResults;
#Creating a reference Folder in the $outputDir
my $refDir = $outputDir."/referenceFiles";
toolbox::makeDir($refDir);
#Transforming in a list of files
my @listOfRefFiles = split /\n/, $refLsResults;
#Performin a ln command per file
while (@listOfRefFiles)
{
  my $currentRefFile = shift @listOfRefFiles;
  my $shortRefFileName = toolbox::extractName($currentRefFile);
  my $refLsCommand = "cp $currentRefFile $refDir/$shortRefFileName";
  ##DEBUG print $refLsCommand,"\n";
  toolbox::run($refLsCommand,"noprint");
}

#Providing the good reference location 
$refFastaFile = $refDir."/".$shortRefFileName;
##DEBUG print $refFastaFile,"\n";

onTheFly::indexCreator($configInfo,$refFastaFile);

#Generate script
my $scriptSingle = "$outputDir/toggleBzz.pl";
my $scriptMultiple = "$outputDir/toggleMultiple.pl";

my $hashOrder=toolbox::extractHashSoft($configInfo,"order"); #Picking up the options for the order of the pipeline
my $hashCleaner=toolbox::extractHashSoft($configInfo,"cleaner"); #Picking up infos for steps to be cleaned / data to be removed all along the pipeline
my $hashCompressor=toolbox::extractHashSoft($configInfo,"compress"); #Picking up infos for steps to be compress 

my ($orderBefore1000,$orderAfter1000,$lastOrderBefore1000);

foreach my $step (sort {$a <=> $b} keys %{$hashOrder}) #Will create two subhash for the order, to launch twice the generateScript
{
    if ($step < 1000)
    {
        $$orderBefore1000{$step}=$$hashOrder{$step};
        $lastOrderBefore1000 = $step;
    }
    else
    {
        $$orderAfter1000{$step}=$$hashOrder{$step};
    }
}


#########################################
# Launching the generated script on all subfolders if steps lower than 1000
#########################################

#Creating global output folder
my $finalDir = $outputDir."/finalResults";
my $intermediateDir = $workingDir."/intermediateResults";

#Graphviz Graphic generator
toolbox::exportLog("#########################################\nINFOS: Generating graphical view of the current pipeline \n#########################################\n",1);
onTheFly::generateGraphviz($hashOrder,$outputDir);


if ($orderBefore1000)
{
    toolbox::exportLog("\n#########################################\nINFOS: Running individual pipeline script \n#########################################\n",1);

    #generate toggleBzzzz.pl
    onTheFly::generateScript($orderBefore1000,$scriptSingle,$hashCleaner,$hashCompressor);
    my $listSamples=toolbox::readDir($workingDir);
    
    #CORRECTING $listSamples if only one individual, ie readDir will provide only the list of files...
    if (scalar @{$listSamples} < 3) #ex: Data/file_1.fastq, Data/file_2.fastq, or a single SAM/BAM/VCF individual
    {
      my @listPath = split /\//, $$listSamples[0];
      pop @listPath;
      my $trueDirName=join ("/",@listPath);
      $trueDirName .= ":"; 
      my @tempList = ($trueDirName);
      $listSamples = \@tempList;
    }
    my $errorList="obiWanKenobi";
    
    #we need those variable for Scheduler launching
    my $jobList="";
    my %jobHash;
        
    foreach my $currentDir(@{$listSamples})
    {
        next unless $currentDir =~ m/:$/; # Will work only on folders
        $currentDir =~ s/:$//;
        my $launcherCommand="$scriptSingle -d $currentDir -c $fileConf -r $refFastaFile";
        $launcherCommand.=" -g $gffFile" if (defined $gffFile);
        
        #Launching through the scheduler launching system  
        my $jobOutput = scheduler::launcher($launcherCommand, "1", $currentDir, $configInfo); #not blocking job, explaining the '1'
        ##DEBUG        toolbox::exportLog("WARNING: $0 : jobID = $jobOutput -- ",2);
        if ($jobOutput == 0)
        {
          #the linear job is not ok, need to pick up the number of jobs
          my $individualName = `basename $currentDir` or warn("\nERROR: $0 : Cannot pick up basename for $currentDir : $!\n");
          chomp $individualName;
          $individualName = $currentDir unless ($individualName); # Basename did not succeed...

          $errorList.="\$\|".$individualName;
          ##DEBUG          print "++$errorList++\n";
          #Need to remove the empty name...
          $errorList =~ s/obiWanKenobi\$\|//;
          ##DEBUG          print "++$errorList++\n";
        }
        next unless ($jobOutput > 1); #1 means the job is Ok and is running in a normal linear way, ie no scheduling
        
        ##DEBUG        toolbox::exportLog("INFOS: $0 : Parallel job",2);
        
        $jobList = $jobList.$jobOutput."|";
        my $baseNameDir=`basename $currentDir` or die("\nERROR : $0 : Cannot pickup the basename for $currentDir: $!\n");
        chomp $baseNameDir;
        $jobHash{$baseNameDir}=$jobOutput;
    }
    
    
    #If parallel mode, we have to wait the end of jobs before populating
    chop $jobList if ($jobList =~ m/\|$/);
    if ($jobList ne "")
    {
      #Have to wait that all jobs are finished
      my $waitOutput = scheduler::waiter($jobList,\%jobHash);
      if ($waitOutput != 1)
      {
        #Creating a chain with the list of individual with an error in the job...
        $errorList=join ("\$\|",@{$waitOutput}); 
      }

    }
    if ($errorList ne "obiWanKenobi")
      {
        #Some errors appears
        #problem somewhere for some individuals, reporting the info
        my $outputErrors = $errorList;
         ##DEBUG        print "==$errorList==\n";
        $outputErrors =~ s/\$\|/,/;
        toolbox::exportLog("\n>>>>>>>>>>>>>>>> WARNINGS: $0 : Some individuals are erroneous and not treated: $outputErrors\n",2);
      }
    
    # Going through the individual tree
    #Transferring unmapped bam generated by tophat from tempory directory into tophat directory
    foreach my $currentDir (@{$listSamples})
    {
      next unless $currentDir =~ m/\//; # Will work only on folders
      next if $currentDir =~ m/$errorList/; # A job in error will not be transfered, to avoid errors.
      my $fileList = toolbox::readDir($currentDir);
      foreach my $file (@{$fileList}) #Copying intermediate data in the intermediate directory
      {
         
        if ($file =~ "tophatTempory")
        {
          $file =~ s/://g;
          my $mvCommand = "mv ".$file."/* ".$currentDir."/*_tophat*/ && rm -rf $file";
          toolbox::run($mvCommand);
        }
      }
        
    }

 
    #Populationg the intermediate directory
    if ($orderAfter1000) #There is a global analysis afterward
    {
        #Creating intermediate directory
        toolbox::makeDir($intermediateDir);
                
        # Going through the individual tree
        foreach my $currentDir (@{$listSamples})
        {
            next unless $currentDir =~ m/\//; # Will work only on folders
            next if $currentDir =~ m/$errorList/; # A job in error will not be transfered, to avoid errors.

            my $lastDir = $currentDir."/".$lastOrderBefore1000."_".$$orderBefore1000{$lastOrderBefore1000};
            $lastDir =~ s/ //g;
            my $fileList = toolbox::readDir($lastDir);
            foreach my $file (@{$fileList}) #Copying intermediate data in the intermediate directory
            {
                my ($basicName)=toolbox::extractPath($file);
                my $lnCommand="ln -s $file $intermediateDir/$basicName";
                toolbox::run($lnCommand,"noprint")         
            }
        }
    }
    else #There is no global analysis afterward
    {
        #Creating final directory
        toolbox::makeDir($finalDir);
                
        # Going through the individual tree
        foreach my $currentDir (@{$listSamples})
        {
            ##DEBUG toolbox::exportLog($currentDir,1);

            next unless $currentDir =~ m/\//; # Will work only on folders
            next if $currentDir =~ m/$errorList/; # A job in error will not be transfered, to avoid errors.

            my $lastDir = $currentDir."/".$lastOrderBefore1000."_".$$orderBefore1000{$lastOrderBefore1000};
            $lastDir =~ s/ //g;
            ##DEBUG toolbox::exportLog($lastDir,1);
            
            my $fileList = toolbox::readDir($lastDir);
            foreach my $file (@{$fileList}) #Copying the final data in the final directory
            {
                next if (not defined $file or $file =~ /^\s*$/);	
		$file =~s/://g;
		my ($basicName)=toolbox::extractPath($file);
                my $cpLnCommand="cp -rf $file $finalDir/$basicName && rm -rf $file && ln -s $finalDir/$basicName $file";
                toolbox::run($cpLnCommand,"noprint")       
            }
        }
    }
}


if ($orderAfter1000)
{
  
    toolbox::exportLog("\n#########################################\n INFOS: Running multiple pipeline script \n#########################################\n",1);

    onTheFly::generateScript($orderAfter1000,$scriptMultiple,$hashCleaner,$hashCompressor);
    
    $workingDir = $intermediateDir if ($orderBefore1000); # Changing the target directory if we have under 1000 steps before.

    my $launcherCommand="$scriptMultiple -d $workingDir -c $fileConf -r $refFastaFile";
    $launcherCommand.=" -g $gffFile" if (defined $gffFile);
    
    my $jobList="";
    my %jobHash;

    #Launching through the scheduler launching system
    my $jobOutput = scheduler::launcher($launcherCommand, "1", "Global analysis", $configInfo); #not blocking job, explaining the '1'
    if ($jobOutput ne 1) #1 means the job is Ok and is running in a normal linear way, ie no scheduling
    {
      $jobList = $jobOutput;
      $jobHash{"global"}=$jobOutput;
    
      #If qsub mode, we have to wait the end of jobs before populating
      chop $jobList if ($jobList =~ m/\|$/);
      if ($jobList ne "")
      {
        #Have to wait that all jobs are finished
        my $waitOutput = scheduler::waiter($jobList,\%jobHash);
        if ($waitOutput != 1)
        {
          toolbox::exportLog("ERROR: $0 : Multiple job is not finished correctly, please check error log.\n",0);
        }
      }
    }
    
    #Creating final directory
    toolbox::makeDir($finalDir);
            
    # Going through the individual tree
    my $lastDir = $workingDir."/".$lastOrder."_".$$orderAfter1000{$lastOrder};
    $lastDir =~ s/ //g;
    ##DEBUG toolbox::exportLog($lastDir,1);
    my $fileList = toolbox::readDir($lastDir);
    foreach my $file (@{$fileList}) #Copying the final data in the final directory
    {
        my ($basicName)=toolbox::extractPath($file);
        my $cpLnCommand="cp -rf $file $finalDir/$basicName && rm -rf $file && ln -s $finalDir/$basicName $file";
        ##DEBUG toolbox::exportLog($cpLnCommand,1);
        toolbox::run($cpLnCommand,"noprint")       
    }
    
}

close F1;

toolbox::exportLog("#########################################\nINFOS: Analysis correctly done. \n#########################################\n",1);
toolbox::exportLog("\nThank you for using TOGGLE!
###########################################################################################################################
#\tCITATION:
#\tTOGGLE: Toolbox for generic NGS analyses.
#\tCécile Monat, Christine Tranchant-Dubreuil, Ayité Kougbeadjo, Cédric Farcy, Enrique Ortega-Abboud,
#\tSouhila Amanzougarene,Sébastien Ravel, Mawussé Agbessi, Julie Orjuela-Bouniol, Maryline Summo and François Sabot.
#\tBMC Bioinformatics 2015, 16:374
###########################################################################################################################",1);  

exit;

=head1 Name

toggleGenerator.pl - Automatic pipeline generator

=head1 Usage


toggleGenerator.pl -d DIR -c FILE -r FILE -o DIR -g FILE

=head1 Required Arguments

      -d DIR    The directory containing initial files
      -c FILE   The configuration file
      -r FILE   The reference sequence (fasta)
      -o DIR    The directory containing output files
      
=head1 For RNAseq analysis

      -g FILE   The gff file containing reference annotations

=head1  Authors

Cecile Monat, Christine Tranchant, Cedric Farcy, Maryline Summo, Julie Orjuela-Bouniol, Sebastien Ravel, Gautier Sarah, and Francois Sabot

Copyright 2014-2015 IRD-CIRAD-INRA-ADNid

=cut
