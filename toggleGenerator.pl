#!/usr/bin/env perl

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

# Perl Modules
use strict;
use warnings 'all';
no warnings 'experimental';
use Data::Dumper;
use Getopt::ArgParse;
# For gz files
use IO::Compress::Gzip qw(gzip $GzipError);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
# TOGGLE Modules
use localConfig;
use pairing;
use toolbox;
use onTheFly;
use scheduler;
use radseq;
use versionSofts;
use checkFormat;

##########################################
# recovery of parameters/arguments given when the program is executed
##########################################
my $version = "Release 0.3.4, 27th of September, 2017";
my @shortVersion = (3,4);

my $url = "toggle.southgreen.fr/install/releaseNotes/index.html";
my $lastRealease = `curl -m 5 --connect-timeout 5 --max-time 5 -s "$url" 2>&1 | grep -m 1 '<li><a href="\#0' | cut -f3 -d'>' | cut -f1 -d'<'`;
chomp($lastRealease);
my $newRelease="";
if ($lastRealease ne $version)
{
    my $shortRelease = $lastRealease;
    $shortRelease =~ s/Release (\d\.\d\.\d), .*/$1/;
    my ($main,$major,$minor) = split /\./, $shortRelease;
    if ($major < $shortVersion[0])
    {
        $newRelease = "\n** NOTE: This TOGGLE version is higher than the production version, you are using a dev version\n\n";
    }
    if ($major == $shortVersion[0] && $minor < $shortVersion[1])
    {
        $newRelease = "\n** NOTE: This TOGGLE version is higher than the production version, you are using a dev version\n\n";
    }
    else
    {
        $newRelease =  "
** NOTE: Latest version of TOGGLE is $lastRealease, and can be obtained at:
    http://toggle.southgreen.fr/\n\n"
    }

}

my $cmd_line=$0." @ARGV"; # for printing in log file

my $parser = Getopt::ArgParse->new_parser(
        prog            => "\n\ntoggleGenerator.pl",
        description     => '',
        epilog          => "
 $newRelease
##########################################################################
# More information:
#\thttp://toggle.southgreen.fr/
#
# Citation:
#\tTOGGLE: Toolbox for generic NGS analyses. Cécile Monat & al..
#\tBMC Bioinformatics 2015, 16:374
###########################################################################\n",
        help            => 'a framework to build quickly NGS pipelines'."\n\n".$version,
        error_prefix    => "\n\tERROR MSG: "
);

$parser->add_args(
                    [
                        '-o','--outputdir',
                        required => 1,
                        type     =>"Scalar",
                        metavar  => "DIR",
                        help     => 'Output folder name (it will be created)',
                        dest     => 'outputdir'
                    ],
                    [
                        '-c','--config',
                        required => 1,
                        type     =>"Scalar",
                        metavar  => "FILE",
                        help     => 'Software configuration file',
                        dest     => 'config'
                    ],
                    [
                        '-d','--directory',
                        required => 1,
                        type     => "Scalar",
                        metavar  => "DIR",
                        help     => 'Directory name with the files to analyse',
                        dest     => 'directory'
                    ],
                    [
                        '-nocheck','--nocheckFastq',
                        required => 0,
                        type     =>"Bool",
                        help     => 'Use if you did not check fastq file',
                        dest     => 'checkFastq'
                    ],
                    [
                        '-v','--version',
                        required => 0,
                        type     =>"Bool",
                        help     => 'Use if you want to know which version of TOGGLE you are using',
                        dest     => 'version'
                    ],
                    [
                        '-g','--gff',
                        required => 0,
                        type     =>"Scalar",
                        metavar  => "FILE",
                        help     => 'gff file name used by topHat for example',
                        dest     => 'gff',
                        default  => "None"
                    ],
                    [
                        '-k','--keyfile',
                        required => 0,
                        type     =>"Scalar",
                        metavar  => "FILE",
                        help     => 'keyfile file name used by radseq (demultiplexing step)',
                        dest     => 'keyfile',
                        default  => "None"
                    ],
                    [
                        '-r','--reference',
                        required => 0,
                        type     =>"Scalar",
                        metavar  => "FILE",
                        help     => 'Reference file name',
                        dest     => 'reference',
                        default  => "None"
                    ]

                );
# for usage if die print help
my $usage = $parser->format_usage();
my $help = join ("\n", @$usage);

$parser->{"error_prefix"} = $help."\n".$parser->{"error_prefix"};

#recovery supplementary arguments undefined by toggle
my @argv= $parser->argv;

if ("-v" ~~ @ARGV or "--version" ~~ @ARGV or "-version" ~~ @ARGV)
{
    print $version.$newRelease;
    exit;
}

my $args = $parser->parse_args();

#Recovery obligatory arguments
my $initialDir = toolbox::relativeToAbsolutePath($parser->namespace->directory, 0);       # recovery of the name of the directory to analyse
my $fileConf = toolbox::relativeToAbsolutePath($parser->namespace->config, 0);            # recovery of the name of the software.configuration.txt file
my $outputDir = toolbox::relativeToAbsolutePath($parser->namespace->outputdir, 0);        # recovery of the output folder

#Recovery optional arguments
my $refFastaFile = toolbox::relativeToAbsolutePath($parser->namespace->reference, 0);   # recovery of the reference file
my $gffFile = toolbox::relativeToAbsolutePath($parser->namespace->gff, 0);              # recovery of the gff file used by topHat and rnaseq analysis
my $keyfile = toolbox::relativeToAbsolutePath($parser->namespace->keyfile, 0);          # recovery of the keyfile used by radseq

#verify if -nocheckfastq arguments exist in args. The fastq format is verified par default if $checkFastq == 0.
# WARNING with the parser : if nocheckfastq argument is add then $checkFastq == 1
my $checkFastq = $parser->namespace->checkFastq;

#stock mandatory files to test if they exist
my @listFilesMandatory=($initialDir, $fileConf);
push (@listFilesMandatory,$refFastaFile) if $refFastaFile !~ m/None$/;
push (@listFilesMandatory,$gffFile) if $gffFile !~ m/None$/;
push (@listFilesMandatory,$keyfile) if $keyfile !~ m/None$/;

##########################################
# Creation of the output folder
##########################################

if (not -d $outputDir) #if the output folder is not existing yet
{
    #creating the output folder
    my $createOutputDirCommand = "mkdir -p $outputDir";
    system ("$createOutputDirCommand") and toolbox::exportLog("\nERROR: $0 : cannot create the output folder $outputDir: $!\nExiting...\n",0);
}

chdir $outputDir;

#Checking if $outputDir is empty

my $lsOutputDir = `ls`;
chomp $lsOutputDir;
if ($lsOutputDir ne "") # The folder is not empty
{
  toolbox::exportLog("\nERROR: $0 : The output directory $outputDir is not empty, TOGGLE will not continue\nPlease provide an empty directory for outputting results.\n\nExiting...\n\n",0);
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
open (F1, ">",$infosFile) or toolbox::exportLog("$0 : open error of $infosFile .... $!\n",0);
print F1 "GLOBAL\n";
print F1 "ANALYSIS_$date\n";

toolbox::exportLog("#########################################\nINFOS: TOGGLE analysis starts \n#########################################\n",1);;
toolbox::exportLog("INFOS: $0 : Command line : $cmd_line\n",1);
toolbox::exportLog("INFOS: Your output folder is $outputDir\n",1);
toolbox::exportLog("INFOS: the current version of TOGGLE is $version\n",1);

# Verify if file arguments exist
foreach my $file (@listFilesMandatory)
{
  toolbox::checkFile($file); #check file
}

##########################################
# Printing software configurations
########################################

toolbox::exportLog("#########################################\nINFOS: Software version/location \n#########################################\n",1);

versionSofts::writeLogVersion($fileConf,$version.$newRelease);

toolbox::exportLog("\n#########################################\nINFOS: Data checking \n#########################################\n",1);
toolbox::checkFile($fileConf);                              # check if this file exists
toolbox::existsDir($initialDir);                            # check if this directory exists

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
my ($firstOrder,$lastOrder) = onTheFly::checkOrder($configInfo,$refFastaFile,$gffFile,$keyfile);


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

if ($previousExtension !~ m/fasta|fastq|vcf|sam|bam/)  # j'ai rajouté fasta pour les besoins de TGICL
{
    toolbox::exportLog("ERROR : $0 : The filetype $previousExtension is not taken in charge by TOGGLE\n",0);
}


#Linking the original data to the output dir

#Creating a specific name for the working directory depending on the type of analysis

my $resultsDir = "output";

my $workingDir = $outputDir."/$resultsDir";
toolbox::makeDir($workingDir);


#########################################
# check if 1=processRadtags in $order
#########################################
my $hashOrder=toolbox::extractHashSoft($configInfo,"order");					#Picking up the options for the order of the pipeline

my @values;
for my $value ( values %{ $hashOrder } )
{
    push(@values,$value);
}

if ("processRadtags" ~~ @values)												# Check if processRadtags in step order
{
    $initialDirContent = radseq::checkOrder($outputDir,$fileConf,$initialDir,$checkFastq,$keyfile);
    $hashOrder = toolbox::rmHashOrder($hashOrder, "processRadtags")
}
#########################################
# END check if 1=processRadtags in $order
#########################################


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
    my $pairsInfos = pairing::pairRecognition($workingDir,$checkFastq);            # from files fasta recognition of paired files
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
if ($refFastaFile ne 'None')
{
    my $referenceShortName = $refFastaFile;
    my $shortRefFileName = toolbox::extractName($refFastaFile); #We have only the reference name, not the full path
    $referenceShortName =~ s/\.\w+$//; #Ref can finish with .fa or .fasta, and we need also .dict file
    my $refLsCommand = " ls ".$referenceShortName.".*";
    my $refLsResults = `$refLsCommand` or toolbox::exportLog("ERROR : $0 : Cannot obtain the list of reference associated files with the command $refLsCommand: $!\n",0);
    chomp $refLsResults;
    #Creating a reference Folder in the $outputDir
    my $refDir = $outputDir."/referenceFiles";
    toolbox::makeDir($refDir);
    #Transforming in a list of files
    my @listOfRefFiles = split /\n/, $refLsResults;

    my $goodFileFasta = ""; #Providing the good reference location
    while (@listOfRefFiles)
    {
      my $currentRefFile = shift @listOfRefFiles;
      $shortRefFileName = toolbox::extractName($currentRefFile);
      if ($shortRefFileName =~ m/\.fa$/ or $shortRefFileName =~ m/\.fna$/ or $shortRefFileName =~ m/\.fasta$/) #Providing the good reference location
      {
        $goodFileFasta = $shortRefFileName; #Providing the good reference location
      }
      my $refLsCommand = "cp $currentRefFile $refDir/$shortRefFileName";
      ##DEBUG print $refLsCommand,"\n";
      toolbox::run($refLsCommand,"noprint");
    }

    #Providing the good reference location
    $refFastaFile = $refDir."/".$goodFileFasta;
    checkFormat::checkFormatFasta($refFastaFile); # checking format fasta
 
    ##DEBUG print $refFastaFile,"\n";
    onTheFly::indexCreator($configInfo,$refFastaFile);
}
#Generate script
my $scriptSingle = "$outputDir/toggleBzz.pl";
my $scriptMultiple = "$outputDir/toggleMultiple.pl";


my $hashCleaner=toolbox::extractHashSoft($configInfo,"cleaner"); #Picking up infos for steps to be cleaned / data to be removed all along the pipeline
my $hashCompressor=toolbox::extractHashSoft($configInfo,"compress"); #Picking up infos for steps to be compress
my $hashmerge=toolbox::extractHashSoft($configInfo,"merge"); #Picking up infos for steps to be merge

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
    onTheFly::generateScript($orderBefore1000,$scriptSingle,$hashCleaner,$hashCompressor,$hashmerge);
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
        my $launcherCommand="$scriptSingle -d $currentDir -c $fileConf ";
        $launcherCommand.=" -r $refFastaFile" if ($refFastaFile ne 'None');
        $launcherCommand.=" -g $gffFile" if ($gffFile ne 'None');
        $launcherCommand.=" -nocheck" if ($checkFastq == 1);

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
        my $baseNameDir=`basename $currentDir` or toolbox::exportLog("\nERROR : $0 : Cannot pickup the basename for $currentDir: $!\n",0);
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

    onTheFly::generateScript($orderAfter1000,$scriptMultiple,$hashCleaner,$hashCompressor,$hashmerge);

    $workingDir = $intermediateDir if ($orderBefore1000); # Changing the target directory if we have under 1000 steps before.

    my $launcherCommand="$scriptMultiple -d $workingDir -c $fileConf ";
    $launcherCommand.=" -r $refFastaFile" if ($refFastaFile ne 'None');
    $launcherCommand.=" -g $gffFile" if ($gffFile ne 'None');
    $launcherCommand.=" -nocheck" if ($checkFastq == 1);


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

toggleGenerator.pl -d DIR -c FILE -o DIR [optional : -r FILE -g FILE -k FILE -noCheckFastq ]
=head1 Required Obligatoy Arguments :

      -d DIR    	The directory containing initial files
      -c FILE   	The configuration file
      -o DIR    	The directory containing output files

=head1 Optional Arguments :

      -r FILE   	The reference sequence (fasta)
      -g FILE   	The gff file containing reference annotations (For RNAseq analysis per exemple)
      -k FILE		The keyFile used to demultiplexing (For stacks analysis)
      -nocheckfastq 	No check format in every fastq file

=head1  Authors

Cecile Monat, Christine Tranchant, Laura Helou, Abdoulaye Diallo, Julie Orjuela-Bouniol, Sebastien Ravel, Gautier Sarah, and Francois Sabot

Copyright 2014-2015 IRD-CIRAD-INRA-ADNid

=cut
