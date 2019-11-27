#!/usr/bin/env perl

###################################################################################################################################
#
# Copyright 2014-2019 IRD-CIRAD-INRA-ADNid
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
# TOGGLe Modules
use localConfig;
use pairing;
use toolbox;
use onTheFly;
use scheduler;
use radseq;
use checkFormat;
use softwareManagement;

##########################################
# recovery of parameters/arguments given when the program is executed
##########################################
my $TOGGLeVersion = "Release 0.3.8, 27th of November, 2019";
my @shortVersion = (3,8);

my $url = "TOGGLe.southgreen.fr/install/releaseNotes/index.html";
my $lastRealease = `curl -m 5 --connect-timeout 5 --max-time 5 -s "$url" 2>&1 | grep -m 1 '<li><a href="\#0' | cut -f3 -d'>' | cut -f1 -d'<'`;
chomp($lastRealease);
my $newRelease="";
if ($lastRealease ne $TOGGLeVersion)
{
	my $shortRelease = $lastRealease;
	$shortRelease =~ s/Release (\d\.\d\.\d), .*/$1/;
	my ($main,$major,$minor) = split /\./, $shortRelease;
	if ($major < $shortVersion[0])
	{
		$newRelease = "\n** NOTE: This TOGGLe version is higher than the production version, you are using a dev version\n\n";
	}
	if ($major == $shortVersion[0] && $minor < $shortVersion[1])
	{
		$newRelease = "\n** NOTE: This TOGGLe version is higher than the production version, you are using a dev version\n\n";
	}
	else
	{
		$newRelease =  "\nNOTE: The Latest version of TOGGLe ($lastRealease) is available at http://TOGGLe.southgreen.fr/\n\n"
	}

}

my $cmd_line=$0." @ARGV"; # for printing in log file

my $parser = Getopt::ArgParse->new_parser(
		prog            => "\n\nTOGGLeGenerator.pl",
		description     => '',
		epilog          => "
 $newRelease
##########################################################################
# More information:
#\thttp://TOGGLe.southgreen.fr/
#\tCITATION:
#\tTOGGLe, a flexible framework for easily building complex workflows and performing robust large-scale NGS analyses.
#\tChristine Tranchant-Dubreuil, Sebastien Ravel, Cecile Monat, Gautier Sarah, Abdoulaye Diallo, Laura Helou, Alexis Dereeper,
#\tNdomassi Tando, Julie Orjuela-Bouniol, Francois Sabot.
#\tbioRxiv, doi: https://doi.org/10.1101/245480
#\thttps://TOGGLe.southgreen.fr/
###########################################################################\n",
		help            => 'a framework to build quickly NGS pipelines'."\n\n".$TOGGLeVersion,
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
						help     => 'Use if you want not to check fastq file',
						dest     => 'checkFastq'
					],
					[
						'-v','--version',
						required => 0,
						type     =>"Bool",
						help     => 'Use if you want to know which version of TOGGLe you are using',
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
					],
					[
						'-report','--report',
						required => 0,
						type     =>"Bool",
						help     => 'Use if you want to generate workflow an analysis reports',
						dest     => 'report'
					],
										[
						'-add','--add',
						required => 0,
						type     =>"Bool",
						help     => 'Use if you want to add new samples to an already run analysis',
						dest     => 'add'
					],
					[
						'-rerun','--rerun',
						required => 0,
						type     =>"Bool",
						help     => 'Use if you want to re-run samples that have encountered error previously',
						dest     => 'rerun'
					]

				);
# for usage if die print help
my $usage = $parser->format_usage();
my $help = join ("\n", @$usage);

$parser->{"error_prefix"} = $help."\n".$parser->{"error_prefix"};

#recovery supplementary arguments undefined by TOGGLe
my @argv= $parser->argv;

if ("-v" ~~ @ARGV or "--version" ~~ @ARGV or "-version" ~~ @ARGV)
{
	print $TOGGLeVersion.$newRelease;
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
my $report = $parser->namespace->report;

# for reentrancy
my $addSample = $parser->namespace->add;
my $rerun = $parser->namespace->rerun;

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
	system ("$createOutputDirCommand") and die "\nERROR: $0 : cannot create the output folder $outputDir: $!\nExiting...\n";
}
chdir $outputDir;

#Checking if $outputDir is empty
my $lsOutputDir = `ls`;
chomp $lsOutputDir;
my $reentrancyLog = "";
if ($lsOutputDir ne "") # The folder is not empty
{
	if ($addSample || $rerun)
	{
		$reentrancyLog = "\nWARN: $0 : cannot create the output folder $outputDir: $!\nProbably you want to add samples or re-run the jobs.\n Every logs and reports will be backuped and new ones will be added (if requested).\n";
		#Backuping pre data
		my @fileList = split /\n/, $lsOutputDir;
		while (@fileList)
		{
			my $currentFile = shift @fileList;
			if ($currentFile =~ m/log\./ || $currentFile =~ m/Report/)
			{
				#Old log and report files
				my $newName = "OLD_".$currentFile."OLD";
				my $mvCom = "mv $outputDir/$currentFile $outputDir/$newName";
				my $mvResult= toolbox::run($mvCom, "noprint");
			}
		}
	 }
	 else
	 {
			die "\nERROR: $0 : The output directory $outputDir is not empty, TOGGLe will not continue\nPlease provide an empty directory for outputting results.\n\nExiting...\n\n";
	 }
}


# Creating log file GLOBAL
my ($sec, $min, $h, $mois_jour, $mois, $an, $sem_jour, $cal_jour, $heure_ete) = localtime(time);
$mois+=1;
$mois = $mois < 10 ? $mois = "0".$mois : $mois;
$mois_jour = $mois_jour < 10 ? $mois_jour = "0".$mois_jour : $mois_jour;
$h = $h < 10 ? $h = "0".$h : $h;
$min = $min < 10 ? $min = "0".$min : $min;
$an+=1900;
my $date="$mois_jour-$mois-$an-$h"."_"."$min";

my $logFile=$outputDir."/TOGGLe_".$date."_log.o";
my $errorFile=$outputDir."/TOGGLe_".$date."_log.e";
system("touch $logFile $errorFile") and die "\nERROR: $0 : cannot create the log files $logFile and $errorFile: $!\nExiting...\n";


toolbox::exportLog("#########################################\nINFOS: TOGGLe analysis starts \n#########################################\n",1);;
toolbox::exportLog("INFOS: $0 : Command line : $cmd_line\n",1);
toolbox::exportLog("INFOS: Your output folder is $outputDir\n",1);
#toolbox::exportLog($reentrancyLog,1) if $reentrancyLog;

### Generate tex file
if ($report)
{
	my $inputFile="commandLine.tex";
	open(my $cmdFh,">", $inputFile) or toolbox::exportLog("$0 : open error of $inputFile .... $!\n",0);
	print $cmdFh "\n $cmd_line \n";
	print $cmdFh "\n $reentrancyLog \n" if $reentrancyLog;
	close $cmdFh;
}

##########################################
# Check directories
##########################################
toolbox::exportLog("\n#########################################\nINFOS: Data checking \n#########################################\n",1);
# Verify if file arguments exist
foreach my $file (@listFilesMandatory)
{
 toolbox::checkFile($file); #check file
}
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
#last true order is the really last step value, even if NA...
my ($firstOrder,$lastOrder,$lastTrueOrder) = onTheFly::checkOrder($configInfo,$refFastaFile,$gffFile,$keyfile);




##########################################
# Transferring data in the output folder and organizing
#########################################

# check that the input dir has no subdirectories
toolbox::checkInitialDirContent($initialDir);

# check that all of the input files are of the same type (homogeneity)
my $initialDirContent=toolbox::readDir($initialDir);
toolbox::controlReadGroup($initialDirContent);

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
		toolbox::exportLog("ERROR : $0 : The compressed file $file format in the initial directory is not taken in charge by TOGGLe.\n",0);
	  }
	}
	#The file is not a compressed one in gz
	if ($extension eq "fq" or $extension eq "fastq") #homogeneisation for fastq extension
	{
		$extension = "fastq";
	}
	if ($extension eq "fa" or $extension eq "fasta" or $extension eq "fna") #homogeneisation for fasta extension
	{
		$extension = "fasta";
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

#Checking if the files are taken in charge by TOGGLe
if ($previousExtension !~ m/fasta|fastq|vcf|sam|bam|ped|gtf|txt/)  # j'ai rajouté fasta pour les besoins de TGICL et moi ped pour SNIPLAY; GTF for stringtie (option merge)
{
	toolbox::exportLog("ERROR : $0 : The filetype $previousExtension is not taken in charge by TOGGLe\n",0);
}

##########################################
# Printing software version
########################################
my $hashOrder=toolbox::extractHashSoft($configInfo,"order");					#Picking up the options for the order of the pipeline

toolbox::exportLog("#########################################\nINFOS: Software version/location \n#########################################\n",1);
softwareManagement::writeLogVersion($hashOrder,$TOGGLeVersion.$newRelease,$outputDir,$report);

#########################################
# check if 1=processRadtags in $order
#########################################
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

# Determine the last order before 1000 to be used with rerun and generateScripts

my ($orderBefore1000,$orderAfter1000,$lastOrderBefore1000,$lastTrueOrderBefore1000);
# lastTrueOrderBefore1000 takes NA steps into account, lastOrderBefore1000 doesn't

#Initializing value of lastOrder in case of only NA OUT pipeline
$lastOrderBefore1000=1;

#Obtaining infos for OUT NA steps
my $hashInOut= softwareManagement::returnSoftInfos();

foreach my $step (sort {$a <=> $b} keys %{$hashOrder}) #Will create two subhash for the order, to launch twice the generateScript
{
	if ($step < 1000)
	{
		$$orderBefore1000{$step}=$$hashOrder{$step};
		$lastTrueOrderBefore1000 = $step;
		$lastOrderBefore1000 = $step unless $hashInOut->{$$hashOrder{$step}}{"OUT"} eq "NA"; # the last step will be everything but a dead-end one.
	}
	else
	{
		$$orderAfter1000{$step}=$$hashOrder{$step};
	}
}

#########################################
# Check working dir
#########################################
#Linking the original data to the output dir
#Creating a specific name for the working directory depending on the type of analysis
my $resultsDir = "output";

my $workingDir = $outputDir."/$resultsDir";

# If we are adding new samples, first list all the existing ones
# TODO For rerun, we might want to check the samples that are already complete here and also add them to the list
my @alreadyRun = ();
if ($addSample || $rerun)
{
	# TODO potential bug : do we check somewhere that $workingDir exists ?
	my $files = `ls $workingDir`;
	chomp $files;
	@alreadyRun = split /\n/, $files;
}
else
{
	toolbox::makeDir($workingDir);
}


# RERUN MODE
# Go through the list of samples which already have a have a working dir and see which ones need to be rerun or run
my @rerunSamples = (); # Samples that need to be lauched with --rerun
if ($rerun)
{
	toolbox::exportLog("#########################################\nINFOS: Rerun mode : Checking the status of individual samples in the output folder \n#########################################\n",1);
	my @newAlreadyRun = (); # alreadyRun will contain only the samples that have already finished
	foreach my $sample (@alreadyRun)
	{
		next if ($sample eq "intermediateResults"); # Skip the global analysis folder, it's not an individual sample
		next if ($sample eq "OLD_intermediateResults");

		# If we are in rerun mode, we have to determine what state this sample was left in
		# We just know that the folder for this sample exists

		# If the 0_initialFiles folder does not exist, toggleBzz has not run at all for this sample
		# If the step file does not exist, not a single step has completed succesfully
		# In both cases, we simply empty the folder and treat it like a new sample

		my $stepFilePath = "$workingDir/$sample/.".$sample."_last_step";
		if (! -d "$workingDir/$sample/0_initialFiles" || ! -f $stepFilePath)
		{
			toolbox::exportLog("$sample : not started correctly (will be run as a new sample)", 1);
			system("rm -r $workingDir/$sample/");
			next;
		}

		open (my $stepFileHandle, '<', $stepFilePath) or die "ERROR : $0 : cannot open file $stepFilePath. $!\nExiting...\n";
		my $lastOkStep = <$stepFileHandle>;

		if ($lastOkStep == $lastTrueOrderBefore1000)
		{
			toolbox::exportLog("$sample : finished succesfully", 1);
			# The individual analysis has finished
			push @newAlreadyRun, $sample;
		}
		else
		{
			toolbox::exportLog("$sample : will rerun", 1);
			# The individual analysis will rerun
			push @rerunSamples, $sample;
		}

	}

	@alreadyRun = @newAlreadyRun;
	toolbox::exportLog("\n", 1);
}
########END of rerun mode###############

my @listOfFiles; #list of files (symbolic links of samples (path to pairing))
my @listSamplesRun; #list of directory samples to run ifadd
my @listAllSamples; # list of all directory samples (already run and added)

foreach my $file (@{$initialDirContent})
{
	my ($shortName)=toolbox::extractPath($file); # name of the file i.e irigin1_1.fastq
	my ($name) = split /_/, $shortName; # i.e irigin1
 if ($name eq $shortName) #to correct list of samples
	{
		($name)= split /\./, $shortName; 	
	}
	
	#if (scalar(keys %{$orderBefore1000})<1)
	#{
	#	$name="globalAnalysis"
	#}
	
	#$$orderBefore1000{$step}=$$hashOrder{$step};
	next if ($name eq "intermediateResults");
	next if ($name eq "OLD_intermediateResults");

	if ($name ~~ @alreadyRun)
	{
		# populating array containing all directory of samples if add or rerun
		push(@listAllSamples, "$workingDir/$name") if (!("$workingDir/$name"  ~~ @listAllSamples));
	}
	else
	{
		# populating array @listOfFiles, @listSamplesRun and @listAllSamples
		# only add files to listOfFiles (which is used for the pairing process) if the sample is new (otherwise the folder is already in place)
		unless ($name ~~ @rerunSamples)
		{
			my $lnCommand = "ln -s $file $workingDir/$shortName";
			toolbox::run($lnCommand,"noprint");
			push(@listOfFiles, "$workingDir/$shortName");
		}
		push(@listSamplesRun, "$workingDir/$name") if (!("$workingDir/$name"  ~~ @listSamplesRun));
		push(@listAllSamples, "$workingDir/$name") if (!("$workingDir/$name"  ~~ @listAllSamples));
	}
}

if ($previousExtension eq "fastq")               # the data are all in FASTQ format
{
	#########################################
	# recognition of pairs of files and create a folder for each pair
	#########################################
	my $pairsInfos = pairing::pairRecognition(\@listOfFiles,$checkFastq);            # from files fasta recognition of paired files
	pairing::createDirPerCouple($pairsInfos,$workingDir);              # from infos of pairs, construction of the pair folder. Also move the files in the folders.

	##DEBUG toolbox::exportLog("INFOS: $0 : toolbox::readDir : $workingDir after create dir per couple: @$listOfFiles\n",1);

}

#Other Data are not always treated singlely, but can work together => check if order hash steps higher than 1000 using the $lastStep value
elsif ($firstOrder<1000) #Other types of data requesting a single treatment
{
	#Create a directory per file
	foreach my $file (@listOfFiles)
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
	toolbox::makeDir($refDir) if (!($addSample || $rerun));
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
	checkFormat::checkFormatFasta($refFastaFile) if (!($addSample || $rerun)); # checking format fasta

	##DEBUG print $refFastaFile,"\n";
	onTheFly::indexCreator($configInfo,$refFastaFile) if (!($addSample || $rerun));
}


############################
#Generate script
############################

my $scriptSingle = "$outputDir/TOGGLeBzz.pl";
my $scriptMultiple = "$outputDir/TOGGLeMultiple.pl";


my $hashCleaner=toolbox::extractHashSoft($configInfo,"cleaner"); #Picking up infos for steps to be cleaned / data to be removed all along the pipeline
my $hashCompressor=toolbox::extractHashSoft($configInfo,"compress"); #Picking up infos for steps to be compress
my $hashmerge=toolbox::extractHashSoft($configInfo,"merge"); #Picking up infos for steps to be merge


#########################################
# Launching the generated script on all subfolders if steps lower than 1000
#########################################

#Creating global output folder
my $finalDir = $outputDir."/finalResults";
my $intermediateDir = $workingDir."/intermediateResults";
my $name="";


#Creating  directory
my $statDir = $outputDir."/statsReport";

toolbox::makeDir($statDir) if $report;

#Graphviz Graphic generator
toolbox::exportLog("#########################################\nINFOS: Generating graphical view of the current pipeline \n#########################################\n",1);
onTheFly::generateGraphviz($hashOrder,$outputDir);

if ($orderBefore1000)
{
	toolbox::exportLog("\n#########################################\nINFOS: Running individual pipeline script \n#########################################\n",1);

	#generate TOGGLeBzzzz.pl
	onTheFly::generateScript($orderBefore1000,$scriptSingle,$hashCleaner,$hashCompressor,$hashmerge) if (!($addSample || $rerun));

	my $errorList="obiWanKenobi";

	#we need those variable for Scheduler launching
	my $jobList="";
	my %jobHash;
	foreach my $currentDir(@listSamplesRun)
	{
				my $individualName = `basename $currentDir` or toolbox::exportLog("\nERROR : $0 : Cannot pickup the basename for $currentDir: $!\n",0);
				chomp $individualName;
		
				my $launcherCommand="$scriptSingle -d $currentDir -c $fileConf ";
				$launcherCommand.=" -r $refFastaFile" if ($refFastaFile ne 'None');
				$launcherCommand.=" -g $gffFile" if ($gffFile ne 'None');
				$launcherCommand.=" -nocheck" if ($checkFastq == 1);
				$launcherCommand.=" -report" if ($report);
				$launcherCommand.=" -rerun" if ($rerun && $individualName ~~ @rerunSamples);
		
				#Launching through the scheduler launching system
				my ($jobOutput, $errorFile) = scheduler::launcher($launcherCommand, "1", $currentDir, $configInfo); #not blocking job, explaining the '1'
				##DEBUG        toolbox::exportLog("WARNING: $0 : jobID = $jobOutput -- \nerrorFile = $errorFile",2);
				if ($jobOutput eq 0)
				{
						#the linear job is not ok, need to pick up the number of jobs
						$errorList.="\$|".$individualName;
						##DEBUG          print "++$errorList++\n";
						#Need to remove the empty name...
						$errorList =~ s/obiWanKenobi\$\|//;
						##DEBUG          print "++$errorList++\n";
				}
				#next unless ($jobOutput > 1); #1 means the job is Ok and is running in a normal linear way, ie no scheduling
		
				##DEBUG        toolbox::exportLog("INFOS: $0 : Parallel job",2);
		
				$jobList = $jobList.$jobOutput."|";
				$jobHash{$individualName}{output}=$jobOutput;
				$jobHash{$individualName}{errorFile}=$errorFile;
	}


	#If parallel mode, we have to wait the end of jobs before populating
	chop $jobList if ($jobList =~ m/\|$/);
	if ($jobList ne "")
	{
	  #Have to wait that all jobs are finished
	  my $waitOutput = scheduler::waiter($jobList,\%jobHash, $outputDir, $report);
	  if ($waitOutput != 1)
	  {
		#Creating a chain with the list of individual with an error in the job...
		$errorList=join ("\$|",@{$waitOutput});
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
	foreach my $currentDir (@listSamplesRun)
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
		# back up the previous directory if it exists
		if (($addSample || $rerun) && -d $intermediateDir)
		{
			my $oldIntermediateDir = "$workingDir/OLD_intermediateResults";
			my $mvCom = "mv $intermediateDir $oldIntermediateDir";
			toolbox::run($mvCom, "noprint");
			toolbox::exportLog("INFOS : Backing up previous intermediateResults directory to $oldIntermediateDir\n", 1);
		}

		#Creating intermediate directory
		toolbox::makeDir($intermediateDir);

		# Going through the individual tree
		foreach my $currentDir (@listAllSamples)
		{
			next unless $currentDir =~ m/\//; # Will work only on folders
			next if $currentDir =~ m/$errorList/; # A job in error will not be transfered, to avoid errors.

			my $lastDir = $currentDir."/".$lastOrderBefore1000."_".$$orderBefore1000{$lastOrderBefore1000};
			$lastDir =~ s/ //g;
			my $fileList = toolbox::readDir($lastDir);
			foreach my $file (@{$fileList}) #Copying intermediate data in the intermediate directory
			{
				my ($basicName)=toolbox::extractPath($file);
				##DEBUG toolbox::exportLog("------>".$basicName."\n",1);

				# Moving stat file into stat directory instead of intermediateDir
				if ($basicName =~ /\.stat$/)
				{
					my $mvCommand="mv $file $statDir/$basicName && rm -f $file";
					toolbox::run($mvCommand,"noprint");
				}
				else
				{
					my $lnCommand="ln -s $file $intermediateDir/$basicName";
					toolbox::run($lnCommand,"noprint");
				}
			}
		}
	}
	else #There is no global analysis afterward
	{
		#Creating final directory
		toolbox::makeDir($finalDir);

		# Going through the individual tree
		foreach my $currentDir (@listSamplesRun)
		{
			##DEBUG toolbox::exportLog($currentDir,1);

			next unless $currentDir =~ m/\//; # Will work only on folders
			next if $currentDir =~ m/$errorList/; # A job in error will not be transfered, to avoid errors.

			my $lastDir = $currentDir."/".$lastTrueOrderBefore1000."_".$$orderBefore1000{$lastTrueOrderBefore1000};
			$lastDir =~ s/ //g;
			##DEBUG toolbox::exportLog($lastDir,1);

			my $fileList = toolbox::readDir($lastDir);
			foreach my $file (@{$fileList}) #Copying the final data in the final directory
			{
				next if (not defined $file or $file =~ /^\s*$/);
				$file =~s/://g;
				my ($basicName)=toolbox::extractPath($file);

				# Moving stat file into stat directory instead of finalDir
				if ($basicName =~ /\.stat$/)
				{
					my $mvCommand="mv $file $statDir/$basicName && rm -f $file";
					toolbox::run($mvCommand,"noprint");
				}
				else
				{
					my $cpLnCommand="cp -rf $file $finalDir/$basicName && rm -rf $file && ln -s $finalDir/$basicName $file";
					toolbox::run($cpLnCommand,"noprint");
				}
			}
		}
	}
}


if ($orderAfter1000)
{

	toolbox::exportLog("\n#########################################\n INFOS: Running multiple pipeline script \n#########################################\n",1);

	onTheFly::generateScript($orderAfter1000,$scriptMultiple,$hashCleaner,$hashCompressor,$hashmerge) unless (!$rerun && $addSample);

if ($firstOrder>=1000)
{
		$intermediateDir = $outputDir."/output/globalAnalysis";
		#my $dirName=$workingDir."/output/globalAnalysis";
		toolbox::makeDir($intermediateDir);
		my $mvCommand = "cd $outputDir/output && mv *.* globalAnalysis/.";
	if (toolbox::run($mvCommand,"noprint") == 1)
	{
		toolbox::exportLog("INFOS : $0 : Transferring all files to $intermediateDir\n",1);
	}
}
	
	$workingDir = $intermediateDir; # if ($orderBefore1000); # Changing the target directory if we have under 1000 steps before.

	my $launcherCommand="$scriptMultiple -d $workingDir -c $fileConf ";
	$launcherCommand.=" -r $refFastaFile" if ($refFastaFile ne 'None');
	$launcherCommand.=" -g $gffFile" if ($gffFile ne 'None');
	$launcherCommand.=" -nocheck" if ($checkFastq == 1);
	$launcherCommand.=" -report" if ($report);


	my $jobList="";
	my %jobHash;

	#Launching through the scheduler launching system
	my ($jobOutput, $errorFile) = scheduler::launcher($launcherCommand, "1", "Global analysis", $configInfo); #not blocking job, explaining the '1'
	if ($jobOutput ne 1) #1 means the job is Ok and is running in a normal linear way, ie no scheduling
	{
	  $jobList = $jobOutput;
	  $jobHash{"global"}{output}=$jobOutput;
	  $errorFile =~ s/intermediateResults_/intermediateResults_/;
	  $jobHash{"global"}{errorFile}=$errorFile;

	  #If qsub mode, we have to wait the end of jobs before populating
	  chop $jobList if ($jobList =~ m/\|$/);
	  if ($jobList ne "")
	  {
		#Have to wait that all jobs are finished
		my $waitOutput = scheduler::waiter($jobList,\%jobHash, $outputDir, $report);
		if ($waitOutput != 1)
		{
		  toolbox::exportLog("ERROR: $0 : Multiple job is not finished correctly, please check error log $errorFile.\n",0);
		}
	  }
	}

	# backup existing finalResults if it exists
	if (($addSample || $rerun) && -d $finalDir)
	{
		my $oldFinalDir = "$outputDir/OLD_finalResults";
		my $mvCom = "mv $finalDir $oldFinalDir";
		toolbox::run($mvCom, "noprint");
		toolbox::exportLog("INFOS : Backing up previous finalResults directory to $oldFinalDir\n", 1);
	}

	#Creating final directory
	toolbox::makeDir($finalDir);

	# Going through the individual tree
	my $lastDir = $workingDir."/".$lastTrueOrder."_".$$orderAfter1000{$lastTrueOrder};
	$lastDir =~ s/ /_/g;
	##DEBUG toolbox::exportLog($lastDir,1);
	my $fileList = toolbox::readDir($lastDir);
	foreach my $file (@{$fileList}) #Copying the final data in the final directory
	{
		my ($basicName)=toolbox::extractPath($file);

						# Moving stat file into stat directory instead of finalDir
		if ($basicName =~ /\.stat$/)
		{
			my $mvCommand="mv $file $statDir/$basicName && rm -f $file";
			toolbox::run($mvCommand,"noprint");
		}
		else
		{
			my $cpLnCommand="cp -rf $file $finalDir/$basicName && rm -rf $file && ln -s $finalDir/$basicName $file";
			##DEBUG toolbox::exportLog($cpLnCommand,1);
			toolbox::run($cpLnCommand,"noprint")
		}
	}

}


toolbox::exportLog("#########################################\nINFOS: Analysis correctly done. \n#########################################\n",1);
toolbox::exportLog("\nThank you for using TOGGLe!
###########################################################################################################################
#\tCITATION:
#\tTOGGLe, a flexible framework for easily building complex workflows and performing robust large-scale NGS analyses.
#\tChristine Tranchant-Dubreuil, Sebastien Ravel, Cecile Monat, Gautier Sarah, Abdoulaye Diallo, Laura Helou, Alexis Dereeper,
#\tNdomassi Tando, Julie Orjuela-Bouniol, Francois Sabot.
#\tbioRxiv, doi: https://doi.org/10.1101/245480
#\thttps://TOGGLe.southgreen.fr/
###########################################################################################################################",1);

onTheFly::generateReports($outputDir, $fileConf) if $report;

exit;

=head1 Name

TOGGLeGenerator.pl - Automatic pipeline generator

=head1 Usage

TOGGLeGenerator.pl -d DIR -c FILE -o DIR [optional : -r FILE -g FILE -k FILE -noCheckFastq ]
=head1 Required Obligatoy Arguments :

	  -d DIR    	The directory containing initial files
	  -c FILE   	The configuration file
	  -o DIR    	The directory containing output files

=head1 Optional Arguments :

	  -r FILE   	The reference sequence (fasta)
	  -g FILE   	The gff file containing reference annotations (For RNAseq analysis per exemple)
	  -k FILE		The keyFile used to demultiplexing (For stacks analysis)
	  -nocheckfastq 	No check format in every fastq file
	  -report    	Generate workflow and analysis reports

=head1  Authors

Christine Tranchant-Dubreuil, Sebastien Ravel, Cecile Monat, Gautier Sarah, Abdoulaye Diallo, Laura Helou, Alexis Dereeper,
Ndomassi Tando, Julie Orjuela-Bouniol, Francois Sabot.

Copyright 2014-2018 IRD-CIRAD-INRA-ADNid

=cut
