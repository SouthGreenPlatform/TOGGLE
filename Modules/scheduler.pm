package scheduler;

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
use Switch;

##################################
#
#LAUNCHING
#
##################################

our ($commandLine, $requirement, $sample, $configInfo, $jobList, %jobHash);

#Requirement means the job needs to be achieved for the next steps or it is not a blocking job?
# A zero value (0) means mandatory
# Any non-zero value means not blocking job

sub launcher { #Global function for launching, will recover the command to be launch and will re-dispatch to normal or other scheduler
    
    ($commandLine,$requirement, $sample, $configInfo) = @_;
    
    #Picking up sample name
    
    $sample=`basename $sample` or die("ERROR : $0 : Cannot pickup the basename for $sample: $!\n");
    chomp $sample;
    
    my $hashCapability = &checkingCapability;
    
    my $runOutput;
       
    switch (1)
    {
        case ($hashCapability->{"sge"} && defined $$configInfo{"sge"}){$runOutput = &sgeRun} #For SGE running
        case ($hashCapability->{"slurm"} && defined $$configInfo{"slurm"}){$runOutput = &slurmRun} #For SLURM running
	case ($hashCapability->{"mprun"} && defined $$configInfo{"mprun"}){$runOutput = &mprunRun} #For mprun running
	case ($hashCapability->{"lsf"} && defined $$configInfo{"lsf"}){$runOutput = &lsfRun} #For lsf running
        
        #If no scheduler available or configurated in the config info file, let's run it in a normal way
        else {$runOutput = &normalRun};
    }
    
    ##DEBUG    toolbox::exportLog("WARNING : scheduler : run output is $runOutput",2);
    
    if ($runOutput == 0 && $requirement == 0)
    {
	#The job has to succeed either it will kill all other jobs
        toolbox::exportLog("ERROR: scheduler::launcher on $sample: ".$commandLine."\nThe job cannot be achieved and is mandatory, thus the whole analysis is stop\n",0);
    
	return 0;
    }
    
    return $runOutput;
}

sub checkingCapability { #Will test the capacity of launching using various schedulers on the current system
    
    my %capabilityValue;
    
    #SGE test
    $capabilityValue{"sge"} = `qsub -help 2>/dev/null | grep usage`; #Will provide a not-empty output if SGE is installed
    chomp $capabilityValue{"sge"};
    
    #SLURM test
    $capabilityValue{"slurm"}=`sbatch -V 2>/dev/null | grep slurm`; #Will provide a not-empty output if SLURM is installed
    chomp $capabilityValue{"slurm"};
    
    #mprun test
    $capabilityValue{"mprun"}=`ccc_mprun -h 2>&1 | grep usage`; #Will provide a not-empty output if mprun is installed
    chomp $capabilityValue{"mprun"};
    
    #lsf test
    $capabilityValue{"lsf"}=`bsub -h 2>&1 | grep -i Synopsis`; #Will provide a not-empty output if lsf is installed
    chomp $capabilityValue{"lsf"};
    
    
    #Returning infos as a reference
    return \%capabilityValue;
}


sub normalRun { #For running in normal mode, ie no scheduler
    
    #BASED ON TOOLBOX::RUN, but will not stop the whole pipeline for an error
    use Capture::Tiny qw(capture);
        
    toolbox::exportLog("INFOS: scheduler::normalRun : $commandLine\n",1);
    
    ##Execute the command
    my ($result,$stderr)=capture {` $commandLine `};
    
    ##Log export according to the error
    if ($?==0) #Success, no error
    {
	##DEBUG toolbox::exportLog("INFOS: scheduler::normalRun on $sample: ".$result."\n--".$stderr."\n",1);
	return 1;
    }
    else  #Error, the job cannot be achieved for any reason
    {
	##DEBUG
	toolbox::exportLog("WARNING: scheduler::normalRun on $sample: ".$result."\n--".$stderr."\nThe $sample data set has provoked and error, and was not analyzed anymore.\n",2);
	return 0;
    }
    
}

sub sgeRun
{ #For SGE cluster, running using qsub
    
    my $sgeOptionsHash=toolbox::extractHashSoft($configInfo,"sge");
    my $sgeOptions=toolbox::extractOptions($sgeOptionsHash);
    
    #Adding qsub options
    my $launcherCommand = "qsub ".$sgeOptions." ".$commandLine;
    $launcherCommand =~ s/ +/ /g; #Replace multiple spaces by a single one, to have a better view...
    $launcherCommand =~ s/ =/=/g; #Replace espace equal to accept -l mem_free=10G SGE option for exemple ...

    my $currentJID = `$launcherCommand`;
    
    if ($!) #There are errors in the launching...
    {
        toolbox::exportLog ("WARNING : $0 : Cannot launch the job for $sample: $!\n",2);
        $currentJID = "";
    }
        
    #Parsing infos and informing logs
    chomp $currentJID;
    ##DEBUG    toolbox::exportLog("INFOS: $0 : $currentJID \n\n",2);
    
    unless ($currentJID) #The job has no output in STDOUT, ie there is a problem...
    {
        return 0; #Returning to launcher subprogram the error type
    }
    
    toolbox::exportLog("INFOS: $0 : Correctly launched for $sample in qsub mode through the command:\n\t$launcherCommand\n",1);
    
    ## DEBUG    toolbox::exportLog("INFOS: $0 : Output for the command is $currentJID\n\n",2);
    
    my @infosList=split /\s/, $currentJID; #the format is such as "Your jobID ("NAME") has been submitted"
    $currentJID = $infosList[2];
    
    ## DEBUG    toolbox::exportLog("INFOS: $0 : jobID is $currentJID\n\n",2);
    
    ##THIS PART IS ONLY FOR DEBUGGING  -  MUST STAY COMMENTED OR MAY PROVOKE TROUBLES!!
    
#    my $runningNodeCommand="qstat | grep $currentJID";
#    my $runningNode="x";
#    my $trying=0;
#    my $stopWaiting = 0;
#    while ($stopWaiting == 0) #If the job is not yet launched or already finished
#    {
#        sleep 3;#Waiting for the job to be launched
#        $runningNode=`$runningNodeCommand`;
#        chomp $runningNode;
#        $runningNode = "x" unless $runningNode; #if empty variable, problem after...
#        if ($runningNode !~ /\s+r\s+/)
#        {# not running yet
#            $trying++;
#	    ## DEBUG	    toolbox::exportLog("WARNING : $trying trys for $currentJID.",2);
#        }
#	if ($trying > 4)
#        {
#            #We already tryed to pick up the node infos 5 times, let's stop
#            $runningNode = "Still unknown (either not running, or already finished)";
#            ## DEBUG toolbox::exportLog("WARNING : $0 : Cannot pickup the running node for the job $currentJID: $!\n",2);
#            $stopWaiting = 1;
#        }
#	next unless $runningNode =~ /\s+r\s+/; # Next if the job is still not running 
#        #my @runningFields = split /\s+/,$runningNode; #To obtain the correct field
#        #$runningNode = $runningFields[7];
#        #$runningNode =~ s/.+@//;#removing queue name providing only node name
#	$stopWaiting = 1;
#    }
    ## DEBUG    toolbox::exportLog("INFOS: $0 : Running node for job $currentJID is $runningNode\n\n",2);
    
    ## DEBUG    toolbox::exportLog("INFOS: $0 : jobID is $currentJID\n\n",2);
    
    return $currentJID;
}

sub slurmRun
{ #for SLURM cluster, running using sbatch
    
    my $slurmOptionsHash=toolbox::extractHashSoft($configInfo,"slurm");
    my $slurmOptions=toolbox::extractOptions($slurmOptionsHash);
    
     #Picking up ENV variable
    my $envOptionsHash=toolbox::extractHashSoft($configInfo,"env");
    my $envOptions=toolbox::extractOptions($envOptionsHash,"=","\n");
    
    #Adding slurm options
    my $launcherCommand = "sbatch ".$slurmOptions;
    my $date =`date +%Y_%m_%d_%H_%M_%S`;
    chomp $date;
    #Creating the bash script for slurm to launch the command
    my $scriptName="~/slurmScript_".$date.".sh";
    my $bashScriptCreationCommand= "echo \"#!/bin/bash\n".$envOptions."\n".$commandLine."\nexit 0;\" | cat - > $scriptName && chmod 777 $scriptName";
    toolbox::run($bashScriptCreationCommand);
    ##DEBUG toolbox::exportLog("INFOS : $0 : Created the slurm bash file",1);
    $launcherCommand.=" $scriptName";
    $launcherCommand =~ s/ +/ /g; #Replace multiple spaces by a single one, to have a better view...
    ##DEBUG toolbox::exportLog($launcherCommand,2);
    my $currentJID = `$launcherCommand`;
    
    if ($!) #There are errors in the launching...
    {
        warn ("WARNING : $0 : Cannot launch the job for $sample: $!\n");
        $currentJID = "";
    }
        
    #Parsing infos and informing logs
    chomp $currentJID;
    
    unless ($currentJID) #The job has no output in STDOUT, ie there is a problem...
    {
        return 0; #Returning to launcher subprogram the error type
    }
    
    toolbox::exportLog("INFOS: $0 : Correctly launched for $sample in sbatch mode through the command:\n\t$launcherCommand\n",1);
    ## DEBUG toolbox::exportLog("INFOS: $0 : Output for the command is $currentJID\n\n",1);
    
    my @infosList=split /\s/, $currentJID; #the format is such as "Submitted batch job 3"
    $currentJID = $infosList[3];

    ##THIS PART IS ONLY FOR DEBUGGING  -  MUST STAY COMMENTED OR MAY PROVOKE TROUBLES!!

    #my $runningNodeCommand="squeue -u \$USER | grep \"  $currentJID \"";
    #my $runningNode="x";
    #my $trying=0;
    #while ($runningNode ne "R") #If the job is not yet launched or already finished
    #{
    #    sleep 3;#Waiting for the job to be launched
    #    $runningNode=`$runningNodeCommand`;
    #    chomp $runningNode;
    #    $runningNode = "x" unless $runningNode; #if empty variable, problem after...
    #    if ($runningNode !~ /\s+R\s+/)
    #    {# not running yet
    #        $trying++;
    #        if ($trying == 5)
    #        {
    #            #We already tryed to pick up the node infos 5 times, let's stop
    #            $runningNode = "still unknown (either not running, or already finished)";
    #            ## DEBUG toolbox::exportLog("WARNING : $0 : Cannot pickup the running node for the job $currentJID: $!\n",2);
    #            last;
    #        }
    #        next;
    #    }
    #    my @runningFields = split /\s+/,$runningNode; #To obtain the correct field
    #    $runningNode = $runningFields[8];
    #}
    ## DEBUG toolbox::exportLog("INFOS: $0 : Running node for job $currentJID is $runningNode\n\n",1);
    
    #Provide the job ID
    return $currentJID;
    
}

sub mprunRun
{ #for SLURM  MPRUN cluster, running using ccc_msub
    
    my $msubOptionsHash=toolbox::extractHashSoft($configInfo,"mprun");
    my $msubOptions=toolbox::extractOptions($msubOptionsHash);
    
    #Picking up ENV variable
    my $envOptionsHash=toolbox::extractHashSoft($configInfo,"env");
    my $envOptions=toolbox::extractOptions($envOptionsHash,"=","\n");
    
    #Adding slurm options
    my $launcherCommand = "ccc_msub ".$msubOptions;
    #Creating the bash script for slurm to launch the command
    my $date =`date +%Y_%m_%d_%H_%M_%S`;
    chomp $date;
    my $scriptName="msubScript_".$date.".sh";
    my $bashScriptCreationCommand= "echo \"#!/bin/bash\n".$envOptions."\n".$commandLine."\nexit 0;\" | cat - > $scriptName && chmod 777 $scriptName";
    toolbox::run($bashScriptCreationCommand);
    ##DEBUG toolbox::exportLog("INFOS : $0 : Created the slurm bash file",1);
    $launcherCommand.=" $scriptName";
    $launcherCommand =~ s/ +/ /g; #Replace multiple spaces by a single one, to have a better view...
    ##DEBUG toolbox::exportLog($launcherCommand,2);
    my $currentJID = `$launcherCommand`;
    
    if ($!) #There are errors in the launching...
    {
        warn ("WARNING : $0 : Cannot launch the job for $sample: $!\n");
        $currentJID = "";
    }
        
    #Parsing infos and informing logs
    chomp $currentJID;
    
    unless ($currentJID) #The job has no output in STDOUT, ie there is a problem...
    {
        return 0; #Returning to launcher subprogram the error type
    }
    
    toolbox::exportLog("INFOS: $0 : Correctly launched for $sample in ccc_msub mode through the command:\n\t$launcherCommand\n",1);
    ## DEBUG toolbox::exportLog("INFOS: $0 : Output for the command is $currentJID\n\n",1);
    
    #Format: Submitted Batch Session 3223145
    my @infosList=split /\s/, $currentJID; #the format is such as "Submitted batch job 3"
    $currentJID = $infosList[3];
    
    ##THIS PART IS ONLY FOR DEBUGGING  -  MUST STAY COMMENTED OR MAY PROVOKE TROUBLES!!

    #my $runningNodeCommand="ccc_mstat -u \$USER| grep \"  $currentJID \"";
    #my $runningNode="x";
    #my $trying=0;
    #while ($runningNode ne "R01") #If the job is not yet launched or already finished
    #{
    #    sleep 3;#Waiting for the job to be launched
    #    $runningNode=`$runningNodeCommand`;
    #    chomp $runningNode;
    #    $runningNode = "x" unless $runningNode; #if empty variable, problem after...
    #    if ($runningNode !~ /\s+R\s+/)
    #    {# not running yet
    #        $trying++;
    #        if ($trying == 5)
    #        {
    #            #We already tryed to pick up the node infos 5 times, let's stop
    #            $runningNode = "still unknown (either not running, or already finished)";
    #            ## DEBUG toolbox::exportLog("WARNING : $0 : Cannot pickup the running node for the job $currentJID: $!\n",2);
    #            last;
    #        }
    #        next;
    #    }
    #    my @runningFields = split /\s+/,$runningNode; #To obtain the correct field
    #    $runningNode = $runningFields[8];
    #}
    ## DEBUG toolbox::exportLog("INFOS: $0 : Running node for job $currentJID is $runningNode\n\n",1);
    
    #Provide the job ID
    return $currentJID;
}

sub lsfRun{ #for LSF cluster, running using bsub
    
    my $bsubOptionsHash=toolbox::extractHashSoft($configInfo,"lsf");
    my $bsubOptions=toolbox::extractOptions($bsubOptionsHash);
    
    #Picking up ENV variable
    my $envOptionsHash=toolbox::extractHashSoft($configInfo,"env");
    my $envOptions=toolbox::extractOptions($envOptionsHash,"=","\n");
    
    #Adding slurm options
    my $launcherCommand = "bsub ".$bsubOptions;
    #Creating the bash script for slurm to launch the command
    my $date =`date +%Y_%m_%d_%H_%M_%S`;
    chomp $date;
    my $scriptName="bsubScript_".$date.".sh";
    my $bashScriptCreationCommand= "echo \"#!/bin/bash\n".$envOptions."\n".$commandLine."\nexit 0;\" | cat - > $scriptName && chmod 777 $scriptName";
    toolbox::run($bashScriptCreationCommand);
    ##DEBUG toolbox::exportLog("INFOS : $0 : Created the slurm bash file",1);
    $launcherCommand.=" $scriptName";
    $launcherCommand =~ s/ +/ /g; #Replace multiple spaces by a single one, to have a better view...
    ##DEBUG toolbox::exportLog($launcherCommand,2);
    my $currentJID = `$launcherCommand`;
    
    if ($!) #There are errors in the launching...
    {
        warn ("WARNING : $0 : Cannot launch the job for $sample: $!\n");
        $currentJID = "";
    }
        
    #Parsing infos and informing logs
    chomp $currentJID;
    
    unless ($currentJID) #The job has no output in STDOUT, ie there is a problem...
    {
        return 0; #Returning to launcher subprogram the error type
    }
    
    toolbox::exportLog("INFOS: $0 : Correctly launched for $sample in bsub mode through the command:\n\t$launcherCommand\n",1);
    ## DEBUG toolbox::exportLog("INFOS: $0 : Output for the command is $currentJID\n\n",1);
    
    #Format: Submitted Batch Session 3223145
    #Changing queue from '' to normal
    #Changing application from '' to user
    #Changing memory request to 4G	
    #Job <21773> is submitted to queue <normal>.
    my @infosList=split /<|>/, $currentJID; 
    $currentJID = $infosList[1];
    
    ##THIS PART IS ONLY FOR DEBUGGING  -  MUST STAY COMMENTED OR MAY PROVOKE TROUBLES!!

    #my $runningNodeCommand="bjobs $currentJID";
    #my $runningNode="x";
    #my $trying=0;
    #while ($runningNode ne "RUN") #If the job is not yet launched or already finished
    #{
    #    sleep 3;#Waiting for the job to be launched
    #    $runningNode=`$runningNodeCommand`;
    #    chomp $runningNode;
    #    $runningNode = "x" unless $runningNode; #if empty variable, problem after...
    #    if ($runningNode !~ /\s+R\s+/)
    #    {# not running yet
    #        $trying++;
    #        if ($trying == 5)
    #        {
    #            #We already tryed to pick up the node infos 5 times, let's stop
    #            $runningNode = "still unknown (either not running, or already finished)";
    #            ## DEBUG toolbox::exportLog("WARNING : $0 : Cannot pickup the running node for the job $currentJID: $!\n",2);
    #            last;
    #        }
    #        next;
    #    }
    #    my @runningFields = split /\s+/,$runningNode; #To obtain the correct field
    #    $runningNode = $runningFields[8];
    #}
    ## DEBUG toolbox::exportLog("INFOS: $0 : Running node for job $currentJID is $runningNode\n\n",1);
    
    #Provide the job ID
    return $currentJID;
}

##################################
#
#WAITING for schedulers ONLY!
#
##################################

sub waiter
{ #Global function for waiting, will recover the jobID to survey and will re-dispatch to scheduler
    
    ($jobList,my $jobsInfos) = @_;
    
    %jobHash = %{$jobsInfos};
        
    my $hashCapability = &checkingCapability;
    my $stopWaiting;
    switch (1)
    {
        case ($hashCapability->{"sge"} && defined $$configInfo{"sge"}){$stopWaiting = &sgeWait} #For SGE running
        case ($hashCapability->{"slurm"} && defined $$configInfo{"slurm"}){$stopWaiting = &slurmWait} #For SLURM running
	case ($hashCapability->{"mprun"} && defined $$configInfo{"mprun"}){$stopWaiting = &mprunWait} #For mprun running
	case ($hashCapability->{"lsf"} && defined $$configInfo{"lsf"}){$stopWaiting = &lsfWait} #For lsf running

    }
    return $stopWaiting;
}


sub sgeWait
{
    
    my $nbRunningJobs = 1;
    my @jobsInError=();
    
    ##Waiting for jobs to finish
    while ($nbRunningJobs)
    {  
      #Picking up the number of currently running jobs
      my $qstatCommand = "qstat | egrep -c \"$jobList\"";
      $nbRunningJobs = `$qstatCommand`;
      chomp $nbRunningJobs;
      sleep 50;
    }
    
    #Compiling infos about sge jobs: jobID, node number, exit status
    sleep 25;#Needed for qacct to register infos...
    toolbox::exportLog("\n#########################################\nJOBS SUMMARY\n#########################################
\n---------------------------------------------------------
Individual\tJobID\tNode\tExitStatus
---------------------------------------------------------",1);
    
    foreach my $individual (sort {$a cmp $b} keys %jobHash)
    {
      my $qacctCommand = "qacct -j ".$jobHash{$individual}." 2>&1";
      my $qacctOutput = `$qacctCommand`;
      my $outputLine;
      
      chomp $qacctOutput;
      if ($qacctOutput =~ "-bash: qacct" or $qacctOutput =~ "installed")
      {
        #IF qacct cannot be run on the node
        $outputLine = "$individual\t$jobHash{$individual}\tNA\tNA\n";
        toolbox::exportLog($outputLine,1);
        next;
      }
      
      my @linesQacct = split /\n/, $qacctOutput;
      $outputLine = $individual."\t".$jobHash{$individual}."\t";
      while (@linesQacct) #Parsing the qacct output
      {
        my $currentLine = shift @linesQacct;
        if ($currentLine =~ m/^hostname/) #Picking up hostname
        {
          $currentLine =~ s/hostname     //;
          $outputLine .= $currentLine."\t";
        }
        elsif ($currentLine =~ m/^exit_status/) #Picking up exit status
        {
          $currentLine =~ s/exit_status  //;
          if ($currentLine == 0) #No errors
	  {
	    $currentLine = "Normal";
	  }
	  else
	  {
	    push @jobsInError, $individual;
	    $currentLine = "Error";
	  }
          $outputLine .= $currentLine;
        }
        else
        {
          next;
        }
        
      }
      toolbox::exportLog($outputLine,1);
      
    }
    toolbox::exportLog("---------------------------------------------------------\n",1);#To have a better table presentation
  
    if (scalar @jobsInError)
    {
	#at least one job has failed
	return \@jobsInError;
    }
    return 1;
}

sub slurmWait
{
    
    my $nbRunningJobs = 1;
    my @jobsInError=();
    
    ##Waiting for jobs to finish
    while ($nbRunningJobs)
    {  
      #Picking up the number of currently running jobs
      my $squeueCommand = "squeue | egrep -c \"$jobList\"";
      $nbRunningJobs = `$squeueCommand`;
      chomp $nbRunningJobs;
      sleep 50;
    }
    
    #Compiling infos about sge jobs: jobID, node number, exit status
    sleep 25;#Needed for qacct to register infos...
    toolbox::exportLog("\n#########################################\nJOBS SUMMARY\n#########################################
\n---------------------------------------------------------
Individual\tJobID\tNode\tExitStatus
---------------------------------------------------------",1);
    
    foreach my $individual (sort {$a cmp $b} keys %jobHash)
    {
      my $sacctCommand = "sacct -j ".$jobHash{$individual};
      my $sacctOutput = `$sacctCommand`;
      my $outputLine;
      chomp $sacctOutput;
      if ($sacctOutput =~ "-bash: qacct" or $sacctOutput =~ "installed")
      {
        #IF sacct cannot be run on the node
        $outputLine = "$individual\t$jobHash{$individual}\tNA\tNA\n";
        toolbox::exportLog($outputLine,1);
        next;
      }
      my @lineSacct=split/\n/, $sacctOutput;
      $outputLine=$individual."\t".$jobHash{$individual}."\t";
      while (@lineSacct)
      {
	my $line = shift @lineSacct;
	
	#Passing the header lines
	next if $line =~ m/JobID/;
	next if $line !~ m/^\d/;
	
	my @fields = split /\s+/, $line;

	if ($fields[5] eq "COMPLETED") #No errors
	{
	    $outputLine .= "Normal";
	}
	else
	{
	    push @jobsInError, $individual;
	    $outputLine .= "Error";
	}
      }
      
      toolbox::exportLog($outputLine,1);
      
    }
    toolbox::exportLog("---------------------------------------------------------\n",1);#To have a better table presentation
  
    if (scalar @jobsInError)
    {
	#at least one job has failed
	return \@jobsInError;
    }
    return 1;
}


sub mprunWait
{
    
    my $nbRunningJobs = 1;
    my @jobsInError=();
    
    ##Waiting for jobs to finish
    while ($nbRunningJobs)
    {  
      #Picking up the number of currently running jobs
      my $squeueCommand = "ccc_mstat -u \$USER | egrep -c \"$jobList\"";
      $nbRunningJobs = `$squeueCommand`;
      chomp $nbRunningJobs;
      sleep 50;
    }
    
    #Compiling infos about sge jobs: jobID, node number, exit status
    sleep 25;#Needed for ccc_macct to register infos...
    toolbox::exportLog("\n#########################################\nJOBS SUMMARY\n#########################################
\n---------------------------------------------------------
Individual\tJobID\tExitStatus
---------------------------------------------------------",1);
    
    foreach my $individual (sort {$a cmp $b} keys %jobHash)
    {
      my $macctCommand = "ccc_macct ".$jobHash{$individual}."|tail -n 1";
      my $macctOutput = `$macctCommand`;
      my $outputLine;
      chomp $macctOutput;
      if ($macctOutput =~ "-bash: ccc_macct" or $macctOutput =~ "installed")
      {
        #IF sacct cannot be run on the node
        $outputLine = "$individual\t$jobHash{$individual}\tNA\tNA\n";
        toolbox::exportLog($outputLine,1);
        next;
      }
      my @fields = split /\s+/, $macctOutput;

      if ($macctOutput=~ m/COMPLETED/) #No errors
      {
	$outputLine = "$individual\t$jobHash{$individual}\tNormal";
      }
      else
      {
	push @jobsInError, $individual;
	$outputLine = "$individual\t$jobHash{$individual}\tError";
      }

      
      toolbox::exportLog($outputLine,1);
      
    }
    toolbox::exportLog("---------------------------------------------------------\n",1);#To have a better table presentation
  
    if (scalar @jobsInError)
    {
	#at least one job has failed
	return \@jobsInError;
    }
    return 1;
}

sub lsfWait
{
    
    my $nbRunningJobs = 1;
    my @jobsInError=();
    
    ##Waiting for jobs to finish
    while ($nbRunningJobs)
    {  
      #Picking up the number of currently running jobs
      my $bjobsQueueCommand = "bjobs -u \$USER | egrep -c \"$jobList\"";
      #Looks like
	#~$ bjobs
	#JOBID USER STAT QUEUE FROM_HOST EXEC_HOST JOB_NAME SUBMIT_TIME
	#21773 rdkassa RUN normal rddor-rdkas rddoridt51 *e_test.sh Apr 7 14:19
      $nbRunningJobs = `$bjobsQueueCommand`;
      chomp $nbRunningJobs;
      sleep 50;
    }
    
    #Compiling infos about sge jobs: jobID, node number, exit status
    sleep 25;#Needed for bacct to register infos...
    toolbox::exportLog("\n#########################################\nJOBS SUMMARY\n#########################################
\n---------------------------------------------------------
Individual\tJobID\tExitStatus
---------------------------------------------------------",1);
    
    foreach my $individual (sort {$a cmp $b} keys %jobHash)
    {
      my $bacctCommand = "bacct ".$jobHash{$individual}."|grep \"Total number of done jobs:\"";
      my $bacctOutput = `$bacctCommand`;
      my $outputLine;
      chomp $bacctOutput;
      if ($bacctOutput =~ "-bash: bacct" or $bacctOutput =~ "installed")
      {
        #IF bacct cannot be run on the node
        $outputLine = "$individual\t$jobHash{$individual}\tNA\tNA\n";
        toolbox::exportLog($outputLine,1);
        next;
      }
      my @fields = split /\s+/, $bacctOutput;

      if ($bacctOutput=~ m/Total number of done jobs:\s*1\s*Total number of exited jobs:\s*0\s*/) #No errors
      {
	$outputLine = "$individual\t$jobHash{$individual}\tNormal";
      }
      else
      {
	push @jobsInError, $individual;
	$outputLine = "$individual\t$jobHash{$individual}\tError";
      }

      
      toolbox::exportLog($outputLine,1);
      
    }
    toolbox::exportLog("---------------------------------------------------------\n",1);#To have a better table presentation
  
    if (scalar @jobsInError)
    {
	#at least one job has failed
	return \@jobsInError;
    }
    return 1;
}
1;

=head1 NAME

    Package I<scheduler> 

=head1 SYNOPSIS

	use scheduler;
    
	scheduler::launcher($launcherCommand, $requirement, $sample, $configInfo);
    
	scheduler::waiter($jobList,$jobInfos);

=head1 DESCRIPTION

    Package scheduler will prepare and launch the different command in a scheduler if possible and requested. Either, it will launch it through the normal way, using a quite identical way to toolbox::run
    !!DIFFERENCE with toolbox::run is that an abnormal job will be a blocking one ONLY if requirement is of a zero value (0). Else, the job error will be warned but it will not block the whole system.

=head2 FUNCTIONS

=head3 scheduler::launcher

This module will prepare and decide to which scheduler sending the command. It takes as arguments the command $launcherCommand, the state of the job (0 mandatory finishing, non-0 not mandatory) $requirement, the sample name $sample and the $configInfo hash for configuration

=head3 scheduler::waiter

This module will allow the job launched through a scheduler to be wait to finish.
It takes as arguments the list of jobs $jobList and the reference of the hash for the informations about the jobs (sample name) $jobInfos



=head1 AUTHORS

Intellectual property belongs to IRD, CIRAD and South Green developpement plateform for all versions also for ADNid for v2 and v3 and INRA for v3
Written by Cecile Monat, Christine Tranchant, Laura Helou, Abdoulaye Diallo, Julie Orjuela-Bouniol, Sebastien Ravel, Gautier Sarah, and Francois Sabot

=head1 SEE ALSO

L<http://www.southgreen.fr/>

=cut
