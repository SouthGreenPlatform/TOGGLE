There are 3 ways to obtain TOGGLE. We recommend the first one for an install in a cluster.

# INSTALL FROM GitHub REPOSITORY

* Create the directory TOGGLE where you want to install TOGGLE and go into this directory

* Get the TOGGLE code Clone the git

````
$git clone https://github.com/SouthGreenPlatform/TOGGLE.git .
````
* Add the Module path to the PERL5LIB environment variable

````
export PERL5LIB=$PERL5LIB:/pathToToggle/Modules
````
* In the same way add the TOGGLE directory to the PATH environment variable

````
export PATH=$PATH:/pathToToggle
````

Rq : you can add this to the ~/.bashrc to make it always available when you log-in.


* Change the permission of the perl files in the directory TOGGLE

```
$cd /pathToToggle
$chmod 755 *pl
```

* Modify the file /pathToToggle/Modules/localConfig.pm
* Modify the shebang of perl in the begining of the script /pathToToggle/toggleGenerator.pl and onTheFly/startBlock.txt 



* Test the pipeline with the test data

Modify the path of adaptator file in the file SNPdiscoveryPaired.config.txt

````
$cutadapt
-O=10
-m=35
-q=20
--overlap=7
-u=8
#If you have a specific adaptator file, please indicate here. 
-adaptatorFile=/path/to/adaptator.txt
````

````
$cd /pathToToggle
$toggleGenerator.pl -c SNPdiscoveryPaired.config.txt -d DATA/testData/fastq/pairedTwoIndividusIrigin/ -r DATA/Bank/referenceIrigin.fasta -o toggleout2/
````

* Check the good running
> > * No error message
> > * toggleOUTPUT has been created in the /tmp folder
> > * the data generated are good

````
tail /tmp/toggleOUTPUT/finalResults/GATKVARIANTFILTRATION.vcf

TO BE CHANGED FOR  THE FINAL VCF USING IRIGIN
````


* Test the RNASeq pipeline with the test data

````
$cd /pathToToggle
toggleGenerator.pl -c RNASeq.config.txt -d DATA/testData/rnaseq/pairedOneIndividu/ -r DATA/Bank/referenceRnaseq.fa -g DATA/Bank/referenceRnaseq.gff3 -o toggleout/
````

* Check the good running
> > * No error message
> > * the directory toggleout has been created
> > * the directory toggleout/output contain all the results
> > * the file toggleout/finalResults/RNASeq.accepted_hits.HTSEQCOUNT.txt has been created
> > * the data generated are good


# INSTALL FROM THE DOCKER IMAGE

A Docker image based on Ubuntu 14.04 is available at http://bioinfo-web.mpl.ird.fr/toggle/toggle.tgz

If Docker is installed on your system, you can use this image:

````
$wget http://bioinfo-web.mpl.ird.fr/toggle/toggle.tgz
$cat toggle.tgz | docker import - toggle
$docker run -i -t toggle /bin/bash
````

You can then download your data in the container (or attach a network/local drive within). TOGGLE can be launch as for a GitHub installation.

# INSTALL FROM THE AUTOMATIC SCRIPT

This script is dedicated for an installation in the user space (no admin rights required).

Just download it, launch it and follow the instructions.

````
$wget http://bioinfo-web.mpl.ird.fr/toggle/TOGGLEinstall.sh
$bash TOGGLEinstall.sh
````

Then you can launch the different scripts in the TOGGLE folder
