There are 3 ways to obtain TOGGLE. We recommend the first one for an install in a cluster.

# INSTALL FROM GitHub REPOSITORY

* Create the directory TOGGLE where you want to install TOGGLE and go into this directory

* Get the TOGGLE code Clone the git

````
$git clone https://github.com/SouthGreenPlatform/TOGGLE.git
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
* Modify the shebang of perl in the begining of the script /pathToToggle/globalAnalysis.pl
* Modify the shebang of perl in the begining of the script /pathToToggle/singleAnalysis.pl
* Modify the shebang of perl in the begining of the script /pathToToggle/pairAnalysis.pl
* Modify the shebang of perl in the begining of the script /pathToToggle/mergeAnalysis.pl

* Run the test script

````
$cd TEST
$sh all_tests.sh
````

* Test the pipiline with the test data

````
$globalAnalysis.pl $PATH_INSTALL/DATA/arcardTest/ $PATH_INSTALL/software.config.txt $PATH_INSTALL/DATA/referenceArcad.fasta
````

* Check the good running
> > * No error message
> > * BamDirectory has been well created into $PATH_INSTALL/DATA
> > * the data generated are good


````
$globalAnalysis.pl $PATH_INSTALL/DATA/iriginTest/ $PATH_INSTALL/software.config.txt $PATH_INSTALL/DATA/referenceIrigin.fasta
````

* Check the good running
> > * No error message
> > * BamDirectory has been well created into $PATH_INSTALL/DATA
> > * the data generated are good

````
tail $PATH_INSTALL/DATA/BamDirectory/GATKVARIANTFILTRATION.vcf

##contig=<ID=2299889,length=39753>
##contig=<ID=2299897,length=19555>
##reference=file:///home/tranchant/TOGGLE-19-02/DATA/referenceIrigin.fasta
#CHROM    POS    ID    REF    ALT    QUAL    FILTER    INFO    FORMAT    irigin1    irigin1Single    irigin2_1    irigin3    irigin3Single
2224477    996    .    TA    T    35.61    PASS    AC=2;AF=1.00;AN=2;DP=2;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;MQ0=0;QD=17.81;SOR=0.693    GT:AD:DP:GQ:PL    ./.    ./.    ./.    1/1:0,2:2:6:69,6,0    ./.
2248321    377    .    C    G    65.65    PASS    AC=2;AF=1.00;AN=2;DP=2;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=29.00;MQ0=0;QD=32.83;SOR=0.693    GT:AD:DP:GQ:PL    ./.    ./.    ./.    1/1:0,2:2:6:90,6,0    ./.
2248321    379    .    C    T    65.65    PASS    AC=2;AF=1.00;AN=2;DP=2;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=29.00;MQ0=0;QD=32.83;SOR=0.693    GT:AD:DP:GQ:PL    ./.    ./.    ./.    1/1:0,2:2:6:90,6,0    ./.
2281178    4213    .    G    A    65.65    PASS    AC=2;AF=1.00;AN=2;DP=2;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;MQ0=0;QD=32.83;SOR=0.693    GT:AD:DP:GQ:PL    ./.    ./.    ./.    1/1:0,2:2:6:90,6,0    ./.
2281178    4214    .    A    G    65.65    PASS    AC=2;AF=1.00;AN=2;DP=2;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;MQ0=0;QD=32.83;SOR=0.693    GT:AD:DP:GQ:PL    ./.    ./.    ./.    1/1:0,2:2:6:90,6,0    ./.
2290182    1013    .    A    G    45.65    PASS    AC=2;AF=1.00;AN=2;DP=2;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=29.00;MQ0=0;QD=22.83;SOR=0.693    GT:AD:DP:GQ:PL    ./.    ./.    ./.    1/1:0,2:2:6:70,6,0    ./.
````
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
$wget http://bioinfo-web.mpl.ird.fr/toggle/installationTOGGLE.sh
$sh installationTOGGLE.sh
````

Then you can launch the different scripts in the TOGGLE folder
