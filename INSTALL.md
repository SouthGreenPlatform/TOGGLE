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
* Modify the shebang of perl in the begining of the script for RNASeq analysis

* Run the test script (BE CAREFUL: will succeed only with specifc versions of the different softwares - see ChangeLog.md)

````
$cd TEST
$sh all_tests.sh
````

* Test the pipeline with the test data

````
$globalAnalysis.pl -d $PATH_INSTALL/DATA/arcardTest/ -c $PATH_INSTALL/software.config.txt -r $PATH_INSTALL/DATA/referenceArcad.fasta
````

* Check the good running
> > * No error message
> > * BamDirectory has been well created into $PATH_INSTALL/DATA
> > * the data generated are good

````
tail $PATH_INSTALL/DATA/BamDirectory/GATKVARIANTFILTRATION.vcf

##contig=<ID=LOC_Os01g62920.1,length=2879>
##contig=<ID=LOC_Os12g32240.1,length=1088>
##reference=file:///home/ravel/TOGGLE/DATA/referenceArcad.fasta
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  arcad1Single    arcad1_1.REPAIRING.BWAALN       arcad2_1        arcad3Single    arcad3_1.REPAIRING.BWAALN
LOC_Os01g44110.1        2316    .       G       T       38.78   PASS    AC=3;AF=0.500;AN=6;BaseQRankSum=-1.422;ClippingRankSum=-1.279;DP=24;FS=4.873;MLEAC=3;MLEAF=0.500;MQ=56.82;MQ0=0;MQRankSum=-0.142;QD=1.62;ReadPosRankSum=0.213;SOR=0.707 GT:AD:DP:GQ:PL  ./.     0/1:8,2:10:19:19,0,200  0/1:3,1:4:30:30,0,113   ./.     0/1:8,2:10:19:19,0,200
LOC_Os01g62920.1        293     .       C       T       74.70   PASS    AC=3;AF=0.500;AN=6;BaseQRankSum=-1.000;ClippingRankSum=0.380;DP=12;FS=0.000;MLEAC=3;MLEAF=0.500;MQ=56.82;MQ0=0;MQRankSum=-0.529;QD=6.23;ReadPosRankSum=-0.689;SOR=0.760 GT:AD:DP:GQ:PL  ./.     0/1:3,2:5:34:34,0,64    0/1:1,1:2:36:36,0,36    ./.     0/1:3,2:5:34:34,0,64
LOC_Os01g62920.1        1602    .       T       C       1050.69 PASS    AC=6;AF=1.00;AN=6;DP=34;FS=0.000;MLEAC=6;MLEAF=1.00;MQ=55.45;MQ0=0;QD=30.90;SOR=1.244   GT:AD:DP:GQ:PL  ./.     1/1:0,13:13:39:371,39,0 1/1:0,8:8:24:334,24,0   ./.     1/1:0,13:13:39:371,39,0
LOC_Os01g62920.1        1721    .       A       G       2197.56 PASS    AC=10;AF=1.00;AN=10;DP=65;FS=0.000;MLEAC=10;MLEAF=1.00;MQ=57.37;MQ0=0;QD=33.81;SOR=0.854        GT:AD:DP:GQ:PL  1/1:0,1:1:3:37,3,0
      1/1:0,28:28:84:928,84,0 1/1:0,7:7:21:292,21,0   1/1:0,1:1:3:37,3,0      1/1:0,28:28:84:928,84,0
LOC_Os01g62920.1        1907    .       T       C       2228.68 PASS    AC=6;AF=1.00;AN=6;DP=60;FS=0.000;MLEAC=6;MLEAF=1.00;MQ=56.82;MQ0=0;QD=33.32;SOR=3.113   GT:AD:DP:GQ:PL  ./.     1/1:0,25:25:75:918,75,0 1/1:0,10:10:30:418,30,0 ./.     1/1:0,25:25:75:918,75,0
LOC_Os01g62920.1        2305    .       A       C       35.92   PASS    AC=1;AF=0.167;AN=6;BaseQRankSum=-1.472;ClippingRankSum=-1.383;DP=28;FS=12.553;MLEAC=1;MLEAF=0.167;MQ=54.43;MQ0=0;MQRankSum=-1.561;QD=4.49;ReadPosRankSum=-1.829;SOR=2.205       GT:AD:DP:GQ:PL  ./.     0/0:10,0:10:30:0,30,735 0/1:6,2:8:66:66,0,392   ./.     0/0:10,0:10:30:0,30,735
LOC_Os01g62920.1        2308    .       T       G       35.92   PASS    AC=1;AF=0.167;AN=6;BaseQRankSum=-0.937;ClippingRankSum=-0.937;DP=28;FS=12.553;MLEAC=1;MLEAF=0.167;MQ=54.43;MQ0=0;MQRankSum=-2.007;QD=4.49;ReadPosRankSum=-2.364;SOR=2.205       GT:AD:DP:GQ:PL  ./.     0/0:10,0:10:30:0,30,735 0/1:6,2:8:66:66,0,392   ./.     0/0:10,0:10:30:0,30,735
LOC_Os12g32240.1        864     .       C       T       974.62  PASS    AC=10;AF=1.00;AN=10;DP=26;FS=0.000;MLEAC=10;MLEAF=1.00;MQ=55.54;MQ0=0;QD=26.53;SOR=4.255        GT:AD:DP:GQ:PL  1/1:0,1:1:3:37,3,0
      1/1:0,10:10:30:379,30,0 1/1:0,4:4:12:167,12,0   1/1:0,1:1:3:37,3,0      1/1:0,10:10:30:379,30,0
````

````
$globalAnalysis.pl -d $PATH_INSTALL/DATA/iriginTest/ -c $PATH_INSTALL/software.config.txt -r $PATH_INSTALL/DATA/referenceIrigin.fasta
````

* Check the good running
> > * No error message
> > * BamDirectory has been well created into $PATH_INSTALL/DATA
> > * the data generated are good

````
tail $PATH_INSTALL/DATA/BamDirectory/GATKVARIANTFILTRATION.vcf

##reference=file:///home/ravel/TOGGLE/DATA/referenceArcad.fasta
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	arcad1Single	arcad1_1.REPAIRING.BWAALN	arcad2_1	arcad3Single	arcad3_1.REPAIRING.BWAALN
LOC_Os01g44110.1	2316	.	G	T	38.78	PASS	AC=3;AF=0.500;AN=6;BaseQRankSum=-1.422;ClippingRankSum=-1.279;DP=24;FS=4.873;MLEAC=3;MLEAF=0.500;MQ=56.82;MQ0=0;MQRankSum=-0.142;QD=1.62;ReadPosRankSum=0.213;SOR=0.707	GT:AD:DP:GQ:PL	./.	0/1:8,2:10:19:19,0,200	0/1:3,1:4:30:30,0,113	./.	0/1:8,2:10:19:19,0,200
LOC_Os01g62920.1	293	.	C	T	74.70	PASS	AC=3;AF=0.500;AN=6;BaseQRankSum=-1.000;ClippingRankSum=0.380;DP=12;FS=0.000;MLEAC=3;MLEAF=0.500;MQ=56.82;MQ0=0;MQRankSum=-0.529;QD=6.23;ReadPosRankSum=-0.689;SOR=0.760	GT:AD:DP:GQ:PL	./.	0/1:3,2:5:34:34,0,64	0/1:1,1:2:36:36,0,36	./.	0/1:3,2:5:34:34,0,64
LOC_Os01g62920.1	1602	.	T	C	1050.69	PASS	AC=6;AF=1.00;AN=6;DP=34;FS=0.000;MLEAC=6;MLEAF=1.00;MQ=55.45;MQ0=0;QD=30.90;SOR=1.244	GT:AD:DP:GQ:PL	./.	1/1:0,13:13:39:371,39,0	1/1:0,8:8:24:334,24,0	./.	1/1:0,13:13:39:371,39,0
LOC_Os01g62920.1	1721	.	A	G	2197.56	PASS	AC=10;AF=1.00;AN=10;DP=65;FS=0.000;MLEAC=10;MLEAF=1.00;MQ=57.37;MQ0=0;QD=33.81;SOR=0.854	GT:AD:DP:GQ:PL	1/1:0,1:1:3:37,3,0	1/1:0,28:28:84:928,84,0	1/1:0,7:7:21:292,21,0	1/1:0,1:1:3:37,3,0	1/1:0,28:28:84:928,84,0
LOC_Os01g62920.1	1907	.	T	C	2228.68	PASS	AC=6;AF=1.00;AN=6;DP=60;FS=0.000;MLEAC=6;MLEAF=1.00;MQ=56.82;MQ0=0;QD=33.32;SOR=3.113	GT:AD:DP:GQ:PL	./.	1/1:0,25:25:75:918,75,0	1/1:0,10:10:30:418,30,0	./.	1/1:0,25:25:75:918,75,0
LOC_Os01g62920.1	2305	.	A	C	35.92	PASS	AC=1;AF=0.167;AN=6;BaseQRankSum=-1.472;ClippingRankSum=-1.383;DP=28;FS=12.553;MLEAC=1;MLEAF=0.167;MQ=54.43;MQ0=0;MQRankSum=-1.561;QD=4.49;ReadPosRankSum=-1.829;SOR=2.205	GT:AD:DP:GQ:PL	./.	0/0:10,0:10:30:0,30,735	0/1:6,2:8:66:66,0,392	./.	0/0:10,0:10:30:0,30,735
LOC_Os01g62920.1	2308	.	T	G	35.92	PASS	AC=1;AF=0.167;AN=6;BaseQRankSum=-0.937;ClippingRankSum=-0.937;DP=28;FS=12.553;MLEAC=1;MLEAF=0.167;MQ=54.43;MQ0=0;MQRankSum=-2.007;QD=4.49;ReadPosRankSum=-2.364;SOR=2.205	GT:AD:DP:GQ:PL	./.	0/0:10,0:10:30:0,30,735	0/1:6,2:8:66:66,0,392	./.	0/0:10,0:10:30:0,30,735
LOC_Os12g32240.1	864	.	C	T	974.62	PASS	AC=10;AF=1.00;AN=10;DP=26;FS=0.000;MLEAC=10;MLEAF=1.00;MQ=55.54;MQ0=0;QD=26.53;SOR=4.255	GT:AD:DP:GQ:PL	1/1:0,1:1:3:37,3,0	1/1:0,10:10:30:379,30,0	1/1:0,4:4:12:167,12,0	1/1:0,1:1:3:37,3,0	1/1:0,10:10:30:379,30,0
````

* Test the RNASeq pipeline with the test data

````
$globalAnalysisRnaSeq.pl $PATH_INSTALL/DATA/rnaseqTest/ $PATH_INSTALL/software.config.txt $PATH_INSTALL/DATA/referenceRnaseq.fa $PATH_INSTALL/DATA/referenceRnaseq.gff3
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
$wget http://bioinfo-web.mpl.ird.fr/toggle/TOGGLEinstall.sh
$bash TOGGLEinstall.sh
````

Then you can launch the different scripts in the TOGGLE folder
