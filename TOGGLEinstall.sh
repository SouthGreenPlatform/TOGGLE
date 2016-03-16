#!/bin/bash

###################################################################################################################################
#
# Copyright 2014-2016 IRD-CIRAD-INRA-ADNid
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

###################################################################################################################################
#
# This script will allow a user-space automatic installation of TOGGLE from the current Master version https://github.com/SouthGreenPlatform/TOGGLE
#
###################################################################################################################################


# USE AN INTEGER TO SIGNAL THE NUMBER OF PROGRAMS YOU NEED TO CHECK
# DEBUG -- FOR DEV: total=5
# DEBUG -- FOR PROD:
total=4

# TEST IF MACHINE IS 64BITS

echo "\n\n############################################"
echo "##\t Checking if your system is a 64 bits"
echo "############################################\n"


MACHINE_TYPE=`uname -m`;
if [ ${MACHINE_TYPE} = "x86_64" ];
then
#	echo "$MACHINE_TYPE";
	echo "Your system is a 64 bits:\t$MACHINE_TYPE";
else
#	echo $MACHINE_TYPE;
	echo "** It seems that your system is not 64bits.**";
	echo "**Please verify that your system is 64 bits**";
	echo "** $MACHINE_TYPE **";
        exit 1;
fi




# TEST IF SOFTWARES ARE INSTALLED, AND PRINT THEIR PATH

echo "\n\n############################################"
echo "##\t Searching for some required software being installed"
echo "############################################\n"


x=0
for i in perl wget git tar

do
#echo "Searching for perl";
	command -v $i > /dev/null 2>&1 && echo "$i is installed in :\t"`command -v $i` && x=$(( $x + 1 )) || { echo >&2 "** TOGGLE requires $i but it is not installed.  Please install it or contact your administrator for installation..."; };
#	x=$(( $x + 1 ));
done



echo "\n\tYou have all required sotfwares already installed"



if [ $x -lt $total ] ;
then
	echo "\nYour system is missing some programs, please install them and try again.";
	exit 1;
elif [ $x -eq $total ]
then
	echo "\nRequired UNIX programs have been found  !!\nYou are one step closer to installing TOGGLE"
else
	echo "\nThis is a bug message from TOGGLE devs.\nThis was not supposed to happen, you might want to check the raw code or contact us through https://github.com/SouthGreenPlatform/TOGGLE/issues"
	exit 1;
fi


# INPUT VARIABLES FOR JAVA

echo "\n\n############################################"
echo "##\t Path to Java "
echo "############################################\n"

echo "\nTOGGLE requires a functional version of Java "
echo "Please input the adress to the executable file\n"

# INPUT JAVA 7 PATH

echo "\nType the absolute path for Java :"
read JAVASEVEN

while true; do

    echo "\n";
    $JAVASEVEN -version;
    echo "\n";

    read -p "IS THIS THE RIGHT JAVA VERSION ?? [Y|N]: " yn
    case $yn in
        [Yy]* ) echo "Java path has been established to: '$JAVASEVEN' "; break;;
        [Nn]* ) echo "Please, Try again:"; read JAVASEVEN ;;
        * ) echo "Please answer 'Y' or 'N':";;
    esac
done

#######################################################
## License infos for all softwares
#######################################################
sleep 2

echo "By installing this software you agree to the Licenses, Copyrights or Copylefts of all the softwares used by TOGGLE, as well as the License of TOGGLE.\n"

sleep 2

echo "Please Visit the individual sites for the details of the corresponding Licenses and the citations details:\n"

sleep 2

echo "
CutAdapt
\tLicense: https://github.com/marcelm/cutadapt/blob/master/LICENSE
\tTo cite: http://journal.embnet.org/index.php/embnetjournal/article/view/200

bwa
\tLicense: http://sourceforge.net/projects/bio-bwa/
\tTo cite: http://www.ncbi.nlm.nih.gov/pubmed/19451168

SAMtools
\tLicense: http://sourceforge.net/projects/samtools/
\tTo cite: http://www.ncbi.nlm.nih.gov/pubmed/19505943

Picard-Tools
\tLicense: No explicit License
\tTo cite: Cite SAMtools paper and their site http://broadinstitute.github.io/picard/

FastQC
\tLicense: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
\tTo cite: http://www.bioinformatics.bbsrc.ac.uk/projects/fastqc/

GATK
\tLicense: https://github.com/broadgsa/gatk-protected/
\tTo cite: http://www.ncbi.nlm.nih.gov/pubmed?term=20644199
\tTo cite: http://www.ncbi.nlm.nih.gov/pubmed?term=21478889
\tTo cite: http://onlinelibrary.wiley.com/doi/10.1002/0471250953.bi1110s43/abstract;jsessionid=D95C25686A6F9F397B710DE983CE10D8.f03t02

TopHat
\tLicense: https://github.com/infphilo/tophat/blob/master/LICENSE
\tTo cite: http://bioinformatics.oxfordjournals.org/content/25/9/1105.abstract

Bowtie2
\tLicense: http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.5/
\tTo cite: http://www.nature.com/nmeth/journal/v9/n4/full/nmeth.1923.html

FASTX-Trimmer
\tLicense: http://hannonlab.cshl.edu/fastx_toolkit/license.html
\tTo cite: http://hannonlab.cshl.edu/fastx_toolkit/

TOGGLE
\tLicense: https://github.com/SouthGreenPlatform/TOGGLE/blob/master/LICENSE
\tTo cite: Monat et al, TOGGLE: toolbox for generic NGS analyses, BMC Bioinformatics, 2015, 16:374
"

sleep 5

echo "!! BE CAREFUL !!

We cannot automatically install HTSeqCount, as it requires PySam to be installed as an administrator..."

sleep 2

#######################################################
## TOGGLE installation per se
#######################################################

echo "\nINSTALLING TOGGLE\n";

echo "\nPlease provide the installation path:"
read INSTALLPATH

#Transforming in an absolute PATH. Using -f option, all composants must exist but the last

#$INSTALLPATH = readlink -f $INSTALLPATH

while true; do
    
    echo "\n\t$INSTALLPATH\n";
    
    read -p "IS THIS THE RIGHT INSTALL PATH ?? [Y|N]: " yn
    case $yn in
        [Yy]* ) echo "Install path is: '$INSTALLPATH' "; break;;
        [Nn]* ) echo "Please, Try again:"; read INSTALLPATH ;;
        * ) echo "Please answer 'Y' or 'N':";;
    esac
done

mkdir $INSTALLPATH

#Cloning current version of TOGGLE

echo "\nCloning the current Git Master Version";
sleep 1

git clone https://github.com/SouthGreenPlatform/TOGGLE.git $INSTALLPATH

cd $INSTALLPATH

#Adding binaries, libraries and a basic localConfig.pm to change

echo "\nDownloading the version to be compiled for CutAdapt, bwa, SAMtools, Picard-Tools, FastQC, GATK, TopHat, Bowtie2 and FASTX-Trimmer"

wget http://bioinfo-web.mpl.ird.fr/toggle/bin.tar.gz

echo "\nDownloading the various Perl modules"

wget http://bioinfo-web.mpl.ird.fr/toggle/perlModules.tar.gz

echo "\nDownloading the localConfig.pm"

wget http://bioinfo-web.mpl.ird.fr/toggle/BAK_localConfig.pm


#######################################################
## Once Toggle has been cloned, this part will be executed:
#######################################################

# DECLARE VARIABLES WITH PATHS 
CURRENTPATH=$INSTALLPATH
TOGGLEPATH=$INSTALLPATH
MODULES=$INSTALLPATH"/Modules"
BINARIES=$INSTALLPATH"/bin"

sleep 1

#COMPILATION
echo "\tDecompressing, compiling and installing softwares"

tar xzvf bin.tar.gz

rm -Rf bin.tar.gz

cd $BINARIES

## compiling bwa
cd bwa

make

cd $BINARIES

##compiling samtools
cd samtools

./configure

make

cd $BINARIES

##compiling picardTools
#NO NEED, JAVA ARCHIVE AVAILABLE

##compiling GATK
#NO NEED, JAVA ARCHIVE AVAILABLE

##compiling cutadapt
cd cutadapt

python setup.py build_ext -i

cd $BINARIES

##compiling fastQC
#NO NEED, JAVA ARCHIVE AVAILABLE

##compiling TopHat
#NO NEED we used a Linux x64 already compiled version from original website

##compiling bowtie and bowtie2
#NO NEED we used a Linux x64 already compiled version from original website

## compiling Fastx_toolkit
#NO NEED we used a Linux x64 already compiled version from original website

sleep 1

echo "\tDecompressing Perl Modules"

cd $INSTALLPATH

tar xvzf perlModules.tar.gz

rm -Rf perlModules.tar.gz

cp -R perlModules/* $MODULES/.

rm -Rf perlModules

sleep 1

echo "\nCONFIGURING YOUR PERSONAL localConfig.pm"

cp BAK_localConfig.pm $MODULES/localConfig.pm
sed -i "s|togglepath|$TOGGLEPATH|g" $MODULES/localConfig.pm
sed -i "s|java7|$JAVASEVEN|g" $MODULES/localConfig.pm
sed -i "s|bwabinary|$BINARIES/bwa/bwa|g" $MODULES/localConfig.pm
sed -i "s|cutadaptbinary|$BINARIES/cutadapt/bin/cutadapt|g" $MODULES/localConfig.pm
sed -i "s|samtoolsbinary|$BINARIES/samtools/samtools|g" $MODULES/localConfig.pm
sed -i "s|picardbinary|$BINARIES/picard-tools|g" $MODULES/localConfig.pm
sed -i "s|fastqcbinary|$BINARIES/FastQC/fastqc|g" $MODULES/localConfig.pm
sed -i "s|GATKbinary|$BINARIES/GenomeAnalysisTK/GenomeAnalysisTK.jar|g" $MODULES/localConfig.pm
sed -i "s|tophat2binary|$BINARIES/tophat/tophat2|g" $MODULES/localConfig.pm
sed -i "s|bowtie2-buildbinary|$BINARIES/bowtie2/bowtie2-build|g" $MODULES/localConfig.pm
sed -i "s|bowtie-buildbinary|$BINARIES/bowtie/bowtie-build|g" $MODULES/localConfig.pm
sed -i "s|fastx_trimmerbinary|$BINARIES/fastx_toolkit/fastx_trimmer|g" $MODULES/localConfig.pm

#Adding toggle in the user PERL5LIB path
sleep 1

echo "\nPERL5LIB=$PERL5LIB:$MODULES\n" | cat - > ~/.bashrc
echo "\nPATH=$PATH:$TOGGLEPATH\n" | cat - > ~/.bashrc

source ~/.bashrc

echo "\nHOORAY !! Configuration finished!\n\nPlease use first the test data as recommanded on the GitHub https://github.com/SouthGreenPlatform/TOGGLE.\n\nThanks for using TOGGLE\n"

exit 0;

