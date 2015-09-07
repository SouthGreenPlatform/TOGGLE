#!/bin/bash

###################################################################################################################################
#
# Copyright 2014-2015 IRD-CIRAD-ADNid
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
# Written by Cecile Monat, Christine Tranchant, Ayite Kougbeadjo, Cedric Farcy, Souhila Amanzougarene, Mawusse Agbessi,
# Enrique Ortega-Abboud, Julie Orjuela-Bouniol, Marilyne Summo, and Francois Sabot
#
###################################################################################################################################


###################################################################################################################################
#
# This script will allow a user-space automatic installation of TOGGLE https://github.com/SouthGreenPlatform/TOGGLE
#
###################################################################################################################################


# USE AN INTEGER TO SIGNAL THE NUMBER OF PROGRAMS YOU NEED TO CHECK
# DEBUG -- FOR DEV: total=5
# DEBUG -- FOR PROD: 
total=4

# TEST IF MACHINE IS 64BITS

echo -e "\n\n############################################"
echo -e "##\t Checking if your system is a 64 bits"
echo -e "############################################\n"


MACHINE_TYPE=`uname -m`;
if [ ${MACHINE_TYPE} = "x86_64" ];
then
#	echo -e "$MACHINE_TYPE";
	echo -e "Your system is a 64 bits:\t$MACHINE_TYPE";
else
#	echo -e $MACHINE_TYPE;
	echo -e "** It seems that your system is not 64bits.**";
	echo -e "**Please verify that your system is 64 bits**";
	echo -e "** $MACHINE_TYPE **";
        exit 1;
fi




# TEST IF SOFTWARES ARE INSTALLED, AND PRINT THEIR PATH

echo -e "\n\n############################################"
echo -e "##\t Searching for some required software being installed"
echo -e "############################################\n"


x=0
for i in perl wget git tar

do
#echo "Searching for perl";
	command -v $i > /dev/null 2>&1 && echo -e "$i is installed in :\t"`command -v $i` && x=$(( $x + 1 )) || { echo >&2 "** TOGGLE requires $i but it is not installed.  Please install it or contact your administrator for installation..."; };
#	x=$(( $x + 1 ));
done



echo -e "\n\tYou have all required sotfwares already installed"



if [ $x -lt $total ] ;
then
	echo -e "\nYour system is missing some programs, please install them and try again.";
	exit 1;
elif [ $x -eq $total ]
then
	echo -e "\nRequired UNIX programs have been found  !!\nYou are one step closer to installing TOGGLE"
else
	echo -e "\nThis is a bug message from TOGGLE devs.\nThis is not supposed to happen, you might want to check the raw code or contact us through https://github.com/SouthGreenPlatform/TOGGLE/issues"
	exit 1;
fi


# INPUT VARIABLES FOR JAVA

echo -e "\n\n############################################"
echo -e "##\t Path to Java 7"
echo -e "############################################\n"

echo -e "\nTOGGLE requires a functional version of Java 7"
echo -e "Please input the adress to the executable file\n"


# INPUT JAVA 6 PATH

#echo -e "Type the absolute path for Java 6:"
#read JAVASIX

# echo -e "\n"
# $JAVASIX -version
# echo -e "\n"

#while true; do
#
#	echo -e "\n";
#	$JAVASIX -version;
#	echo -e "\n";
#
#    read -p "IS THIS THE RIGHT JAVA 6 VERSION ?? [Y|N]: " yn
#    case $yn in
#		[Yy]* ) echo "Java 6 path has been establised to: '$JAVASIX' "; break;;
#		[Nn]* ) echo "Please, Try again:"; read JAVASIX ;;
#        * ) echo "Please answer yes or no.";;
#    esac
#done



# INPUT JAVA 7 PATH

echo -e "\nType the absolute path for Java 7:"
read JAVASEVEN


while true; do

	#
	#echo -e "\n";
	#$JAVASEVEN -version;
	#echo -e "\n";

    read -p "IS THIS THE RIGHT JAVA 7 VERSION ?? [Y|N]: " yn
    case $yn in
        [Yy]* ) echo "Java 7 path has been established to: '$JAVASEVEN' "; break;;
        [Nn]* ) echo "Please, Try again:"; read JAVASEVEN ;;
        * ) echo "Please answer 'Y' or 'N':";;
    esac
done

#######################################################
## License infos for all softwares
#######################################################

echo -e "By installing this software you agree to the Licenses, Copyrights or Copylefts of all the softwares used by TOGGLE, as well as the License of TOGGLE.

Please Visit the individual sites for the details of the corresponding Licenses and the citations details:

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
\tTo cite: *to come*
"


#######################################################
## TOGGLE installation per se
#######################################################

echo -e "\nINSTALLING TOGGLE\n";

echo -e "\nPlease provide the installation path:"
read INSTALLPATH


while true; do
    
    echo -e "\n";
	$INSTALLPATH;
    echo -e "\n";
    
    read -p "IS THIS THE RIGHT INSTALL PATH ?? [Y|N]: " yn
    case $yn in
        [Yy]* ) echo "Install path is: '$INSTALLPATH' "; break;;
        [Nn]* ) echo "Please, Try again:"; read INSTALLPATH ;;
        * ) echo "Please answer 'Y' or 'N':";;
    esac
done

mkdir $INSTALLPATH

echo -e "\nCloning the current Git Master Version";

git clone https://github.com/SouthGreenPlatform/TOGGLE.git $INSTALLPATH

cd $INSTALLPATH

echo -e "\nDownloading the compiled version for CutAdapt, bwa, SAMtools, Picard-Tools, FastQC, GATK, TopHat, Bowtie2 and FASTX-Trimmer"

wget http://bioinfo-web.mpl.ird.fr/toggle/bin.zip

echo -e "\nDownloading the various Perl modules"

wget http://bioinfo-web.mpl.ird.fr/toggle/perlModules.zip

echo -e "\nDownloading the localConfig.pm"

wget http://bioinfo-web.mpl.ird.fr/toggle/BAK_localConfig.pm


#######################################################
## Once Toggle has been cloned, this part wil be executed:
#######################################################

# DECLARE VARIABLES WITH PATHS 
CURRENTPATH=$INSTALLPATH
TOGGLEPATH=$INSTALLPATH
MODULES=$INSTALLPATH"/Modules"
BINARIES=$INSTALLPATH"/bin"

echo -e "\nUnzipping compiled versions and Perl modules"

unzip bin.zip

rm -Rf bin.zip

unzip perlModules.zip

rm -Rf perlModules.zip

cp -R perlModules/* $MODULES/.

rm -Rf perlModules

echo -e "\nCONFIGURING YOUR PERSONAL localConfig.pm"

cp BAK_localConfig.pm $MODULES/localConfig.pm
sed -i "s|togglepath|$TOGGLEPATH|g" $MODULES/localConfig.pm
sed -i "s|java7|$JAVASEVEN|g" $MODULES/localConfig.pm
sed -i "s|bwabinary|$BINARIES/bwa\.kit/bwa|g" $MODULES/localConfig.pm
sed -i "s|cutadaptbinary|$BINARIES/cutadapt/bin/cutadapt|g" $MODULES/localConfig.pm
sed -i "s|samtoolsbinary|$BINARIES/samtools/samtools|g" $MODULES/localConfig.pm
sed -i "s|picardbinary|$BINARIES/picard/dist/picard.jar|g" $MODULES/localConfig.pm
sed -i "s|fastqcbinary|$BINARIES/FastQC/fastqc|g" $MODULES/localConfig.pm
sed -i "s|GATKbinary|$BINARIES/GATK/GenomeAnalysisTK.jar|g" $MODULES/localConfig.pm
sed -i "s|tophat2binary|$BINARIES/tophat/tophat2|g" $MODULES/localConfig.pm
sed -i "s|bowtie2-buildbinary|$BINARIES/bowtie2/bowtie2-build|g" $MODULES/localConfig.pm
sed -i "s|bowtie-buildbinary|$BINARIES/bowtie/bowtie-build|g" $MODULES/localConfig.pm
sed -i "s|fastx_trimmerbinary|$BINARIES/fastx_toolkit/bin/fastx_trimmer|g" $MODULES/localConfig.pm


echo -e "\nHOORAY !! Configuration finished!\n\nPlease launch the TEST/all_tests.sh script to be sure the whole installation works properly; for that type directly here in the terminal:\n\ncd /path/to/TOGGLE/TEST && sh all_tests.sh.\n\n You can also use the test data as recommanded on the GitHub https://github.com/SouthGreenPlatform/TOGGLE.\n\nThanks for using TOGGLE\n"

exit 0;

#######################################################
#######################################################



# x=1
# while [ $x -le 5 ]
# do
# 	echo "Welcome $x times"
# 	x=$(( $x + 1 ))
# done


# while true; do
#     read -p "\nDo you wish to install this program?" yn
#     case $yn in
#         [Yy]* ) echo -e "You said YES"; break;;
#         [Nn]* ) echo -e "You said NO";
#         * ) echo "Please answer yes or no.";;
#     esac
# done



##########################################################################################


###### Autres commandes pour tester si un soft est installé.

# command -v foo >/dev/null 2>&1 || { echo >&2 "I require foo but it's not installed.  Aborting."; exit 1; }
# type foo >/dev/null 2>&1 || { echo >&2 "I require foo but it's not installed.  Aborting."; exit 1; }
# hash foo 2>/dev/null || { echo >&2 "I require foo but it's not installed.  Aborting."; exit 1; }

