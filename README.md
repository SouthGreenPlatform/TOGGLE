TOGGLE : Toolbox for generic NGS analyses
===========

Licencied under CeCill-C (http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html) and GPLv3 
Intellectual property belongs to IRD, CIRAD and SouthGreen developpement plateform 
Written by Cecile Monat, Christine Tranchant, Ayite Kougbeadjo, Cedric Farcy, Mawusse Agbessi, Marilyne Summo, and Francois Sabot

This suite of modules is aimed to offer a set of wrappers for a lot of tools in NGS analysis, in order to be able to prototype really fastly and efficiently a pipeline for NGS data.

All the tools are provided with a test script companion, to check a correct behaviour and to ensure a production level script.

Examples of pipelines are also provided.

All the modules options are based on a file called "software.config.txt", structured as follows:

        #Config file
        # True comments with "#"
        #Program name code by $NAME
        #under the program name (empty lines do not count), the option are symbolized by key=value (ex -n=5) or key alone (ex m)
        
        #first program config
        $cutadapt
        -e=0.1
        -O=10
        -m=35
        -q=20
        --overlap=7
        
        #2nd program config
        $samtools view single
        -h
        -b
        -F=0x04

The toolbox module will parse this file to send in the different modules the options.
THIS IS THE ONLY MANDATORY MODULE.


-------

Versions Notes


Release 0.2, 14st of March, 2015

Second version release