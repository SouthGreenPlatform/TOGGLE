#! /bin/bash

rm DATA/referenceIrigin.dict DATA/referenceIrigin.fasta.amb DATA/referenceIrigin.fasta.ann DATA/referenceIrigin.fasta.bwt DATA/referenceIrigin.fasta.fai DATA/referenceIrigin.fasta.pac DATA/referenceIrigin.fasta.sa

rm -rf DATA/BamDirectory

rm -rf DATA/iriginTest/*
cp DATAbak/iriginTest/* DATA/iriginTest/

