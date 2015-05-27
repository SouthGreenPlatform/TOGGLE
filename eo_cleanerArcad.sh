#! /bin/bash

rm DATA/referenceArcad.dict DATA/referenceArcad.fasta.amb DATA/referenceArcad.fasta.ann DATA/referenceArcad.fasta.bwt DATA/referenceArcad.fasta.fai DATA/referenceArcad.fasta.pac DATA/referenceArcad.fasta.sa

rm -rf DATA/BamDirectory

rm -rf DATA/arcadTest/*
cp DATAbak/arcadTest/* DATA/arcadTest/

