#!/bin/bash

NORMAL=$1
TUMOR=$2


qsub -v DIR=${PWD},NORMALSAMPLE=${NORMAL},TUMORSAMPLE=${TUMOR},VCF=${NORMAL}-${TUMOR}.HC_All.snpEff.vcf,NORMALDAT=${NORMAL}.proj.md.jr.bam.clc.cln.dat,TUMORDAT=${TUMOR}.proj.md.jr.bam.clc.cln.dat,OFILE=${TUMOR} /home/jaldrich/scripts/cna2015/ngs_cna2015.pbs
