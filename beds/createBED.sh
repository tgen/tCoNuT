#!/bin/bash
#/*************************************************************************
## *  [2010] - [2016] Translational Genomics Research Institute (TGen)
## *  All Rights Reserved.
## *
## * Major Contributor(s):
##    Jessica Aldrich, Austin Christofferson, David Craig, Jonathan Keats
## * Minor Contributor(s):
##
##  Generate BED file for tCoNuT to be used with exome or panel sequencing. The BED file provides locations
##  of regions to be removed not targeted by sequencing. Segments are in 100 base segments. Fourth column of BED file
##  is an integer indicating if the segment was covered by a target [0 == not covered; >0 == covered].  BEDtools intersect
##  is used to find the overlap between targets file provided with exome sequencing kits and a BED file (patient.B.bed) provided 
##  on the tCoNuT GitHub page. Consider doing an intersect between the assay targets BED file and a BED file containing all the 
##  known exons and using that as the TARGETSFILE. Agilent provides a "padded" targets BED which is the preferred one. 
##
##  
##  Example:
##  CHR StartPosition EndPosition Integer
##  1	0	100	0
##  1	101	200	1
##  ..
## 
##  INPUT:
##	TARGETSFILE is the BED file containing the targets for the exome/panel sequencing (i.e. S07604514_Covered.bed or S07604514_Padded.bed)
##	OUTFILE is name to be used for the output BED file to be used by tCoNuT
##
##  Example:
## 	./createBED.sh S07604514_Padded.bed AgilentV6.cna.bed where S07604514_Padded.bed was downloaded from Agilent.
##
##  Requirements:  BEDtools http://bedtools.readthedocs.io/en/latest/
##
##/

TARGETSFILE=$1
OUTFILE=$2

# pad each side of a target with 100 bases
cat ${TARGETSFILE} | awk '{print $1"\t"$2-100"\t"$3+100}' > padded_targets.bed
bedtools intersect -a patient.B.bed -b padded_targets.bed -c > ${OUTFILE}
