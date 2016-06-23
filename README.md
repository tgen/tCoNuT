# tCoNuT
TGen Copy Number Tool

# About
tCoNuT is a read depth based comparative copy number tool designed for whole genome, exome and panel NGS data. In addition, tCoNuT pipeline provides scripts for calculating B-allele frequencies.  B-allele frequencies can be used in conjuction with copy number to determine regions of LOH and provide additional evidence of a copy number change. The tool requires a control sample which can be matched or unmatched to the affected or tumor sample.

Please see tCoNuT wiki for a diagram of the tCoNuT workflow.

#Requirements
A local installation of Perl and R are required.  

R needs the DNAcopy package found at
https://bioconductor.org/packages/release/bioc/html/DNAcopy.html.

validateCNAVariantsVCF.pl requires Statistics::R module
http://search.cpan.org/~gmpassos/Statistics-R-0.02/lib/Statistics/R.pm

tCoNuT requires the MATLAB Runtime (MCR) v9. Link and instructions for installation are found at http://www.mathworks.com/products/compiler/mcr/. Compiled MATLAB code does not require a MATLAB license just requires the MCR.

tCoNuT pipeline was developed and compiled (specifically MATLAB code) on Linux 64 systems. Most of the scripts are platform independent. Uncompiled MATLAB code (*.m) found in tCoNuT/tCoNuT folder is not platform dependent but would require a license of MATLAB to run.

#Usage
Please refer to tCoNuT workflow for overview and ngs_cna2015.pbs for examples on how to call each script.

<b>Step 1 (prior tCoNuT):</b> Align paired-end sequencing data (BAMs) for each control and affected/tumor samples. Currently, tCoNuT can only be used on human data. Run HaploType Caller(HC) on BAMs then annotate with SnpEff/SnpSift.

<b>Step 2:</b> Create DAT files using tgen_CloneCov.v0092.pl for each BAM.

<b>Step 3:</b> Run parseMergeVCF.pl on HC VCF to get baf.txt and merged.vcf.txt

```
${tCoNuTdir}/parseMergeVCF.pl ${HCVCF} ${CONTROLSAMPLENAME} ${AFFECTEDSAMPLENAME}
```

<b>Step 4:</b> Run tCoNuT with DAT and merged.vcf.txt files

```
##
## Parameters for tCoNuT.  Currently set for EXOME data.
##
TARGETSFILE=Agilent_Clinical_Research_Exome_hs37d5.cna.bed # Copy number BED file of Agilent Clinic Research exome targets
HETFILE=merged.vcf.txt  #   from parseMergeVCF.pl

smWin=6                 #   <<<< THIS CAN BE ADJUSTED - smoothing window size >>>>
fcThresh=0.75           #   <<<< THIS CAN BE ADJUSTED - fold-change threshold - plot >>>>
res=2                   #   <<<< THIS CAN BE ADJUSTED - min resolution >>>>
maxGap=1000             #   <<<< THIS CAN BE ADJUSTED - max distance between consecutive postions >>>>

##
## Average Coverage calculated Picard hsMetrics can be used for minimum depth (${NORMX} and ${TUMORX})
##
hetDepthN=${NORMX}      #   <<<< THIS CAN BE ADJUSTED - min depth of diploid het position >>>>
hetDepthT=${TUMORX}     #   <<<< THIS CAN BE ADJUSTED - min depth of diploid het position >>>>
hetDev=0.025            #   <<<< THIS CAN BE ADJUSTED - allowable deviation from ref allele frequency of 0.5+/-0.025

readDepth=$( echo "$hetDepthN * 3" | bc )  # Set max read depth to 3 x control sample's average target coverage (from Picard HS metrics)

${tCoNuTdir}/tCoNuT/run_tCoNuT.sh ${MCRPATH} ${NORMALDAT} ${TUMORDAT} ${OFILE} ${HETFILE} ${smWin} ${fcThresh} ${assayID} ${res} ${readDepth} ${maxGap} ${hetDepthN} ${hetDepthT} ${hetDev} ${TARGETSFILE}
```

