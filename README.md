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
Step 1 (prior tCoNuT): Aligned paired-end sequencing files (BAMs) for each control and affected/tumor sample. Currently, tCoNuT can only be used on human data. Run HaploType Caller(HC) on BAMs then annotate with SnpEff/SnpSift.
Step 2: Create DAT files using tgen_CloneCov.v0092.pl for each BAM.
Step 3: Run parseMergeVCF.pl on HC VCF to get baf.txt and merged.vcf.txt
