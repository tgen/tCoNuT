# tCoNuT
TGen Copy Number Tool

# About
tCoNuT is a read depth based comparative copy number tool designed for whole genome, exome and panel NGS data. In addition, tCoNuT pipeline provides scripts for calculating B-allele frequencies.  B-allele frequencies can be used in conjuction with copy number to determine regions of LOH and provide additional evidence of a copy number change. The tool requires a control sample which can be matched or unmatched to the affected or tumor sample.

Please see tCoNuT wiki for a diagram of the tCoNuT workflow.

#Requirements
An installation of Perl and R are required.  R needs the DNAcopy package found at https://bioconductor.org/packages/release/bioc/html/DNAcopy.html.

validateCNAVariantsVCF.pl requires Statistics::R module http://search.cpan.org/~gmpassos/Statistics-R-0.02/lib/Statistics/R.pm

tCoNuT requires the MATLAB Runtime (MCR) v9. Link and instructions for installation are found at http://www.mathworks.com/products/compiler/mcr/. 
