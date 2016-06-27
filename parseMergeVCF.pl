#!/packages/perl/5.14.2/bin/perl
##
## This script parses a SnpEff/SnpSift annotated HaplotypeCaller VCF run on
## control and affected/tumor samples together. It identifies heterzygous SNPs 
## in the control sample and grabs the corresponding read counts and allele frequencies 
## from both control and affected/tumor samples at the SNP position. In order to filter 
## for the high quality common SNPs, the VCF needs to be annotated with dbSNP GMAF/CAF 
## which can be done using SnpSift. This script was written specifically for HaplotypeCaller/SnpEff
## but may work on VCFs from other variant callers.
##
## gatk's HaploType Caller https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php
## SnpEff/SnpSift http://snpeff.sourceforge.net/ 
##
## parseMergeVCF.pl normal-tumor.hc.snpEff.vcf normalSampleName tumorSampleName readDepth
##
## ARGUMENTS:
##	[0] HaploType Caller SnpEff annotated VCF
##	[1] control sample name used in HaplotypeCaller
##      [2] affected/tumor sample name used in HaplotypeCaller
##	[3] read depth is the min depth of SNP considered for BAF
##
## OUTPUT:
##	normal-tumor.baf.txt is text file with header: Chromosome Position NormalReadDepth NormalRefAllele NormalAltAllele TumorReadDepth TumorRefAllele TumorAltAllele CONTROLBAF CONTROLMAF BAF MAF ABSBAFDEV
##      merged.vcf.txt is a text file with header: Chromosome Position NormalReadDepth NormalRefAllele NormalAltAllele TumorReadDepth TumorRefAllele TumorAltAllele
##
## *  [2010] - [2016] Translational Genomics Research Institute (TGen)
## *  All Rights Reserved.
## *
## * Major Contributor(s):
##    Jessica Aldrich
## * Minor Contributor(s):

use List::Util qw( sum min );
use List::MoreUtils 'first_index';

open(VCF,"$ARGV[0]");
open(OFILE,">","merged.vcf.txt");

$normalName=$ARGV[1];
$tumorName=$ARGV[2];
$readDepth=$ARGV[3]

open(OFILE2,">","$normalName-$tumorName.baf.txt");

print OFILE "Chromosome\tPosition\tNormalReadDepth\tNormalRefAllele\tNormalAltAllele\tTumorReadDepth\tTumorRefAllele\tTumorAltAllele\n";
print OFILE2 "Chromosome\tPosition\tNormalReadDepth\tNormalRefAllele\tNormalAltAllele\tTumorReadDepth\tTumorRefAllele\tTumorAltAllele\tCONTROLBAF\tCONTROLMAF\tBAF\tMAF\tABSBAFDEV\n";


LOOP:while ($line=<VCF>){

	chomp($line);

	if ($line=~/^##/){next LOOP};
	if ($line=~/^#CHROM/){	
		@tmp=split(/\t+/,$line);
		$szTMP=@tmp;
                for ( my $i=9; $i <= $szTMP; $i++){
                if ($tmp[$i] eq $normalName){
                        $norm=$i;
                }
                if ($tmp[$i] eq $tumorName){
                        $tumor=$i;
                }
                }

	}
	
	@temp=split(/\t+/,$line);

	if ($temp[9] eq '.' || $temp[10] eq '.'){next LOOP};
	
	@COL8=split(/:/,$temp[8]);	
	
	$GTind = first_index { /GT/ } @COL8;
	$ADind = first_index { /AD/ } @COL8;
	$DPind = first_index { /DP/ } @COL8;
	if ($DPind == -1){next LOOP};	

	@control = split(/:/,$temp[$norm]);
	@affected = split(/:/,$temp[$tumor]);
	
	@con_allele=split(/,/,$control[$ADind]);
	@aff_allele=split(/,/,$affected[$ADind]);	

	if ( $temp[7] =~ /GMAF=/ || $temp[7] =~ /CAF=/ ){
	
		if ($temp[7] =~ /GMAF=/) {
			$temp[7] =~ /GMAF=(.*?);/;
			$GMAF = $1;
		}

		if ($temp[7] =~ /CAF=/ ){
			$temp[7] =~ /CAF=(.*?);/;
			@CAF = split(/,/,$1);
			@sortCAF = sort { $b <=> $a } @CAF;
			$GMAF = $sortCAF[1];
		}

		$charR = length($temp[3]);
		$charA = length($temp[4]);
	
		# Write out baf.txt
		if ( $GMAF > 0.05 && $control[0] eq '0/1' && $charR < 2 && $charA < 2 && sum(@con_allele) > $readDepth && sum(@aff_allele) > $readDepth && $con_allele[1]/sum(@con_allele) > 0.05) {
	
			$baf = $aff_allele[1]/sum(@aff_allele);
			$maf = min(@aff_allele)/sum(@aff_allele);
			$control_baf = $con_allele[1]/sum(@con_allele);
			$control_maf = min(@con_allele)/sum(@con_allele);
			$bafdev = abs(0.5-$baf);		

			if($temp[0] eq "X"){$temp[0]=23};
        		if($temp[0] eq "Y"){$temp[0]=24};

                	print OFILE2 "$temp[0]\t$temp[1]\t$control[$DPind]\t$con_allele[0]\t$con_allele[1]\t$affected[$DPind]\t$aff_allele[0]\t$aff_allele[1]\t$control_baf\t$control_maf\t$baf\t$maf\t$bafdev\n";

        	}		
		
		# Write out merged.vcf.txt
		if ( $GMAF > 0.05 && $temp[2] =~ /^rs.*/ && $control[0] eq '0/1' && sum(@con_allele) > 10 && sum(@aff_allele) > 10) {
		
		        print OFILE "$temp[0]\t$temp[1]\t$control[$DPind]\t$con_allele[0]\t$con_allele[1]\t$affected[$DPind]\t$aff_allele[0]\t$aff_allele[1]\n";
	
	        }
	}
}

close(VCF);
close(OFILE);
close(OFILE2);
