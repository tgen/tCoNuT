#!/usr/bin/perl 
#
# Script to merge two snpEff annotated mpileup VCFs.  Generates table for CNA median centering (merged.vcf.txt).  
# Annotates TUMOR mpileup VCF with LOH information.  Allele shift, test of two proportions, and p-value for test statistics.
#
# ./alleleShiftMpileupVCF.pl ztable.txt CONTROL.mpileup.snpEff.vcf TUMOR.mpileup.snpEff.vcf
#
# Output: alleleshift.vcf
#
# Jessica Aldrich
# TGen
# Sept 18, 2014 
# revised: Dec 2, 2014 - limited hets to qual of 225 and 1000 genomes (KGPROD=true).


#open(ZFILE,"/Users/jaldrich/Data/strexomeCNA/ztable.txt");

$ztable=$ARGV[0];
open(ZFILE,"$ztable");
LOOP: while (<ZFILE>){
	$line=$_;
	chomp $line;
	if ($line=~/^Zscore/){next LOOP};
	@fields=split(/\t/,$line);   
	$ztable{$fields[0]}=$fields[3];	
}
close (ZFILE);

$file1=$ARGV[1];
$file2=$ARGV[2];
open(OFILE,">","alleleshift.vcf");

open(OFILE2,">","alleleshift.txt");
print OFILE2 "CHR\tPOS\tCONTROLDEPTH\tCONTROLALT\tCONTROLALTFREQ\tTUMORDEPTH\tTUMORALT\tBAF\tALLELESHIFT\tZVALUE\tPVALUE\n";

open (FILE1, "$file1") or die "Can't find $file1\n";
LOOP:  while (<FILE1>) {
	$line=$_;
  	if ($line=~/^#/){next LOOP};
  	@fields=split(/\t/,$line);
  	$chr=$fields[0];
  	$pos=$fields[1];
  	$dbsnp=$fields[2];
	$qual=$fields[5];
  	if ($dbsnp=~/^rs/ && $fields[9]=~/^0\/1/ && $qual==225 && $fields[7]=~/KGPROD/) {
    		$lookup{$dbsnp}=1;
    	if ($line=~/;DP4=(.*?);/) {
      		$temp=$1;
      		@DP4=split(/\,/,$temp);
      		$control{$dbsnp}{'ref'}=$DP4[0]+$DP4[1];
      		$control{$dbsnp}{'alt'}=$DP4[2]+$DP4[3];
    		}
  	}
}
close (FILE1);

open (FILE2, "$file2") or die "Can't find $file2\n";
LOOP:  while (<FILE2>) {
	$line=$_;
 
  	if ($line=~/^#/){
		print OFILE "$line";
		next LOOP;
	}
  	@fields=split(/\t/,$line);
  	$chr=$fields[0];
  	$pos=$fields[1];
  	$dbsnp=$fields[2];
  	if ($fields[9]=~/^0\/1/) {
    		if (exists($lookup{$dbsnp})) {
      			if ($line=~/;DP4=(.*?);/) {
        			$temp=$1;
        			@DP4=split(/\,/,$temp);
        			$tumor{$dbsnp}{'ref'}=$DP4[0]+$DP4[1];
        			$tumor{$dbsnp}{'alt'}=$DP4[2]+$DP4[3];
      			}
      			$DP1=$control{$dbsnp}{'ref'}+$control{$dbsnp}{'alt'};
      			$DP2=$tumor{$dbsnp}{'ref'}+$tumor{$dbsnp}{'alt'};

      			$altFreqControl = $control{$dbsnp}{'alt'}/$DP1;
      			$altFreqTumor = $tumor{$dbsnp}{'alt'}/$DP2;	      

     			$altShift = abs($altFreqControl-$altFreqTumor);     

      			$phat = ($control{$dbsnp}{'alt'} + $tumor{$dbsnp}{'alt'})/($DP1+$DP2);
      			$zstat = ($altFreqControl - $altFreqTumor)/sqrt($phat*(1-$phat)*((1/$DP1)+(1/$DP2)));	
			
			if (abs($zstat) > 5){
				$pvalue = "3.13741E-07";
			}else{
				$pvalue = $ztable{sprintf('%.3f',abs($zstat))};			
   			}

      			print OFILE "$fields[0]\t$fields[1]\t$fields[2]\t$fields[3]\t$fields[4]\t$fields[5]\t$fields[6]\t$fields[7];";
      			print OFILE "DPCONTROL=$DP1;DPTUMOR=$DP2;CONTROLALT=$control{$dbsnp}{'alt'};TUMORALT=$tumor{$dbsnp}{'alt'};";
			print OFILE "CONTROLALTFREQ=$altFreqControl;BAF=$altFreqTumor;ALLELESHIFT=$altShift;ZVALUE=$zstat;PVALUE=$pvalue\t";
      			print OFILE "$fields[8]\t$fields[9]";
			
			print OFILE2 "$fields[0]\t$fields[1]\t$DP1\t$control{$dbsnp}{'alt'}\t$altFreqControl\t$DP2\t$tumor{$dbsnp}{'alt'}\t$altFreqTumor\t$altShift\t$zstat\t$pvalue\n"

    		}
    		else{
			print OFILE "$line";
    		}

  	}
  	else{
      		print OFILE "$line";	
	}
}

close(FILE2);
close(OFILE);
close(OFILE2);
