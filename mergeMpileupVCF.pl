#!/usr/bin/perl 
#
#Script to merge two snpEff annotated mpileup VCFs.  Generates table for CNA median centering (merged.vcf.txt).  
#
#./mergeVCF.pl CONTROL.vcf TUMOR.vcf
#
#creates merged.vcf.txt which is a table with columns: chromosome position controlDP controlREFcount controlALTcount tumorDP tumorREFcount tumorALTcount
#



$file1=$ARGV[0];
$file2=$ARGV[1];
open(OFILE,">","merged.vcf.txt");

open (FILE1, "$file1") or die "Can't find $file1\n";
LOOP:  while (<FILE1>) {
  $line=$_;
  if ($line=~/^#/){next LOOP};
  @fields=split(/\t/,$line);
  $chr=$fields[0];
  $pos=$fields[1];
  $dbsnp=$fields[2];
  if ($dbsnp=~/^rs/ && $fields[9]=~/^0\/1/) {
    $lookup{$dbsnp}=1;
    if ($line=~/;DP4=(.*?);/) {
      $temp=$1;
      @DP4=split(/\,/,$temp);
      $control{$dbsnp}{'ref'}=$DP4[0]+$DP4[1];
      #$control{$dbsnp}{'ref'}=$DP4[1];
      $control{$dbsnp}{'alt'}=$DP4[2]+$DP4[3];
      #$control{$dbsnp}{'altR'}=$temp[3];
    }
  }
}
close (FILE1);

open (FILE2, "$file2") or die "Can't find $file2\n";
LOOP:  while (<FILE2>) {
  $line=$_;
  if ($line=~/^#/){next LOOP};
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
        #$tumor{$dbsnp}{'refR'}=$DP4[1];
        $tumor{$dbsnp}{'alt'}=$DP4[2]+$DP4[3];
        #$tumor{$dbsnp}{'altR'}=$DP4[3];
      }
      $DP1=$control{$dbsnp}{'ref'}+$control{$dbsnp}{'alt'};
      $DP2=$tumor{$dbsnp}{'ref'}+$tumor{$dbsnp}{'alt'};
      print OFILE "$chr\t$pos\t$DP1\t$control{$dbsnp}{'ref'}\t$control{$dbsnp}{'alt'}\t";
      print OFILE "$DP2\t$tumor{$dbsnp}{'ref'}\t$tumor{$dbsnp}{'alt'}\n";	
    }
  }
}
close(FILE2);
close(OFILE);
