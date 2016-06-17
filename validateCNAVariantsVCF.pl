#!/packages/perl/5.14.2/bin/perl
##!/usr/bin/perl
##
##  This script incorporates BAF information for each segment in copy number data if available.
##  It calculates the median(MEDIANBAF), median minor allele frequency (min(alt,ref)) for both 
##  control (MEDIANCONTROLMAF) and affected/tumor (MEDIANMAF), median abs(0.5-BAF) (MEDIANBAFDEV), and number
##  of BAFs (NBAF) across a segment. In addition, a Wilcoxon Rank Sum Test (RANKSUMVALUE) and a test of two 
##  proportions (ZSTAT and ZPVALUE) are performed between the control and affected/tumor alternate allele 
##  frequencies. This additional information can be used to support a shift from diploidy in copy number.
##
##  ARGUMENTS:
##	ARGV[0] is SAMPLE.cna.seg.vcf from annotSeg.pl
##	ARGV[1] is baf.txt from parseMergeVCF.pl
##	ARGV[2] is ztable.txt which contains z-scores for ztest.
##
##  OUTPUT:
##	ARGV[0].val.vcf is a VCF file with addition fields in INFO for BAF and tests for each segment.
##
##  validateCNAVariantsVCF.pl cna.seg.vcf baf.txt ztable.txt
##
## *  [2010] - [2016] Translational Genomics Research Institute (TGen)
## *  All Rights Reserved.
## *
## * Major Contributor(s):
##    Jessica Aldrich
## * Minor Contributor(s):


$VERS=1.0;

use List::Util qw( sum min max );
use Statistics::R

$segVCF=$ARGV[0];
$alleleTXT=$ARGV[1];

open(SEGFILE,$segVCF) or die "Can't find $segVCF\n";
open(ALLELEFILE,$alleleTXT) or die "Can't find $alleleTXT\n";
open(OFILE,">","$segVCF.val.vcf");


$ztable=$ARGV[2];
open(ZFILE,"$ztable");
LOOP: while (<ZFILE>){
        $line=$_;
        chomp $line;
        if ($line=~/^Zscore/){next LOOP};
        @fields=split(/\t/,$line);
        $ztable{$fields[0]}=$fields[3];
}
close (ZFILE);


LOOP: while (<ALLELEFILE>){

	$line=$_;
	chomp $line;
	
	if ($line=~/^CHROM/){next LOOP};

	@fields = split(/\t/,$line);

	if($fields[0] eq "X"){$fields[0]=23};
        if($fields[0] eq "Y"){$fields[0]=24};

	$CONTROLMAF = $fields[9];
	$BAF = $fields[10];
	$MAF = $fields[11];
	$BAFDev = $fields[12];

	$TUMOR_MIN = min(($fields[6],$fields[7]));
	$TUMOR_MAX = max(($fields[6],$fields[7]));

	$CONTROL_MIN = min(($fields[3],$fields[4]));
	$CONTROL_MAX = max(($fields[3],$fields[4]));

	$snp{$fields[0]}{$fields[1]}{"CONTROLMAF"} = sprintf('%.4f',$CONTROLMAF);
        $snp{$fields[0]}{$fields[1]}{"BAF"} = sprintf('%.4f',$BAF);
        $snp{$fields[0]}{$fields[1]}{"MAF"} = sprintf('%.4f',$MAF);
        $snp{$fields[0]}{$fields[1]}{"BAFDev"} = sprintf('%.4f',$BAFDev);
	
	$snp{$fields[0]}{$fields[1]}{"TMIN"} = $TUMOR_MIN;
	$snp{$fields[0]}{$fields[1]}{"TMAX"} = $TUMOR_MAX;
	$snp{$fields[0]}{$fields[1]}{"CMIN"} = $CONTROL_MIN;
	$snp{$fields[0]}{$fields[1]}{"CMAX"} = $CONTROL_MAX;
}

close(ALLELEFILE);

LOOP: while (<SEGFILE>){
	
	$line=$_;
	chomp $line;
	
	#EXTRACT DLR from seg.vcf
	if ($line=~/^##DLR/){@dlr=split('=',$line);};
	
	if ($line=~/Deletion/){

		print OFILE "##INFO=<ID=MEDIANBAF,Number=.,Type=Float,Description=\"Median BAF value from across segment\">\n";
		print OFILE "##INFO=<ID=MEDIANMAF,Number=.,Type=Float,Description=\"Median MAF value from across segment\">\n";
		print OFILE "##INFO=<ID=MEDIANCONTROLMAF,Number=.,Type=Float,Description=\"Median CONTROLMAF value from across segment\">\n";
		print OFILE "##INFO=<ID=MEDIANBAFDEV,Number=.,Type=Float,Description=\"Median value of abs(0.5-BAF)\">\n";
		print OFILE "##INFO=<ID=NBAF,Number=.,Type=Integer,Description=\"Number of BAF values from a collection of hets across an identified segment\">\n";
		print OFILE "##INFO=<ID=RANKSUMPVALUE,Number=.,Type=Float,Description=\"P-value from a Wilcoxon Rank Sum Test between MAF and Control MAF\">\n";
		print OFILE "##INFO=<ID=ZSTAT,Number=.,Type=Float,Description=\"Test statistics from a test of 2 proportions between proportion of tumor alternate allele frequencies and normal alternate allele frequencies across a segment\">\n";
		print OFILE "##INFO=<ID=ZPVALUE,Number=.,Type=Float,Description=\"P-value from a test of 2 proportions between proportion of tumor alternate allele frequencies and normal alternate allele frequencies across a segment\">\n";
		

		print OFILE "$line\n";
		next LOOP;

	}

	if ($line=~/^##source/){
		
		print OFILE "$line\n";
		print OFILE "##source=\"validateCNAVariantsVCF.pl v$VERS\"\n";
		next LOOP;
	}

	if ($line=~/^#/){print OFILE "$line\n"; next LOOP};
	
	@fields = split(/\t/,$line);


	$START = $fields[1];
        $END = $fields[2];

        for (my $i=$START; $i <= $END; $i++){
        	if ( exists($snp{$fields[0]}{$i})){
			push(@controlmaf,$snp{$fields[0]}{$i}{"CONTROLMAF"});
	                push(@baf,$snp{$fields[0]}{$i}{"BAF"});
                        push(@maf,$snp{$fields[0]}{$i}{"MAF"});
			push(@bafdev,$snp{$fields[0]}{$i}{"BAFDev"});
	               	push(@tmin,$snp{$fields[0]}{$i}{"TMIN"});
			push(@tmax,$snp{$fields[0]}{$i}{"TMAX"});
			push(@cmin,$snp{$fields[0]}{$i}{"CMIN"});
			push(@cmax,$snp{$fields[0]}{$i}{"CMAX"});
			}
                }
                
		#if(exists($baf[0])){
		if(exists($baf[0]) & @baf > 1){
			$szBAF = @baf;
			$mdControlMAF = sprintf('%.3f',median(@controlmaf));
			$mdMAF = sprintf('%.3f',median(@maf));
                        $mdBAF = sprintf('%.3f',median(@baf));
			$mdBAFDEV = sprintf('%.3f',median(@bafdev));
			$pvalue = sprintf('%.3f',testStat(\@controlmaf,\@maf));
			$zstat = sprintf('%.3f',testOfProp(\@tmin,\@tmax,\@cmin,\@cmax));
			
			if (abs($zstat) > 5){
                                $zpvalue = "3.13741E-07";
                        }else{
                                $zpvalue = sprintf('%.3f',$ztable{sprintf('%.3f',abs($zstat))});
                        }


                        print OFILE "$line";
                        print OFILE ";MEDIANBAF=$mdBAF;MEDIANMAF=$mdMAF;MEDIANCONTROLMAF=$mdControlMAF;MEDIANBAFDEV=$mdBAFDEV;NBAF=$szBAF;RANKSUMPVALUE=$pvalue;ZSTAT=$zstat;ZPVALUE=$zpvalue\n";

                }else{
                        print OFILE "$line\n";
                }
                undef @baf;
                undef @bafdev;
		undef @maf;
		undef @controlmaf;
		undef @tmin;
		undef @tmax;
		undef @cmin;
		undef @cmax;	
}



sub mean {
    return sum(@_)/@_;
}

sub median
{
    my @vals = sort {$a <=> $b} @_;
    my $len = @vals;
    if($len%2) #odd?
    {
        return $vals[int($len/2)];
    }
    else #even
    {
        return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
    }
}

sub testStat
{
	my ($data1,$data2) = @_;
	my @data1 = @{$data1};
	my @data2 = @{$data2};
	
	my $R1 = Statistics::R -> new ();
	$R1 -> set('data1',\@data1);
	$R1 -> set('data2',\@data2);
	$R1-> run(q`w<-wilcox.test(data1,data2)`);
	my $pvalue = $R1->get('w$p.value'); 
	$R1 -> stop();	
	return $pvalue;



}

sub testOfProp
{

        my ($tmin,$tmax,$cmin,$cmax) = @_;
        my @tmin = @{$tmin};
        my @tmax = @{$tmax};
        my @cmin = @{$cmin};
        my @cmax = @{$cmax};

	$ControlAlt = sum(@cmin);
	$TumorAlt = sum(@tmin);

        $ControlDP=sum(@cmin)+sum(@cmax);
        $TumorDP=sum(@tmin)+sum(@tmax);

        $ControlAltFreq = $ControlAlt/$ControlDP;
        $TumorAltFreq = $TumorAlt/$TumorDP;


	$num = $ControlAltFreq - $TumorAltFreq;
	$pCommon = ($ControlAlt + $ TumorAlt)/($ControlDP + $TumorDP);
	$denom = sqrt($pCommon * (1-$pCommon) * ((1/$ControlDP) + (1/$TumorDP)));
	$zstat = $num/$denom;
	return $zstat; 	
}
