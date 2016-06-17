#!/usr/bin/perl
#
# This script adds DLR values to the *seg.vcf

$stats=$ARGV[0];
$vcf=$ARGV[1];

open(STATS,$stats) or die "Can't find $dlrs\n";
LOOP: while(<STATS>){

        $line=$_;
        chomp $line;

        if ($line=~/^TumorPhysCoverage/){next LOOP};

        @tmp = split(/\t/,$line);
        $dlr = $tmp[1];
}
close(STATS);

open (FILE, $vcf);
$tmpfile="$vcf.dlr.vcf";
open (OFILE,">",$tmpfile);

LOOP: while (<FILE>){

	$line=$_;
	chomp $line;

	if ($line=~/^#/){print OFILE "$line\n";next LOOP}

	print OFILE "$line;DLR=$dlr\n";	

}
close(FILE);
close(OFILE);



