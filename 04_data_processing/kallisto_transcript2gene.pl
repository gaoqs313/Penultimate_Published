#!/usr/bin/perl -w

open(IN, "$ARGV[0]");
open(OUT, ">$ARGV[1]");
<IN>;

while(<IN>)
{
	chomp;
	@l = split(/\t/,);
	@a = split(/\_\_/, $l[0]);
	$hash{$a[0]} += $l[4];
	$chr{$a[0]} = $a[3];
}

print OUT "GeneID\tEnsemblID\tChr\tTPM\n";
open(REF, "/scratch_space/qgao/Penultimate/1.reference/mouse/GTF/Mus_musculus.GRCm38.98.geneid.txt");
<REF>;

while(<REF>)
{
	chomp;
	@l = split(/\t/,);
	print OUT "$l[1]\t$l[0]\t$chr{$l[0]}\t$hash{$l[0]}\n";
}
