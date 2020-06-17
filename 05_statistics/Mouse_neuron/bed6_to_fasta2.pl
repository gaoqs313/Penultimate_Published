#!/usr/bin/perl -w

open(REF, "/scratch_space/qgao/Penultimate/1.reference/mouse/Genome/mm10.fa");
while(<REF>)
{
	chomp;
	if($_=~m/^>/)
	{
		@a = split(/\s+/, substr($_, 1));
		$id = $a[0];
	}else
	{
		$hash{$id}.=uc($_);
	}
}

open(IN, "neuron_same.tsv");
<IN>;
open(OUT, ">neuron_same.fa");
open(OUT1, ">id_match2.txt");
$i=1;
while(<IN>)
{
	chomp;
	@l = split(/\t/,);
	@m = split(/\:/, $l[0]);
	$chr = $m[0];
	$start = $m[4] - 51;
	$len = $m[5] - $m[4] + 101;
	$strand = $m[9];
	$frag = substr($hash{$chr}, $start, $len);
	if($strand eq "-")
	{
		$frag=~tr/{A,T,C,G}/{T,A,G,C}/;
		$frag=reverse($frag);
	}
	print OUT "\>BID$i\n$frag\n";
	print OUT1 "BID$i\t$l[0]\n";
	$i++;
}
