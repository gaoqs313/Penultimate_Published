#!/usr/bin/perl -w

$dir="/scratch_space/qgao/Penultimate/5.miso/mouse";
open(OUT, ">neuron_psi.txt");

my @samples;
open(SAMP, "/scratch_space/qgao/Penultimate/2.data/neuron_SRR_Acc_List.txt");
while(<SAMP>)
{
	chomp;
	@line = split(/\t/,);
	$sample = $line[0];
	push @samples, $sample;
	open(IN, "$dir/$sample/summary/$sample.miso_summary");
	<IN>;
	while(<IN>)
	{
		chomp;
		@l = split(/\t/,);
		@m = split(/\,/, $l[6]);
		$num = 0;
		for($i=0;$i<=$#m;$i++)
		{
			$num += substr($m[$i], 2);
		}
		if($num >= 20)
		{
			$psi{$l[0]}{$sample} = $l[1] * 100;
		}
	}
}
print OUT "Event\t".join("\t", @samples)."\n";

open(REF, "/scratch_space/qgao/Penultimate/1.reference/mouse/B_MISO/penultimate.novelty.txt");
while(<REF>)
{
        chomp;
        @l = split(/\t/,);
	$id = $l[0];
	print OUT "$id";
	foreach $s (@samples)
	{
		$p = (exists $psi{$id}{$s}) ? $psi{$id}{$s} : NA;
		print OUT "\t$p";
	}
	print OUT "\n";
}


