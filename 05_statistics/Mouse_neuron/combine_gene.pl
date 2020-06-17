#!/usr/bin/perl -w

$dir="/research/rgs01/scratch/users/qgao/Penultimate/10.kallisto";
open(OUT, ">neuron_expr.txt");

my @samples;
open(SAMP, "/scratch_space/qgao/Penultimate/2.data/neuron_SRR_Acc_List.txt");
while(<SAMP>)
{
	chomp;
	@line = split(/\t/,);
	$sample = $line[0];
	push @samples, $sample;
	open(IN, "$dir/$sample/$sample.tsv");
	<IN>;
	while(<IN>)
	{
		chomp;
		@l = split(/\t/,);
		next if ($l[2] eq "Y" || $l[2] eq "MT");
		$hash{$l[0]}{$sample} = $l[3];
	}
}
print OUT "Gene\t".join("\t", @samples)."\n";

foreach $g (sort keys %hash)
{
	print OUT "$g";
	foreach $s (@samples)
	{
		print OUT "\t$hash{$g}{$s}";
	}
	print OUT "\n";
}

