#!/usr/bin/perl -w

open(IN, "/research/rgs01/scratch/users/qgao/Penultimate/11.revision/A_Mouse_Penultimate/mouse_support_info.txt");
<IN>;
while(<IN>)
{
	chomp;
	@l = split(/\t/,);
	$hash{$l[0]} = $l[3];
}

open(DA, "liftOver.results");
open(OUT, ">separated.liftOver.results");

while(<DA>)
{
	chomp;
	@l = split(/\t/,);
	print OUT "$_\t$hash{$l[0]}\n";
}
