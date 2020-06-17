#!/usr/bin/perl -w

open(IN, "liftOver.results");
open(OUT, ">format.liftOver.results.txt");

print OUT "Mouse_ExonID\tMouse_TranscriptID\tMouse_GeneID\tMouse_SupportedByRNASeq\tHuman_ExonID\tHuman_TranscriptID\tHuman_GeneID\tHuman_SupportedByRNASeq\n";

while(<IN>)
{
	chomp;
	@l = split(/\t/,);
	$mouse = ($l[3] == 0) ? "NO" : "YES";
	$human = ($l[7] == 0) ? "NO" : "YES";
	print OUT "$l[0]\t$l[1]\t$l[2]\t$mouse\t$l[4]\t$l[5]\t$l[6]\t$human\n";
}
