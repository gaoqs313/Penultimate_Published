#!/usr/bin/perl -w

open(IN, "../A_Annotation/AA_content.tsv");
while(<IN>)
{
	chomp;
	@l = split(/\t/,);
	next if $l[1] ne "P";
	if($l[2] eq "fl")
	{
		$old{$l[0]} = $l[3];
	}elsif($l[2] eq "fs")
	{
		$new{$l[0]} = $l[3];
	}
}

open(IN, "../A_Annotation/StopCodonInPenultimateExon.Mus_musculus.GRCm38.98.chr.gtf.uniq.filtered.aa");
while(<IN>)
{
	chomp;
	@l = split(/\_\_/,);
	$id = substr($l[0], 1);
	$next = <IN>;
	chomp $next;
	$next = substr($next, 0, length($next) -1);
        <IN>;
        $nxt = <IN>;
        chomp $nxt;
        $nxt = substr($nxt, 0, length($nxt) -1);
        <IN>;<IN>;<IN>;<IN>;<IN>;<IN>;<IN>;<IN>;
        $oldaa{$id} = $next;
        $newaa{$id} = $nxt;
}

open(IN, "../A_Annotation/StopCodonInLastExon.Mus_musculus.GRCm38.98.chr.gtf.uniq.filtered.aa");
while(<IN>)
{
        chomp;
        @l = split(/\_\_/,);
        $id = substr($l[0], 1);
        $next = <IN>;
        chomp $next;
        $next = substr($next, 0, length($next) -1);
	<IN>;
	$nxt = <IN>;
        chomp $nxt;
        $nxt = substr($nxt, 0, length($nxt) -1);
	<IN>;<IN>;<IN>;<IN>;<IN>;<IN>;<IN>;<IN>;
        $oldaa{$id} = $next;
        $newaa{$id} = $nxt;
}

open(REF, "/scratch_space/qgao/Penultimate/1.reference/mouse/GTF/Mus_musculus.GRCm38.98.chr.gtf");
while(<REF>)
{
	chomp;
	next if $_=~m/^#/;
	@l = split(/\t/,);
	next if $l[2] ne "transcript";
	@a = split(/\"/, $l[8]);
	$trans{$a[5]} = $a[9];
}

open(RNA, "/scratch_space/qgao/Penultimate/7.figures/mouse_support_info.txt");
<RNA>;
while(<RNA>)
{
	chomp;
	@l = split(/\t/,);
	$rna{$l[0]} = ($l[4] == 0 ) ? "NO" : "YES";
}

open(IN, "../B_MISO/penultimate.novelty.txt");
open(OUT, ">Mouse_supp.txt");
print OUT "ExonID_in_MISO_format\tTranscriptID\tGeneID\tFulllengthAA\tFrameshiftAA\tFulllengthProlineContent\tFrameshiftProlineContent\tAnnotated\tStopCodonPosition\tSupportedByRNASeq\n";
while(<IN>)
{
        chomp;
        @l = split(/\t/,);
	$anno= ($l[1] ne "Novel") ? "YES" : "NO";
	$position = ($l[2] eq "Penultimate") ? "Penultimate" : "Ultimate";
	print OUT "$l[0]\t$l[3]\t$trans{$l[3]}\t$oldaa{$l[3]}\t$newaa{$l[3]}\t$old{$l[3]}\t$new{$l[3]}\t$anno\t$position\t$rna{$l[0]}\n";
}
