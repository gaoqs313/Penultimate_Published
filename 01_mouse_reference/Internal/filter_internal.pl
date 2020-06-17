#!/usr/bin/perl -w

# only keep the internal exons after start codon and with length not divisible by 3

open(IN, "/scratch_space/qgao/Penultimate/1.reference/mouse/GTF/Mus_musculus.GRCm38.98.chr.gtf");
while(<IN>)
{
        chomp;
        next if $_=~m/^#/;
        @l = split(/\t/,);
        next if ($l[2] ne "start_codon" && $l[2] ne "stop_codon" && $l[2] ne "exon");
        @m = split(/\"/, $l[8]);
        for($i=0; $i<$#m; $i++)
        {
                if($m[$i]=~m/gene_biotype/)
                {
                        $biotype_idx = $i + 1;
                }elsif($m[$i]=~m/transcript_id/)
                {
                        $transcript_idx = $i + 1;
                }elsif($m[$i]=~m/transcript_biotype/)
                {
                        $transcript_biotype = $i + 1;
                }
        }
        next if $m[$transcript_biotype] ne "protein_coding";
        next if $m[$biotype_idx] ne "protein_coding"; # only consider protein coding transcripts
        if($l[2] eq "start_codon")
        {
                $leftStartCodon{$m[$transcript_idx]} = $l[3];
                $rightStartCodon{$m[$transcript_idx]} = $l[4];
        }
}

open(DA, "internal.novelty.txt");
open(OUT, ">filtered.internal.novelty.txt");
open(OUT2, ">filtered.internal.frame_preserving.novelty.txt");

while(<DA>)
{
	chomp;
	@l = split(/\t/,);
	@m = split(/\:/, $l[0]);
	$strand=substr($l[0], -1);
	if($strand eq "+")
	{
		next if $m[4] <= $rightStartCodon{$l[2]};
	}else
	{
		next if $m[5] >= $leftStartCodon{$l[2]};
	}
	$len=$m[5]-$m[4]+1;
	if($len % 3)
	{
		print OUT "$_\n";
	}else
	{
		print OUT2 "$_\n";
	}
}
