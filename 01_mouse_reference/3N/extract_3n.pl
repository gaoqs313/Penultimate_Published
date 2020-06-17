#!/usr/bin/perl -w

$GTF=$ARGV[0];
open(IN, "$GTF");
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
		}elsif($m[$i]=~m/exon_number/)
		{
			$exon_idx = $i + 1;
		}elsif($m[$i]=~m/transcript_biotype/)
		{
			$transcript_biotype = $i + 1;
		}
	}
	next if $m[$transcript_biotype] ne "protein_coding";
	next if $m[$biotype_idx] ne "protein_coding"; # only consider protein coding transcripts
	$strand{$m[$transcript_idx]} = $l[6];
	if($l[2] eq "stop_codon")
	{
		$leftStopCodon{$m[$transcript_idx]} = $l[3];
		$rightStopCodon{$m[$transcript_idx]} = $l[4];
	}elsif($l[2] eq "start_codon")
	{
                $leftStartCodon{$m[$transcript_idx]} = $l[3];
                $rightStartCodon{$m[$transcript_idx]} = $l[4];
		$numberStartCodon{$m[$transcript_idx]} = $m[$exon_idx];
	}else
	{
		$startExon{$m[$transcript_idx]}{$m[$exon_idx]} = $l[3];
		$endExon{$m[$transcript_idx]}{$m[$exon_idx]} = $l[4];
	}
}

my @G = split(/\//, $GTF);
$GTFN = $G[-1];
open(OUT1, ">StopCodonInLastExon.final.3N.$GTFN");

foreach $t (keys %leftStopCodon)
{
	$n = keys %{$startExon{$t}};
	next if $n < 3; # at least 3 exons
	next if !exists $leftStartCodon{$t}; # start codon annotated
	$p = $n - 1;
	$q = $n - 2;
	next if ($numberStartCodon{$t} == $n || $numberStartCodon{$t} == $p); # start codon in last two exons
	if($startExon{$t}{$n} <= $leftStopCodon{$t} && $endExon{$t}{$n} >= $rightStopCodon{$t})
	{ # stop codon in last exon 
		if(!(($endExon{$t}{$p}-$startExon{$t}{$p}+1)%3))
		{ # penultimate exon divisible by 3
			if($strand{$t} eq "+")
			{ # start codon not in penultimate exon
				if($rightStartCodon{$t} < $startExon{$t}{$p})
				{
					$lastID{$t} = join("\t", $endExon{$t}{$q}, $startExon{$t}{$p}, $endExon{$t}{$p}, $startExon{$t}{$n}, $leftStopCodon{$t}); # upstream stop, penultimate start, penultimate end, downstream start, stop codon
				}
			}else
			{
				if($leftStartCodon{$t} > $endExon{$t}{$p})
				{
					$lastID{$t} = join("\t", $startExon{$t}{$q}, $endExon{$t}{$p}, $startExon{$t}{$p}, $endExon{$t}{$n}, $rightStopCodon{$t});
				}
			}

		}
	}
}

open(IN, "$GTF");
while(<IN>)
{
        chomp;
        next if $_=~m/^#/;
        @l = split(/\t/,);
        next if ($l[2] ne "start_codon" && $l[2] ne "stop_codon" && $l[2] ne "exon");
        @m = split(/\"/, $l[8]);
        for($i=0; $i<$#m; $i++)
        {
                if($m[$i]=~m/transcript_id/)
                {
                        $transcript_idx = $i + 1;
                }
        }	
	next if $_=~m/ENSMUSG00000024411/;
	next if $m[$transcript_idx] eq "ENSMUST00000234466";
	next if $m[$transcript_idx] eq "ENSMUST00000079081";
        next if $m[$transcript_idx] eq "ENSMUST00000105358";
        next if $m[$transcript_idx] eq "ENSMUST00000182896";
        next if $m[$transcript_idx] eq "ENSMUST00000233503";
        next if $m[$transcript_idx] eq "ENSMUST00000238953";
	if(exists $lastID{$m[$transcript_idx]})
	{
		if(exists $hash{$m[$transcript_idx]})
		{ # only output the transcript with the same location once
			print OUT1 "$_\n";
		}else
		{
			if(!exists $mark{$lastID{$m[$transcript_idx]}})
			{
				print OUT1 "$_\n";
				$mark{$lastID{$m[$transcript_idx]}} = "";
				$hash{$m[$transcript_idx]} = "";
			}
		}
	}
}

