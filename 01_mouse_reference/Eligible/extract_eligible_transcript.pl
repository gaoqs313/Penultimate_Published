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
	$chr{$m[$transcript_idx]} = $l[0];
	if($l[2] eq "stop_codon")
	{
		$leftStopCodon{$m[$transcript_idx]} = $l[3];
		$rightStopCodon{$m[$transcript_idx]} = $l[4];
		$numberStopCodon{$m[$transcript_idx]} = $m[$exon_idx];
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
open(OUT1, ">StopCodonInLastExon.$GTFN");
open(OUT2, ">StopCodonInPenultimateExon.$GTFN");

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
		if(($endExon{$t}{$p}-$startExon{$t}{$p}+1)%3)
		{ # penultimate exon not divisible by 3
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
	}elsif($startExon{$t}{$p} <= $leftStopCodon{$t} && $endExon{$t}{$p} >= $rightStopCodon{$t})
	{ # stop codon in penultimate exon
		if($strand{$t} eq "+")
		{
			$penultimateID{$t} = join("\t", $endExon{$t}{$q}, $startExon{$t}{$p}, $endExon{$t}{$p}, $startExon{$t}{$n}, $leftStopCodon{$t});
		}else
		{
			$penultimateID{$t} = join("\t", $startExon{$t}{$q}, $endExon{$t}{$p}, $startExon{$t}{$p}, $endExon{$t}{$n}, $rightStopCodon{$t});
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
	}elsif(exists $penultimateID{$m[$transcript_idx]})
	{
                if(exists $hash2{$m[$transcript_idx]})
                {
                        print OUT2 "$_\n";
                }else
                {
                        if(!exists $mark{$penultimateID{$m[$transcript_idx]}})
                        {
                                print OUT2 "$_\n";
                                $mark{$penultimateID{$m[$transcript_idx]}} = "";
				$hash2{$m[$transcript_idx]} = "";
                        }
                }
	}
}

open(FASTA, "$ARGV[1]");
while(<FASTA>)
{
        chomp;
        if($_=~m/^>/)
        {
                @l = split(/\s+/, substr($_, 1));
                $id = $l[0];
        }else
        {
                $fasta{$id}.= uc($_);
        }
}

open(OUT3, ">StopCodonInLastExon.$GTFN.fa");
open(OUT4, ">StopCodonInPenultimateExon.$GTFN.fa");
open(OUT5, ">StopCodonInLastExon.$GTFN.aa");
open(OUT6, ">StopCodonInPenultimateExon.$GTFN.aa");
open(OUT7, ">StopCodonInLastExon.$GTFN.uniq.aa");
open(OUT8, ">StopCodonInPenultimateExon.$GTFN.uniq.aa");
open(OUT9, ">StopCodonInLastExon.$GTFN.uniq.filtered.aa");
open(OUT10, ">StopCodonInPenultimateExon.$GTFN.uniq.filtered.aa");

sub codon2aa{
    my ($codon)=@_;
    $codon= uc $codon;
    my %genetic_code = (
        'TTT'=>'F', #Phenylalanine
        'TTC'=>'F', #Phenylalanine
        'TTA'=>'L', #Leucine
        'TTG'=>'L', #Leucine
        'TCA'=>'S', #Serine
        'TCC'=>'S', #Serine
        'TCG'=>'S',  #Serine
        'TCT'=>'S', #Serine
        'TAC'=>'Y', #Tyrosine
        'TAT'=>'Y', #Tyrosine
        'TAA'=>'_', #Stop
        'TAG'=>'_', #Stop  
        'TGC'=>'C', #Cysteine
        'TGT'=>'C', #Cysteine
        'TGA'=>'_', #Stop
        'TGG'=>'W', #Tryptophan
        
        'CTA'=>'L', #Leucine
        'CTC'=>'L', #Leucine
        'CTG'=>'L', #Leucine
        'CTT'=>'L', #Leucine
        'CCA'=>'P', #Proline
        'CCC'=>'P', #Proline
        'CCG'=>'P', #Proline
        'CCT'=>'P', #Proline
        'CAT'=>'H', #Histidine
        'CAC'=>'H', #Histidine
        'CAA'=>'Q', #Glutamine
        'CAG'=>'Q', #Glutamine
        'CGA'=>'R', #Arginine
        'CGC'=>'R', #Arginine
        'CGG'=>'R', #Arginine
        'CGT'=>'R', #Arginine
        
        'ATA'=>'I', #Isoleucine
        'ATC'=>'I', #Isoleucine
        'ATT'=>'I', #Isoleucine
        'ATG'=>'M', #Methionine
        'ACA'=>'T', #Threonine
        'ACC'=>'T', #Threonine
        'ACG'=>'T', #Threonine
        'ACT'=>'T', #Threonine
        'AAC'=>'N', #Asparagine
        'AAT'=>'N', #Asparagine
        'AAA'=>'K', #Lysine
        'AAG'=>'K', #Lysine
        'AGC'=>'S', #Serine#Valine
        'AGT'=>'S', #Serine
        'AGA'=>'R', #Arginine
        'AGG'=>'R', #Arginine
        
        'GTA'=>'V', #Valine
        'GTC'=>'V', #Valine
        'GTG'=>'V', #Valine
        'GTT'=>'V', #Valine
        'GCA'=>'A', #Alanine
        'GCC'=>'A', #Alanine
        'GCG'=>'A', #Alanine
        'GCT'=>'A', #Alanine
        'GAC'=>'D', #Aspartic Acid
        'GAT'=>'D', #Aspartic Acid
        'GAA'=>'E', #Glutamic Acid
        'GAG'=>'E', #Glutamic Acid
        'GGA'=>'G', #Glycine
        'GGC'=>'G', #Glycine
        'GGG'=>'G', #Glycine
        'GGT'=>'G', #Glycine
    );
    if(exists $genetic_code{$codon}){
        return $genetic_code{$codon};
    }
    else{
    		return "*";
    }
}

sub getProtein{
    my ($dna)=@_;
    my $protein='';
    my $codon3;
    for(my $i=0; $i<(length($dna)-2); $i+=3){			
        $codon3=substr($dna,$i,3);
        $protein.=codon2aa($codon3);
    }
    return $protein;
}

foreach $t (keys %hash)
{
        $n = keys %{$startExon{$t}};
	my (@starts, @ends, $anno, $skip);
	my ($anno_all, $anno_last, $anno_penult, $anno_utr);
	if($strand{$t} eq "+")
	{
		for($i=$numberStartCodon{$t}; $i<=$numberStopCodon{$t}; $i++)
		{
			if($i == $numberStartCodon{$t})
			{
				$anno.= substr($fasta{$chr{$t}}, $leftStartCodon{$t}-1, $endExon{$t}{$i}-$leftStartCodon{$t}+1);
				$skip.= substr($fasta{$chr{$t}}, $leftStartCodon{$t}-1, $endExon{$t}{$i}-$leftStartCodon{$t}+1);
				$anno_all+= $endExon{$t}{$i}-$leftStartCodon{$t}+1;
			}elsif($i == $numberStopCodon{$t})
			{
				$anno.= substr($fasta{$chr{$t}}, $startExon{$t}{$i}-1, $rightStopCodon{$t}-$startExon{$t}{$i}+1);
				$skip.= substr($fasta{$chr{$t}}, $startExon{$t}{$i}-1, $rightStopCodon{$t}-$startExon{$t}{$i}+3001);
				$anno_all+= $rightStopCodon{$t}-$startExon{$t}{$i}+1;
				$anno_last+= $rightStopCodon{$t}-$startExon{$t}{$i}+1;
				$anno_utr = substr($fasta{$chr{$t}}, $rightStopCodon{$t}, 302);
			}elsif($i == $numberStopCodon{$t}-1)
			{
				$anno.= substr($fasta{$chr{$t}}, $startExon{$t}{$i}-1, $endExon{$t}{$i}-$startExon{$t}{$i}+1);
                                $anno_all+= $endExon{$t}{$i}-$startExon{$t}{$i}+1;
				$anno_penult+= $endExon{$t}{$i}-$startExon{$t}{$i}+1;
			}else
			{
				$anno.= substr($fasta{$chr{$t}}, $startExon{$t}{$i}-1, $endExon{$t}{$i}-$startExon{$t}{$i}+1);
				$skip.= substr($fasta{$chr{$t}}, $startExon{$t}{$i}-1, $endExon{$t}{$i}-$startExon{$t}{$i}+1);
				$anno_all+= $endExon{$t}{$i}-$startExon{$t}{$i}+1;
			}
		}
	}else
	{
                for($i=$numberStopCodon{$t}; $i>=$numberStartCodon{$t}; $i--)
                {
                        if($i == $numberStopCodon{$t})
                        {
                                $anno.= substr($fasta{$chr{$t}}, $leftStopCodon{$t}-1, $endExon{$t}{$i}-$leftStopCodon{$t}+1);
                                $skip.= substr($fasta{$chr{$t}}, $leftStopCodon{$t}-3001, $endExon{$t}{$i}-$leftStopCodon{$t}+3001);
				$anno_all+= $endExon{$t}{$i}-$leftStopCodon{$t}+1;
				$anno_last+= $endExon{$t}{$i}-$leftStopCodon{$t}+1;
				$anno_utr = substr($fasta{$chr{$t}}, $leftStopCodon{$t}-303, 302);
                        }elsif($i == $numberStartCodon{$t})
                        {
                                $anno.= substr($fasta{$chr{$t}}, $startExon{$t}{$i}-1, $rightStartCodon{$t}-$startExon{$t}{$i}+1);
                                $skip.= substr($fasta{$chr{$t}}, $startExon{$t}{$i}-1, $rightStartCodon{$t}-$startExon{$t}{$i}+1);
				$anno_all+= $rightStartCodon{$t}-$startExon{$t}{$i}+1;
                        }elsif($i == $numberStopCodon{$t}-1)
                        {
                                $anno.= substr($fasta{$chr{$t}}, $startExon{$t}{$i}-1, $endExon{$t}{$i}-$startExon{$t}{$i}+1);
				$anno_all+= $endExon{$t}{$i}-$startExon{$t}{$i}+1;
				$anno_penult+= $endExon{$t}{$i}-$startExon{$t}{$i}+1;
                        }else
                        {
                                $anno.= substr($fasta{$chr{$t}}, $startExon{$t}{$i}-1, $endExon{$t}{$i}-$startExon{$t}{$i}+1);
                                $skip.= substr($fasta{$chr{$t}}, $startExon{$t}{$i}-1, $endExon{$t}{$i}-$startExon{$t}{$i}+1);
				$anno_all+= $endExon{$t}{$i}-$startExon{$t}{$i}+1;
                        }
                }
                $anno=~tr/[AGCT]/[TCGA]/;
                $anno=reverse($anno);
                $skip=~tr/[AGCT]/[TCGA]/;
                $skip=reverse($skip);
		$anno_utr=~tr/[AGCT]/[TCGA]/;
		$anno_utr=reverse($anno_utr);
	}
	next if length($anno) % 3; # wrong annotation: coding length not divisible by 3
	$anno_aa=getProtein($anno);
	$skip_aa=getProtein($skip);
	$anno_aa_stop=substr($anno_aa, 0,  index($anno_aa, "_") + 1);
	$skip_aa_stop=substr($skip_aa, 0,  index($skip_aa, "_") + 1);
	next if $anno_aa ne $anno_aa_stop; # wrong annotation: premature stop codone
        print OUT3 "\>$t\_\_$strand{$t}\_\_Penultimate_splice_in\n$anno\n\>$t\_\_$strand{$t}\_\_Penultimate_splice_out\n$skip\n";
	print OUT5 "\>$t\_\_$strand{$t}\_\_Penultimate_splice_in\n$anno_aa\n\>$t\_\_$strand{$t}\_\_Penultimate_splice_out\n$skip_aa_stop\n";
        $anno_aa_notwo_len = int(($anno_all - $anno_penult - $anno_last) / 3);
        $skip_aa_uniq = substr($skip_aa_stop, $anno_aa_notwo_len);
        $anno_aa_uniq = substr($anno_aa, $anno_aa_notwo_len);
	$anno_aa_cds0 = (length($anno_aa) >=100) ? substr($anno_aa, length($anno_aa)-100) : $anno_aa; # last 100 AA 
	$anno_aa_utr0 = getProtein(substr($anno_utr,0,300));
        $anno_aa_utr1 = getProtein(substr($anno_utr,1,300));
        $anno_aa_utr2 = getProtein(substr($anno_utr,2,300));
	$anno_aa_last_len = int(($anno_last + 2)/3);
	print OUT7 "\>$t\_\_$strand{$t}\_\_Penultimate_splice_in\_\_$anno_aa_last_len\n$anno_aa_uniq\n\>$t\_\_$strand{$t}\_\_Penultimate_splice_out\_\_$anno_aa_last_len\n$skip_aa_uniq\n";
        print OUT7 "\>$t\_\_$strand{$t}\_\_Penultimate_splice_in_tail100\n$anno_aa_cds0\n\>$t\_\_$strand{$t}\_\_Penultimate_splice_out_utr0\n$anno_aa_utr0\n\>$t\_\_$strand{$t}\_\_Penultimate_splice_out_utr1\n$anno_aa_utr1\n\>$t\_\_$strand{$t}\_\_Penultimate_splice_out_utr2\n$anno_aa_utr2\n";
	if(length($skip_aa_uniq)>20 && length($skip_aa_uniq)-$anno_aa_last_len>=10)
	{
		$final{$t} = "";
		print OUT9 "\>$t\_\_$strand{$t}\_\_Penultimate_splice_in\_\_$anno_aa_last_len\n$anno_aa_uniq\n\>$t\_\_$strand{$t}\_\_Penultimate_splice_out\_\_$anno_aa_last_len\n$skip_aa_uniq\n";
	        print OUT9 "\>$t\_\_$strand{$t}\_\_Penultimate_splice_in_tail100\n$anno_aa_cds0\n\>$t\_\_$strand{$t}\_\_Penultimate_splice_out_utr0\n$anno_aa_utr0\n\>$t\_\_$strand{$t}\_\_Penultimate_splice_out_utr1\n$anno_aa_utr1\n\>$t\_\_$strand{$t}\_\_Penultimate_splice_out_utr2\n$anno_aa_utr2\n";
	}
}

foreach $t (keys %hash2)
{
        $n = keys %{$startExon{$t}};
        my (@starts, @ends, $anno, $skip);
        my ($anno_all, $anno_last, $anno_utr);
        if($strand{$t} eq "+")
        {
                for($i=$numberStartCodon{$t}; $i<=$numberStopCodon{$t}+1; $i++)
                {
                        if($i == $numberStartCodon{$t})
                        {
                                $anno.= substr($fasta{$chr{$t}}, $leftStartCodon{$t}-1, $endExon{$t}{$i}-$leftStartCodon{$t}+1);
                                $skip.= substr($fasta{$chr{$t}}, $leftStartCodon{$t}-1, $endExon{$t}{$i}-$leftStartCodon{$t}+1);
				$anno_all+= $endExon{$t}{$i}-$leftStartCodon{$t}+1;
                        }elsif($i == $numberStopCodon{$t})
                        {
                                $anno.= substr($fasta{$chr{$t}}, $startExon{$t}{$i}-1, $rightStopCodon{$t}-$startExon{$t}{$i}+1);
				$anno_all+= $rightStopCodon{$t}-$startExon{$t}{$i}+1;
				$anno_last+= $rightStopCodon{$t}-$startExon{$t}{$i}+1;
				$anno_utr.= substr($fasta{$chr{$t}}, $rightStopCodon{$t}, $endExon{$t}{$i}-$rightStopCodon{$t});
                        }elsif($i == $numberStopCodon{$t}+1)
                        {
                                $skip.= substr($fasta{$chr{$t}}, $startExon{$t}{$i}-1, $endExon{$t}{$i}-$startExon{$t}{$i}+3001);
				$anno_utr.= substr($fasta{$chr{$t}}, $startExon{$t}{$i}-1, $endExon{$t}{$i}-$startExon{$t}{$i}+302);
                        }else
                        {
                                $anno.= substr($fasta{$chr{$t}}, $startExon{$t}{$i}-1, $endExon{$t}{$i}-$startExon{$t}{$i}+1);
                                $skip.= substr($fasta{$chr{$t}}, $startExon{$t}{$i}-1, $endExon{$t}{$i}-$startExon{$t}{$i}+1);
				$anno_all+= $endExon{$t}{$i}-$startExon{$t}{$i}+1;
                        }
                }
        }else
        {
                for($i=$numberStopCodon{$t}+1; $i>=$numberStartCodon{$t}; $i--)
                {
                        if($i == $numberStopCodon{$t})
                        {
                                $anno.= substr($fasta{$chr{$t}}, $leftStopCodon{$t}-1, $endExon{$t}{$i}-$leftStopCodon{$t}+1);
				$anno_all+= $endExon{$t}{$i}-$leftStopCodon{$t}+1;
				$anno_last+= $endExon{$t}{$i}-$leftStopCodon{$t}+1;
				$anno_utr.= substr($fasta{$chr{$t}}, $startExon{$t}{$i}, $leftStopCodon{$t}-$startExon{$t}{$i});
                        }elsif($i == $numberStartCodon{$t})
                        {
                                $anno.= substr($fasta{$chr{$t}}, $startExon{$t}{$i}-1, $rightStartCodon{$t}-$startExon{$t}{$i}+1);
                                $skip.= substr($fasta{$chr{$t}}, $startExon{$t}{$i}-1, $rightStartCodon{$t}-$startExon{$t}{$i}+1);
				$anno_all+= $rightStartCodon{$t}-$startExon{$t}{$i}+1;
                        }elsif($i == $numberStopCodon{$t}+1)
                        {
                                $skip.= substr($fasta{$chr{$t}}, $startExon{$t}{$i}-3001, $endExon{$t}{$i}-$startExon{$t}{$i}+3001);
                                $anno_utr.= substr($fasta{$chr{$t}}, $startExon{$t}{$i}-303, $endExon{$t}{$i}-$startExon{$t}{$i}+302);		
                        }else
                        {
                                $anno.= substr($fasta{$chr{$t}}, $startExon{$t}{$i}-1, $endExon{$t}{$i}-$startExon{$t}{$i}+1);
                                $skip.= substr($fasta{$chr{$t}}, $startExon{$t}{$i}-1, $endExon{$t}{$i}-$startExon{$t}{$i}+1);
				$anno_all+= $endExon{$t}{$i}-$startExon{$t}{$i}+1;
                        }
                }
                $anno=~tr/[AGCT]/[TCGA]/;
                $anno=reverse($anno);
                $skip=~tr/[AGCT]/[TCGA]/;
                $skip=reverse($skip);
                $anno_utr=~tr/[AGCT]/[TCGA]/;
                $anno_utr=reverse($anno_utr);
        }
        next if length($anno) % 3; # wrong annotation: coding length not divisible by 3
        $anno_aa=getProtein($anno);
        $skip_aa=getProtein($skip);
        $anno_aa_stop=substr($anno_aa, 0,  index($anno_aa, "_") + 1);
        $skip_aa_stop=substr($skip_aa, 0,  index($skip_aa, "_") + 1);
        next if $anno_aa ne $anno_aa_stop; # wrong annotation: premature stop codon
        print OUT4 "\>$t\_\_$strand{$t}\_\_Penultimate_splice_in\n$anno\n\>$t\_\_$strand{$t}\_\_Penultimate_splice_out\n$skip\n";
        print OUT6 "\>$t\_\_$strand{$t}\_\_Penultimate_splice_in\n$anno_aa\n\>$t\_\_$strand{$t}\_\_Penultimate_splice_out\n$skip_aa_stop\n";
        $anno_aa_nolast_len = int(($anno_all - $anno_last) / 3);
        $skip_aa_uniq = substr($skip_aa_stop, $anno_aa_nolast_len);
        $anno_aa_uniq = substr($anno_aa, $anno_aa_nolast_len);
        $anno_aa_cds0 = (length($anno_aa) >=100) ? substr($anno_aa, length($anno_aa)-100) : $anno_aa; # last 100 AA
        $anno_aa_utr0 = getProtein(substr($anno_utr,0,300));
        $anno_aa_utr1 = getProtein(substr($anno_utr,1,300));
        $anno_aa_utr2 = getProtein(substr($anno_utr,2,300));
        print OUT8 "\>$t\_\_$strand{$t}\_\_Penultimate_splice_in\_\_0\n$anno_aa_uniq\n\>$t\_\_$strand{$t}\_\_Penultimate_splice_out\_\_0\n$skip_aa_uniq\n";
        print OUT8 "\>$t\_\_$strand{$t}\_\_Penultimate_splice_in_tail100\n$anno_aa_cds0\n\>$t\_\_$strand{$t}\_\_Penultimate_splice_out_utr0\n$anno_aa_utr0\n\>$t\_\_$strand{$t}\_\_Penultimate_splice_out_utr1\n$anno_aa_utr1\n\>$t\_\_$strand{$t}\_\_Penultimate_splice_out_utr2\n$anno_aa_utr2\n";
        if(length($skip_aa_uniq)>20)
        {
		$final{$t} = "";
        	print OUT10 "\>$t\_\_$strand{$t}\_\_Penultimate_splice_in\_\_0\n$anno_aa_uniq\n\>$t\_\_$strand{$t}\_\_Penultimate_splice_out\_\_0\n$skip_aa_uniq\n";
        	print OUT10 "\>$t\_\_$strand{$t}\_\_Penultimate_splice_in_tail100\n$anno_aa_cds0\n\>$t\_\_$strand{$t}\_\_Penultimate_splice_out_utr0\n$anno_aa_utr0\n\>$t\_\_$strand{$t}\_\_Penultimate_splice_out_utr1\n$anno_aa_utr1\n\>$t\_\_$strand{$t}\_\_Penultimate_splice_out_utr2\n$anno_aa_utr2\n";
        }
}


open(IN1, "StopCodonInLastExon.$GTFN");
open(IN2, "StopCodonInPenultimateExon.$GTFN");
open(OUT11, ">StopCodonInLastExon.final.$GTFN");
open(OUT12, ">StopCodonInPenultimateExon.final.$GTFN");

while(<IN1>)
{
        chomp;
        @l = split(/\t/,);
        next if $l[2] ne "exon";
        @m = split(/\"/, $l[8]);
        for($i=0; $i<$#m; $i++)
        {
                if($m[$i]=~m/transcript_id/)
                {
                        $transcript_idx = $i + 1;
                }
        }
	if(exists $final{$m[$transcript_idx]})
	{
		print OUT11 "$_\n";
	}
}

while(<IN2>)
{
        chomp;
        @l = split(/\t/,);
        next if $l[2] ne "exon";
        @m = split(/\"/, $l[8]);
        for($i=0; $i<$#m; $i++)
        {
                if($m[$i]=~m/transcript_id/)
                {
                        $transcript_idx = $i + 1;
                }
        }
        if(exists $final{$m[$transcript_idx]})
        {
                print OUT12 "$_\n";
        }
}
