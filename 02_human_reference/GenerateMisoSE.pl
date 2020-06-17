#!/usr/bin/perl -w

open(IN, "/scratch_space/qgao/Penultimate/1.reference/human/GTF/Homo_sapiens.GRCh38.98.chr.gtf");
while(<IN>)
{
        chomp;
        next if $_=~m/^#/;
        @l = split(/\t/,);
        next if $l[2] ne "exon";
        @m = split(/\"/, $l[8]);
        for($i=0; $i<$#m; $i++)
        {
                if($m[$i]=~m/transcript_id/)
                {
                        $transcript_idx = $i + 1;
                }elsif($m[$i]=~m/exon_number/)
                {
                        $exon_idx = $i + 1;
                }
        }
        $chrs{$m[$transcript_idx]} = $l[0];
	$strands{$m[$transcript_idx]} = $l[6];
        $startExon{$m[$transcript_idx]}{$m[$exon_idx]} = $l[3];
        $endExon{$m[$transcript_idx]}{$m[$exon_idx]} = $l[4];
}

foreach $t (keys %startExon)
{
        $n = keys %{$startExon{$t}};
        next if $n < 2;
        for($e=1; $e<$n; $e++)
        {
		$f = $e + 1;
		if($strands{$t} eq "+")
		{
			$junc=join("\t", $chrs{$t}, $endExon{$t}{$e}, $startExon{$t}{$f});
		}else
		{
			$junc=join("\t", $chrs{$t}, $startExon{$t}{$e}, $endExon{$t}{$f});
		}
		$annotate{$junc} = "";
	}
}

open(OUT,">hg38.SE.gff3");
open(OUT2, ">hg38.SE.intron.bed");
open(INT, ">hg38.SE.internal.gff3");
open(OUT3, ">penultimate.novelty.txt");
open(OUT4, ">internal.novelty.txt");

open(IN, "/scratch_space/qgao/Penultimate/1.reference/human/A_Annotation/StopCodonInLastExon.final.Homo_sapiens.GRCh38.98.chr.gtf");
while(<IN>)
{
	chomp;
	@l=split(/\t/,);
	next if $l[2] ne "exon";
        @m = split(/\"/, $l[8]);
        for($i=0; $i<$#m; $i++)
        {
                if($m[$i]=~m/transcript_id/)
                {
                        $transcript_idx = $i + 1;
                }elsif($m[$i]=~m/exon_number/)
                {
                        $exon_idx = $i + 1;
                }elsif($m[$i]=~m/gene_name/)
		{
			$gene_idx = $i + 1;
		}
        }
	$strand{$m[$transcript_idx]}=$l[6];
	$chr{$m[$transcript_idx]}=$l[0];
	$gene{$m[$transcript_idx]}=$m[$gene_idx];
	$start{$m[$transcript_idx]}{$m[$exon_idx]}=$l[3];
	$end{$m[$transcript_idx]}{$m[$exon_idx]}=$l[4];
	$type{$m[$transcript_idx]}="Last";
}

open(IN, "/scratch_space/qgao/Penultimate/1.reference/human/A_Annotation/StopCodonInPenultimateExon.final.Homo_sapiens.GRCh38.98.chr.gtf");
while(<IN>)
{
        chomp;
        @l=split(/\t/,);
        next if $l[2] ne "exon";
        @m = split(/\"/, $l[8]);
        for($i=0; $i<$#m; $i++)
        {
                if($m[$i]=~m/transcript_id/)
                {
                        $transcript_idx = $i + 1;
                }elsif($m[$i]=~m/exon_number/)
                {
                        $exon_idx = $i + 1;
                }elsif($m[$i]=~m/gene_name/)
                {
                        $gene_idx = $i + 1;
                }
        }
        $strand{$m[$transcript_idx]}=$l[6];
        $chr{$m[$transcript_idx]}=$l[0];
        $gene{$m[$transcript_idx]}=$m[$gene_idx];
        $start{$m[$transcript_idx]}{$m[$exon_idx]}=$l[3];
        $end{$m[$transcript_idx]}{$m[$exon_idx]}=$l[4];
        $type{$m[$transcript_idx]}="Penultimate";
}

foreach $t (keys %start)
{
	$k = keys %{$start{$t}};
	next if $k < 3;
	for($e=2; $e<$k; $e++)
	{
		$id=$chr{$t}.":".$start{$t}{$e-1}.":".$end{$t}{$e-1}.":".$strand{$t}."@".$chr{$t}.":".$start{$t}{$e}.":".$end{$t}{$e}.":".$strand{$t}."@".$chr{$t}.":".$start{$t}{$e+1}.":".$end{$t}{$e+1}.":".$strand{$t};
		$a= ($start{$t}{$e-1} > $start{$t}{$e+1}) ? $start{$t}{$e+1} : $start{$t}{$e-1};
		$b= ($end{$t}{$e+1} > $end{$t}{$e-1}) ? $end{$t}{$e+1} : $end{$t}{$e-1};
		if($strand{$t} eq "+")
		{
			$juncid= join("\t", $chr{$t}, $end{$t}{$e-1}, $start{$t}{$e+1});
		}else
		{
			$juncid= join("\t", $chr{$t}, $start{$t}{$e-1}, $end{$t}{$e+1});
		}
		$novel = (exists $annotate{$juncid}) ? "Annotated" : "Novel";
		if($e == $k -1)
		{
			next if exists $mark{$id};
			$mark{$id} = "";
			print OUT "$chr{$t}\tSE\tgene\t$a\t$b\t\.\t$strand{$t}\t\.\tName\=$id\;gid=$id\;ID=$id\n";
         		print OUT "$chr{$t}\tSE\tmRNA\t$a\t$b\t\.\t$strand{$t}\t\.\tgid=$id\;ID=$id\.A\;Parent=$id\n";
                	print OUT "$chr{$t}\tSE\texon\t".$start{$t}{$e-1}."\t".$end{$t}{$e-1}."\t"."\.\t$strand{$t}\t\.\tgid=$id\;ID=$id\.A.up\;Parent=$id\.A\n";
                	print OUT "$chr{$t}\tSE\texon\t".$start{$t}{$e}."\t".$end{$t}{$e}."\t"."\.\t$strand{$t}\t\.\tgid=$id\;ID=$id\.A.se\;Parent=$id\.A\n";
                	print OUT "$chr{$t}\tSE\texon\t".$start{$t}{$e+1}."\t".$end{$t}{$e+1}."\t"."\.\t$strand{$t}\t\.\tgid=$id\;ID=$id\.A.dn\;Parent=$id\.A\n";
                	print OUT "$chr{$t}\tSE\tmRNA\t$a\t$b\t\.\t$strand{$t}\t\.\tgid=$id\;ID=$id\.B\;Parent=$id\n";
                	print OUT "$chr{$t}\tSE\texon\t".$start{$t}{$e-1}."\t".$end{$t}{$e-1}."\t"."\.\t$strand{$t}\t\.\tgid=$id\;ID=$id\.B.up\;Parent=$id\.B\n";
                	print OUT "$chr{$t}\tSE\texon\t".$start{$t}{$e+1}."\t".$end{$t}{$e+1}."\t"."\.\t$strand{$t}\t\.\tgid=$id\;ID=$id\.B.dn\;Parent=$id\.B\n";	
                	if($strand{$t} eq "+")
                	{
                	        print OUT2 "$chr{$t}\t$end{$t}{$e-1}\t".($start{$t}{$e}-1)."\t$t\_\_$gene{$t}\_\_up\t$id\t$strand{$t}\n";
                	        print OUT2 "$chr{$t}\t$end{$t}{$e}\t".($start{$t}{$e+1}-1)."\t$t\_\_$gene{$t}\_\_dn\t$id\t$strand{$t}\n";
                	}else
                	{
                	        print OUT2 "$chr{$t}\t$end{$t}{$e}\t".($start{$t}{$e-1}-1)."\t$t\_\_$gene{$t}\_\_up\t$id\t$strand{$t}\n";
                	        print OUT2 "$chr{$t}\t$end{$t}{$e+1}\t".($start{$t}{$e}-1)."\t$t\_\_$gene{$t}\_\_dn\t$id\t$strand{$t}\n";
                	}
			print OUT3 "$id\t$novel\t$type{$t}\t$t\n";
		}else
		{
                        next if (exists $mark{$id} || exists $mark2{$id}); 
                        $mark2{$id} = "";
                        print INT "$chr{$t}\tSE\tgene\t$a\t$b\t\.\t$strand{$t}\t\.\tName\=$id\;gid=$id\;ID=$id\n";
                        print INT "$chr{$t}\tSE\tmRNA\t$a\t$b\t\.\t$strand{$t}\t\.\tgid=$id\;ID=$id\.A\;Parent=$id\n";
                        print INT "$chr{$t}\tSE\texon\t".$start{$t}{$e-1}."\t".$end{$t}{$e-1}."\t"."\.\t$strand{$t}\t\.\tgid=$id\;ID=$id\.A.up\;Parent=$id\.A\n";
                        print INT "$chr{$t}\tSE\texon\t".$start{$t}{$e}."\t".$end{$t}{$e}."\t"."\.\t$strand{$t}\t\.\tgid=$id\;ID=$id\.A.se\;Parent=$id\.A\n";
                        print INT "$chr{$t}\tSE\texon\t".$start{$t}{$e+1}."\t".$end{$t}{$e+1}."\t"."\.\t$strand{$t}\t\.\tgid=$id\;ID=$id\.A.dn\;Parent=$id\.A\n";
                        print INT "$chr{$t}\tSE\tmRNA\t$a\t$b\t\.\t$strand{$t}\t\.\tgid=$id\;ID=$id\.B\;Parent=$id\n";
                        print INT "$chr{$t}\tSE\texon\t".$start{$t}{$e-1}."\t".$end{$t}{$e-1}."\t"."\.\t$strand{$t}\t\.\tgid=$id\;ID=$id\.B.up\;Parent=$id\.B\n";
                        print INT "$chr{$t}\tSE\texon\t".$start{$t}{$e+1}."\t".$end{$t}{$e+1}."\t"."\.\t$strand{$t}\t\.\tgid=$id\;ID=$id\.B.dn\;Parent=$id\.B\n";
			print OUT4 "$id\t$novel\t$t\n";
		}
	}
}
