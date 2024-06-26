#!/usr/bin/perl -w
use strict;
die "Usage: $0  \"unprocessed gtf\"  \"genome_location\"  \"output basename\" \"\(optional\) extend N bases\" " if (@ARGV < 3);
# store postion of exons and orders.
# output mRNA sequences together with exon sequences
# put extended sequences in lower case
# new in v4, treat transcript without strandness (mostly single-exonic) as from "+"
my $anno=$ARGV[0];
my $genome=$ARGV[1];
my $fileout=$ARGV[2];
my $Extend=0;
if (scalar(@ARGV) > 3) {$Extend=$ARGV[3];}
my %uniq;   # store all exons according to chr
my %Exons;  # store the number of exons of each transcript
open IN,$anno;
open OUT,">".$fileout.".exon.fa";
open OUT1,">".$fileout.".gene.fa";
#open OUT2,">".$fileout.".tmp";
while(<IN>) {
    chomp;
    my @a=split("\t",$_);
    if ($a[0]=~m/^chromosome/i) {$a[0]=~s/chromosome//i;}
    if ($a[0]=~m/^chr/i) {$a[0]=~s/chr//i;}
    if ((length($a[0]) <= 3) and ($a[2] eq "exon")) {
        my @b=split(/\"/,$a[8]);
        ### Cel  gene_id "MTCE.3"; transcript_id "MTCE.3"; exon_number "1"; gene_name "NU6M_CAEEL"; transcript_name "NU6M_CAEEL";
        ### Hsa  gene_id "ENSG00000223972"; transcript_id "ENST00000456328"; exon_number "1"; gene_biotype "pseudogene"; gene_name "DDX11L1"; transcript_name "DDX11L1-002"; tss_id "TSS26614";
        ### Mmu  gene_id "Xkr4"; gene_name "Xkr4"; p_id "P15240"; transcript_id "NM_001011874.1"; tss_id "TSS13146";
        ### Dme  gene_id "FBgn0031208"; transcript_id "FBtr0300689"; exon_number "1"; gene_name "CG11023"; p_id "P13746"; transcript_name "CG11023-RB"; tss_id "TSS8369";
        my $len=scalar(@b);
        my $gene_id="";
        my $transcript_id="";
        my $gene_name="";
        my $exon_number=0;
        my $strand="+";
        if ($a[6] eq "-") {$strand="-";}
        for (my $i=0; $i<$len; $i=$i+2) {
            if ($b[$i]=~m/gene\_id/) {$gene_id=$b[$i+1];}
            elsif ($b[$i]=~m/transcript\_id/) {$transcript_id=$b[$i+1];}
            elsif ($b[$i]=~m/gene\_name/) {$gene_name=$b[$i+1];}
            elsif ($b[$i]=~m/exon\_number/) {$exon_number=$b[$i+1];}
        }
        #                                           chr  strand   start end                           biotype
        #                   0        1              2     3       4     5     6          7            8
        my $info=join("__",$gene_id,$transcript_id,$a[0],$strand,$a[3],$a[4],$gene_name,$exon_number,$a[1]);
        $uniq{$a[0]}{$gene_id."__".$transcript_id}{$exon_number}=$info;
        $Exons{join("__",$a[0],$gene_id,$transcript_id)}++;
        #print OUT2 $info,"\n";
    }
}
close IN;
my %reported;
my $f=0; 
foreach my $chr (sort keys %uniq) {
    #if (length($chr) > 3) {next}
    if ($f eq 1) {print "\t=======\tyes\n";}
    elsif ($f eq -1) {print "\t=======\tno\n";}
    print "searching chromosome : ".$chr.".fa";
    $f=-1;
    if (exists $reported{$chr}) {next;}
    open (IN1, $genome."/".$chr.".fa") or next;
    <IN1>;
    my $seq;
    $f=1;
    $reported{$chr}=1;
    while(<IN1>) {chomp; $seq=$seq.$_;}
    close IN1;
    foreach my $gene (sort keys %{$uniq{$chr}}) {
        my $mRNA="";
        my $left=9999999999;
        my $right=-1;
        my $biotype="";
        foreach my $exon (sort{$a <=> $b} keys %{$uniq{$chr}{$gene}}) {
            my @a=split("__",$uniq{$chr}{$gene}{$exon});
            if ($left > $a[4]) {$left=$a[4]}
            if ($right < $a[5]) {$right=$a[5]}
            $biotype=$a[8];
            #my $exon_seq=substr($seq,($a[4]-1),($a[5]-$a[4]+1));
            my $exon_seq="";
            if ($a[3] eq "+") {
                if (($exon == 1) and ($Exons{$chr."__".$gene} == 1)) {
                    my $up=lc substr($seq,($a[4]-1-$Extend),$Extend);
                    my $down=lc substr($seq,$a[5],$Extend);
                    my $body=substr($seq,($a[4]-1),($a[5]-$a[4]+1));
                    $exon_seq=$up.$body.$down;
                    #$exon_seq=substr($seq,($a[4]-1-$Extend),($a[5]-$a[4]+1+$Extend+$Extend));
                }
                elsif ($exon == 1) {
                    my $up=lc substr($seq,($a[4]-1-$Extend),$Extend);
                    my $body=substr($seq,($a[4]-1),($a[5]-$a[4]+1));
                    $exon_seq=$up.$body;
                    #$exon_seq=substr($seq,($a[4]-1-$Extend),($a[5]-$a[4]+1+$Extend));
                }
                elsif ($exon == $Exons{$chr."__".$gene}) {
                    my $down=lc substr($seq,$a[5],$Extend);
                    my $body=substr($seq,($a[4]-1),($a[5]-$a[4]+1));
                    $exon_seq=$body.$down;
                    #$exon_seq=substr($seq,($a[4]-1),($a[5]-$a[4]+1+$Extend));
                }
                else { $exon_seq=substr($seq,($a[4]-1),($a[5]-$a[4]+1)); }
                print OUT ">".$uniq{$chr}{$gene}{$exon},"\n",$exon_seq,"\n";
                if ($mRNA eq "") {$mRNA=$exon_seq;}
                else {$mRNA=$mRNA.$exon_seq;}
            }
            else {
                if (($exon == 1) and ($Exons{$chr."__".$gene} == 1)) {
                    my $up=lc substr($seq,($a[4]-1-$Extend),$Extend);
                    my $down=lc substr($seq,$a[5],$Extend);
                    my $body=substr($seq,($a[4]-1),($a[5]-$a[4]+1));
                    $exon_seq=$up.$body.$down;
                    #$exon_seq=substr($seq,($a[4]-1-$Extend),($a[5]-$a[4]+1+$Extend+$Extend));
                }
                elsif ($exon == 1) {
                    my $up=lc substr($seq,$a[5],$Extend);
                    my $body=substr($seq,($a[4]-1),($a[5]-$a[4]+1));
                    $exon_seq=$body.$up;
                    #$exon_seq=substr($seq,($a[4]-1),($a[5]-$a[4]+1+$Extend));
                }
                elsif ($exon == $Exons{$chr."__".$gene}) {
                    my $down=lc substr($seq,($a[4]-1-$Extend),$Extend);
                    my $body=substr($seq,($a[4]-1),($a[5]-$a[4]+1));
                    $exon_seq=$down.$body;
                    #$exon_seq=substr($seq,($a[4]-1-$Extend),($a[5]-$a[4]+1+$Extend));
                }
                else { $exon_seq=substr($seq,($a[4]-1),($a[5]-$a[4]+1)); }
                $exon_seq=~tr/[atcgATCG]/[tagcTAGC]/;
                my $tmp=reverse scalar $exon_seq;
                print OUT ">".$uniq{$chr}{$gene}{$exon},"\n",$tmp,"\n";
                if ($mRNA eq "") {$mRNA=$tmp;}
                else {$mRNA=$mRNA.$tmp;}
            }
        }
        print OUT1 ">".$gene."__".$biotype."__".$chr,"\n",$mRNA,"\n";
    }
}
foreach my $chr (sort keys %uniq) {
    #if (length($chr) > 3) {next}
    if ($f eq 1) {print "\t=======\tyes\n";}
    elsif ($f eq -1) {print "\t=======\tno\n";}
    print "searching chromosome : chr".$chr.".fa";
    $f=-1;
    if (exists $reported{$chr}) {next;}
    open (IN1, $genome."/chr".$chr.".fa") or next;
    <IN1>;
    my $seq;
    $f=1;
    $reported{$chr}=1;
    while(<IN1>) {chomp; $seq=$seq.$_;}
    close IN1;
    foreach my $gene (sort keys %{$uniq{$chr}}) {
        my $mRNA="";
        my $left=9999999999;
        my $right=-1;
        my $biotype="";
        foreach my $exon (sort{$a <=> $b} keys %{$uniq{$chr}{$gene}}) {
            my @a=split("__",$uniq{$chr}{$gene}{$exon});
            if ($left > $a[4]) {$left=$a[4]}
            if ($right < $a[5]) {$right=$a[5]}
            $biotype=$a[8];
            #my $exon_seq=substr($seq,($a[4]-1),($a[5]-$a[4]+1));
            my $exon_seq="";
            if ($a[3] eq "+") {
                if (($exon == 1) and ($Exons{$chr."__".$gene} == 1)) {
                    my $up=lc substr($seq,($a[4]-1-$Extend),$Extend);
                    my $down=lc substr($seq,$a[5],$Extend);
                    my $body=substr($seq,($a[4]-1),($a[5]-$a[4]+1));
                    $exon_seq=$up.$body.$down;
                    #$exon_seq=substr($seq,($a[4]-1-$Extend),($a[5]-$a[4]+1+$Extend+$Extend));
                }
                elsif ($exon == 1) {
                    my $up=lc substr($seq,($a[4]-1-$Extend),$Extend);
                    my $body=substr($seq,($a[4]-1),($a[5]-$a[4]+1));
                    $exon_seq=$up.$body;
                    #$exon_seq=substr($seq,($a[4]-1-$Extend),($a[5]-$a[4]+1+$Extend));
                }
                elsif ($exon == $Exons{$chr."__".$gene}) {
                    my $down=lc substr($seq,$a[5],$Extend);
                    my $body=substr($seq,($a[4]-1),($a[5]-$a[4]+1));
                    $exon_seq=$body.$down;
                    #$exon_seq=substr($seq,($a[4]-1),($a[5]-$a[4]+1+$Extend));
                }
                else { $exon_seq=substr($seq,($a[4]-1),($a[5]-$a[4]+1)); }
                print OUT ">".$uniq{$chr}{$gene}{$exon},"\n",$exon_seq,"\n";
                if ($mRNA eq "") {$mRNA=$exon_seq;}
                else {$mRNA=$mRNA.$exon_seq;}
            }
            else {
                if (($exon == 1) and ($Exons{$chr."__".$gene} == 1)) {
                    my $up=lc substr($seq,($a[4]-1-$Extend),$Extend);
                    my $down=lc substr($seq,$a[5],$Extend);
                    my $body=substr($seq,($a[4]-1),($a[5]-$a[4]+1));
                    $exon_seq=$up.$body.$down;
                    #$exon_seq=substr($seq,($a[4]-1-$Extend),($a[5]-$a[4]+1+$Extend+$Extend));
                }
                elsif ($exon == 1) {
                    my $up=lc substr($seq,($a[4]-1-$Extend),$Extend);
                    my $body=substr($seq,($a[4]-1),($a[5]-$a[4]+1));
                    $exon_seq=$up.$body;
                    #$exon_seq=substr($seq,($a[4]-1),($a[5]-$a[4]+1+$Extend));
                }
                elsif ($exon == $Exons{$chr."__".$gene}) {
                    my $down=lc substr($seq,$a[5],$Extend);
                    my $body=substr($seq,($a[4]-1),($a[5]-$a[4]+1));
                    $exon_seq=$body.$down;
                    #$exon_seq=substr($seq,($a[4]-1-$Extend),($a[5]-$a[4]+1+$Extend));
                }
                else { $exon_seq=substr($seq,($a[4]-1),($a[5]-$a[4]+1)); }
                $exon_seq=~tr/[atcgATCG]/[tagcTAGC]/;
                my $tmp=reverse scalar $exon_seq;
                print OUT ">".$uniq{$chr}{$gene}{$exon},"\n",$tmp,"\n";
                if ($mRNA eq "") {$mRNA=$tmp;}
                else {$mRNA=$mRNA.$tmp;}
            }
        }
        print OUT1 ">".$gene."__".$biotype."__".$chr,"\n",$mRNA,"\n";
    }
}
foreach my $chr (sort keys %uniq) {
    #if (length($chr) > 3) {next}
    if ($f eq 1) {print "\t=======\tyes\n";}
    elsif ($f eq -1) {print "\t=======\tno\n";}
    print "searching chromosome : chromosome".$chr.".fa";
    $f=-1;
    if (exists $reported{$chr}) {next;}
    open (IN1, $genome."/chromosome".$chr.".fa") or next;
    <IN1>;
    my $seq;
    $f=1;
    $reported{$chr}=1;
    while(<IN1>) {chomp; $seq=$seq.$_;}
    close IN1;
    foreach my $gene (sort keys %{$uniq{$chr}}) {
        my $mRNA="";
        my $left=9999999999;
        my $right=-1;
        my $biotype="";
        foreach my $exon (sort{$a <=> $b} keys %{$uniq{$chr}{$gene}}) {
            my @a=split("__",$uniq{$chr}{$gene}{$exon});
            if ($left > $a[4]) {$left=$a[4]}
            if ($right < $a[5]) {$right=$a[5]}
            $biotype=$a[8];
            #my $exon_seq=substr($seq,($a[4]-1),($a[5]-$a[4]+1));
            my $exon_seq="";
            if ($a[3] eq "+") {
                if (($exon == 1) and ($Exons{$chr."__".$gene} == 1)) {
                    my $up=lc substr($seq,($a[4]-1-$Extend),$Extend);
                    my $down=lc substr($seq,$a[5],$Extend);
                    my $body=substr($seq,($a[4]-1),($a[5]-$a[4]+1));
                    $exon_seq=$up.$body.$down;
                    #$exon_seq=substr($seq,($a[4]-1-$Extend),($a[5]-$a[4]+1+$Extend+$Extend));
                }
                elsif ($exon == 1) {
                    my $up=lc substr($seq,($a[4]-1-$Extend),$Extend);
                    my $body=substr($seq,($a[4]-1),($a[5]-$a[4]+1));
                    $exon_seq=$up.$body;
                    #$exon_seq=substr($seq,($a[4]-1-$Extend),($a[5]-$a[4]+1+$Extend));
                }
                elsif ($exon == $Exons{$chr."__".$gene}) {
                    my $down=lc substr($seq,$a[5],$Extend);
                    my $body=substr($seq,($a[4]-1),($a[5]-$a[4]+1));
                    $exon_seq=$body.$down;
                    #$exon_seq=substr($seq,($a[4]-1),($a[5]-$a[4]+1+$Extend));
                }
                else { $exon_seq=substr($seq,($a[4]-1),($a[5]-$a[4]+1)); }
                print OUT ">".$uniq{$chr}{$gene}{$exon},"\n",$exon_seq,"\n";
                if ($mRNA eq "") {$mRNA=$exon_seq;}
                else {$mRNA=$mRNA.$exon_seq;}
            }
            else {
                if (($exon == 1) and ($Exons{$chr."__".$gene} == 1)) {
                    my $up=lc substr($seq,($a[4]-1-$Extend),$Extend);
                    my $down=lc substr($seq,$a[5],$Extend);
                    my $body=substr($seq,($a[4]-1),($a[5]-$a[4]+1));
                    $exon_seq=$up.$body.$down;
                    #$exon_seq=substr($seq,($a[4]-1-$Extend),($a[5]-$a[4]+1+$Extend+$Extend));
                }
                elsif ($exon == 1) {
                    my $up=lc substr($seq,($a[4]-1-$Extend),$Extend);
                    my $body=substr($seq,($a[4]-1),($a[5]-$a[4]+1));
                    $exon_seq=$up.$body;
                    #$exon_seq=substr($seq,($a[4]-1),($a[5]-$a[4]+1+$Extend));
                }
                elsif ($exon == $Exons{$chr."__".$gene}) {
                    my $down=lc substr($seq,$a[5],$Extend);
                    my $body=substr($seq,($a[4]-1),($a[5]-$a[4]+1));
                    $exon_seq=$body.$down;
                    #$exon_seq=substr($seq,($a[4]-1-$Extend),($a[5]-$a[4]+1+$Extend));
                }
                else { $exon_seq=substr($seq,($a[4]-1),($a[5]-$a[4]+1)); }
                $exon_seq=~tr/[atcgATCG]/[tagcTAGC]/;
                my $tmp=reverse scalar $exon_seq;
                print OUT ">".$uniq{$chr}{$gene}{$exon},"\n",$tmp,"\n";
                if ($mRNA eq "") {$mRNA=$tmp;}
                else {$mRNA=$mRNA.$tmp;}
            }
        }
        print OUT1 ">".$gene."__".$biotype."__".$chr,"\n",$mRNA,"\n";
    }
}
if ($f eq 1) {print "\t=======\tyes\n";}
elsif ($f eq -1) {print "\t=======\tno\n";}
close OUT;
