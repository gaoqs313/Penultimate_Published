#!/usr/bin/perl -w

open(OUT, ">AA_content.tsv");

@AA = ('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y');

open(IN, "StopCodonInLastExon.Homo_sapiens.GRCh38.98.chr.gtf.uniq.filtered.aa");
while(<IN>)
{
	$id = substr($_, 1, 15);
	$fl = <IN>;
	chomp $fl;
	$fl=~s/\_//g;
	foreach $a (@AA)
	{
		$num = () = $fl =~ /$a/g;
		$perc = $num/length($fl);
		print OUT "$id\t$a\tfl\t$perc\n";
	}
	<IN>;
	$fs = <IN>;
	chomp $fs;
	$fs=~s/\_//g;
        foreach $a (@AA)
        {
                $num = () = $fs =~ /$a/g;
                $perc = $num/length($fs);
                print OUT "$id\t$a\tfs\t$perc\n";
        }
        <IN>;
        $cds = <IN>;
        chomp $cds;
	$cds=~s/\_//g;
        foreach $a (@AA)
        {
                $num = () = $cds =~ /$a/g;
                $perc = $num/length($cds);
                print OUT "$id\t$a\tcds\t$perc\n";
        }
        <IN>;
        $utr1 = <IN>;
        chomp $utr1;
	$utr1=~s/\_//g;
        <IN>;
        $utr2 = <IN>;
        chomp $utr2;
        $utr2=~s/\_//g;
        <IN>;
        $utr3 = <IN>;
        chomp $utr3;
        $utr3=~s/\_//g;
	$utr=$utr1.$utr2.$utr3;
        foreach $a (@AA)
        {
                $num = () = $utr =~ /$a/g;
                $perc = $num/length($utr);
                print OUT "$id\t$a\tutr\t$perc\n";
        }
}

open(IN, "StopCodonInPenultimateExon.Homo_sapiens.GRCh38.98.chr.gtf.uniq.filtered.aa");
while(<IN>)
{
        $id = substr($_, 1, 15);
        $fl = <IN>;
        chomp $fl;
        $fl=~s/\_//g;
	if(length($fl) == 0)
	{
		<IN>;<IN>;<IN>;<IN>;<IN>;<IN>;<IN>;<IN>;<IN>;<IN>;next;
	}
        foreach $a (@AA)
        {
                $num = () = $fl =~ /$a/g;
                $perc = $num/length($fl);
                print OUT "$id\t$a\tfl\t$perc\n";
        }
        <IN>;
        $fs = <IN>;
        chomp $fs;
        $fs=~s/\_//g;
        foreach $a (@AA)
        {
                $num = () = $fs =~ /$a/g;
                $perc = $num/length($fs);
                print OUT "$id\t$a\tfs\t$perc\n";
        }
        <IN>;
        $cds = <IN>;
        chomp $cds;
        $cds=~s/\_//g;
        foreach $a (@AA)
        {
                $num = () = $cds =~ /$a/g;
                $perc = $num/length($cds);
                print OUT "$id\t$a\tcds\t$perc\n";
        }
        <IN>;
        $utr1 = <IN>;
        chomp $utr1;
        $utr1=~s/\_//g;
        <IN>;
        $utr2 = <IN>;
        chomp $utr2;
        $utr2=~s/\_//g;
        <IN>;
        $utr3 = <IN>;
        chomp $utr3;
        $utr3=~s/\_//g;
        $utr=$utr1.$utr2.$utr3;
        foreach $a (@AA)
        {
                $num = () = $utr =~ /$a/g;
                $perc = $num/length($utr);
                print OUT "$id\t$a\tutr\t$perc\n";
        }
}
