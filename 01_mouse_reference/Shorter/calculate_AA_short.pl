#!/usr/bin/perl -w

open(OUT, ">short_AA_content.tsv");

@AA = ('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y');

open(IN, "StopCodonInLastExon.shorter.Mus_musculus.GRCm38.98.chr.gtf.uniq.filtered.aa");
while(<IN>)
{
	if($_=~m/ENSMUST00000238488/)
	{
		<IN>;<IN>;<IN>;<IN>;<IN>;<IN>;<IN>;<IN>;<IN>;<IN>;<IN>;
		next;
	}
	$id = substr($_, 1, 18);
	$fl = <IN>;
	chomp $fl;
	$fl=~s/\_//g;
	<IN>;
	$fs = <IN>;
	chomp $fs;
	$fs=~s/\_//g;
	$lfs = length($fs);
	foreach $a (@AA)
	{
		$num = () = $fl =~ /$a/g;
		$perc = $num/length($fl);
		print OUT "$id\t$a\tfl\t$perc\t$lfs\n";
	}
        foreach $a (@AA)
        {
                $num = () = $fs =~ /$a/g;
                $perc = ($lfs > 0) ? $num/length($fs) : 0;
                print OUT "$id\t$a\tfs\t$perc\t$lfs\n";
        }
        <IN>;
        $cds = <IN>;
        chomp $cds;
	$cds=~s/\_//g;
        foreach $a (@AA)
        {
                $num = () = $cds =~ /$a/g;
                $perc = $num/length($cds);
                print OUT "$id\t$a\tcds\t$perc\t$lfs\n";
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
                print OUT "$id\t$a\tutr\t$perc\t$lfs\n";
        }
}

