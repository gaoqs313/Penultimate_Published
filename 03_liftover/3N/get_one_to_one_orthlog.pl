#!/usr/bin/perl -w

open(IN, "mm10.SE.intron.lift2human.bed");
while(<IN>)
{
	chomp;
	@l = split(/\t/,);
	@m = split(/\_\_/, $l[3]);
	$id = join("\t", @m[0..1]);
	$m_m{$id}{$m[2]} = join("\t", @m[3..5]);
	$m_h{$id}{$m[2]} = join("\t", @l[0..2]);
	$mmiso{$id} = $l[4];
}

open(IN, "hg38.SE.intron.lift2mouse.bed");
while(<IN>)
{
        chomp;
        @l = split(/\t/,);
        @m = split(/\_\_/, $l[3]);
        $id = join("\t", @m[0..1]);
        $h_h{$id}{$m[2]} = join("\t", @m[3..5]);
	$h_m{$id}{$m[2]} = join("\t", @l[0..2]);
	$hmiso{$id} = $l[4];
}

foreach $i (keys %m_m)
{
	$key = keys %{$m_m{$i}};
	next if $key == 1;
	@up = split(/\t/, $m_m{$i}{'up'});
	@dn = split(/\t/, $m_m{$i}{'dn'});
	$mpos = join("\t", $up[0], $up[1], $up[2], $dn[1], $dn[2]);
	@upl = split(/\t/, $m_h{$i}{'up'});
	@dnl = split(/\t/, $m_h{$i}{'dn'});
	$hpos = join("\t", $upl[0], $upl[1], $upl[2], $dnl[1], $dnl[2]);
	$mouse2human{$mpos} = $hpos;
	$mousegene{$mpos} = $i;
	$mmiso{$mpos} = $mmiso{$i};
}

foreach $i (keys %h_h)
{
        $key = keys %{$h_h{$i}};
        next if $key == 1;
        @up = split(/\t/, $h_h{$i}{'up'});
        @dn = split(/\t/, $h_h{$i}{'dn'});
        $hpos = join("\t", $up[0], $up[1], $up[2], $dn[1], $dn[2]);
        @upl = split(/\t/, $h_m{$i}{'up'});
        @dnl = split(/\t/, $h_m{$i}{'dn'});
        $mpos = join("\t", $upl[0], $upl[1], $upl[2], $dnl[1], $dnl[2]);
	$human2mouse{$hpos} = $mpos;
	$humangene{$hpos} = $i;
	$hmiso{$hpos} = $hmiso{$i};
}

open(MM, "/scratch_space/qgao/Penultimate/11.revision/B_Mouse_Penultimate_vs_Internal/mouse_3N_support_info.txt");
<MM>;
while(<MM>)
{
	chomp;
	@l = split(/\t/,);
	$mrna{$l[0]} = $l[4];
}

open(HH, "/scratch_space/qgao/Penultimate/11.revision/C_Human_Penultimate/human_3N_support_info.txt");
<HH>;
while(<HH>)
{
        chomp;
        @l = split(/\t/,);
        $hrna{$l[0]} = $l[4];
}

open(OUT, ">liftOver.results");
foreach $k (keys %mouse2human)
{
	if(exists $human2mouse{$mouse2human{$k}})
	{
		next if ($k ne $human2mouse{$mouse2human{$k}});
		print OUT "$mmiso{$k}\t$mousegene{$k}\t$mrna{$mmiso{$k}}\t$hmiso{$mouse2human{$k}}\t$humangene{$mouse2human{$k}}\t$hrna{$hmiso{$mouse2human{$k}}}\n"; 
	}
}


