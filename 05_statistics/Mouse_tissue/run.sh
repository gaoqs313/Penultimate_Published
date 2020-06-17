
## combine penultimate longer exon
perl combine_mouse_tissue_penultimate.pl
#3233
Rscript mouse_tissue_penultimate_stat.R

## combine internal frameshift exons
perl combine_internal_fs.pl
#87093
Rscript internal_fs.R

## combine internal preserving exons
perl combine_internal_fp.pl
#60007
Rscript internal_fp.R

## combine penultimate 3N exons
perl combine_pen_3n.pl
# 8512
Rscript pen_3n.R

## combine penultimate shorter exons
perl combine_pen_short.pl
# 9791
Rscript pen_short.R

## combine plot
Rscript combine_plot.R

