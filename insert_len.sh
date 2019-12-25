
SAMPLE=$1
BAM=$2
GFF=$3

source activate miso
pe_utils --compute-insert-len $BAM $GFF --output-dir $SAMPLE/

