#!/bin/bash
#module load bedtools
end_longer=$1
start_longer=$2
ref_end_longer=$3
ref_start_longer=$4
full_bed=$5
#first make the distinct plus case
set_diff () {
  bedtools intersect -v -a $1 -b $2 > /tmp/one.bed
  bedtools intersect -v -a $2 -b $1 > /tmp/two.bed
  cat /tmp/one.bed /tmp/two.bed > $3

}
ref='/tmp/all_ref_exons.bed'
grow='/tmp/grow.bed'
ref_ss='/tmp/ref_ss.bed'
set_diff $end_longer $start_longer $grow
set_diff $ref_end_longer $ref_start_longer $ref
RANDOM=896
shuf -n $(cat $grow| wc -l) $ref > $ref_ss
set_diff $grow $ref_ss $full_bed

#
#
# bedtools intersect -v -a $end_longer -b $start_longer > /tmp/el_nol.bed
# bedtools intersect -v -a $start_longer -b $end_longer > /tmp/sl_nol.bed
# cat /tmp/el_nol.bed /tmp/sl_nol.bed > /tmp/all_growing_exons.bed
#
# bedtools intersect -v -a $ref_end_longer -b $ref_start_longer
# bedtools intersect -v -a $ref_start_longer -b $ref_end_longer; } | cat   |\
#      head -n $(wc -l $grow) > $ref
#
# { bedtools intersect -v -a $grow -b $ref &\
#     bedtools intersect -v -a $ref -b $grow; } | cat > $full_bed
