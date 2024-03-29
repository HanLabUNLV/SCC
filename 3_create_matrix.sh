#!/bin/bash

mkdir labelled
rm class_meta.txt

for bedfile in state_bedfiles/scored/*bed; do
	
	id=${bedfile%_scored.bed}
	id=${id##*/}
	echo $id >> class_meta.txt

	awk -v id=$id '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"id}' $bedfile > labelled/${id}.bed


done

cat labelled/* > labelled_all.bed


bedtools sort -i labelled_all.bed | bedtools merge -i - -d -200 -o collapse -c 4,5,6,7,8 > merged.bed

# get original labels
# -F 100% of the region should be in the ChromHMM region
bedtools intersect -a chromHMM_total.bed -b merged.bed -wo -F 1.0 | bedtools sort -i - > final_merged.bed

./3_1_bed_to_matrix.py
