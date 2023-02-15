#!/bin/bash

####### LOOPS FILE, 7th COLUMN NEEDS TO BE SCORE
###	chr_contact_window,pos1_contact_window,pos2_contact_window  chr_base,pos1_base,pos2_base,contact_score label_base

# remove / in "ZNF/Rpts"
zcat all_loops.bedpe.gz | awk '{ print gensub(/ZNF\/Rpts/, "ZNF-Rpts", "g");}' -  > all_loops.bedpe

loops=$1
#loops=all_loops.bedpe

### Prep chromHMM bedfile for intersect
chromHMM=$2
#chromHMM=chromHMM_epimap_calls.bed.gz

zcat $chromHMM | awk '{ print gensub(/ZNF\/Rpts/, "ZNF-Rpts", "g");}' -  > x
awk '{ print gensub(/(chr[0-9XY]+)/, "\\1_h1", "g");}' x > chromHMM_total.bed

cut -f1-4 chromHMM_total.bed > chromHMM_coordonly.bed
cut -f4 chromHMM_coordonly.bed | sort | uniq > class_meta.txt

echo "prepped loops and ChromHMM"


### STEP 1: Intersect contacts (micro-c) with chromHMM labels
# don't have to do other side because it's accounted for in step 1 get_contact_scores

echo "intersecting loops and chromhmm"
bedtools intersect -a $loops -b chromHMM_coordonly.bed -wo > intersect_loops_all.bed

### STEP 2: merge duplicate elements overlapping more than one Micro-C 1kb window
# average Micro-C score across all windows that chromHMM label is found in
# chr_base,pos1_base,pos2_base,score,label_contact
./get_mean_score_overlap_windows.py intersect_loops_all.bed contact_states.bed


#### STEP 3: SEPARATE AND SCORE BY STATE
mkdir state_bedfiles
mkdir state_bedfiles/scored

echo "summing states"
while read -r state; do

	grep -P "\t$state$" contact_states.bed > state_bedfiles/${state}.bed

	bedtools sort -i state_bedfiles/${state}.bed > z
	bedtools merge -o sum,collapse,collapse,collapse -c 4,5,6,7 -i z -d -200 > state_bedfiles/scored/${state}_scored.bed


done < class_meta.txt

