#!/bin/bash

## STEP 0: preprocess labels for multiple cell types

# add labels to different cell types micro-c scores
zcat ../1_get_contact_scores/endoderm_2M.scored.bed.gz_noquies.gz | awk '{ print gensub(/(chr[0-9XY]+)/, "\\1_endo", "g");}' > labelled_endoderm.loops
zcat ../1_get_contact_scores/H1_2M.scored.bed.gz_noquies.gz | awk '{ print gensub(/(chr[0-9XY]+)/, "\\1_h1", "g");}' > labelled_h1.loops
zcat ../1_get_contact_scores/HFF_2M.scored.bed.gz_noquies.gz | awk '{ print gensub(/(chr[0-9XY]+)/, "\\1_hff", "g");}' > labelled_hff.loops
zcat ../1_get_contact_scores/HeLa-S3_2M.scored.bed.gz_noquies.gz | awk '{ print gensub(/(chr[0-9XY]+)/, "\\1_hela", "g");}' > labelled_hela.loops

cat labelled_endoderm.loops labelled_h1.loops labelled_hff.loops labelled_hela.loops | awk '{ print gensub(/ZNF\/Rpts/, "ZNF-Rpts", "g");}' -  > all_loops.bedpe

# add labels to different cell types chromHMM labels
zcat ../chromHMM_epimap_calls/H1.bed.gz | awk '{ print gensub(/ZNF\/Rpts/, "ZNF-Rpts", "g");}' -  > x
awk '{ print gensub(/(chr[0-9XY]+)/, "\\1_h1", "g");}' x > chromHMM_total.bed

zcat ../chromHMM_epimap_calls/endoderm.bed.gz | awk '{ print gensub(/ZNF\/Rpts/, "ZNF-Rpts", "g");}' -  > x
awk '{ print gensub(/(chr[0-9XY]+)/, "\\1_endo", "g");}' x >> chromHMM_total.bed

zcat ../chromHMM_epimap_calls/HFF.bed.gz | awk '{ print gensub(/ZNF\/Rpts/, "ZNF-Rpts", "g");}' -  > x
awk '{ print gensub(/(chr[0-9XY]+)/, "\\1_hff", "g");}' x >> chromHMM_total.bed

zcat ../chromHMM_epimap_calls/HeLa-S3.bed.gz | awk '{ print gensub(/ZNF\/Rpts/, "ZNF-Rpts", "g");}' -  > x
awk '{ print gensub(/(chr[0-9XY]+)/, "\\1_hela", "g");}' x >> chromHMM_total.bed


### Swap all_loops.bedpe order
###	chr_cont,pos1_cont,pos2_cont  chr_base,pos1_base,pos2_base,label_base  contact_score
### base label in $9

####### LOOPS FILE, 7th COLUMN NEEDS TO BE SCORE
awk '{print $4"\t"$5"\t"$6"\t"$1"\t"$2"\t"$3"\t"$7}' all_loops.bedpe > x
mv x all_loops.bedpe
