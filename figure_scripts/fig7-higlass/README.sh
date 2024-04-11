


grep 'endoderm' GeneHancer_Tissues_v5.12.txt | cut -f1 | sort | uniq > endo_enhancers.txt
grep -E 'H1-h?ESC' -i GeneHancer_Tissues_v5.12.txt | cut -f1 | sort | uniq > h1_enhancers.txt

# add = to front of .txt files

./find_enhancers.py

bedtools intersect -a h1_enhancers.gff -b ../create_bed_network_files/network_h1.bed -wb > h1_net_enhancers.bed
bedtools intersect -a endo_enhancers.gff -b ../create_bed_network_files/network_endo.bed -wb > endo_net_enhancers.bed

# vim %s/--:--/\t

cut -f1,2,3,4,5,6,10,11,12,13,14,15,16,18,19,20,21,22,23 h1_net_enhancers.bed > trim_h1_net_enhancers.bed
cut -f1,2,3,4,5,6,10,11,12,13,14,15,16,18,19,20,21,22,23 endo_net_enhancers.bed > trim_endo_net_enhancers.bed


cut -f1,2,3,4,5,6,10,11,12,13,14,15,16,18,19,20,21,22,23 upreg_endo_net_enhancers.bed > trim_upreg_endo_net_enhancers.bed



# how to get bedfile of endoderm enhancers:

./make_enh_bed_ghid.py
clodius aggregate bedfile --chromsizes-filename colored_beds/chromsizes.hg38 GeneHancer_v5.12_endoderm.bed

higlass-manage ingest --filetype beddb --assembly hg38  --datatype bedlike GeneHancer_v5.12_endoderm.bed.beddb --name "Endo CCS"
