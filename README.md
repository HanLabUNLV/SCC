# SOCCS
Sum Of Chromatin State by Contact

## Step 1

Get contact scores for every region with a chromHMM label.

```
# python 1_get_scores_mcool.py <mcool file with res> <ChromHMM bedfile> <outfile> [margin size]

python 1_get_scores_mcool.py 4DNFI9E222YJ.mcool::/resolutions/1000 \
							 endoderm_ChromHMM.bed.gz \
							 all_loops.bedpe.gz
```

## Step 2

Get labels of every contact region found in step 1. In this 
analysis we also merge cell types in this step, see `preprocess_scripts`.

```
# chromHMM_total.bed
#  chr1    10000    177200    Quies    0    .    10000    177200    255,255,255
#  chr1    257849    297849    Quies    0    .    257849    297849    255,255,255
#  chr1    586020    777820    Quies    0    .    586020    777820    255,255,255

# all_loops.bedpe.gz
#  chr1    778000    779000    chr1_endo    777820    778420    0.17871180991657024
#  chr1    779000    780000    chr1_endo    777820    778420    0.23094778295675608
#  chr1    787000    788000    chr1_endo    777820    778420    0.05381787728341296

./2_get_contact_labels.sh all_loops.bedpe.gz chromHMM_total.bed
```

## Step 3

Create a summed matrix of contacts from previous step

```
# include path to scored bedfiles within script (created in step 2)
# 	state_bedfiles/scored/*bed

./3_create_matrix.sh
```

The resulting `matrix_connections.tsv` is ready for clustering.





# Supplemental Data for the manuscript

[Supplemental_data_S1](https://drive.google.com/file/d/15KZXQA44aqJyhd_Yzg-EkTM-0C3UHrTO/view?usp=share_link): SCC matrix resulting from the sum of chromatin state by contact 

[Supplemental_data_S2](https://drive.google.com/file/d/1Khxh88MKg6TZDbATef04-Df4MKYfA5ke/view?usp=share_link): cluster membership for each row resulting from k-means clustering on SCC matrix

[Supplemental_table_S1](https://docs.google.com/spreadsheets/d/1Fe3m6UaxT8Dmn32-xCZBX0l2AJsxaUT3/edit?usp=share_link&ouid=116370094857147559521&rtpof=true&sd=true): differential expression analysis 


