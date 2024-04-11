
# format h1_labelled.bed to have colors:

setwd("~/Documents/UNLV/Year5/new_kmeans_majority_chromhmm_in_window/fig7-higlass/colored_beds//")


# chromHMM <- read_tsv("../chromHMM_epimap/hff.bed",col_names = c("chr","pos1","pos2","label",
                                                                # "score","strand","x","y","color"))
# colors <- chromHMM %>% select(label,color) %>% unique()

# subcategories:
# c(
#   7 = NULL
#   6 low contact = "Low_cont"
#   5 het contact = "Het_strong_cont","Het_cont"
#   4 biv contact = "Biv_cont","ReprPC_cont"
#   3 weak enahncer contact =  "EnhWk_TxWk_cont","EnhWk_cont"
#   2 enh contact = "EnhA_EnhG_cont","EnhA_EnhWk_cont"
#   1 tss contact = "Tss_Enh_strong_cont","Tss_TssFlnkD_cont","Tss_noEnh_cont","Tss_TssFlnkU_cont"
#   0 tx contact = "Tx_cont","TxWk_EnhG_cont",
# )

# ord <- rev(c("cont_Tss_Enh","cont_Tss_EnhG","cont_Tss_EnhWk","cont_Tss_noEnh",
#              "cont_EnhA","cont_EnhA_EnhWk","cont_Tx_EnhG1","cont_Tx_EnhG2",
#              "cont_TxWk_Enh","cont_TxWk_EnhWk","cont_Low_TssBiv","cont_Low_EnhWk",
#              "cont_Biv","cont_ReprPC_Biv","cont_ReprPC",
#              "cont_Het_Rpts","cont_Het_Rpts_strong","cont_Het_Rpts_strongest"))

clusters <- c("cont_Low_TssBiv" = "220,220,220",
              "cont_Low_EnhWk" = "220,220,220",
              
              "cont_Tss_noEnh" = "255,69,0",
              "cont_Tss_Enh" = "255,0,0",
              "cont_Tss_EnhWk" = "255,0,0",
              "cont_Tss_EnhG" = "255,0,0",
              
              "cont_Het_Rpts" = "128,128,128",
              "cont_Het_Rpts_strong" = "128,128,128",
              "cont_Het_Rpts_strongest" = "128,128,128",
              
              
              "cont_EnhA" = "255,195,77",
              "cont_EnhA_EnhWk" = "255,195,77",
              
              "cont_Tx_EnhG1" = "0,128,0",
              "cont_Tx_EnhG2" = "0,128,0",
              "cont_TxWk_Enh" = "255,255,0",
              "cont_TxWk_EnhWk" = "255,255,0",
              
              
              "cont_Biv" = "189,183,107",
              "cont_ReprPC_Biv" = "189,183,107",
              "cont_ReprPC" = "189,183,107")

hff <- read_tsv("../../fig5-all_sankey/create_bed_network_files/hff_labelled.bed",col_names = F)
hff$color <- clusters[hff$X5]
hff$score <- '0'
hff$strand <- '.'
hff$dup1 <- hff$X2
hff$dup2 <- hff$X3

write_tsv(hff %>% select(X1,X2,X3,X5,score,strand,dup1,dup2,color),"colored_hff.bed",col_names = F)

#clodius aggregate bedfile --chromsizes-filename ../example/chromsizes.hg38 hff.bed
#higlass-manage ingest --filetype beddb --assembly hg38  --datatype bedlike hff.bed.beddb


endo <- read_tsv("../../fig5-all_sankey/create_bed_network_files/endo_labelled.bed",col_names = F)
endo$color <- clusters[endo$X5]
endo$score <- '0'
endo$strand <- '.'
endo$dup1 <- endo$X2
endo$dup2 <- endo$X3

write_tsv(endo %>% select(X1,X2,X3,X5,score,strand,dup1,dup2,color),"colored_endo.bed",col_names = F)


h1 <- read_tsv("../../fig5-all_sankey/create_bed_network_files/h1_labelled.bed",col_names = F)
h1$color <- clusters[h1$X5]
h1$score <- '0'
h1$strand <- '.'
h1$dup1 <- h1$X2
h1$dup2 <- h1$X3

write_tsv(h1 %>% select(X1,X2,X3,X5,score,strand,dup1,dup2,color),"colored_h1.bed",col_names = F)


hela <- read_tsv("../../fig5-all_sankey/create_bed_network_files/hela_labelled.bed",col_names = F)
hela$color <- clusters[hela$X5]
hela$score <- '0'
hela$strand <- '.'
hela$dup1 <- hela$X2
hela$dup2 <- hela$X3

write_tsv(hela %>% select(X1,X2,X3,X5,score,strand,dup1,dup2,color),"colored_hela.bed",col_names = F)

