
# Fig 5 heatmap of TSS labels by gene
setwd("~/Documents/UNLV/Year5/new_kmeans_majority_chromhmm_in_window/fig6-gene_heatmap/")

library(readr)
library(ComplexHeatmap)
library(dplyr)
library(tidyr)

# subcategories:
# c(
#   7 = NULL
#   6 low contact = "cont_Low_TssBiv", "cont_Low_EnhWk"
#   5 het contact = "cont_Het_Rpts", "cont_Het_Rpts_strong", "cont_Het_Rpts_strongest"
#   4 ReprPC contact =  "cont_ReprPC_Biv","cont_ReprPC"
#   3 biv contact = "cont_Biv"
#   2 enh contact = "cont_EnhA","cont_EnhA_EnhWk"
#   1 tss contact = "cont_Tss_Enh","cont_Tss_EnhG","cont_Tss_noEnh","cont_Tss_EnhWk"
#   0 tx contact = "cont_Tx_EnhG1","cont_Tx_EnhG2","cont_TxWk_Enh","cont_TxWk_EnhWk"
# )

ccs_colors = c("cont_Tss_Enh" ="#FFFF00",
           "cont_Tss_EnhG"="#FFFF00",
           "cont_Tss_noEnh"="#FFFF00",
           "cont_Tss_EnhWk"="#FFFF00",
           
           "cont_EnhA"="#FF0000",
           "cont_EnhA_EnhWk"="#FF0000",
           
           "cont_ReprPC_Biv"= "#964B00",
           "cont_ReprPC" = "#964B00",
           
           "cont_Tx_EnhG1" = "#00CD00",
           "cont_Tx_EnhG2" = "#00CD00",
           "cont_TxWk_EnhWk" = "#00CD00",
           "cont_TxWk_Enh" = "#00CD00",
           
           "cont_Biv" = "#FFC1C1",

           "cont_Het_Rpts"= "#A020F0",
           "cont_Het_Rpts_strong" = "#A020F0",
           "cont_Het_Rpts_strongest" = "#A020F0",
           
           "cont_Low_TssBiv" = "#75816B",
           "cont_Low_EnhWk" = "#d3d3d3",
           
           "NA" = 'white'
           )

chrom_colors = c("TssA" ="#FFFF00",
               "TssFlnkU"="#FFFF00",
               "TssFlnkD"="#FFFF00",
               "TssFlnk"="#FFFF00",
               
               "EnhA1"="#FF0000",
               "EnhA2"="#FF0000",
               "EnhG1"= "#FF0000",
               "EnhG2" = "#FF0000",
               "EnhWk" = "#FF0000",
               
               
               "Tx" = "#00CD00",
               "TxWk" = "#00CD00",
               
               "EnhBiv" = "#FFC1C1",
               "TssBiv" = "#FFFF9F",
               "ReprPCWk"= "#964B00",
               "ReprPC"= "#964B00",
               
               "Het"= "#A020F0",
               "ZNF-Rpts" = "#A020F0",
               
               "Low_cont" = "#000000",
               "NA" = 'white'
)

# # Gene upreg
# 
# h1_genes <- read_tsv("../fig8-enrichment/expressed_genes/h1_pone.0192625.s007.tsv",col_types = 'ccccccccccccccccccccc') %>% 
#         select(HUES1_TPM,ensembl_id,hgnc_symbol) %>% subset(!is.na(hgnc_symbol))
# colnames(h1_genes) <- c('expr','ensembl_id','hgnc_symbol')
# h1_genes$expr <- as.numeric(h1_genes$expr)
# gene_order_h1 <- h1_genes$hgnc_symbol[order(h1_genes$expr,decreasing = T)]
# 
# 
# endo_genes <- read_tsv("../fig8-enrichment/expressed_genes/endo_Table_S8C-Table_1.tsv",col_types = 'cccccccc') %>%
#         select(AveExpr,Gene,`Gene Name`)
# colnames(endo_genes) <- c('expr','ensembl_id','hgnc_symbol')
# endo_genes <- endo_genes %>% subset(!is.na(hgnc_symbol))
# endo_genes$expr <- as.numeric(endo_genes$expr)
# gene_order_endo <- endo_genes$hgnc_symbol[order(endo_genes$expr,decreasing = T)]
# 
# 
# hff_genes <- read_tsv("../fig8-enrichment/expressed_genes/hff_GSM2448852_0h_3-aligned.tsv",col_types = 'cccccc') %>% 
#         select(`0h_3_FPKM`,gene_id,gene_name)
# colnames(hff_genes) <- c('expr','ensembl_id','hgnc_symbol')
# hff_genes <- hff_genes %>% subset(!is.na(hgnc_symbol))
# hff_genes$expr <- as.numeric(hff_genes$expr)
# gene_order_hff <- hff_genes$hgnc_symbol[order(hff_genes$expr,decreasing = T)]
# 
# 
# gene_exp_endo_h1 <- read_tsv("gene_exp/dec_h1_de.tsv")
# gene_exp_endo_hff <- read_tsv("gene_exp/dec_hff_de.tsv")
# gene_exp_hff_h1 <- read_tsv("gene_exp/hff_h1_de.tsv") %>% subset(complete.cases(.))


gene_exp_hff_h1 <- read_tsv("gene_exp/hff_h1_de.tsv")

sig_hff_h1 <- gene_exp_hff_h1[gene_exp_hff_h1$log2FoldChange > 2 & gene_exp_hff_h1$padj < 1e-10,]

hff_genes <- sig_hff_h1$gene

# all_rows <- apply(ss_h1_ccs,1,table)


# "TX","TSS","Enh","WkEnh","Biv","Het","Low"]) 
# .range(["rgba(93,93,93,.25)",
# "#00CD00","#FFFF00","#FF0000","#FFC1C1","#A020F0","#000000"


# ALL CCS


# 
# h1_ccs <- read_tsv("h1_tss_overlap.ccs_500bp.bed",col_names = F,na = character())
# hela_ccs <- read_tsv("hela_tss_overlap.ccs_500bp.bed",col_names = F)
# hff_ccs <- read_tsv("hff_tss_overlap.ccs_500bp.bed",col_names = F)
# endo_ccs <- read_tsv("endo_tss_overlap.ccs_500bp.bed",col_names = F)
# 
# hff_chrom <- read_tsv("hff_tss_overlap.chrom_500bp.bed",col_names = F)
# h1_chrom <- read_tsv("h1_tss_overlap.chrom_500bp.bed",col_names = F)
# hela_chrom <- read_tsv("hela_tss_overlap.chrom_500bp.bed",col_names = F)
# endo_chrom <- read_tsv("endo_tss_overlap.chrom_500bp.bed",col_names = F)


h1_ccs_mat <- as.matrix(h1_ccs[,-1])
rownames(h1_ccs_mat) <- h1_ccs$X1
h1_chrom_mat <- as.matrix(h1_chrom[,-1])
rownames(h1_chrom_mat) <- h1_chrom$X1


hff_ccs_mat <- as.matrix(hff_ccs[,-1])
rownames(hff_ccs_mat) <- hff_ccs$X1
hff_chrom_mat <- as.matrix(hff_chrom[,-1])
rownames(hff_chrom_mat) <- hff_chrom$X1


hela_ccs_mat <- as.matrix(hela_ccs[,-1])
rownames(hela_ccs_mat) <- hela_ccs$X1
hela_chrom_mat <- as.matrix(hela_chrom[,-1])
rownames(hela_chrom_mat) <- hela_chrom$X1

endo_ccs_mat <- as.matrix(endo_ccs[,-1])
rownames(endo_ccs_mat) <- endo_ccs$X1
endo_chrom_mat <- as.matrix(endo_chrom[,-1])
rownames(endo_chrom_mat) <- endo_chrom$X1


get_percents <- function(mat, lev_list,label) {
        
        all_list <- apply(mat,1,FUN = function(x) {data.frame(table(x) / 500)})
        
        
        V = tibble(name = rownames(mat),all_list)
        all_percentages <- V %>% unnest_longer(col = all_list)
        all_percentages <-  data.frame(all_percentages)
        colnames(all_percentages) <- c("gene",label)
        
        all_percentages[[label]]$x <- factor(all_percentages[[label]]$x)
        levels(all_percentages[[label]]$x) <- lev_list
        
        all_percentages
        
}


chromhmm_levs <- list("TSS" = c("TssA","TssFlnk","TssFlnkD","TssFlnkU"),
                      "Enhancer" = c("EnhA1","EnhA2","EnhG1","EnhG2"),
                      "Weak Enhancer" = c("EnhWk"),
                      "Tx" = c("Tx","TxWk"),
                      "Bivalent" = c("TssBiv","EnhBiv"),
                      "Repressed" = c("ReprPC","ReprPCWk","ZNF-Rpts","Het"),
                      "No Label" = c("NA"))
cis_levs <- list("Tss Contact"=c("cont_Tss_Enh","cont_Tss_noEnh","cont_Tss_EnhG",
                                 "cont_Tss_EnhWk"),
                 "Enh Contact"=c("cont_EnhA","cont_EnhA_EnhWk"),
                 "Tx Contact"=c("cont_Tx_EnhG2","cont_Tx_EnhG1","cont_TxWk_Enh","cont_TxWk_EnhWk"),
                 "Bivalent Contact" = c("cont_Biv"),
                 "ReprPC Contact" = c("cont_ReprPC_Biv","cont_ReprPC"),
                 "Het Contact" = c("cont_Het_Rpts","cont_Het_Rpts_strong","cont_Het_Rpts_strongest"),
                 "Low Contact"= c("cont_Low_TssBiv", "cont_Low_EnhWk"),
                 "No Label" = "NA")

h1_cis_perc <- get_percents(h1_ccs_mat,cis_levs,"H1_CIS")
h1_chr_perc <- get_percents(h1_chrom_mat,chromhmm_levs,"H1_Chromhmm")

hff_cis_perc <- get_percents(hff_ccs_mat,cis_levs,"HFF_CIS")
hff_chr_perc <- get_percents(hff_chrom_mat,chromhmm_levs,"HFF_Chromhmm")

df_list <- list(h1_cis_perc,h1_chr_perc,hff_cis_perc,hff_chr_perc)
all_merged <- df_list %>% purrr::reduce(full_join, by='gene',)

all_percentages_ord <- all_merged[order(all_merged$HFF_CIS$x,all_merged$HFF_CIS$Freq,
                                        all_merged$HFF_Chromhmm$x,all_merged$HFF_Chromhmm$Freq,
                                        all_merged$H1_CIS$x,all_merged$H1_CIS$Freq,
                                        all_merged$H1_Chromhmm$x,all_merged$H1_Chromhmm$Freq,
                                        
                                        decreasing = c(F,T,F,T,
                                                       F,T,F,T)),]

saveRDS(all_percentages_ord,"all_percentages_ord_all_labs.rds")
gene_order <- unique(all_percentages_ord$gene)

hff_genes <-gene_order[gene_order %in% hff_genes]
hff_genes <-gene_order


# H1_ccs



missing_hff_genes <- hff_genes[!(hff_genes %in% rownames(h1_ccs_mat))]
na_mat <- matrix('NA',nrow = length(missing_hff_genes),ncol = 500)
rownames(na_mat) <- missing_hff_genes

h1_ccs_mat <- rbind2(h1_ccs_mat,na_mat)
ss_h1_ccs <- h1_ccs_mat[hff_genes,]

# png("h1_ccs_tss_heatmap.png",height = 1500,width = 900,res = 200)
ccs_h1_hm <- Heatmap(ss_h1_ccs,use_raster = F,cluster_columns = F,cluster_rows = F,col = ccs_colors,
        show_heatmap_legend = T,show_row_names = F,show_column_names = F,row_title = "Genes",
        row_order = hff_genes, name = "CIS label",
        column_title = "H1 CIS")
# dev.off()




missing_hff_genes <- hff_genes[!(hff_genes %in% rownames(h1_chrom_mat))]
na_mat <- matrix('NA',nrow = length(missing_hff_genes),ncol = 500)
rownames(na_mat) <- missing_hff_genes

ss_h1_chrom <- rbind2(h1_chrom_mat,na_mat)
ss_h1_chrom <- ss_h1_chrom[hff_genes,]
# png("h1_chrom_tss_heatmap.png",height = 1500,width = 900,res = 200)
chrom_h1_hm <- Heatmap(ss_h1_chrom,use_raster = F,cluster_columns = F,cluster_rows = T, col =chrom_colors,
        show_heatmap_legend = T,show_row_names = F,show_column_names = F,row_title = "Genes",
        row_order = hff_genes, name = "ChromHMM label",
        column_title = "H1 ChromHMM")
# dev.off()

########## HFF

# H1_ccs


missing_hff_genes <- hff_genes[!(hff_genes %in% rownames(hff_ccs_mat))]
na_mat <- matrix('NA',nrow = length(missing_hff_genes),ncol = 500)
rownames(na_mat) <- missing_hff_genes

hff_ccs_mat <- rbind2(hff_ccs_mat,na_mat)
ss_hff_ccs <- hff_ccs_mat[hff_genes,]

# png("hff_ccs_tss_heatmap.png",height = 1500,width = 900,res = 200)
ccs_hff_hm <- Heatmap(ss_hff_ccs,use_raster = F,cluster_columns = F,cluster_rows = T, col = ccs_colors,
        show_heatmap_legend = F,show_row_names = F,show_column_names = F,row_title = "Genes",
        row_order = hff_genes, 
        column_title = "HFF CIS")
# dev.off()



missing_hff_genes <- hff_genes[!(hff_genes %in% rownames(hff_chrom_mat))]
na_mat <- matrix('NA',nrow = length(missing_hff_genes),ncol = 500)
rownames(na_mat) <- missing_hff_genes

ss_hff_chrom <- rbind2(hff_chrom_mat,na_mat)
ss_hff_chrom <- ss_hff_chrom[hff_genes,]

# png("hff_chrom_tss_heatmap.png",height = 1500,width = 900,res = 200)
chrom_hff_hm <- Heatmap(ss_hff_chrom,use_raster = F,cluster_columns = F,cluster_rows = T, col = chrom_colors,
        show_heatmap_legend = F,show_row_names = F,show_column_names = F,row_title = "Genes",
        row_order = hff_genes, 
        column_title = "HFF ChromHMM")
# dev.off()



#########



missing_hff_genes <- hff_genes[!(hff_genes %in% rownames(hela_ccs_mat))]
na_mat <- matrix('NA',nrow = length(missing_hff_genes),ncol = 500)
rownames(na_mat) <- missing_hff_genes

hela_ccs_mat <- rbind2(hela_ccs_mat,na_mat)
ss_hela_ccs <- hela_ccs_mat[hff_genes,]


# png("hela_ccs_tss_heatmap.png",height = 1500,width = 900,res = 200)
ccs_hela_hm <- Heatmap(ss_hela_ccs,use_raster = F,cluster_columns = F,cluster_rows = T,col = ccs_colors,
        show_heatmap_legend = F,show_row_names = F,show_column_names = F,row_title = "Genes",
        row_order = hff_genes, 
        column_title = "HeLa CIS")
# dev.off()



missing_hff_genes <- hff_genes[!(hff_genes %in% rownames(hela_chrom_mat))]
na_mat <- matrix('NA',nrow = length(missing_hff_genes),ncol = 500)
rownames(na_mat) <- missing_hff_genes

hela_chrom_mat <- rbind2(hela_chrom_mat,na_mat)
ss_hela_chrom <- hela_chrom_mat[hff_genes,]

# png("hela_chrom_tss_heatmap.png",height = 1500,width = 900,res = 200)
chrom_hela_hm <- Heatmap(ss_hela_chrom,use_raster = F,cluster_columns = F,cluster_rows = T,col = chrom_colors,
        show_heatmap_legend = F,show_row_names = F,show_column_names = F,row_title = "Genes",
        row_order = hff_genes, 
        column_title = "HeLa ChromHMM")
# dev.off()



########




missing_hff_genes <- hff_genes[!(hff_genes %in% rownames(endo_ccs_mat))]
na_mat <- matrix('NA',nrow = length(missing_hff_genes),ncol = 500)
rownames(na_mat) <- missing_hff_genes

endo_ccs_mat <- rbind2(endo_ccs_mat,na_mat)
ss_endo_ccs <- endo_ccs_mat[hff_genes,]

# png("endo_ccs_tss_heatmap.png",height = 1500,width = 900,res = 200)
ccs_endo_hm <- Heatmap(ss_endo_ccs,use_raster = F,cluster_columns = F,cluster_rows = T,col = ccs_colors,
        show_heatmap_legend = F,show_row_names = F,show_column_names = F,row_title = "Genes",
        row_order = hff_genes, 
        column_title = "Endoderm CIS")
# dev.off()



missing_hff_genes <- hff_genes[!(hff_genes %in% rownames(endo_chrom_mat))]
na_mat <- matrix('NA',nrow = length(missing_hff_genes),ncol = 500)
rownames(na_mat) <- missing_hff_genes

endo_chrom_mat <- rbind2(endo_chrom_mat,na_mat)
ss_endo_chrom <- endo_chrom_mat[hff_genes,]

# png("endo_chrom_tss_heatmap.png",height = 1500,width = 900,res = 200)
chrom_endo_hm <- Heatmap(ss_endo_chrom,use_raster = F,cluster_columns = F,cluster_rows = T,col = chrom_colors,
        show_heatmap_legend = F,show_row_names = F,show_column_names = F,row_title = "Genes",
        row_order = hff_genes,
        column_title = "Endoderm ChromHMM")
# dev.off()

# 
# png("new_figs/ChromHMM_tss_heatmap_HFF_H1_upreg.png",height = 700,width = 1050,res = 200)
# draw(chrom_h1_hm + chrom_hff_hm, column_title="ChromHMM Labels - significant HFF genes >2 log2FC")
# dev.off()

png("new_figs/tss_heatmap_ALL.png",height = 700,width = 2100,res = 200)
draw(chrom_hff_hm + ccs_hff_hm + chrom_h1_hm + ccs_h1_hm, column_title="All TSS Regions (HFF ordered)")
dev.off()
# 
# png("new_figs/endo_H1_upreg.png",height = 700,width = 2100,res = 200)
# draw(ccs_h1_hm + ccs_endo_hm + chrom_h1_hm + chrom_endo_hm, column_title="Significant Endoderm genes >2 log2FC")
# dev.off()
# 
# png("new_figs/endo_HFF_downreg.png",height = 700,width = 2100,res = 200)
# draw(ccs_hff_hm + ccs_endo_hm + chrom_hff_hm + chrom_endo_hm, column_title="Significant Endoderm genes <-2 log2FC")
# dev.off()


 # png("ChromHMM_tss_heatmap1.png",height = 700,width = 1500,res = 200)
# draw(chrom_h1_hm + chrom_endo_hm + chrom_hff_hm + chrom_hela_hm, column_title="ChromHMM Labels")
# dev.off()
# 
# png("CCS_tss_heatmap.png1",height = 700,width = 1500,res = 200)
# draw(ccs_h1_hm + ccs_endo_hm + ccs_hff_hm + ccs_hela_hm, column_title="CCS Labels")
# dev.off()


