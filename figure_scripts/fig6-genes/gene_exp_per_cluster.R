# gene expression per cluster


setwd("~/Documents/UNLV/Year5/new_kmeans_majority_chromhmm_in_window/fig6-gene_heatmap/gene_exp/")

library(tidyverse)
library(ggsignif)

# all_cells <- read_tsv("../mat_raw.bed",col_names = c("chr","pos1","pos2","cell_type","region"))

all_cells <- readRDS("../../data/mat_raw_mean.rds") %>% separate(region,
                                                    into = c("chr","cell_type","pos1","pos2"),
                                                    remove = F)

# RAW COUNT NORMALIZATION
exp_data_raw <- read_csv("GSE75748_bulk_cell_type_ec.csv")
exp_data <- round(exp_data_raw[,c(2,3,4,5,9,10,14,15,16)])

meta <- cbind.data.frame(cell_type = c(rep("H1",4),"DEC","DEC",rep("HFF",3)))
rownames(meta) <- colnames(exp_data)
rownames(exp_data) <- exp_data_raw$gene

library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData = exp_data, colData = meta, design = ~ cell_type)
dds <- estimateSizeFactors(dds)

dds <- DESeq(dds)

# gene_expression_H1_DEC <- data.frame(results(dds,contrast = c("cell_type","H1","DEC")))
gene_expression_DEC_H1 <- data.frame(results(dds,contrast = c("cell_type","DEC","H1")))

gene_expression_DEC_HFF <- data.frame(results(dds,contrast = c("cell_type","DEC","HFF")))
gene_expression_HFF_H1 <- data.frame(results(dds,contrast = c("cell_type","HFF","H1")))

# gene_expression_H1_DEC$gene <- rownames(gene_expression_H1_DEC)
gene_expression_DEC_HFF$gene <- rownames(gene_expression_DEC_HFF)
gene_expression_DEC_H1$gene <- rownames(gene_expression_DEC_H1)
gene_expression_HFF_H1$gene <- rownames(gene_expression_HFF_H1)



# write_tsv(gene_expression_H1_DEC,"h1_dec_de.tsv")
write_tsv(gene_expression_DEC_H1,"dec_h1_de.tsv")
write_tsv(gene_expression_DEC_HFF,"dec_hff_de.tsv")
write_tsv(gene_expression_HFF_H1,"hff_h1_de.tsv")

# 
# normalized_counts <- data.frame(counts(dds, normalized=TRUE))
# 
# normalized_counts$gene <- exp_data_raw$gene
# normalized_counts <- normalized_counts %>% rowwise() %>% 
#   mutate(mean_H1 = mean(c(H1_rep1,H1_rep2,H1_rep3,H1_rep4)),
#          mean_DEC =mean(c(DEC_rep1,DEC_rep2)),
#          mean_HFF = mean(c(HFF_rep1,HFF_rep2,HFF_rep3)))
# 
# write_tsv(normalized_counts,"gene_expression_normalized.tsv")
# 
# 
## INCLUDE GENE OVERLAPS

write_tsv(all_cells[,c(2,4,5,3,6,1)],"mat_raw.bed",col_names = F)
# bedtools intersect -a RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.TSS500bp.bed -b mat_raw.bed -f .51 -wo > all_tss_overlap_mat.bed

overlaps <-read_tsv("all_tss_overlap_mat.bed",
                    col_names = c('tss_chr','tss_pos1','tss_pos2','gene','score','strand',
                                  'chr','pos1','pos2','cell_type',
                                  'base_region','region','overlap')) %>% separate_rows(gene,sep = ',')
annot_regions <- merge(overlaps,all_cells,
      by = 'region')

h1_reg <- annot_regions %>% subset(cell_type.x == 'h1')
h1_labelled <- read_tsv("../../fig4-all_sankey/create_bed_network_files/h1_labelled.bed", 
                     col_names = c("chr","pos1","pos2","base_region","ccs_h1","cell_type"))
h1_all <- merge(h1_labelled,h1_reg,by.x = c("chr","pos1","pos2"), by.y =c("chr.x","pos1.x","pos2.x")) %>%
  select(chr,pos1,pos2,base_region = base_region.x,ccs_h1,gene = gene) %>% unique()



hff_reg <- annot_regions %>% subset(cell_type.x == 'hff')
hff_labelled <- read_tsv("../../fig4-all_sankey/create_bed_network_files/hff_labelled.bed", 
                        col_names = c("chr","pos1","pos2","base_region","ccs_hff","cell_type"))
hff_all <- merge(hff_labelled,hff_reg,by.x = c("chr","pos1","pos2"), by.y =c("chr.x","pos1.x","pos2.x")) %>%
  select(chr,pos1,pos2,base_region = base_region.x,ccs_hff,gene = gene) %>% unique()



endo_reg <- annot_regions %>% subset(cell_type.x == 'endo')
endo_labelled <- read_tsv("../../fig4-all_sankey/create_bed_network_files/endo_labelled.bed", 
                        col_names = c("chr","pos1","pos2","base_region","ccs_endo","cell_type"))
endo_all <- merge(endo_labelled,endo_reg,by.x = c("chr","pos1","pos2"), by.y =c("chr.x","pos1.x","pos2.x")) %>%
  select(chr,pos1,pos2,base_region = base_region.x,ccs_endo,gene = gene) %>% unique()


m <- merge(endo_all,h1_all,by = 'gene')
mm <- merge(m,gene_expression_DEC_H1,by = 'gene')

saveRDS(mm,"dec_h1_all_merged.rds")

all <- mm %>% select(gene,ccs_endo,ccs_h1) %>% unique()
all_mat <- table(all$ccs_endo,all$ccs_h1) + 1e-5
all_mat <- all_mat / sum(all_mat)



m <- merge(endo_all,hff_all,by = 'gene')
mm <- merge(m,gene_expression_DEC_HFF,by = 'gene')
saveRDS(mm,"dec_hff_all_merged.rds")

all <- mm %>% select(gene,ccs_endo,ccs_hff) %>% unique()
all_mat <- table(all$ccs_endo,all$ccs_hff) + 1e-5
all_mat <- all_mat / sum(all_mat)

# all <- mm %>% select(gene,base_region.x,base_region.y) %>% unique()
# all_mat <- table(all$base_region.x,all$base_region.y) + 1e-5

x <- mm %>% subset(padj < 1e-5 & log2FoldChange > 2) %>% select(gene,ccs_endo,ccs_hff) %>% unique()
# x <- mm %>% subset(padj < 1e-5 & log2FoldChange > 3) %>% select(gene,base_region.x,base_region.y) %>% unique()


# lab_order <- c("Het_cont","Biv_cont1","Biv_cont2","Low_cont",
#                "EnhWk_cont","EnhWk_cont2",
#                "Tx_cont","TxWk_cont","EnhG_cont","Enh_cont",
#                "EnhA2_cont","TSS_cont","TSS_Enh_cont",
#                "TSS_noEnh_cont","TSSFlnkD_cont")
# 
# 
# mat <- table(x$ccs_endo,x$ccs_hff)
# # mat <- table(x$base_region.x,x$base_region.y)
# 
# all_mat_ss <- all_mat[match(rownames(mat),rownames(all_mat)),match(colnames(mat),colnames(all_mat))]
# all_mat_ss <- all_mat_ss * sum(mat)
# scaled_mat = mat/all_mat_ss
# # scaled_mat = mat
# 
# scaled_mat=scaled_mat[,match(lab_order,colnames(scaled_mat))]
# scaled_mat=scaled_mat[match(lab_order,rownames(scaled_mat)),]
# colnames(scaled_mat) <- lab_order
# rownames(scaled_mat) <- lab_order
# 
# 
# row.subsections <- c(3,3,3,2,4)
# row_split = data.frame(rep(c("A_Repressed", "B_Weak", "C_Transcription", "D_Enhancer", "E_TSS"), row.subsections))
# 
# 
# col_fun = colorRamp2(c(0, 1, 2), c("blue","gray95","red"))
# # col_fun = colorRamp2(c(0, 1,50), c("black","gray95","red"))
# greater = Heatmap(scaled_mat,cluster_rows = F,cluster_columns = F,
#         col = col_fun,
#         # name = '% of total regions\nwith that change',
#         name = 'ratio obs/exp',
#         cluster_row_slices = FALSE,
#         cluster_column_slices = FALSE,
#         # row_order = lab_order, column_order = lab_order,
#         split=row_split,
#         row_title = "Endoderm CCS",gap = unit(0.25,'cm'),
#         column_title = "H1 CCS",
#         row_title_side = 'left',
#         column_title_side = "bottom",
#         show_heatmap_legend = F,
#         cell_fun = function(j, i, x, y, width, height, fill) {
#           grid.text(sprintf("%.1f", scaled_mat[i, j]), x, y, gp = gpar(fontsize = 9))
#         }
#         )
# 
# 
# 
# x <- mm %>% subset(padj < 1e-5 & log2FoldChange < 2) %>% select(gene,ccs_endo,ccs_h1) %>% unique()
# 
# 
# lab_order <- c("Het_cont","Biv_cont1","Biv_cont2","Low_cont",
#                "EnhWk_cont","EnhWk_cont2",
#                "Tx_cont","TxWk_cont","EnhG_cont","Enh_cont",
#                "EnhA2_cont","TSS_cont","TSS_Enh_cont",
#                "TSS_noEnh_cont","TSSFlnkD_cont")
# 
# 
# mat <- table(x$ccs_endo,x$ccs_h1)
# 
# all_mat_ss <- all_mat[match(rownames(mat),rownames(all_mat)),match(colnames(mat),colnames(all_mat))]
# all_mat_ss <- all_mat_ss * sum(mat)
# scaled_mat = mat/all_mat_ss
# # scaled_mat = mat
# 
# scaled_mat=scaled_mat[,match(lab_order,colnames(scaled_mat))]
# scaled_mat=scaled_mat[match(lab_order,rownames(scaled_mat)),]
# colnames(scaled_mat) <- lab_order
# rownames(scaled_mat) <- lab_order
# 
# row.subsections <- c(3,3,3,2,4)
# row_split = data.frame(rep(c("A_Repressed", "B_Weak", "C_Transcription", "D_Enhancer", "E_TSS"), row.subsections))
# 
# 
# col_fun = colorRamp2(c(0, 1, 2), c("blue","gray95","red"))
# # col_fun = colorRamp2(c(0, 1,50), c("black","gray95","red"))
# 
# less = Heatmap(t(scaled_mat),cluster_rows = F,cluster_columns = F,
#         col = col_fun,
#         name = 'less',
#         cluster_row_slices = FALSE,
#         cluster_column_slices = FALSE,
#         row_order = lab_order, column_order = lab_order,
#         row_title = "H1 CCS",gap = unit(0.25,'cm'),
#         column_title = "Endoderm CCS",
#         row_title_side = 'left',
#         column_title_side = "bottom",
#         split=row_split,
#         show_heatmap_legend = F,
#         cell_fun = function(j, i, x, y, width, height, fill) {
#           grid.text(sprintf("%.1f", t(scaled_mat)[i, j]), x, y, gp = gpar(fontsize = 9))
#         }
# )
# 
# 
# 
# 
# x <- mm %>% subset(padj > 0.01 & log2FoldChange > -1 & log2FoldChange < 1) %>% select(gene,ccs_endo,ccs_hff) %>% unique()
# # x <- mm %>% select(gene,ccs_endo,ccs_hff) %>% unique()
# 
# lab_order <- c("Het_cont","Biv_cont1","Biv_cont2","Low_cont",
#                "EnhWk_cont","EnhWk_cont2",
#                "Tx_cont","TxWk_cont","EnhG_cont","Enh_cont",
#                "EnhA2_cont","TSS_cont","TSS_Enh_cont",
#                "TSS_noEnh_cont","TSSFlnkD_cont")
# 
# 
# mat <- table(x$ccs_endo,x$ccs_hff)
# 
# all_mat_ss <- all_mat[match(rownames(mat),rownames(all_mat)),match(colnames(mat),colnames(all_mat))]
# all_mat_ss <- all_mat_ss * sum(mat)
# scaled_mat = mat/all_mat_ss
# # scaled_mat = mat
# 
# scaled_mat=scaled_mat[,match(lab_order,colnames(scaled_mat))]
# scaled_mat=scaled_mat[match(lab_order,rownames(scaled_mat)),]
# colnames(scaled_mat) <- lab_order
# rownames(scaled_mat) <- lab_order
# 
# row.subsections <- c(3,3,3,2,4)
# row_split = data.frame(rep(c("A_Repressed", "B_Weak", "C_Transcription", "D_Enhancer", "E_TSS"), row.subsections))
# 
# 
# col_fun = colorRamp2(c(0, 1, 2), c("blue","gray95","red"))
# # col_fun = colorRamp2(c(0, 1,50), c("black","gray95","red"))
# 
# null = Heatmap(scaled_mat,cluster_rows = F,cluster_columns = F,
#                col = col_fun,
#                name = 'null',show_heatmap_legend = F,
#                cluster_row_slices = FALSE,
#                cluster_column_slices = FALSE,
#                column_title = "not significant",
#                row_order = lab_order, column_order = lab_order,
#                split=row_split,
#                row_title = " ",gap = unit(0.25,'cm'),
#                row_title_side = 'left',
#                cell_fun = function(j, i, x, y, width, height, fill) {
#                  grid.text(sprintf("%.1f", scaled_mat[i, j]), x, y, gp = gpar(fontsize = 9))
#                }
# )
# 
# 
# draw(greater + less,column_title = "HFF CCS (denominator)",row_title = "Endoderm CCS (numerator)")
# 
# library("gridExtra")
# library("cowplot")
# 
# gr_h1 = grid.grabExpr(draw(greater,column_title = "Endo/H1 >2 log2 fold change"))
# le_h1 = grid.grabExpr(draw(less,column_title ="H1/Endo >2 log2 fold change"))
# 
# plot_grid(gr_hff, le_hff,
#           gr_h1, le_h1,
#           labels = "AUTO", ncol = 2, nrow = 2)
# 
# 
# ######################
# x <- mm %>% select(gene,ccs_endo,
#                    ccs_hff,
#                    base_endo = base_region.x,
#                    base_hff = base_region.y,
#                    padj,
#                    log2FoldChange) %>% unique()
#                                                                                      
# upreg <- x %>% subset(padj < 1e-5 &log2FoldChange < -2)
# 
# row_ann <- HeatmapAnnotation("Endo Base Label" = upreg$`base_endo`,
#                              "HFF Base Label" = upreg$`base_hff`,
#                              col = list("HFF Base Label" = c("EnhA1"="red","EnhA2" = "red4",
#                                                             "EnhBiv" = "purple", "EnhG1" = "palevioletred2",
#                                                             "EnhG2" = "palevioletred", "EnhWk" = "rosybrown1",
#                                                             "Het" = "gray90", "Quies" = "gray48",
#                                                             "ReprPC" = "gray", "ReprPCWk" = "gray60",
#                                                             "TssA" = "yellow", "TssBiv" = "plum", "TssFlnk" = "orange",
#                                                             "TssFlnkD" ="lightyellow","TssFlnkU" = "khaki", 
#                                                             "Tx" = "green3", "TxWk" = "lightgreen", "ZNF-Rpts" = "black"),
#                                         "Endo Base Label" = c("EnhA1"="red","EnhA2" = "red4",
#                                                               "EnhBiv" = "purple", "EnhG1" = "palevioletred2",
#                                                               "EnhG2" = "palevioletred", "EnhWk" = "rosybrown1",
#                                                               "Het" = "gray90", "Quies" = "gray48",
#                                                               "ReprPC" = "gray", "ReprPCWk" = "gray60",
#                                                               "TssA" = "yellow", "TssBiv" = "plum", "TssFlnk" = "orange",
#                                                               "TssFlnkD" ="lightyellow","TssFlnkU" = "khaki", 
#                                                               "Tx" = "green3", "TxWk" = "lightgreen", "ZNF-Rpts" = "black")
#                              ),
#                              annotation_legend_param = list(nrow = 9),which= 'row')
# 
# library(RColorBrewer)
# pal = brewer.pal(n=8,"Dark2")
# pal = c("purple","orange","rosybrown1","red","lightgreen",
#         "yellow",'green3','gray90')
# 
# colors = c("Biv_cont1" = pal[1], "Biv_cont2" = pal[1],
#            "Het_cont" = 'black', "Low_cont" = pal[8],
#            
#            "EnhWk_cont" = pal[3],
#            "EnhWk_cont2" = pal[3],
#            
#            "Enh_cont" = pal[4], "EnhA2_cont" = pal[4],
#            
#            "EnhG_cont" = pal[5],
#            "Tx_cont" = pal[5],
#            
#            "TxWk_cont" = pal[7],
#            
#            "TSS_cont" = pal[6],
#            "TSS_Enh_cont" = pal[6],
#            
#            "TSSFlnkD_cont" = pal[2],
#            "TSS_noEnh_cont" = pal[2]
#            
# )  
# 
# 
# png("h1_dec/Endo_HFF_de-genes_downreg_sort2.png",height = 1600,width = 1000,res = 200)
# Heatmap(as.matrix(cbind.data.frame(upreg$ccs_endo,upreg$ccs_hff)), col = colors,
#         left_annotation = row_ann,
#         # heatmap_legend_param = list(
#         #   nrow = 8
#         # ),
#         row_order = order(
#           upreg$ccs_hff,
#           upreg$ccs_endo,
#           upreg$base_endo,
#           upreg$base_hff),
#         name = "CCS",
#         column_labels = c("Endo CCS","HFF CCS"),
#         column_names_rot = 0,
#         column_names_centered = T,
#         row_title = "Base Labeled regions",
#         column_title = "Endo / HFF <-2 log2 fold change, sort hff css, endo ccs")
# dev.off()
# 

