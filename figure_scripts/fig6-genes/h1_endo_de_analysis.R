
setwd("~/Documents/UNLV/Year5/create_cluster_matrix/gene_exp/h1_dec")

library(readr)

library(GenomicRanges)
library(AnnotationHub)
library(ggplot2)


# OVERLAP EACH SEPARATELY WITH TSS (H1, Endo) Some sort of allowance (75% of TSS)
# Then merge by TSS and look at fold change.



total_int <- read_tsv(paste0("../../create_bed_network_files/network_h1.bed"), 
                      col_names = F, col_types = 'cccccccccccccc')

colnames(total_int) <- c("1_chr","1_pos1","1_pos2","1_base_region","1_cluster","1_cell_type","match",
                         "2_chr","2_pos1","2_pos2","2_base_region","2_cluster","2_cell_type","overlap")
total_int$chr_type <- paste(total_int$`1_chr`,total_int$`1_cell_type`,sep='_')
h1_endo_int <- total_int %>% subset(`2_cell_type` == "endo") 

h1_endo_int <- unique(h1_endo_int[,c(1,2,3,4,5,11,12)])
# h1_endo_int_m <- total_int %>% subset(`1_chr` == `2_chr` &
#                                                                           `1_pos1` == `2_pos1` &
#                                                                           `1_pos2` == `2_pos2`)


h1_contacts_annotated <- read_tsv("../DEC_H1_de_overlaps.tsv")
h1_contacts_annotated <- h1_contacts_annotated[complete.cases(h1_contacts_annotated),] %>%
  subset(cell_type == 'h1')

logfoldchange_h1_endo <- merge(h1_endo_int,h1_contacts_annotated,
                                by.x = c("1_chr","1_pos1","1_pos2"),
                                by.y = c('chr','pos1','pos2'))


# quantify changes of CCS in DE genes, and changes in chromatin state

sig_h1_conts <- logfoldchange_h1_endo %>% subset(padj < 1e-5)

# volcano plot
ggplot(sig_h1_conts,aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point()

active_ccs <- c("Enh_cont","EnhA2_cont","EnhG_cont","EnhWk_cont","EnhWk_cont2",
                "TSS_cont","TSS_Enh_cont","TSS_noEnh_cont","TSSFlnkD_cont",
                "Tx_cont","TxWk_cont")
inactive_ccs <- c("Biv_cont1", "Biv_cont2","Low_cont", "Het_cont","EnhWk_cont","EnhWk_cont2")


active_chromhmm <- c("EnhA1", "EnhA2", "EnhG1","EnhG2","EnhWk","TssA","TssFlnk",
                     "TssFlnkD", "TssFlnkU","Tx","TxWk")
inactive_chromhmm <- c("EnhWk","TxWk","EnhBiv", "TssBiv", "Het", "ReprPC", "ReprPCWk",  "ZNF-Rpts")


upreg <- sig_h1_conts %>% subset(log2FoldChange > 2)

upreg_ss <- upreg[upreg$`2_cluster` %in%  active_ccs &
                    upreg$`1_cluster` %in% inactive_ccs &
                    upreg$`2_cluster` != upreg$`1_cluster`,]

upreg_chromhmm_ss <- upreg[upreg$`2_base_region` %in% active_chromhmm &
                             upreg$`1_base_region` %in% inactive_chromhmm &
                             upreg$`1_base_region` != upreg$`2_base_region`,]


upreg_both <- upreg[upreg$`2_base_region` %in% active_chromhmm &
                      upreg$`1_base_region` %in% inactive_chromhmm &
                      upreg$`2_cluster` %in%  active_ccs &
                      upreg$`1_cluster` %in% inactive_ccs,]

length(unique(upreg_ss$genes))/length(unique(upreg$genes))
length(unique(upreg_chromhmm_ss$genes))/length(unique(upreg$genes))
length(unique(c(upreg_chromhmm_ss$genes,upreg_ss$genes)))/length(unique(upreg$genes))
length(unique(upreg_both$genes))/length(unique(upreg$genes))


# upreg <- upreg[upreg$`1_cluster` != upreg$`2_cluster` &
                 # upreg$`1_base_region` == upreg$`2_base_region`,]


row_ann <- HeatmapAnnotation("H1 Base Label" = upreg$`1_base_region`,
                             "Endo Base Label" = upreg$`2_base_region`,
                             col = list("H1 Base Label" = c("EnhA1"="red","EnhA2" = "red4",
                                                              "EnhBiv" = "purple", "EnhG1" = "palevioletred2",
                                                              "EnhG2" = "palevioletred", "EnhWk" = "rosybrown1",
                                                              "Het" = "gray90", "Quies" = "gray48",
                                                              "ReprPC" = "gray", "ReprPCWk" = "gray60",
                                                              "TssA" = "yellow", "TssBiv" = "plum", "TssFlnk" = "orange",
                                                              "TssFlnkD" ="lightyellow","TssFlnkU" = "khaki", 
                                                              "Tx" = "green3", "TxWk" = "lightgreen", "ZNF-Rpts" = "black"),
                                        "Endo Base Label" = c("EnhA1"="red","EnhA2" = "red4",
                                                             "EnhBiv" = "purple", "EnhG1" = "palevioletred2",
                                                             "EnhG2" = "palevioletred", "EnhWk" = "rosybrown1",
                                                             "Het" = "gray90", "Quies" = "gray48",
                                                             "ReprPC" = "gray", "ReprPCWk" = "gray60",
                                                             "TssA" = "yellow", "TssBiv" = "plum", "TssFlnk" = "orange",
                                                             "TssFlnkD" ="lightyellow","TssFlnkU" = "khaki", 
                                                             "Tx" = "green3", "TxWk" = "lightgreen", "ZNF-Rpts" = "black")
                             ),
                             annotation_legend_param = list(nrow = 9),which= 'row')

library(RColorBrewer)
pal = brewer.pal(n=8,"Dark2")
pal = c("purple","orange","rosybrown1","red","lightgreen",
        "yellow",'green3','gray90')

colors = c("Biv_cont1" = pal[1], "Biv_cont2" = pal[1],
           "Het_cont" = 'black', "Low_cont" = pal[8],
           
           "EnhWk_cont" = pal[3],
           "EnhWk_cont2" = pal[3],
           
           "Enh_cont" = pal[4], "EnhA2_cont" = pal[4],
           
           "EnhG_cont" = pal[5],
           "Tx_cont" = pal[5],
           
           "TxWk_cont" = pal[7],
           
           "TSS_cont" = pal[6],
           "TSS_Enh_cont" = pal[6],
           
           "TSSFlnkD_cont" = pal[2],
           "TSS_noEnh_cont" = pal[2]
           
)  


# png("Endo_H1_de-genes_upreg.png",height = 1600,width = 1000,res = 200)
Heatmap(as.matrix(cbind.data.frame(upreg$`1_cluster`,upreg$`2_cluster`)), col = colors,
        left_annotation = row_ann,
        # heatmap_legend_param = list(
        #   nrow = 8
        # ),
        row_order = order(
          upreg$`1_cluster`,
          upreg$`2_cluster`),
        name = "CCS",
        column_labels = c("H1 CCS","Endo CCS"),
        column_names_rot = 0,
        column_names_centered = T,
        row_title = "Base Labeled regions",
        column_title = "Endo / H1 > 1 log2 fold change")
# dev.off()

downreg <- sig_h1_conts %>% subset(log2FoldChange < -2)

row_ann <- HeatmapAnnotation("H1 Base Label" = downreg$`1_base_region`,
                             "Endo Base Label" = downreg$`2_base_region`,
                             col = list("H1 Base Label" = c("EnhA1"="red","EnhA2" = "red4",
                                                              "EnhBiv" = "purple", "EnhG1" = "palevioletred2",
                                                              "EnhG2" = "palevioletred", "EnhWk" = "rosybrown1",
                                                              "Het" = "gray90", "Quies" = "gray48",
                                                              "ReprPC" = "gray", "ReprPCWk" = "gray60",
                                                              "TssA" = "yellow", "TssBiv" = "plum", "TssFlnk" = "orange",
                                                              "TssFlnkD" ="lightyellow","TssFlnkU" = "khaki", 
                                                              "Tx" = "green3", "TxWk" = "lightgreen", "ZNF-Rpts" = "black"),
                                        "Endo Base Label" = c("EnhA1"="red","EnhA2" = "red4",
                                                             "EnhBiv" = "purple", "EnhG1" = "palevioletred2",
                                                             "EnhG2" = "palevioletred", "EnhWk" = "rosybrown1",
                                                             "Het" = "gray90", "Quies" = "gray48",
                                                             "ReprPC" = "gray", "ReprPCWk" = "gray60",
                                                             "TssA" = "yellow", "TssBiv" = "plum", "TssFlnk" = "orange",
                                                             "TssFlnkD" ="lightyellow","TssFlnkU" = "khaki", 
                                                             "Tx" = "green3", "TxWk" = "lightgreen", "ZNF-Rpts" = "black")
                             ),
                             annotation_legend_param = list(nrow = 9),which= 'row')


png("Endo_H1_de-genes_downreg.png",height = 1600,width = 1000,res = 200)
Heatmap(as.matrix(cbind.data.frame(downreg$`1_cluster`,downreg$`2_cluster`)), col = colors,
        left_annotation = row_ann,
        # heatmap_legend_param = list(
        #   nrow = 8
        # ),
        row_order = order(
          downreg$`2_cluster`,
          downreg$`1_cluster`),
        name = "CCS",
        column_labels = c("H1 CCS","Endo CCS"),
        column_names_rot = 0,
        column_names_centered = T,
        row_title = "Base Labeled regions",
        column_title = "Endo / H1 < -2 log2 fold change")
dev.off()

draw(upreg_hm + downreg_hm)



# endo_gene_exp <- endo_gene_exp %>% subset(`2_cell_type` != 'hela' & !(`1_cluster` %in% c(2,4,11,3))) %>%
#                   subset(`1_base_region` == `2_base_region` & `1_cluster` != `2_cluster`) %>%
#                   select(`1_chr`,`1_pos1`,`1_pos2`,genes,`1_cell_type`,`2_cell_type`,
#                                                                `1_cluster`,`2_cluster`,`1_base_region`,`2_base_region`,mean_H1,mean_DEC)
# 
# View(endo_gene_exp)
# 
# 
# 
# #grep 'chr1_.*:225941700' matrix_connections.tsv
# # chr1_endo:225941700-225943100	EnhWk	0.4795382018	0.007187582384	0.1770544833	0	0	0.4387969679	0	0	0	0.04444331483	0.004254884231	0.3605454616	0.01064881107	0.1403565823	0.3464735644	0.004733154195	0.01346294607	0
# # chr1_h1:225941700-225943100	EnhWk	0.01652899671	0.01616454175	0.3650429063	0.009862764189	0	0.1289330546	0	0.004824239066	0.536966902	0.7620327212	0.042133789	0.08659379923	0.01600825793	0.1571331361	0.01249728576	0.006648006743	0.02883481405	0.01310408587
# endo_values <- c(0.4795382018,	0.007187582384,	0.1770544833,	0,	0,	0.4387969679,	0,	0,	0,	0.04444331483,
#           0.004254884231,	0.3605454616,	0.01064881107,	0.1403565823,	0.3464735644,	0.004733154195,	0.01346294607,	0)
# 
# h1_values <- c(0.01652899671,	0.01616454175,	0.3650429063,	0.009862764189,	0,
#                  0.1289330546,	0,	0.004824239066,	0.536966902,	0.7620327212,	0.042133789,	0.08659379923,
#                  0.01600825793,	0.1571331361,	0.01249728576,	0.00664800674,	0.02883481405,	0.01310408587)
# diffs <- data.frame(difference = c(endo_values - h1_values))
# diffs$contact_type <- factor(c("EnhA1",	"EnhA2",	"EnhBiv",	"EnhG1",
#                         "EnhG2",	"EnhWk",	"Het",	"Quies",
#                         "ReprPC",	"ReprPCWk",	"TssA",	"TssBiv",
#                         "TssFlnkD",	"TssFlnk",	"TssFlnkU",	"Tx",
#                         "TxWk",	"ZNF/Rpts"),levels = c(
#                           "EnhA1",	"EnhA2","EnhWk","TxWk","TssA",
#                           "TssFlnkU","TssFlnk","TssFlnkD","EnhG1",
#                           "EnhG2","Tx",	"EnhBiv",	"TssBiv","ReprPC",	"ReprPCWk",
#                           "Het",	"Quies","ZNF/Rpts"
#                         ))
# 
# hot = 'orange2'
# cold = 'skyblue3'
# 
# png("LEFTY1_anecdote.png",height = 600,width = 900,res = 150)
# ggplot(diffs,aes(x = contact_type,y = difference,fill = contact_type)) +
#   geom_bar(stat = "identity") +   
#   scale_fill_manual("Contact", values = c("EnhA1" = hot,	"EnhA2" = hot,	"EnhBiv" = cold,	"EnhG1" = hot,
#                                           "EnhG2" = hot,	"EnhWk" = hot,	"Het" = cold,	"Quies" = cold,
#                                           "ReprPC"= cold,	"ReprPCWk"= cold,	"TssA" = hot,	"TssBiv" = cold,
#                                           "TssFlnkD" = hot,	"TssFlnk" = hot,	"TssFlnkU" = hot,	"Tx" = hot,
#                                           "TxWk"= hot,	"ZNF/Rpts" = cold)) +
#   labs(title = "Endoderm contacts - H1 contacts (chr1:225941700-225943100)",
#        subtitle = "region in contact with LEFTY1 (endoderm gene)",
#        y = "Contact Strength Difference",
#        x = "Contact Region Type") + theme(axis.text.x=element_text(angle = 70, hjust = 1),legend.position = "none")
# dev.off()
# 
