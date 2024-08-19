
setwd("~/Documents/UNLV/Year5/new_kmeans_majority_chromhmm_in_window/fig9-enrichment/")
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

usage <- read_tsv("F5.hg38.enhancers.expression.usage.matrix")
enh_bed <- read_tsv("F5.hg38.enhancers.bed",col_names = c("chr",'pos1','pos2','enh_id','score',"strand",
                                                          "thickStart","thickEnd",'x','x2','x3','x4'))

samples <- read_tsv("samples.txt",col_names = c('cell_type','sample','id'))

small_usage <- usage[,colnames(usage) %in% samples$id] 
small_usage$id <- usage$id
small_usage <- small_usage[rowSums(small_usage[,1:6]) > 0,]

small_usage <- pivot_longer(small_usage,cols = 1:6,names_to = "sample",values_to = "enhancer_pa")

enhancer_cell_type <- merge(samples[,c(1,3)],small_usage,by.x = 'id',by.y = 'sample')

all_enhancers <- enhancer_cell_type[enhancer_cell_type$enhancer_pa != 0,]
all_enhancers <- merge(all_enhancers,enh_bed,by.x = 'id.y',by.y = 'enh_id')

# # http://reftss.clst.riken.jp/datafiles/current/human/gene_annotation/refTSS_v3.3_human_coordinate.ann.hg38.rds
# all_tss <- readRDS("promoters/refTSS_v3.3_human_coordinate.ann.hg38.rds")
# 
# all_tss_margin <- GRanges(seqnames = seqnames(all_tss),
#                           range = ranges(all_tss) + 200,
#                           strand = strand(all_tss),
#                           mcols = data.frame(all_tss)[,6:14])
# saveRDS(all_tss_margin,file = "promoters/refTSS_v3.3_human_coordinate.ann.hg38_200margin.rds")

tss <- read_tsv("../fig6-gene_heatmap/gene_exp/RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.TSS500bp.bed",
         col_names = c("chr",'pos1','pos2','gene','score','strand'))

library(GenomicRanges)
all_tss <- GRanges(tss$chr,
                   IRanges(as.numeric(tss$pos1),
                           as.numeric(tss$pos2)),
                   mcols = data.frame(Gene_symbol = tss$gene))

h1_genes_hi <- read_tsv("expressed_genes/h1_genes_high.tsv")$hgnc_symbol
h1_genes_low <- read_tsv("expressed_genes/h1_genes_low.tsv")$hgnc_symbol
h1_genes_zero <- read_tsv("expressed_genes/h1_genes_zero.tsv")$hgnc_symbol

endo_genes_hi <- read_tsv("expressed_genes/endo_genes_high.tsv")$hgnc_symbol
endo_genes_low <- read_tsv("expressed_genes/endo_genes_low.tsv")$hgnc_symbol
endo_genes_zero <- read_tsv("expressed_genes/endo_genes_zero.tsv")$hgnc_symbol


hff_genes_hi <- read_tsv("expressed_genes/hff_genes_high.tsv")$hgnc_symbol
hff_genes_low <- read_tsv("expressed_genes/hff_genes_low.tsv")$hgnc_symbol
hff_genes_zero <- read_tsv("expressed_genes/hff_genes_zero.tsv")$hgnc_symbol


hela_genes_hi <- read_tsv("expressed_genes/hela_genes_high.tsv")$hgnc_symbol
hela_genes_low <- read_tsv("expressed_genes/hela_genes_low.tsv")$hgnc_symbol
hela_genes_zero <- read_tsv("expressed_genes/hela_genes_zero.tsv")$hgnc_symbol

# GET INTERSECT OF CLUSTERS BY ENHANCERS

library(GenomicRanges)
library(AnnotationHub)
library(rtracklayer)


clusters_hela <- read_tsv("../fig5-all_sankey/create_bed_network_files/hela_labelled.bed",
                          col_names = c("chr","pos1","pos2","label","cluster","cell_type"))
clusters_hff <- read_tsv("../fig5-all_sankey/create_bed_network_files/hff_labelled.bed",
                          col_names = c("chr","pos1","pos2","label","cluster","cell_type"))
clusters_h1 <- read_tsv("../fig5-all_sankey/create_bed_network_files/h1_labelled.bed",
                          col_names = c("chr","pos1","pos2","label","cluster","cell_type"))
clusters_endo <- read_tsv("../fig5-all_sankey/create_bed_network_files/endo_labelled.bed",
                          col_names = c("chr","pos1","pos2","label","cluster","cell_type"))
 
get_enrichment <- function(clusters, query_regions,query_type,cell_type) {
  
  # if (!exists("gc_db")) {
  #   
  #   ah <- AnnotationHub()
  #   # query(ah, c("Gencode", "gff", "human","GRCh38","basic"))
  #   # gc <- ah[["AH75120"]]
  #   # gc <- ah[["AH49556"]]
  #   
  #   gc_db <<- ah[["AH75120"]]
  # }
  
  cluster_regions <- GRanges(clusters$chr,
                             IRanges(as.numeric(clusters$pos1),
                                     as.numeric(clusters$pos2)))
  
  querys <- clusters[cluster_regions %over% query_regions,]

  query_counts <- table(querys$cluster)
  clus_counts <- table(clusters$cluster)
  if (length(query_counts) < 18) {
    missing_category <- names(clus_counts)[!(names(clus_counts) %in% names(query_counts))]
    
    for (i in 1:length(missing_category)) {
      query_counts[[missing_category[i]]] = 0  # for hela
      query_counts <- query_counts[match(names(clus_counts),names(query_counts))]
    }
    
  }
  

  query_list <- list()
  for (clus_k in 1:18) {
    # Initialize variables
    # Enhancer test
    m <- clus_counts[clus_k]          # num IN base cluster
    n <- sum(clus_counts) - m         # num NOT IN base cluster
    k <- sum(query_counts)              # total num IN enhancer clusters
    x <- query_counts[clus_k]           # total num IN enhancer and cluster
    
    query_counts2 <- as.table(rbind(query_counts,clus_counts - query_counts))
    dimnames(query_counts2) <- list(enhancer = c("Enh","NonEnh"),
                              cluster = names(query_counts))
    
    query_mat_small <- cbind(`1` = query_counts2[,clus_k], `else` = c(sum(query_counts2[1,setdiff(1:18,clus_k)]),
                                                                   sum(query_counts2[2,setdiff(1:18,clus_k)])))
    
    num_tests <- 13*18
    query_list[[clus_k]] <- c(names(query_counts)[clus_k],query_type,chisq.test(query_mat_small)$stdres[1,1],
                            chisq.test(query_mat_small)$p.value*num_tests,cell_type)
    
    
  }
  a <- data.frame(do.call(rbind,query_list))
  
  a
}

tss_h1_low <- get_enrichment(clusters_h1,
                         all_tss[all_tss$mcols.Gene_symbol %in% h1_genes_low],
                         "low\nexpr\nTSS","H1")
tss_h1_hi <- get_enrichment(clusters_h1,
                             all_tss[all_tss$mcols.Gene_symbol %in% h1_genes_hi],
                             "high\nexpr\nTSS","H1")
tss_h1_zero <- get_enrichment(clusters_h1,
                            all_tss[all_tss$mcols.Gene_symbol %in% h1_genes_zero],
                            "zero\nexpr\nTSS","H1")


tss_hff_low <- get_enrichment(clusters_hff,
                             all_tss[all_tss$mcols.Gene_symbol %in% hff_genes_low],
                             "low\nexpr\nTSS","HFF")
tss_hff_hi <- get_enrichment(clusters_hff,
                            all_tss[all_tss$mcols.Gene_symbol %in% hff_genes_hi],
                            "high\nexpr\nTSS","HFF")
tss_hff_zero <- get_enrichment(clusters_hff,
                                all_tss[all_tss$mcols.Gene_symbol %in% hff_genes_zero],
                                "zero\nexpr\nTSS","HFF")

tss_endo_low <- get_enrichment(clusters_endo,
                             all_tss[all_tss$mcols.Gene_symbol %in% endo_genes_low],
                             "low\nexpr\nTSS","Endoderm")
tss_endo_hi <- get_enrichment(clusters_endo,
                            all_tss[all_tss$mcols.Gene_symbol %in% endo_genes_hi],
                            "high\nexpr\nTSS","Endoderm")
tss_endo_zero <- get_enrichment(clusters_endo,
                               all_tss[all_tss$mcols.Gene_symbol %in% endo_genes_zero],
                               "zero\nexpr\nTSS","Endoderm")

tss_hela_low <- get_enrichment(clusters_hela,
                               all_tss[all_tss$mcols.Gene_symbol %in% endo_genes_low],
                               "low\nexpr\nTSS","HeLa")
tss_hela_hi <- get_enrichment(clusters_hela,
                              all_tss[all_tss$mcols.Gene_symbol %in% endo_genes_hi],
                              "high\nexpr\nTSS","HeLa")
tss_hela_zero <- get_enrichment(clusters_hela,
                                all_tss[all_tss$mcols.Gene_symbol %in% endo_genes_zero],
                                "zero\nexpr\nTSS","HeLa")

# tss_hela <- get_enrichment(clusters_hela,all_tss,"all TSS","HeLa")

h1_enhancers <- all_enhancers %>% subset(cell_type == "H1")
enh_h1 <- get_enrichment(clusters_h1,
                         GRanges(h1_enhancers$chr,
                                 IRanges(as.numeric(h1_enhancers$pos1),
                                         as.numeric(h1_enhancers$pos2)))
                         ,"active\nenhancer","H1")

hela_enhancers <- all_enhancers %>% subset(cell_type == "HeLa")
enh_hela <- get_enrichment(clusters_hela,
                           GRanges(hela_enhancers$chr,
                                   IRanges(as.numeric(hela_enhancers$pos1),
                                           as.numeric(hela_enhancers$pos2)))
                         ,"active\nenhancer","HeLa")


h1_se <- read_tsv("SEdb/h1_SE.bed")
se_h1 <- get_enrichment(clusters_h1,
                        GRanges(h1_se$se_chr,
                                            IRanges(as.numeric(h1_se$se_start),
                                                    as.numeric(h1_se$se_end)))
                        ,"SE","H1")

hff_se <- read_tsv("SEdb/hff_SE.bed")
se_hff <- get_enrichment(clusters_hff,
                        GRanges(hff_se$se_chr,
                                IRanges(as.numeric(hff_se$se_start),
                                        as.numeric(hff_se$se_end)))
                        ,"SE","HFF")

endo_se <- read_tsv("SEdb/endoderm_SE.bed")
se_endo <- get_enrichment(clusters_endo,
                        GRanges(endo_se$se_chr,
                                IRanges(as.numeric(endo_se$se_start),
                                        as.numeric(endo_se$se_end)))
                        ,"SE","Endoderm")

hela_se <- read_tsv("SEdb/hela_SE.bed")
se_hela <- get_enrichment(clusters_hela,
                        GRanges(hela_se$se_chr,
                                IRanges(as.numeric(hela_se$se_start),
                                        as.numeric(hela_se$se_end)))
                        ,"SE","HeLa")



z <- bind_rows(se_hela,se_endo,se_hff,se_h1,
               enh_hela,enh_h1,
               tss_endo_hi,tss_endo_low,tss_endo_zero,
               tss_hela_hi,tss_hela_low,tss_hela_zero,
               tss_h1_hi,tss_h1_low,tss_h1_zero,
               tss_hff_hi,tss_hff_low,tss_hff_zero)
colnames(z) <- c("cluster","type","stdres","p.value","cell_type")


z$p.value <- as.numeric(z$p.value)
z$stdres <- as.numeric(z$stdres)

ord <- rev(c("cont_Tss_Enh","cont_Tss_EnhG","cont_Tss_EnhWk","cont_Tss_noEnh",
             "cont_EnhA","cont_EnhA_EnhWk","cont_Tx_EnhG1","cont_Tx_EnhG2",
             "cont_TxWk_Enh","cont_TxWk_EnhWk","cont_Low_TssBiv","cont_Low_EnhWk",
             "cont_Biv","cont_ReprPC_Biv","cont_ReprPC",
             "cont_Het_Rpts","cont_Het_Rpts_strong","cont_Het_Rpts_strongest"))


z$cluster <- factor( z$cluster, levels = ord)

z$type <- factor( z$type, levels = c("active\nenhancer", "SE", "high\nexpr\nTSS", "low\nexpr\nTSS",  "zero\nexpr\nTSS"))
# z$stdres <- ifelse(z$stdres > 50, 50, z$stdres)


 png("enrichment_enh_tss_updated.png",height = 1700,width =2700,res = 300)
 # png("enrichment_enh_tss_large.png",height = 1800,width =4000,res = 500)
 # png("enrichment_enh_tss_poster.png",height = 800,width =1600,res = 200)

 ggplot(z, aes(y = cluster,x = type,size = -log10(p.value + 1e-250), color = stdres)) +
   geom_point() + 
   scale_color_gradient2(midpoint = 0,low = 'darkblue',mid = 'lightgray',high='red',
                         limits = c(-50,155)) +
   labs(title = "Enrichment Chi Squared Test",x = "Element type",
        size = '-log(pval)',color = "Std Residuals") +
   facet_grid(~cell_type,space = 'free',scales = 'free')
  dev.off()


 
