
setwd("~/Documents/UNLV/Year5/create_cluster_matrix/enhancers/expressed_genes/")

library(readr)
library(dplyr)
library(tidyr)


h1_genes <- read_tsv("h1_pone.0192625.s007.tsv",col_types = 'ccccccccccccccccccccc') %>% 
  select(HUES1_TPM,ensembl_id,hgnc_symbol) %>% subset(!is.na(hgnc_symbol)) %>% subset(HUES1_TPM != 0)
colnames(h1_genes) <- c('expr','ensembl_id','hgnc_symbol')
h1_genes$expr <- as.numeric(h1_genes$expr)

endo_genes <- read_tsv("endo_Table_S8C-Table_1.tsv",col_types = 'cccccccc') %>% select(AveExpr,Gene,`Gene Name`)
colnames(endo_genes) <- c('expr','ensembl_id','hgnc_symbol')
endo_genes <- endo_genes %>% subset(!is.na(hgnc_symbol))
endo_genes$expr <- as.numeric(endo_genes$expr)

hff_genes <- read_tsv("hff_GSM2448852_0h_3-aligned.tsv",col_types = 'cccccc') %>% select(`0h_3_FPKM`,gene_id,gene_name) %>% 
  subset(`0h_3_FPKM` != 0)
colnames(hff_genes) <- c('expr','ensembl_id','hgnc_symbol')
hff_genes <- hff_genes %>% subset(!is.na(hgnc_symbol))
hff_genes$expr <- as.numeric(hff_genes$expr)

all_genes <- unique(c(hff_genes$hgnc_symbol,endo_genes$hgnc_symbol,h1_genes$hgnc_symbol))


h1_genes_hi <- h1_genes %>% filter(expr > quantile(expr, 0.75)) # gives you the highest 25%
h1_genes_low <- h1_genes %>% filter(expr < quantile(expr, 0.25)) %>% subset(expr != 0) # gives you the lowest 25%
h1_genes_zero <- h1_genes %>% subset(expr == 0) # zero expression

missing_genes <- all_genes[!(all_genes %in% h1_genes$hgnc_symbol)]
h1_genes_zero <- rbind(h1_genes_zero,data.frame(expr = rep(0,length(missing_genes)),ensembl_id = rep(NA,length(missing_genes)),hgnc_symbol = missing_genes))

endo_genes_hi <- endo_genes %>% filter(expr > quantile(expr, 0.75))
endo_genes_low <- endo_genes %>% filter(expr < quantile(expr, 0.25)) %>% subset(expr != 0)

missing_genes <- all_genes[!(all_genes %in% endo_genes$hgnc_symbol)]
endo_genes_zero <- data.frame(expr = rep(0,length(missing_genes)),ensembl_id = rep(NA,length(missing_genes)),hgnc_symbol = missing_genes)


hff_genes_hi <- hff_genes %>% filter(expr > quantile(expr, 0.75))
hff_genes_low <- hff_genes %>% filter(expr < quantile(expr, 0.25)) %>% subset(expr != 0)
hff_genes_zero <- read_tsv("hff_GSM2448852_0h_3-aligned.tsv") %>% select(expr = `0h_3_FPKM`,ensembl_id = gene_id,hgnc_symbol = gene_name) %>% subset(expr == 0) # zero expression
missing_genes <- all_genes[!(all_genes %in% hff_genes$hgnc_symbol)]
hff_genes_zero <- rbind(hff_genes_zero,data.frame(expr = rep(0,length(missing_genes)),ensembl_id = rep(NA,length(missing_genes)),hgnc_symbol = missing_genes))


write_tsv(h1_genes_hi,"h1_genes_high.tsv")
write_tsv(h1_genes_low,"h1_genes_low.tsv")
write_tsv(h1_genes_zero,"h1_genes_zero.tsv")

write_tsv(endo_genes_hi,"endo_genes_high.tsv")
write_tsv(endo_genes_low,"endo_genes_low.tsv")
write_tsv(endo_genes_zero,"endo_genes_zero.tsv")

write_tsv(hff_genes_hi,"hff_genes_high.tsv")
write_tsv(hff_genes_low,"hff_genes_low.tsv",)
write_tsv(hff_genes_zero,"hff_genes_zero.tsv")

