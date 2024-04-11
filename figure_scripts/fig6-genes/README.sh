
bedtools intersect -a gene_exp/RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.TSS500bp.bed -b ../fig5-all_sankey/create_bed_network_files/h1_labelled.bed -wo > h1_tss_overlap.bed
bedtools intersect -a gene_exp/RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.TSS500bp.bed -b ../fig5-all_sankey/create_bed_network_files/hff_labelled.bed -wo > hff_tss_overlap.bed
bedtools intersect -a gene_exp/RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.TSS500bp.bed -b ../fig5-all_sankey/create_bed_network_files/hela_labelled.bed -wo > hela_tss_overlap.bed
bedtools intersect -a gene_exp/RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.TSS500bp.bed -b ../fig5-all_sankey/create_bed_network_files/endo_labelled.bed -wo > endo_tss_overlap.bed
