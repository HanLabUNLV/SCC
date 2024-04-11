#!/usr/bin/env python

import sys

def get_genes(line):
	genes = set()
	ll = line.strip().split()
	li = ll[8].split(';')
	for i in li:
		if i.startswith("connected_gene="):
			genes.add(i.split('=')[1])
	return genes

gh_ids = set()
with open("enh_ids_endo_upreg.tsv",'r') as infile:
	for line in infile:
		gh_ids.add(line.strip())

outfile = open("candidate_target_genes_de_hff.tsv",'w')

total_geneset = set()
with open("upreg_endo_net_enhancers.bed",'r') as infile:
	for line in infile:
		enh = line[line.index('GH'):line.index('GH') + 11]
		if enh in gh_ids:
			genes = get_genes(line)
			print(line.strip() + '\t' + ','.join(genes))
			total_geneset.update(genes)


outfile.write('\n'.join(total_geneset))

