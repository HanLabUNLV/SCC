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
with open("endo_ghids.txt",'r') as infile:
	for line in infile:
		gh_ids.add(line.strip())

outfile = open("GeneHancer_v5.12_endoderm.bed",'w')

total_geneset = set()
with open("GeneHancer_v5.12.bed",'r') as infile:
	for line in infile:
		enh_id = line.strip().split()[3]
		if enh_id in gh_ids:
			outfile.write(line)

