#!/usr/bin/env python

import sys

#https://stackoverflow.com/questions/61029281/check-if-two-python-sets-have-at-least-one-common-element
def intersect(a, b):
    if len(a) > len(b):
        a, b = b, a
    for c in a:
        if c in b:
            return True



candidate_genes = open(sys.argv[1],'r')

gene_set = set()
for line in candidate_genes:
	gene_set.add(line.strip())

with open("endo_enhancers.gff",'r') as infile:
	infile.readline()
	for line in infile:
		line_genes = set()

		ll = line.strip().split()
		li = ll[8].split(';')
		for i in li:
			if i.startswith("connected_gene="):
				line_genes.add(i.split('=')[1])
		
		if intersect(gene_set,line_genes):
			print(line.strip('\n') + '--:--' + ','.join(line_genes))
		
