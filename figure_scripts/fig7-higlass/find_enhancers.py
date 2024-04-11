#!/usr/bin/env python

import sys

#h1_enhancer_list = open("h1_enhancers.txt",'r')
endo_enhancer_list = open("endo_enhancers.txt",'r')

h1_enhancers = set()
endo_enhancers = set()

for line in h1_enhancer_list:
	h1_enhancers.add(line.strip())

for line in endo_enhancer_list:
	endo_enhancers.add(line.strip())

out_h1_enhancer_list = open("h1_enhancers.gff",'w')
out_endo_enhancer_list = open("endo_enhancers.gff",'w')

with open("GeneHancer_v5.12.gff",'r') as infile:
	infile.readline()
	for line in infile:
		enh = line[line.index('=GH'):line.index('=GH') + 12]
		
		if enh in h1_enhancers:
			out_h1_enhancer_list.write(line.strip() + '--:--' + enh[1:] + '\n')
		if enh in endo_enhancers:
			out_endo_enhancer_list.write(line.strip() + '--:--' + enh[1:] + '\n')
			

out_endo_enhancer_list.close()
out_h1_enhancer_list.close()
h1_enhancer_list.close()
endo_enhancer_list.close()


