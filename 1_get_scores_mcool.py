#!/usr/bin/env python

## https://cooltools.readthedocs.io/en/latest/notebooks/viz.html#Inspecting-C-data

import numpy as np
import pandas as pd
import cooltools
import cooler
import gzip

import sys

usage = '''Usage: python 1_get_scores_mcool.py <mcool file with res> <ChromHMM bedfile> <outfile> [margin size]'''

# get resolutions
#cooler.fileops.list_coolers('./4DNFI9GMP2J8.mcool')

if len(sys.argv) < 4:
	print(usage)
	sys.exit(1)
else:
	mcool_file = sys.argv[1]
	chromHMM_path = sys.argv[2]
	outfile_path = sys.argv[3]
	if len(sys.argv) == 5:
		margin_size = int(sys.argv[4])
	else: margin_size = 2_000_000

clr = cooler.Cooler(mcool_file)

chromHMM_file = gzip.open(chromHMM_path, 'rt')
outfile = gzip.open(outfile_path,'wt')

chrom_size = clr.chromsizes

for line in chromHMM_file:
	ll = line.strip().split()

	region1 = (ll[0],ll[1],ll[2])

	start = int(ll[1]) - margin_size
	end = int(ll[2]) + margin_size

	if start < 0: start = 0
	if end > chrom_size[ll[0]]: end = chrom_size[ll[0]]
	region2 = (ll[0],start,end)

	normalized_score_mat = clr.matrix(balance=True).fetch(region1,region2)
	bins_1 = clr.bins().fetch((region1))
	bins_2 = clr.bins().fetch((region2))

	chroms= list(bins_2['chrom'])
	pos1= list(bins_2['start'])
	pos2= list(bins_2['end'])

	row,column = normalized_score_mat.shape
	for r in range(row):
		for c in range(column):
			z = normalized_score_mat[r,c]
			if np.isnan(z) or z == 0:
				continue
			else:
				outfile.write('\t'.join([chroms[c],str(pos1[c]),str(pos2[c]),ll[0],ll[1],ll[2],str(z),ll[3]]) + '\n')
		
outfile.close()

