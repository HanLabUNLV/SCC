#!/usr/bin/env python

import sys
from collections import defaultdict
# take avg score of all windows in micro c for states that overlap multiple windows

score_dict = defaultdict(float)
curr_n_dict = defaultdict(int)

final_line_dict = defaultdict(str)

intersect_file = sys.argv[1]
output_file = sys.argv[2]

with open(intersect_file,'r') as infile:
	for line in infile:
		ll = line.strip().split()

		lab = '\t'.join(ll[3:6] + ll[7:9])

		score = float(ll[6])
		curr_n_dict[lab] += 1

		score_dict[lab] += (score - score_dict[lab]) / curr_n_dict[lab]
		final_line_dict[lab] = line

print("Finished loading Dict")

outfile = open(output_file,'w')

#chr1_endo    777820    778420    0.17871180991657024    chr1_endo    777820    778420    TssFlnk    420

for lab in final_line_dict:
	outline = final_line_dict[lab].strip().split()
	outfile.write('\t'.join(outline[3:6] + [str(score_dict[lab])] + outline[7:11]) + '\n')

outfile.close()
