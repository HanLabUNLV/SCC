#!/usr/bin/env python3

import sys

#chr1	GeneHancer	Promoter/Enhancer	777600	781201	1.54	.	.	genehancer_id=GH01J000777;connected_gene=LINC01409;score=278.84;connected_gene=ENSG00000228327;score=250.67;connected_gene=LOC100288069;score=250.67;connected_gene=piR-60253-001;score=250.67;connected_gene=SLC35E2A;score=26.93;connected_gene=INTS11;score=24.52;connected_gene=ENSG00000227775;score=23.31;connected_gene=WASH7P;score=23.10;connected_gene=ENSG00000240731;score=21.86;connected_gene=ATAD3B;score=19.00;connected_gene=ENSG00000260179;score=14.40;connected_gene=LINC00115;score=13.89;connected_gene=CCNL2;score=13.71;connected_gene=CDK11A;score=13.39;connected_gene=ENSG00000269737;score=12.47;connected_gene=SLC35E2B;score=12.44;connected_gene=KLHL17;score=11.97;connected_gene=ENSG00000230092;score=11.51;connected_gene=CICP3;score=11.11;connected_gene=NOC2L;score=11.09;connected_gene=PERM1;score=10.88;connected_gene=UBE2J2;score=10.10;connected_gene=PUSL1;score=9.90;connected_gene=ENSG00000285268;score=0.67;connected_gene=lnc-SAMD11-8;score=0.41;connected_gene=LOC101928626;score=0.26;connected_gene=OR4F16;score=0.23;connected_gene=SAMD11;score=0.11--:--GH01J000777

genes = set()
with open(sys.argv[1],'r') as infile:
	for line in infile:
		ll = line.strip().split()
		li = ll[8].split(';')
		for i in li:
			if i.startswith("connected_gene="):
				genes.add(i.split('=')[1])
		
for g in genes:
	print(g)

