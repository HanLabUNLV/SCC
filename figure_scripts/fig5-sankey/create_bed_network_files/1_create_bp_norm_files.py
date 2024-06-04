#!/usr/bin/env python


#  ../create_bed_network_files/h1_labelled.bed
#chr10	48460	51260	Het	Low_cont	h1
#chr10	73460	73660	TssFlnk	Low_cont	h1
#chr10	73660	74060	TssBiv	Tss_noEnh_cont	h1
#chr10	74060	74260	TssFlnk	Low_cont	h1
#chr10	74260	75660	ReprPCWk	Tss_noEnh_cont	h1
#chr10	75660	76060	ReprPC	Biv_cont	h1
#chr10	76060	76660	TssBiv	Low_cont	h1

# OLD
#category = {
#   "Low_cont":6,
#   "EnhWk_cont1":6,
#   "Het_cont1":5,
#   "Het_cont2":5,
#   "Biv_cont":4,
#   "ReprPC_cont":4,
#   "EnhWk_cont2":3,
#   "Enh_cont":2,
#   "EnhA1_cont":2,
#   "Tss_noEnh_cont":1,"Tss_Enh_cont":1,"TssFlnkD_cont":1,"TssFlnkU_cont":1,
#   "Tx_cont1":0,"Tx_cont2":0
# }

for f in ["../create_bed_network_files/h1_labelled.bed","../create_bed_network_files/hff_labelled.bed","../create_bed_network_files/hela_labelled.bed","../create_bed_network_files/endo_labelled.bed"]:
	outfile = open(f[:-4] + "_200bp_window.bed",'w')
	with open(f,'r') as infile:
		for line in infile:
			ll = line.strip().split()
			pos1 = int(ll[1])
			pos2 = int(ll[2])
	
			if pos2-pos1 != 200:
				num_times = (pos2-pos1) / 200
				first_pos = pos1
				for i in range(int(num_times)):
					outfile.write(ll[0] + '\t' + str(first_pos) + '\t' + str(first_pos + 200) + '\t' + '\t'.join(ll[3:]) + '\n')
					first_pos += 200
			else:
				outfile.write(line)
	outfile.close()


#cut -f 1,2,3 *200bp_window.bed | sort | uniq > all_regions.txt
