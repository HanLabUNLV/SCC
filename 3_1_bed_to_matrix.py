#!/usr/bin/env python

classDict = dict()
counter = 0

base = []
contact = []

for line in open('class_meta.txt'):
	classDict[line.strip()] = counter
	base.append(line.strip() + '_base')
	contact.append(line.strip() + '_contact')
	counter += 1

outfile = open("matrix_connections.tsv",'w')
#outfile.write("region\t" + '\t'.join(base) + '\t' + '\t'.join(contact) + '\n')
outfile.write("region\tbase_region\t" + '\t'.join(contact) + '\n')

bedfile = open('final_merged.bed','r')

#chr10_endo    73660    74260    TssFlnk    0    .    73660    74260    255,69,0    chr10_endo    73660    74260    0.09434047415,0.05729790147,0.05729790147,0.0148445316    chr10_endo,chr10_endo,chr10_endo,chr10_endo    75660,76260,76060,95060    76060,76660,76260,95260	ReprPC,TssFlnk,TssBiv,Het    600

#chr10_endo    510660    511060    EnhBiv    0    .    510660    511060    189,183,107    chr10_endo    510000    515000    0.01357866557,0.01357866557,0.01357866557,0.01357866557,0.01357866557,0.01357866557,0.02715733115    ReprPCWk,EnhA1,ReprPC,EnhBiv,TssBiv,TssFlnkU,EnhWk    400
#chr10_endo    511060    511460    EnhA2    0    .    511060    511460    255,195,77    chr10_endo    510000    515000    0.01357866557,0.01357866557,0.01357866557,0.01357866557,0.01357866557,0.01357866557,0.02715733115    ReprPCWk,EnhA1,ReprPC,EnhBiv,TssBiv,TssFlnkU,EnhWk    400
#chr10_endo    511460    511660    EnhWk    0    .    511460    511660    255,255,0    chr10_endo    510000    515000    0.01357866557,0.01357866557,0.01357866557,0.01357866557,0.01357866557,0.01357866557,0.02715733115    ReprPCWk,EnhA1,ReprPC,EnhBiv,TssBiv,TssFlnkU,EnhWk    200
#chr10_endo    785260    785460    EnhBiv    0    .    785260    785460    189,183,107    chr10_endo    785000    790000    0.01357866557,0.01357866557,0.01357866557,0.01357866557    EnhA2,EnhBiv,TxWk,EnhWk    200

#chr10_endo    73660    74260    TssFlnk    0    .    73660    74260    255,69,0    chr10_endo    73660    74260    0.05729790147,0.09434047415,0.02239986117,0.05729790147,0.0148445316    TssBiv,ReprPC,ReprPC,TssFlnk,Het    600
#chr10_endo    74260    75660    ReprPCWk    0    .    74260    75660    192,192,192    chr10_endo    74260    75660    0.1710774916,0.1917425648,0.2114261375,0.0148445316    TssFlnk,TssBiv,ReprPC,Het    1400


class_len = len(base)

for line in bedfile:

	contact_outlist = ['0'] * class_len
	base_outlist = ['0'] * class_len

	ll = line.strip().split()

	scores = ll[12].split(',')
	labs = ll[16].split(',')

	for i in range(len(scores)):
		contact_outlist[classDict[labs[i]]] = str(scores[i])

	# one-hot encoding base region
	#base_outlist[classDict[ll[3]]] = '1'
	#outfile.write(ll[0] + ':' + ll[1] + '-' + ll[2] + '\t' + '\t'.join(base_outlist) + '\t' + '\t'.join(contact_outlist) +'\n')

	outfile.write(ll[0] + ':' + ll[1] + '-' + ll[2] + '\t' + ll[3] + '\t' + '\t'.join(contact_outlist) +'\n')

		

bedfile.close()
outfile.close()
