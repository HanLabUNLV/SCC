#!/usr/bin/env python


import sys

#  h1_tss_overlap.bed

#chr1	778438	778938	LOC100288069	0	-	chr1	777820	779020	TssFlnk	Tss_noEnh_cont	h1	500
#chr1	827272	827772	LINC00115	0	-	chr1	826620	827820	TssFlnk	Low_cont	h1	500
#chr1	827340	827840	LINC01128	0	+	chr1	826620	827820	TssFlnk	Low_cont	h1	480
#chr1	876552	877052	FAM41C	0	-	chr1	870620	877420	ReprPCWk	EnhWk_cont1	h1	500
#chr1	919442	919942	LOC100130417	0	-	chr1	917620	920420	ReprPC	Biv_cont	h1	500
#chr1	925490	925990	SAMD11	0	+	chr1	923420	926220	TssBiv	Biv_cont	h1	500
#chr1	959049	959549	NOC2L	0	-	chr1	958420	959420	TssFlnk	Biv_cont	h1	371



for f in ["h1_tss_overlap.bed","hff_tss_overlap.bed","hela_tss_overlap.bed","endo_tss_overlap.bed"]:
	outfile_ccs_name = f[:-3] + 'ccs_500bp.bed'
	outfile_chrom_name = f[:-3] + 'chrom_500bp.bed'
	print(outfile_ccs_name)
	print(outfile_chrom_name)
	with open(f,'r') as infile:
		outfile_ccs = open(outfile_ccs_name,'w')
		outfile_chrom = open(outfile_chrom_name,'w')
		ccs_labs = []
		chrom_labs = []
		bases = []
		prev_gene = ''
		for line in infile:	
			ll = line.strip().split()

			gene = ll[3]
			if prev_gene == '':  # First Time
				prev_gene = gene

			print(bases)
			if prev_gene != gene:
				if sum(bases) != 500:
					if strand == '+':
						bases.append(500 - sum(bases))
						ccs_labs.append('NA')
						chrom_labs.append('NA')
					else:
						bases.insert(0,500 - sum(bases))
						ccs_labs.insert(0,'NA')
						chrom_labs.insert(0,'NA')
			
					
				outfile_ccs.write(prev_gene + '\t')
				outfile_chrom.write(prev_gene + '\t')
				out_ccs_string = ''
				out_chrom_string = ''
				for i in range(len(bases)):
					out_ccs_string += '\t'.join([ccs_labs[i]] * bases[i]) + '\t'
					out_chrom_string += '\t'.join([chrom_labs[i]] * bases[i]) + '\t'
				outfile_ccs.write(out_ccs_string.strip() + '\n')
				outfile_chrom.write(out_chrom_string.strip() + '\n')

				ccs_labs.clear()
				chrom_labs.clear()
				bases.clear()
				print("printed " + prev_gene)
			else:
				print(gene)

			strand = ll[5]
			cell_type = ll[-2]
			ccs_label = ll[-3]
			chrom_label = ll[-4]
			
			if strand == '+':
				bases.append(int(ll[-1]))
				ccs_labs.append(ccs_label)
				chrom_labs.append(chrom_label)
			else:
				bases.insert(0,int(ll[-1]))
				ccs_labs.insert(0,ccs_label)
				chrom_labs.insert(0,chrom_label)

			
			prev_gene = gene

		# LAST GENE
		if sum(bases) != 500:
			if strand == '+':
				bases.append(500 - sum(bases))
				ccs_labs.append('NA')
				chrom_labs.append('NA')
			else:
				bases.insert(0,500 - sum(bases))
				ccs_labs.insert(0,'NA')
				chrom_labs.insert(0,'NA')
	
			
			outfile_ccs.write(prev_gene + '\t')
			outfile_chrom.write(prev_gene + '\t')
			out_ccs_string = ''
			out_chrom_string = ''
			for i in range(len(bases)):
				out_ccs_string += '\t'.join([ccs_labs[i]] * bases[i]) + '\t'
				out_chrom_string += '\t'.join([chrom_labs[i]] * bases[i]) + '\t'
			outfile_ccs.write(out_ccs_string.strip() + '\n')
			outfile_chrom.write(out_chrom_string.strip() + '\n')

			ccs_labs.clear()
			chrom_labs.clear()
			bases.clear()
			print("printed " + prev_gene)


		outfile_ccs.close()
		outfile_chrom.close()



