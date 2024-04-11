#!/bin/bash

#clodius aggregate bedfile --chromsizes-filename chromsizes.hg38 colored_endo.bed
#clodius aggregate bedfile --chromsizes-filename chromsizes.hg38 colored_h1.bed
#clodius aggregate bedfile --chromsizes-filename chromsizes.hg38 colored_hff.bed
#clodius aggregate bedfile --chromsizes-filename chromsizes.hg38 colored_hela.bed


higlass-manage ingest --filetype beddb --assembly hg38  --datatype bedlike colored_endo.bed.beddb --name "test Endo CIS"
#higlass-manage ingest --filetype beddb --assembly hg38  --datatype bedlike colored_hela.bed.beddb --name "HeLa CIS"
#higlass-manage ingest --filetype beddb --assembly hg38  --datatype bedlike colored_h1.bed.beddb --name "H1 CIS"
#higlass-manage ingest --filetype beddb --assembly hg38  --datatype bedlike colored_hff.bed.beddb --name "HFF CIS"
#

#clodius aggregate bedfile --chromsizes-filename chromsizes.hg38 colored_hff.bed
#higlass-manage ingest --filetype beddb --assembly hg38  --datatype bedlike colored_hff.bed.beddb --name "HeLa Chromatin State"
