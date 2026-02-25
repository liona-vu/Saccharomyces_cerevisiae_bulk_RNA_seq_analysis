#!/bin/bash
# Script used to quantify reads with Salmon
# Based on tutorial: https://combine-lab.github.io/salmon/getting_started/

# Build index with Salmon first
salmon index -t ncbi_dataset/ncbi_dataset/data/GCF_000146045.2/rna.fna -i scerevisiae_index

#then run this loop to quantify 
salmon index -t ncbi_dataset/ncbi_dataset/data/GCF_000146045.2/rna.fna -i scerevisiae_index
for fn in raw_data/SRR105516{57..65};
do
samp=`basename ${fn}`
#echo "***Processing sample ${samp}***"
salmon quant -i scerevisiae_index -l A \
	-r ${fn}/${samp}.fastq.gz \
	--gcBias \ #recommended by DESeq2 vignette
	-p 6 --validateMappings -o quants/${samp}
done


