#!/bin/bash
# Script used to quantify reads with Salmon
# Based on tutorial: https://combine-lab.github.io/salmon/getting_started/

for fn in raw_data/SRR105516{57..65};
do
samp=`basename ${fn}`
#echo "***Processing sample ${samp}***"
salmon quant -i scerevisiae_index -l A \
	-r ${fn}/${samp}.fastq.gz \
	--gcBias \
	-p 6 --validateMappings -o quants/${samp}
done
