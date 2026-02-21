#!/bin/bash
# Check S. cerevisiae data with FASTQC

for i in raw_data/SRR1055{57..65};
do
  fastqc {i}.fastq.gz -o fastqc_stats/
  cd ..;
done
