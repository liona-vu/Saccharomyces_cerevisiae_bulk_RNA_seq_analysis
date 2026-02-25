#!/bin/bash
# Based on the tutorial below:
# https://combine-lab.github.io/salmon/getting_started/
# Run the script to download all S. cerevisiae data

mkdir raw_data
cd raw_data

for i in `seq 57 65`; 
do
  mkdir SRR105516${i};
  cd SRR105516${i};
  wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/0${i}/SRR105516${i}/SRR105516${i}.fastq.gz
  cd ..;
done

cd ..
