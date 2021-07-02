#!/bin/bash

set -e

EXP_DIR="/media/large/volume3/processed_RNAseq_s/2020-12-14"
SAMPLES_LIST="/home/mraevsky/LTS_vs_STS_manuscript/Slovenia_RNAseq_samples.txt"
DOWNLOAD_DIR="/home/mraevsky/LTS_vs_STS_manuscript/exp_data"
mkdir -p $DOWNLOAD_DIR

while IFS="" read -r sample; do
    echo "$(date) - Processing $sample"
    find $EXP_DIR -type f -iname "${sample}*" 2>/dev/null |
        xargs -I{} cp {} $DOWNLOAD_DIR
done <"$SAMPLES_LIST"

## From local machine
# scp -r oncobox:/home/mraevsky/LTS_vs_STS_manuscript/exp_data .

rm -r $DOWNLOAD_DIR


