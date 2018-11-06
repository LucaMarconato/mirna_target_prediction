#!/bin/bash

if [ ! -f ./data/sorted.bam ]; then
    samtools sort -@ 4 -o ./data/sorted.bam ./data/GSM759888_hela.bam
fi
if [ ! -f ./data/counts.txt ]; then
    featureCounts -t exon -g gene_id -a ./data/Homo_sapiens.GRCh37.64.gtf -o ./data/counts.txt ./data/sorted.bam
fi
if [ ! -f ./data/gene_expression_profile.tsv ]; then
    printf "gene_id\treads\n" > ./data/gene_expression_profile.tsv
    tail -n +3 ./data/counts.txt | cut -f 1,6 | tr " " "\t" >> ./data/gene_expression_profile.tsv
fi
