#!/bin/bash

raw_files=$(find . -name '*\.htseq\.counts')
for raw_file in $raw_files; do
    echo "processing ${raw_file}"
    processed_file=$(sed 's|raw.*|processed/gene_expression_profile.tsv|g' <<< $raw_file)
    echo "into ${processed_file}"
    lines_raw=$(wc -l $raw_file | awk '{print$1}')
    printf "gene_id\treads\n" > $processed_file
    cat $raw_file | grep 'ENSG' >> $processed_file
    lines_processed=$(wc -l $processed_file | awk '{print$1}')
    # in the processed file 5 lines are dropped
    # (starting respectively with __no_feature, __ambiguous, __too_low_aQual, __not_aligned, __alignment_not_unique)
    # and the header is added
    if (( $(echo "${lines_raw} != ${lines_processed} + 4" | bc -l) != 0 )); then
        echo "error: lines_raw = ${lines_raw}, lines_processed = ${lines_processed}"
        exit 1
    fi
done
echo "done"
