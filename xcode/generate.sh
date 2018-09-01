#!/bin/bash
# generate an Xcode project into the folder _build
cmake -G Xcode -H.. -B_build
mkdir -p _build/Debug
mkdir -p _build/Debug/data
mkdir -p _build/Debug/data/processed
rsync -avh ../data/processed/mirnas_with_scored_interactions.tsv _build/Debug/data/processed/mirnas_with_scored_interactions.tsv
rsync -avh ../data/processed/sites_with_scored_interactions.tsv _build/Debug/data/processed/sites_with_scored_interactions.tsv
rsync -avh ../data/processed/scored_interactions_processed.tsv _build/Debug/data/processed/scored_interactions_processed.tsv
rsync -avh ../data/patients _build/Debug/data/
