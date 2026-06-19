#!/bin/bash

# Scripts to run after parallel mapping scripts are completed

# Make trimming and mapping QC plots
Rscript src/checkQC.R >> logs/postMapping.log 2>&1

# combine samples for each tissue
Rscript src/mergeSamples.R >> logs/postMapping.log 2>&1

# Identify transcripts around ZFY
Rscript src/assembleTranscripts.R >> logs/postMapping.log 2>&1

# Plot splice variation
Rscript src/plotSashimi.R >> logs/postMapping.log 2>&1

# Read gene expression data from featureCounts, extract key genes TPM
Rscript src/readFeatureCounts.R >> logs/postMapping.log 2>&1

# Tar the figures
tar czf report.tar.gz report/*
tar czf stringtie.tar.gz data/stringtie/*
echo "Post mapping done" >> logs/postMapping.log 2>&1