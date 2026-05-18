#!/bin/bash

# Scripts to run after parallel mapping scripts are completed

# Make trimming and mapping QC plots
Rscript src/checkQC.R >> logs/checkQC.log 2>&1

# combine samples for each tissue
Rscript src/mergeSamples.R >> logs/mergeSamples.log 2>&1

# Identify transcripts around ZFY
Rscript src/assembleTranscripts.R >> logs/assembleTranscripts.log 2>&1

# Plot splice variation
Rscript src/plotSashimi.R >> logs/plotSashimi.log 2>&1

# Tar the figures
tar czf report.tar.gz report/*