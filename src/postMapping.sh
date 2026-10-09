#!/bin/bash

# If a mapping script is starting up, wait for lock to release
# Otherwise we may conflict with the metadata files being written
while [ -e data/preMapping.lck ]; do
		sleep 30
done

# Scripts to run after parallel mapping scripts are completed
# Create summary reports and get gene bam files for later visualisation
echo "`date '+%Y-%m-%d %X'` Beginning post mapping" > logs/postMapping.log 2>&1

# Details of the read depths and samples analysed
Rscript src/plotSampleSummary.R >> logs/postMapping.log 2>&1

# Make trimming and mapping QC plots
echo "`date '+%Y-%m-%d %X'` Checking QC" > logs/postMapping.log 2>&1
Rscript src/checkQC.R >> logs/postMapping.log 2>&1

# combine samples for each tissue and run StringTie
echo "`date '+%Y-%m-%d %X'` Merging samples" > logs/postMapping.log 2>&1
Rscript src/mergeSamples.R >> logs/postMapping.log 2>&1

# plot StringTie transcripts from merged samples
echo "`date '+%Y-%m-%d %X'` Plotting transcripts" > logs/postMapping.log 2>&1
Rscript src/plotNovelTranscripts.R >> logs/postMapping.log 2>&1

# Count strand ratios across samples
echo "`date '+%Y-%m-%d %X'` Counting strand ratios" >> logs/postMapping.log 2>&1
rm report/strand_ratios.txt
for f in data/*/[SDE]RR*.*.bam; do
  # Get the read counts from the file with strand info
  FORWARD=$(samtools view --count --tag XS:+ $f)
  REVERSE=$(samtools view --count --tag XS:- $f)
  
  # Count all reads, whether they have strand tags or not
  TOTAL=$(samtools view --count $f)

  echo "$f ${FORWARD} ${REVERSE} ${TOTAL}" >> report/strand_ratios.txt
done

# Plot splice variation
echo "`date '+%Y-%m-%d %X'` Plotting splicing" >> logs/postMapping.log 2>&1
Rscript src/plotSashimi.R >> logs/postMapping.log 2>&1

# Combine individual plots for easier comparison across tissues and species
echo "`date '+%Y-%m-%d %X'` Aggregating sashimi plots" >> logs/postMapping.log 2>&1
Rscript src/aggregateSashimiPlots.R >> logs/postMapping.log 2>&1

# Find mate pairs that suggest E2 is missing
echo "`date '+%Y-%m-%d %X'` Looking for spanning reads" >> logs/postMapping.log 2>&1
Rscript src/findE2SpanningPairs.R >> logs/postMapping.log 2>&1

# Read gene expression data from featureCounts, extract key genes and calc TPM
echo "`date '+%Y-%m-%d %X'` Calculating TPMs" >> logs/postMapping.log 2>&1
Rscript src/readFeatureCounts.R >> logs/postMapping.log 2>&1

# Tar the figures and output data
echo "`date '+%Y-%m-%d %X'` Making tarballs" >> logs/postMapping.log 2>&1
tar czf report.tar.gz report/*
tar czf reads.tar.gz data/*/*.*.bam* # only those with gene name included
tar czf featureCounts.tar.gz data/*/*.counts.txt

echo "`date '+%Y-%m-%d %X'` Post mapping done" >> logs/postMapping.log 2>&1