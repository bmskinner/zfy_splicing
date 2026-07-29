#!/bin/bash

# Scripts to run after parallel mapping scripts are completed
echo "Beginning post mapping" > logs/postMapping.log 2>&1

# Details of the read depths and samples analysed
Rscript src/plotSampleSummary.R >> logs/postMapping.log 2>&1

# Make trimming and mapping QC plots
Rscript src/checkQC.R >> logs/postMapping.log 2>&1

# combine samples for each tissue
# Rscript src/mergeSamples.R >> logs/postMapping.log 2>&1

# Count strand ratios across samples
rm report/strand_ratios.txt
for f in data/*/[SDR]RR*.*.bam; do
  src/countStrandRatio.sh ${f}
done


# Plot splice variation
Rscript src/plotSashimi.R >> logs/postMapping.log 2>&1

# Read gene expression data from featureCounts, extract key genes TPM
Rscript src/readFeatureCounts.R >> logs/postMapping.log 2>&1

# Tar the figures and output data
tar czf report.tar.gz report/*
tar czf stringtie.tar.gz data/stringtie/*
tar czf featureCounts.tar.gz data/*/*.counts.txt

echo "Post mapping done" >> logs/postMapping.log 2>&1