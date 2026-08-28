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
for f in data/*/[SDE]RR*.*.bam; do
  # Get the read counts from the file with strand info
  FORWARD=$(samtools view --count --tag XS:+ $f)
  REVERSE=$(samtools view --count --tag XS:- $f)
  
  # Count all reads, whether they have strand tags or not
  TOTAL=$(samtools view --count $f)

  echo "$f ${FORWARD} ${REVERSE} ${TOTAL}" >> report/strand_ratios.txt
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