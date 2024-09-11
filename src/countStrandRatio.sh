#!/bin/bash
# Calculate ratio of + to - reads

BAMFILE=$1

# Get the read counts from the file with strand info
FORWARD=$(samtools view --count --tag XS:+ $BAMFILE)
REVERSE=$(samtools view --count --tag XS:- $BAMFILE)

# Count all reads, whether they have strand tags or not
TOTAL=$(samtools view --count $BAMFILE)

echo "${BAMFILE} ${FORWARD} ${REVERSE} ${TOTAL}" >> data/stringtie/ratios.txt