#!/bin/bash
# Calculate ratio of + to - reads

BAMFILE=$1

FORWARD=$(samtools view --count --tag XS:+ $BAMFILE)
REVERSE=$(samtools view --count --tag XS:- $BAMFILE)

TOTAL=$(echo " ${FORWARD} + ${REVERSE} " | bc)
FRACTION=$(echo "scale=2; ${FORWARD} / ${TOTAL}" | bc)

echo "${BAMFILE} ${FORWARD} ${REVERSE} ${TOTAL} ${FRACTION}" >> data/stringtie/ratios.txt