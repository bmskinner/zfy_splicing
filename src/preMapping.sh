#!/bin/bash

# Scripts to be run before parallel mapping scripts are submitted

# Ensure all genome and annotations are present
bash src/makeIndexedGenomes.sh >> logs/makeIndexedGenomes.log 2>&1

# Select samples to map from metadata
Rscript src/selectSamples.R >> logs/selectSamples.log 2>&1