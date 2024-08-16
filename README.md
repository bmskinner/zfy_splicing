# ZFY splicing analysis 

Data from Cardosa-Moreira 2019 for 7 species across multiple tissues and timepoints.

Here we remap all samples to the latest genome assemblies, allowing multiple read mapping (since Zfx/Zfy are very similar).

We exclude rabbit, since there is no Y assembly.

Individual sample bam files are merged by species and tissue type.

Splice junctions for Zfx and Zfy are extracted from the merged bam files.

## How to run

```{bash}
# Assign timepoints and select samples for analysis
Rscript ./src/selectSamples.R

# Build genome indexes
./src/makeIndexedGenomes.sh

# Map to the latest genome assemblies
./src/mapSamples.sh

# Combine bam files and extract splice sites
Rscript ./src/mergeSamples.R

# Plot splicing patterns
Rscript ./src/plotSashimi.R
```
