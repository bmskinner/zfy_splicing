# ZFY splicing analysis 

Data from RNA-seq analysis of whole tissues from

- PRJEB26695 Cardosa-Moreira 2019 (7 species across multiple tissues and timepoints)
- PRJEB33381 (adult tissues)

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

# check the trimming, fastqc and mapping
Rscript ./src/checkQC.R

# Combine bam files and extract splice sites
Rscript ./src/mergeSamples.R

# Plot splicing patterns
Rscript ./src/plotSashimi.R
```
