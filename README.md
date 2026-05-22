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
# Get genome and sample info
# Map to the latest genome assemblies
./src/mapSamples.sh

# Perform QC, merge bams, extract ZFX splicing and plot
./src/postMapping.sh

```
