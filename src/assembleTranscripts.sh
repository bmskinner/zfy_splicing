#!/bin/bash

# Assemble transcripts on the sex chromosomes using StringTie
# StringTie expected on the PATH

# - x ignore these chromosomes
# -C output a file with the given name with all transcripts in the provided reference file that are fully covered by reads
# -G reference annotations
# -f minimum isoform abundance of the predicted transcripts as a fraction of the most abundant transcript assembled at a given locus
# -l label prefix to give the transcript names
# -o output GTF of assembled transcripts
# stringtie [-o <output.gtf>] -G <ref_ann.gff> -l <label> -f 0.01  -C <cov_refs.gtf>  -x <seqid_list> <read_alignments.bam>

# Use StringTie to assemble the novel chicken ZFX transcript


# Extract just the region containing chicken ZFX
# 1:118200000-118400000

mkdir -p data/stringtie

samtools view -o data/stringtie/chicken.testis.zfx.reads.bam data/merged/chicken.testis.adult.bam "1:118296000-118319000"

~/bin/stringtie-2.2.3.Linux_x86_64/stringtie -o data/stringtie/chicken.testis.zfx.gtf -p 1 -l chicken -f 0.01 data/stringtie/chicken.testis.zfx.reads.bam