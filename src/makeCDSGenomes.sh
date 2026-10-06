#!/bin/bash

mkdir -p cds/genomes
mkdir -p cds/fasta
mkdir -p cds/data

# Create indexed genomes for any mRNA or CDS ZFX/Y sequences
echo "`date '+%Y-%m-%d %X'` Making CDS genomes" >> logs/preMapping.log 2>&1
cd cds/genomes
for f in ../fasta/*.fa; do
  GENOME=$(echo $f | sed -e 's/.fa//')
  if [ ! -e ${GENOME}.1.ht2 ]; then
    hisat2-build ${f} ${GENOME}
  fi
done


