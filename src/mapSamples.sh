#!/bin/bash

# Remap samples and generate bam files for splice junction detection

CHICKEN_SAMPLES=$(cat metadata/chicken.filt.csv | cut -f 1 -d , | tail -n +2)
OPOSSUM_SAMPLES=$(cat metadata/opossum.filt.csv | cut -f 1 -d , | tail -n +2)
MOUSE_SAMPLES=$(cat metadata/mouse.filt.csv | cut -f 1 -d , | tail -n +2)
HUMAN_SAMPLES=$(cat metadata/human.filt.csv | cut -f 1 -d , | tail -n +2)
# RABBIT_SAMPLES=$(cat metadata/rabbit.filt.csv | cut -f 1 -d , | tail -n +2)
RAT_SAMPLES=$(cat metadata/rat.filt.csv | cut -f 1 -d , | tail -n +2)
MACAQUE_SAMPLES=$(cat metadata/macaque.filt.csv | cut -f 1 -d , | tail -n +2)

# $1 species e.g chicken - should match a folder name in ./data
# $2 ERR id e.g ERR2576379
# $3 genome build id e.g. GRCg7b - should match a .ht2 prefix in ./genomes
map_sample () {
	SPECIES=$(echo $1 | tr -d '"')
	ERR=$(echo $2 | tr -d '"')
	GENOME=$(echo $3 | tr -d '"')

	mkdir -p data/${SPECIES}
	
	# Create if missing
	if [ ! -e data/${SPECIES}/${ERR}.bam ]; then

		# Don't work on a sample already being processed
		if [ ! -e data/${SPECIES}/${ERR}.lck ]; then

			touch data/${SPECIES}/${ERR}.lck
			echo "${ERR}: bam not found"

			# Check for partial downloads before running fasterq-dump
			if [ ! -e data/${SPECIES}/${ERR}.fastq.gz ]; then
				# Fetch data
				if [ ! -e data/${SPECIES}/${ERR}.fastq ]; then
					echo -n "${ERR}: downloading fastq"
					fasterq-dump -o data/${SPECIES}/${ERR}.fastq ${ERR}
				fi
				gzip data/${SPECIES}/${ERR}.fastq
			fi

			# Trim
			if [ ! -e data/${SPECIES}/${ERR}_trimmed.fq.gz ]; then
				echo "${ERR}: trimming"
				trim_galore -o data/${SPECIES} --suppress_warn --fastqc --fastqc_args "-t 8 --outdir report --nogroup --extract" data/${SPECIES}/${ERR}.fastq.gz
			fi

			# map
			echo "${ERR}: mapping"
			# -k controls number of multimapping locations (default 5 for linear index)
			# --downstream-transcriptome-assembly forces longer anchors at novel splice sites (more rigorous)
			# Analysis run with and without --dta to compare effects; we don't need to assemble transcripts, just see if there is greater splicing in testis
			# hisat2 -x genomes/${GENOME} -p 8 -U data/${SPECIES}/${ERR}_trimmed.fq.gz --new-summary --downstream-transcriptome-assembly -S data/${SPECIES}/${ERR}.sam
			hisat2 -x genomes/${GENOME} -p 8 -U data/${SPECIES}/${ERR}_trimmed.fq.gz --new-summary -S data/${SPECIES}/${ERR}.sam
			samtools sort -T data/${SPECIES}/${ERR} -@ 8 -o data/${SPECIES}/${ERR}.bam data/${SPECIES}/${ERR}.sam
			samtools index -c -@ 7 data/${SPECIES}/${ERR}.bam # index with csi due to long chromosomes in opossum
			rm data/${SPECIES}/${ERR}.sam

			# Remove lock file
			rm data/${SPECIES}/${ERR}.lck
		fi
	fi
}

for ERR in ${CHICKEN_SAMPLES}; do
	map_sample chicken ${ERR} GRCg7b
done

for ERR in ${OPOSSUM_SAMPLES}; do
	map_sample opossum ${ERR} ASM229v1
done

for ERR in ${MOUSE_SAMPLES}; do
	map_sample mouse ${ERR} GRCm39
done

for ERR in ${HUMAN_SAMPLES}; do
	map_sample human ${ERR} GRCh38
done

# for ERR in ${RABBIT_SAMPLES}; do
# 	map_sample rabbit ${ERR} OryCun2.0
# done

for ERR in ${RAT_SAMPLES}; do
	map_sample rat ${ERR} mRatBN7.2
done

for ERR in ${MACAQUE_SAMPLES}; do
	map_sample macaque ${ERR} Mmul_10
done
