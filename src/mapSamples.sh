#!/bin/bash

mkdir -p report/FASTQC

# Remap samples and generate bam files for splice junction detection
# Select the Run, Library type, species and genome columns
SE_SAMPLES=$(cat metadata/*.filt.csv | cut -f 1,3,4,5 -d , | grep -e '[S|E|D]RR' | grep -e 'SINGLE')
PE_SAMPLES=$(cat metadata/*.filt.csv | cut -f 1,3,4,5 -d , | grep -e '[S|E|D]RR' | grep -e 'PAIRED')

# Map a single end sample
# $1 species e.g chicken - should match a folder name in ./data
# $2 ERR id e.g ERR2576379
# $3 genome build id e.g. GRCg7b - should match a .ht2 prefix in ./genomes
map_se_sample () {
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

			# Check for existing downloads before running fasterq-dump
			if [ ! -e data/${SPECIES}/${ERR}.fastq.gz ] && [ ! -e data/${SPECIES}/${ERR}_trimmed.fq.gz ]; then
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
				trim_galore -o data/${SPECIES} --suppress_warn --fastqc --fastqc_args "-t 8 --outdir report/FASTQC --nogroup --extract" data/${SPECIES}/${ERR}.fastq.gz
			fi

			# map
			echo "${ERR}: mapping"
			# -k controls number of multimapping locations (default 5 for linear index)
			# --downstream-transcriptome-assembly forces longer anchors at novel splice sites (more rigorous)
			# --dta-cfflinks does this and also looks for novel splice sites, stored in tag XS:A:[+-]
			# Analysis run with and without --dta to compare effects; we don't need to assemble transcripts, just see if there is greater splicing in testis
			# hisat2 -x genomes/${GENOME} -p 8 -U data/${SPECIES}/${ERR}_trimmed.fq.gz --new-summary --downstream-transcriptome-assembly -S data/${SPECIES}/${ERR}.sam
			hisat2 -x genomes/${GENOME} -p 8 -U data/${SPECIES}/${ERR}_trimmed.fq.gz --new-summary --dta-cufflinks -S data/${SPECIES}/${ERR}.sam
			samtools sort -T data/${SPECIES}/${ERR} -@ 8 -o data/${SPECIES}/${ERR}.bam data/${SPECIES}/${ERR}.sam
			samtools index -c -@ 7 data/${SPECIES}/${ERR}.bam # index with csi due to long chromosomes in opossum
			rm data/${SPECIES}/${ERR}.sam
			
			# Remove original FASTQ, we have the trimmed reads still
			if [ -e data/${SPECIES}/${ERR}.fastq.gz ]; then
			  rm data/${SPECIES}/${ERR}.fastq.gz
			fi

			# Remove lock file
			rm data/${SPECIES}/${ERR}.lck
		fi
	fi
}


# Map a paired end sample
# $1 species e.g chicken - should match a folder name in ./data
# $2 ERR id e.g ERR2576379
# $3 genome build id e.g. GRCg7b - should match a .ht2 prefix in ./genomes
map_pe_sample () {
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

			# Check for existing downloads before running fasterq-dump
			if [ ! -e data/${SPECIES}/${ERR}_2.fastq.gz ] && [ ! -e data/${SPECIES}/${ERR}_2_val_2.fq.gz ]; then
				# Fetch data
				if [ ! -e data/${SPECIES}/${ERR}_2.fastq ]; then
					echo -n "${ERR}: downloading fastq"
					fasterq-dump -O data/${SPECIES}/ ${ERR}
				fi
				gzip data/${SPECIES}/${ERR}_1.fastq
				gzip data/${SPECIES}/${ERR}_2.fastq
			fi

			# Trim
			if [ ! -e data/${SPECIES}/${ERR}_1_val_1.fq.gz ]; then
				echo "${ERR}: trimming"
				trim_galore -o data/${SPECIES} --paired --suppress_warn --fastqc --fastqc_args "-t 8 --outdir report/FASTQC --nogroup --extract" data/${SPECIES}/${ERR}_1.fastq.gz data/${SPECIES}/${ERR}_2.fastq.gz
			fi

			# map
			echo "${ERR}: mapping"
			# -k controls number of multimapping locations (default 5 for linear index)
			hisat2 -x genomes/${GENOME} -p 8 -1 data/${SPECIES}/${ERR}_1_val_1.fq.gz -2 data/${SPECIES}/${ERR}_2_val_2.fq.gz --new-summary -S data/${SPECIES}/${ERR}.sam
			samtools sort -T data/${SPECIES}/${ERR} -@ 8 -o data/${SPECIES}/${ERR}.bam data/${SPECIES}/${ERR}.sam
			samtools index -c -@ 7 data/${SPECIES}/${ERR}.bam # index with csi due to long chromosomes in opossum
			rm data/${SPECIES}/${ERR}.sam
			
			# Remove original FASTQ, we have the trimmed reads still
			if [ -e data/${SPECIES}/${ERR}_2.fastq.gz ]; then
			  rm data/${SPECIES}/${ERR}_*.fastq.gz
			fi

			# Remove lock file
			rm data/${SPECIES}/${ERR}.lck
		fi
	fi
}

echo "Processing single end samples"
for LINE in ${SE_SAMPLES}; do

	ERR=$(echo ${LINE} | cut -f 1 -d , )
	SPECIES=$(echo ${LINE} | cut -f 3 -d , )
	GENOME=$(echo ${LINE} | cut -f 4 -d , )
	map_se_sample ${SPECIES} ${ERR} ${GENOME}
done

echo "Processing paired end samples"
for LINE in ${PE_SAMPLES}; do

	ERR=$(echo ${LINE} | cut -f 1 -d , )
	SPECIES=$(echo ${LINE} | cut -f 3 -d , )
	GENOME=$(echo ${LINE} | cut -f 4 -d , )
	map_pe_sample ${SPECIES} ${ERR} ${GENOME}
done
