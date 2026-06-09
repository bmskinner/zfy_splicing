#!/bin/bash

mkdir -p report/FASTQC

# Remap samples and generate bam files for splice junction detection
# Select the Run, Library type, species and genome columns
SE_SAMPLES=$(cat metadata/*.filt.csv | cut -f 1,3,4,5,6 -d , | grep -e '[S|E|D]RR' | grep -e 'SINGLE')
PE_SAMPLES=$(cat metadata/*.filt.csv | cut -f 1,3,4,5,6 -d , | grep -e '[S|E|D]RR' | grep -e 'PAIRED')

# Map a single end sample
# $1 species e.g chicken - should match a folder name in ./data
# $2 ERR id e.g ERR2576379
# $3 genome build id e.g. GRCg7b - should match a .ht2 prefix in ./genomes
# $4 gtf file corresponding to the genome build e.g. ./genomes/example.gtf
map_se_sample () {
	SPECIES=$(echo $1 | tr -d '"')
	ERR=$(echo $2 | tr -d '"')
	GENOME=$(echo $3 | tr -d '"')
	GTF_FILE=$(echo $4 | tr -d '"')
	
	echo "${ERR}: beginning mapping to ${SPECIES} against ${GENOME} and ${GTF_FILE}" >> logs/mapping.log 2>&1

	mkdir -p data/${SPECIES}
	
	if [ ! -e ${GTF_FILE} ]; then
	  echo "${ERR}: Could not find GTF file ${GTF_FILE}, skipping" >> logs/mapping.log 2>&1
	  continue
	fi
	
	# Create if missing
	if [ ! -e data/${SPECIES}/${ERR}.bam ]; then

		# Don't work on a sample already being processed
		if [ ! -e data/${SPECIES}/${ERR}.lck ]; then

			touch data/${SPECIES}/${ERR}.lck
			echo "${ERR}: bam not found" >> logs/${ERR}.mapping.log 2>&1

			# Check for existing downloads before running fasterq-dump
			if [ ! -e data/${SPECIES}/${ERR}.fastq.gz ] && [ ! -e data/${SPECIES}/${ERR}_trimmed.fq.gz ]; then
				# Fetch data
				if [ ! -e data/${SPECIES}/${ERR}.fastq ]; then
					echo -n "${ERR}: downloading fastq" >> logs/${ERR}.mapping.log 2>&1
					fasterq-dump -o data/${SPECIES}/${ERR}.fastq ${ERR}
				fi
				gzip data/${SPECIES}/${ERR}.fastq
			fi

			# Trim
			if [ ! -e data/${SPECIES}/${ERR}_trimmed.fq.gz ]; then
				echo "${ERR}: trimming" >> logs/${ERR}.mapping.log 2>&1
				trim_galore -o data/${SPECIES} --suppress_warn --fastqc --fastqc_args "-t 8 --outdir report/FASTQC --nogroup --extract" data/${SPECIES}/${ERR}.fastq.gz
			fi

			# map
			echo "${ERR}: mapping" >> logs/${ERR}.mapping.log 2>&1
			# -k controls number of multimapping locations (default 5 for linear index)
			# --downstream-transcriptome-assembly forces longer anchors at novel splice sites (more rigorous)
			# --dta-cfflinks does this and also looks for novel splice sites, stored in tag XS:A:[+-]
			# Note that a stranded library may be with respect to forward or reverse strand depending on prep method;
			# Use --rna-strandness R to specify single-end RNA-seq data is reverse stranded, F for forward strand.
			# Analysis run with and without --dta to compare effects; we don't need to assemble transcripts, just see if there is greater splicing in testis
			# hisat2 -x genomes/${GENOME} -p 8 -U data/${SPECIES}/${ERR}_trimmed.fq.gz --new-summary --downstream-transcriptome-assembly -S data/${SPECIES}/${ERR}.sam
			hisat2 -x genomes/${GENOME} -p 8 -U data/${SPECIES}/${ERR}_trimmed.fq.gz --new-summary --dta-cufflinks -S data/${SPECIES}/${ERR}.sam >> logs/${ERR}.mapping.log 2>&1
			samtools sort -T data/${SPECIES}/${ERR} -@ 8 -o data/${SPECIES}/${ERR}.bam data/${SPECIES}/${ERR}.sam
			samtools index -c -@ 7 data/${SPECIES}/${ERR}.bam # index with csi due to long chromosomes in opossum
			rm data/${SPECIES}/${ERR}.sam
			
			# Run feature counts for expression quantification
			if [ ! -e data/${SPECIES}/${ERR}.counts.txt ]; then
		  	featureCounts -t exon -g gene_id -a ${GTF_FILE} -o data/${SPECIES}/${ERR}.counts.txt data/${SPECIES}/${ERR}.bam
		  fi
			
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
# $4 gtf file corresponding to the genome build e.g. ./genomes/example.gtf
map_pe_sample () {
	SPECIES=$(echo $1 | tr -d '"')
	ERR=$(echo $2 | tr -d '"')
	GENOME=$(echo $3 | tr -d '"')
	GTF_FILE=$(echo $4 | tr -d '"')
	
	echo "${ERR}: beginning mapping to ${SPECIES} against ${GENOME} and ${GTF_FILE}" >> logs/mapping.log 2>&1

	mkdir -p data/${SPECIES}
	
	if [ ! -e ${GTF_FILE} ]; then
	  echo "${ERR}: Could not find GTF file ${GTF_FILE}, skipping" >> logs/mapping.log 2>&1
	  continue
	fi
	
	# Create if missing
	if [ ! -e data/${SPECIES}/${ERR}.bam ]; then

		# Don't work on a sample already being processed
		if [ ! -e data/${SPECIES}/${ERR}.lck ]; then

			touch data/${SPECIES}/${ERR}.lck
			echo "${ERR}: bam not found" >> logs/${ERR}.mapping.log 2>&1

			# Check for existing downloads before running fasterq-dump
			if [ ! -e data/${SPECIES}/${ERR}_2.fastq.gz ] && [ ! -e data/${SPECIES}/${ERR}_2_val_2.fq.gz ]; then
				# Fetch data
				if [ ! -e data/${SPECIES}/${ERR}_2.fastq ]; then
					echo -n "${ERR}: downloading fastq" >> logs/${ERR}.mapping.log 2>&1
					fasterq-dump -O data/${SPECIES}/ ${ERR}
				fi
				gzip data/${SPECIES}/${ERR}_1.fastq
				gzip data/${SPECIES}/${ERR}_2.fastq
			fi

			# Trim
			if [ ! -e data/${SPECIES}/${ERR}_1_val_1.fq.gz ]; then
				echo "${ERR}: trimming" >> logs/${ERR}.mapping.log 2>&1
				trim_galore -o data/${SPECIES} --paired --suppress_warn --fastqc --fastqc_args "-t 8 --outdir report/FASTQC --nogroup --extract" data/${SPECIES}/${ERR}_1.fastq.gz data/${SPECIES}/${ERR}_2.fastq.gz
			fi

			# map
			echo "${ERR}: mapping" >> logs/${ERR}.mapping.log 2>&1
			# -k controls number of multimapping locations (default 5 for linear index)
			hisat2 -x genomes/${GENOME} -p 8 -1 data/${SPECIES}/${ERR}_1_val_1.fq.gz -2 data/${SPECIES}/${ERR}_2_val_2.fq.gz --new-summary -S data/${SPECIES}/${ERR}.sam >> logs/${ERR}.mapping.log 2>&1
			samtools sort -T data/${SPECIES}/${ERR} -@ 8 -o data/${SPECIES}/${ERR}.bam data/${SPECIES}/${ERR}.sam
			samtools index -c -@ 7 data/${SPECIES}/${ERR}.bam # index with csi due to long chromosomes in opossum
			rm data/${SPECIES}/${ERR}.sam
			
			# Run feature counts for expression quantification
			if [ ! -e data/${SPECIES}/${ERR}.counts.txt ]; then
			  featureCounts -p --countReadPairs -t exon -g gene_id -a ${GTF_FILE} -o data/${SPECIES}/${ERR}.counts.txt data/${SPECIES}/${ERR}.bam
			fi
			
			# Remove original FASTQ, we have the trimmed reads still
			if [ -e data/${SPECIES}/${ERR}_2.fastq.gz ]; then
			  rm data/${SPECIES}/${ERR}_*.fastq.gz
			fi

			# Remove lock file
			rm data/${SPECIES}/${ERR}.lck
		fi
	fi
}


# Ensure samples selected and genome indexes available.
# We may have parallel scripts running, so ensure only one runs this step
if [ ! -e data/preMapping.lck ]; then
	touch data/preMapping.lck
	echo "Running genome indexing and sample selection" > logs/preMapping.log 2>&1
	# Ensure all genome and annotations are present
	bash src/makeIndexedGenomes.sh >> logs/preMapping.log 2>&1
	if [ $? -ne 0 ]; then
		echo "Error making genome indexes, exiting" >> logs/preMapping.log 2>&1
		rm data/preMapping.lck
		exit 1
	fi
	# Select samples to map from metadata
	Rscript src/selectSamples.R >> logs/preMapping.log 2>&1
	if [ $? -ne 0 ]; then
		echo "Error running sample selection, exiting" >> logs/preMapping.log 2>&1
		rm data/preMapping.lck
		exit 1
	fi
	rm data/preMapping.lck
else
	# Other instances of this script wait for lock to release
	while [ -e data/preMapping.lck ]; do
		sleep 30
	done
fi

echo "Processing single end samples" >> logs/mapping.log 2>&1
for LINE in ${SE_SAMPLES}; do
	ERR=$(echo ${LINE} | cut -f 1 -d , )
	SPECIES=$(echo ${LINE} | cut -f 3 -d , )
	GENOME=$(echo ${LINE} | cut -f 4 -d , )
	GTF_FILE=$(echo ${LINE} | cut -f 5 -d , )
	echo ${LINE} >> logs/mapping.log 2>&1
	map_se_sample ${SPECIES} ${ERR} ${GENOME} ${GTF_FILE}
done

echo "Processing paired end samples" >> logs/mapping.log 2>&1
for LINE in ${PE_SAMPLES}; do
	ERR=$(echo ${LINE} | cut -f 1 -d , )
	SPECIES=$(echo ${LINE} | cut -f 3 -d , )
	GENOME=$(echo ${LINE} | cut -f 4 -d , )
	GTF_FILE=$(echo ${LINE} | cut -f 5 -d , )
	echo ${LINE} >> logs/mapping.log 2>&1
	map_pe_sample ${SPECIES} ${ERR} ${GENOME} ${GTF_FILE}
done
