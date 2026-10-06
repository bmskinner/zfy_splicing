#!/bin/bash

# Map single end sample against a CDS genome - no annotation available, and 
# unaligned reads are discarded.
map_se_sample () {
  ERR=$(echo $1 | tr -d '"')
	SPECIES=$(echo $2 | tr -d '"')
	GENOME=$(echo $3 | tr -d '"')

	mkdir -p cds/data/${SPECIES}
	
	# Run if final output is missing
	if [ ! -e cds/data/${SPECIES}/${ERR}.bam ]; then
	
		# Don't work on a sample already being processed
		if [ ! -e cds/data/${SPECIES}/${ERR}.lck ]; then

			touch cds/data/${SPECIES}/${ERR}.lck
			echo "`date '+%Y-%m-%d %X'` ${ERR}: bam not found" >> logs/${ERR}.mapping.log 2>&1

			# Check for existing downloads before running fasterq-dump
			if [ ! -e cds/data/${SPECIES}/${ERR}.fastq.gz ] && [ ! -e cds/data/${SPECIES}/${ERR}_trimmed.fq.gz ]; then
				# Fetch data
				if [ ! -e cds/data/${SPECIES}/${ERR}.fastq ]; then
					echo "`date '+%Y-%m-%d %X'` ${ERR}: downloading fastq" >> logs/${ERR}.mapping.log 2>&1
					fasterq-dump -o cds/data/${SPECIES}/${ERR}.fastq ${ERR}
				fi
				gzip cds/data/${SPECIES}/${ERR}.fastq
			fi

			# Trim
			if [ ! -e cds/data/${SPECIES}/${ERR}_trimmed.fq.gz ]; then
				echo "`date '+%Y-%m-%d %X'` ${ERR}: trimming" >> logs/${ERR}.mapping.log 2>&1
				trim_galore -o data/${SPECIES} --suppress_warn --fastqc --fastqc_args "-t 8 --outdir report/FASTQC --nogroup --extract" cds/data/${SPECIES}/${ERR}.fastq.gz
			fi

			# map
			if [ ! -e cds/data/${SPECIES}/${ERR}.bam.csi ]; then
  			echo "`date '+%Y-%m-%d %X'` ${ERR}: mapping" >> logs/${ERR}.mapping.log 2>&1

      	hisat2 -x cds/genomes/${GENOME} --threads 8 -U cds/data/${SPECIES}/${ERR}_trimmed.fq.gz --new-summary --no-unal -S cds/data/${SPECIES}/${ERR}.sam >> logs/${ERR}.mapping.log 2>&1
  			samtools sort -T cds/data/${SPECIES}/${ERR} -@ 8 -o cds/data/${SPECIES}/${ERR}.bam cds/data/${SPECIES}/${ERR}.sam
  			samtools index -c -@ 7 cds/data/${SPECIES}/${ERR}.bam # index with csi due to long chromosomes in opossum
  			rm cds/data/${SPECIES}/${ERR}.sam
  		fi

			# Remove original FASTQ, we have the trimmed reads still
			if [ -e cds/data/${SPECIES}/${ERR}.fastq.gz ]; then
			  rm cds/data/${SPECIES}/${ERR}.fastq.gz
			fi

			# Remove lock file
			rm cds/data/${SPECIES}/${ERR}.lck
		fi
	fi
}

# Map paired end sample against a CDS genome - no annotation available, and 
# unaligned reads are discarded.
map_pe_sample () {
  ERR=$(echo $1 | tr -d '"')
	SPECIES=$(echo $2 | tr -d '"')
	GENOME=$(echo $3 | tr -d '"')
	
	mkdir -p cds/data/${SPECIES}
	
	# Run if final output is missing
	if [ ! -e cds/data/${SPECIES}/${ERR}.bam ]; then
	
		# Don't work on a sample already being processed
		if [ ! -e cds/data/${SPECIES}/${ERR}.lck ]; then

			touch cds/data/${SPECIES}/${ERR}.lck
			echo "`date '+%Y-%m-%d %X'` ${ERR}: bam not found" >> logs/${ERR}.mapping.log 2>&1

			# Check for existing downloads before running fasterq-dump
			if [ ! -e cds/data/${SPECIES}/${ERR}_2.fastq.gz ] && [ ! -e data/${SPECIES}/${ERR}_2_val_2.fq.gz ]; then
				# Fetch data
				if [ ! -e cds/data/${SPECIES}/${ERR}_2.fastq ]; then
					echo "`date '+%Y-%m-%d %X'` ${ERR}: downloading fastq" >> logs/${ERR}.mapping.log 2>&1
					fasterq-dump -O cds/data/${SPECIES}/ ${ERR}
				fi
				gzip cds/data/${SPECIES}/${ERR}_1.fastq
				gzip cds/data/${SPECIES}/${ERR}_2.fastq
			fi

			# Trim
			if [ ! -e cds/data/${SPECIES}/${ERR}_1_val_1.fq.gz ]; then
				echo "`date '+%Y-%m-%d %X'` ${ERR}: trimming" >> logs/${ERR}.mapping.log 2>&1
				trim_galore -o cds/data/${SPECIES} --paired --suppress_warn --fastqc --fastqc_args "-t 8 --outdir report/FASTQC --nogroup --extract" cds/data/${SPECIES}/${ERR}_1.fastq.gz cds/data/${SPECIES}/${ERR}_2.fastq.gz
			fi

			# map
			if [ ! -e cds/data/${SPECIES}/${ERR}.bam.csi ]; then
  			echo "`date '+%Y-%m-%d %X'` ${ERR}: mapping against ${GENOME}" >> logs/${ERR}.mapping.log 2>&1
  			# -k controls number of multimapping locations (default 5 for linear index)
  			hisat2 -x cds/genomes/${GENOME} --threads 8 -1 cds/data/${SPECIES}/${ERR}_1_val_1.fq.gz -2 cds/data/${SPECIES}/${ERR}_2_val_2.fq.gz --new-summary --no-unal -S cds/data/${SPECIES}/${ERR}.sam >> logs/${ERR}.mapping.log 2>&1
  			samtools sort -T cds/data/${SPECIES}/${ERR} -@ 8 -o cds/data/${SPECIES}/${ERR}.bam cds/data/${SPECIES}/${ERR}.sam
  			samtools index -c -@ 7 cds/data/${SPECIES}/${ERR}.bam # index with csi due to long chromosomes in opossum
  			rm cds/data/${SPECIES}/${ERR}.sam
			fi

			# Remove original FASTQ, we have the trimmed reads still
			if [ -e cds/data/${SPECIES}/${ERR}_2.fastq.gz ]; then
			  rm cds/data/${SPECIES}/${ERR}_*.fastq.gz
			fi

			# Remove lock file
			rm cds/data/${SPECIES}/${ERR}.lck
		fi
	fi
}

SE_SAMPLES=$(cat cds/cds.samples.csv | cut -f 1,2,3 -d , | grep -e '[S|E|D]RR' | grep -e 'SINGLE')
PE_SAMPLES=$(cat cds/cds.samples.csv | cut -f 1,2,3 -d , | grep -e '[S|E|D]RR' | grep -e 'PAIRED')

echo "`date '+%Y-%m-%d %X'` Processing single end CDS samples" >> logs/mapping.log 2>&1
for LINE in ${SE_SAMPLES}; do
	ERR=$(echo ${LINE} | cut -f 1 -d , )
	SPECIES=$(echo ${LINE} | cut -f 2 -d , )
	GENOME=$(echo ${LINE} | cut -f 3 -d , )
	map_se_sample ${ERR} ${SPECIES} ${GENOME}
done

echo "`date '+%Y-%m-%d %X'` Processing paired end CDS samples" >> logs/mapping.log 2>&1
for LINE in ${PE_SAMPLES}; do
	ERR=$(echo ${LINE} | cut -f 1 -d , )
	SPECIES=$(echo ${LINE} | cut -f 2 -d , )
	GENOME=$(echo ${LINE} | cut -f 3 -d , )
	map_pe_sample ${ERR} ${SPECIES} ${GENOME}
done
