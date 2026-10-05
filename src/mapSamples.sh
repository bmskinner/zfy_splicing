#!/bin/bash

mkdir -p report/FASTQC

# Map a single end sample
# $1 species e.g chicken - should match a folder name in ./data
# $2 ERR id e.g ERR2576379
# $3 genome build id e.g. GRCg7b - should match a .ht2 prefix in ./genomes
# $4 gtf file corresponding to the genome build e.g. ./genomes/example.gtf
# $5 Gene ID region to be exracted from the GTF e.g. ZFX
# $6 Location of the gene in the GTF e.g. 1:123-456
map_se_sample () {
	SPECIES=$(echo $1 | tr -d '"')
	ERR=$(echo $2 | tr -d '"')
	GENOME=$(echo $3 | tr -d '"')
	GTF_FILE=$(echo $4 | tr -d '"')
	GENE_ID=$(echo $5 | tr -d '"')
	LOCATION=$(echo $6 | tr -d '"')
	
	mkdir -p data/${SPECIES}
	
	if [ ! -e ${GTF_FILE} ]; then
	  echo "`date '+%Y-%m-%d %X'` ${ERR}: Could not find GTF file ${GTF_FILE}, skipping" >> logs/mapping.log 2>&1
	  continue
	fi
	
	# Run if final output is missing
	if [ ! -e data/${SPECIES}/${ERR}.${GENE_ID}.gtf ]; then
	
		echo "`date '+%Y-%m-%d %X'` ${ERR}: finding reads covering ${GENE_ID} at ${LOCATION} in ${SPECIES}" >> logs/mapping.log 2>&1

		# Don't work on a sample already being processed
		if [ ! -e data/${SPECIES}/${ERR}.lck ]; then

			touch data/${SPECIES}/${ERR}.lck
			echo "`date '+%Y-%m-%d %X'` ${ERR}: bam not found" >> logs/${ERR}.mapping.log 2>&1

			# Check for existing downloads before running fasterq-dump
			if [ ! -e data/${SPECIES}/${ERR}.fastq.gz ] && [ ! -e data/${SPECIES}/${ERR}_trimmed.fq.gz ]; then
				# Fetch data
				if [ ! -e data/${SPECIES}/${ERR}.fastq ]; then
					echo "`date '+%Y-%m-%d %X'` ${ERR}: downloading fastq" >> logs/${ERR}.mapping.log 2>&1
					fasterq-dump -o data/${SPECIES}/${ERR}.fastq ${ERR}
				fi
				gzip data/${SPECIES}/${ERR}.fastq
			fi

			# Trim
			if [ ! -e data/${SPECIES}/${ERR}_trimmed.fq.gz ]; then
				echo "`date '+%Y-%m-%d %X'` ${ERR}: trimming" >> logs/${ERR}.mapping.log 2>&1
				trim_galore -o data/${SPECIES} --suppress_warn --fastqc --fastqc_args "-t 8 --outdir report/FASTQC --nogroup --extract" data/${SPECIES}/${ERR}.fastq.gz
			fi

			# map
			if [ ! -e data/${SPECIES}/${ERR}.bam.csi ]; then
  			echo "`date '+%Y-%m-%d %X'` ${ERR}: mapping" >> logs/${ERR}.mapping.log 2>&1
  			# -k controls number of multimapping locations (default 5 for linear index)
  			# --downstream-transcriptome-assembly forces longer anchors at novel splice sites (more rigorous)
  			# --dta-cfflinks does this and also looks for novel splice sites, stored in tag XS:A:[+-]
  			# Note that a stranded library may be with respect to forward or reverse strand depending on prep method;
  			# Use --rna-strandness R to specify single-end RNA-seq data is reverse stranded, F for forward strand.
  			# Analysis run with and without --dta to compare effects; we don't need to assemble transcripts, just see if there is greater splicing in testis
  			# hisat2 -x genomes/${GENOME} -p 8 -U data/${SPECIES}/${ERR}_trimmed.fq.gz --new-summary --downstream-transcriptome-assembly -S data/${SPECIES}/${ERR}.sam
  			hisat2 -x genomes/${GENOME} --threads 8 -U data/${SPECIES}/${ERR}_trimmed.fq.gz --new-summary --dta-cufflinks -S data/${SPECIES}/${ERR}.sam >> logs/${ERR}.mapping.log 2>&1
  			samtools sort -T data/${SPECIES}/${ERR} -@ 8 -o data/${SPECIES}/${ERR}.bam data/${SPECIES}/${ERR}.sam
  			samtools index -c -@ 7 data/${SPECIES}/${ERR}.bam # index with csi due to long chromosomes in opossum
  			rm data/${SPECIES}/${ERR}.sam
  		fi
			
			# Run feature counts for expression quantification
			if [ ! -e data/${SPECIES}/${ERR}.counts.txt ]; then
		  	featureCounts -t exon -g gene_id -a ${GTF_FILE} -o data/${SPECIES}/${ERR}.counts.txt data/${SPECIES}/${ERR}.bam
		  fi
		  
		  # Extract the region of interest only and index
		  samtools view -@ 7 -o data/${SPECIES}/${ERR}.${GENE_ID}.bam data/${SPECIES}/${ERR}.bam ${LOCATION}
		  samtools index -c -@ 7 data/${SPECIES}/${ERR}.${GENE_ID}.bam
		  
		  # Assemble transcripts with Stringtie
		  stringtie -o data/${SPECIES}/${ERR}.${GENE_ID}.gtf -p 1 -l ${SPECIES} -G ${GTF_FILE} -f 0.01 data/${SPECIES}/${ERR}.${GENE_ID}.bam
		  
		  # If stringtie fails, the usual reason is failure to parse the gtf file
		  # - most common reason is rows with an empty transcript_id field. Create a
		  # new GTF if needed and try again.			
			if [ $? -ne 0 ]; then
		    if [ ! -e  ${GTF_FILE}.no.gene.gtf ]; then
		      awk '$3 != "gene" ' ${GTF_FILE} > ${GTF_FILE}.no.gene.gtf
		    fi
		    stringtie -o data/${SPECIES}/${ERR}.${GENE_ID}.gtf -p 1 -l ${SPECIES} -G ${GTF_FILE}.no.gene.gtf -f 0.01 data/${SPECIES}/${ERR}.${GENE_ID}.bam >> logs/${ERR}.mapping.log 2>&1
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
# $5 Gene ID region to be exracted from the GTF e.g. ZFX
# $6 Location of the gene in the GTF e.g. 1:123-456
map_pe_sample () {
	SPECIES=$(echo $1 | tr -d '"')
	ERR=$(echo $2 | tr -d '"')
	GENOME=$(echo $3 | tr -d '"')
	GTF_FILE=$(echo $4 | tr -d '"')
	GENE_ID=$(echo $5 | tr -d '"')
	LOCATION=$(echo $6 | tr -d '"')
	
	mkdir -p data/${SPECIES}
	
	if [ ! -e ${GTF_FILE} ]; then
	  echo "`date '+%Y-%m-%d %X'` ${ERR}: Could not find GTF file ${GTF_FILE}, skipping" >> logs/mapping.log 2>&1
	  continue
	fi
	
	# Run if final output is missing
	if [ ! -e data/${SPECIES}/${ERR}.${GENE_ID}.gtf ]; then
	
	  echo "`date '+%Y-%m-%d %X'` ${ERR}: finding reads covering ${GENE_ID} at ${LOCATION} in ${SPECIES}" >> logs/mapping.log 2>&1

		# Don't work on a sample already being processed
		if [ ! -e data/${SPECIES}/${ERR}.lck ]; then

			touch data/${SPECIES}/${ERR}.lck
			echo "`date '+%Y-%m-%d %X'` ${ERR}: Stringtie GTF output not found" >> logs/${ERR}.mapping.log 2>&1

			# Check for existing downloads before running fasterq-dump
			if [ ! -e data/${SPECIES}/${ERR}_2.fastq.gz ] && [ ! -e data/${SPECIES}/${ERR}_2_val_2.fq.gz ]; then
				# Fetch data
				if [ ! -e data/${SPECIES}/${ERR}_2.fastq ]; then
					echo "`date '+%Y-%m-%d %X'` ${ERR}: downloading fastq" >> logs/${ERR}.mapping.log 2>&1
					fasterq-dump -O data/${SPECIES}/ ${ERR}
				fi
				gzip data/${SPECIES}/${ERR}_1.fastq
				gzip data/${SPECIES}/${ERR}_2.fastq
			fi

			# Trim
			if [ ! -e data/${SPECIES}/${ERR}_1_val_1.fq.gz ]; then
				echo "`date '+%Y-%m-%d %X'` ${ERR}: trimming" >> logs/${ERR}.mapping.log 2>&1
				trim_galore -o data/${SPECIES} --paired --suppress_warn --fastqc --fastqc_args "-t 8 --outdir report/FASTQC --nogroup --extract" data/${SPECIES}/${ERR}_1.fastq.gz data/${SPECIES}/${ERR}_2.fastq.gz
			fi

			# map
			if [ ! -e data/${SPECIES}/${ERR}.bam.csi ]; then
  			echo "`date '+%Y-%m-%d %X'` ${ERR}: mapping against ${GENOME}" >> logs/${ERR}.mapping.log 2>&1
  			# -k controls number of multimapping locations (default 5 for linear index)
  			hisat2 -x genomes/${GENOME} --threads 8 -1 data/${SPECIES}/${ERR}_1_val_1.fq.gz -2 data/${SPECIES}/${ERR}_2_val_2.fq.gz --new-summary --dta-cufflinks -S data/${SPECIES}/${ERR}.sam >> logs/${ERR}.mapping.log 2>&1
  			samtools sort -T data/${SPECIES}/${ERR} -@ 8 -o data/${SPECIES}/${ERR}.bam data/${SPECIES}/${ERR}.sam
  			samtools index -c -@ 7 data/${SPECIES}/${ERR}.bam # index with csi due to long chromosomes in opossum
  			rm data/${SPECIES}/${ERR}.sam
			fi
			
			# Run feature counts for expression quantification
			if [ ! -e data/${SPECIES}/${ERR}.counts.txt ]; then
			  featureCounts -p --countReadPairs -t exon -g gene_id -a ${GTF_FILE} -o data/${SPECIES}/${ERR}.counts.txt data/${SPECIES}/${ERR}.bam
			fi
			
			# Extract the region of interest only and index
		  samtools view -@ 7 -o data/${SPECIES}/${ERR}.${GENE_ID}.bam data/${SPECIES}/${ERR}.bam ${LOCATION}
		  samtools index -c -@ 7 data/${SPECIES}/${ERR}.${GENE_ID}.bam
		  
		  # Assemble transcripts with Stringtie
		  stringtie -o data/${SPECIES}/${ERR}.${GENE_ID}.gtf -p 1 -l ${SPECIES} -G ${GTF_FILE} -f 0.01 data/${SPECIES}/${ERR}.${GENE_ID}.bam >> logs/${ERR}.mapping.log 2>&1
		  
		  # If stringtie fails, the usual reason is failure to parse the gtf file
		  # - most common reason is rows with an empty transcript_id field. Create a
		  # new GTF if needed and try again.			
			if [ $? -ne 0 ]; then
		    if [ ! -e  ${GTF_FILE}.no.gene.gtf ]; then
		      awk '$3 != "gene" ' ${GTF_FILE} > ${GTF_FILE}.no.gene.gtf
		    fi
		    stringtie -o data/${SPECIES}/${ERR}.${GENE_ID}.gtf -p 1 -l ${SPECIES} -G ${GTF_FILE}.no.gene.gtf -f 0.01 data/${SPECIES}/${ERR}.${GENE_ID}.bam >> logs/${ERR}.mapping.log 2>&1
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

# Map a long read sample (PacBio or Nanopore)
# $1 species e.g chicken - should match a folder name in ./data
# $2 ERR id e.g ERR2576379
# $3 genome build id e.g. GRCg7b - should match a .ht2 prefix in ./genomes
# $4 gtf file corresponding to the genome build e.g. ./genomes/example.gtf
# $5 Gene ID region to be exracted from the GTF e.g. ZFX
# $6 Location of the gene in the GTF e.g. 1:123-456
map_long_read_sample () {
	SPECIES=$(echo $1 | tr -d '"')
	ERR=$(echo $2 | tr -d '"')
	GENOME=$(echo $3 | tr -d '"')
	GTF_FILE=$(echo $4 | tr -d '"')
	GENE_ID=$(echo $5 | tr -d '"')
	LOCATION=$(echo $6 | tr -d '"')

	mkdir -p data/${SPECIES}
	
	if [ ! -e ${GTF_FILE} ]; then
	  echo "`date '+%Y-%m-%d %X'` ${ERR}: Could not find GTF file ${GTF_FILE}, skipping" >> logs/mapping.log 2>&1
	  continue
	fi
	
	# Run if final output is missing
	if [ ! -e data/${SPECIES}/${ERR}.${GENE_ID}.gtf ]; then

		# Don't work on a sample already being processed
		if [ ! -e data/${SPECIES}/${ERR}.lck ]; then

			touch data/${SPECIES}/${ERR}.lck
			echo "`date '+%Y-%m-%d %X'` ${ERR}: Stringtie GTF output not found" >> logs/${ERR}.mapping.log 2>&1

			# Check for existing downloads before running fasterq-dump
			if [ ! -e data/${SPECIES}/${ERR}.fastq.gz ]; then
				# Fetch data
				if [ ! -e data/${SPECIES}/${ERR}_2.fastq ]; then
					echo "${ERR}: downloading fastq" >> logs/${ERR}.mapping.log 2>&1
					fasterq-dump -O data/${SPECIES}/ ${ERR}
				fi
				gzip data/${SPECIES}/${ERR}.fastq
			fi

				# Trim and QC
			if [ ! -e data/${SPECIES}/${ERR}_trimmed.fq.gz ]; then
				echo "`date '+%Y-%m-%d %X'` ${ERR}: QC and trimming" >> logs/${ERR}.mapping.log 2>&1
				fastplong --in data/${SPECIES}/${ERR}.fastq.gz --out data/${SPECIES}/${ERR}_trimmed.fq.gz --thread 6 --qualified_quality_phred 9 --json report/QC/fastp/${ERR}.json --report_title "${ERR} ${SPECIES}"
				if [ $? -ne 0 ]; then
			    echo "`date '+%Y-%m-%d %X'` Error running fastplong, exiting" >> logs/${ERR}.mapping.log 2>&1
			    exit 1
			  fi
			fi

			# mapping
			if [ ! -e data/${SPECIES}/${ERR}.bam.csi ]; then
  			echo "`date '+%Y-%m-%d %X'` ${ERR}: mapping against ${GENOME}" >> logs/${ERR}.mapping.log 2>&1
  			minimap2 -a genomes/${GENOME}.mmi -x splice:hq -u b -t 7 data/${SPECIES}/${ERR}_trimmed.fq.gz > data/${SPECIES}/${ERR}.sam
  			samtools sort -T data/${SPECIES}/${ERR} -@ 8 -o data/${SPECIES}/${ERR}.bam data/${SPECIES}/${ERR}.sam
  			samtools index -c -@ 7 data/${SPECIES}/${ERR}.bam # index with csi incase of long chromosomes
  			rm data/${SPECIES}/${ERR}.sam
			fi
			
			# Run feature counts for expression quantification
			if [ ! -e data/${SPECIES}/${ERR}.counts.txt ]; then
			  featureCounts -p --countReadPairs -t exon -g gene_id -a ${GTF_FILE} -o data/${SPECIES}/${ERR}.counts.txt data/${SPECIES}/${ERR}.bam
			fi
			
			# Extract the region of interest only and index
		  samtools view -@ 7 -o data/${SPECIES}/${ERR}.${GENE_ID}.bam data/${SPECIES}/${ERR}.bam ${LOCATION}
		  samtools index -c -@ 7 data/${SPECIES}/${ERR}.${GENE_ID}.bam
		  
		  # Assemble transcripts with Stringtie
		  stringtie -o data/${SPECIES}/${ERR}.${GENE_ID}.gtf -p 1 -l ${SPECIES} -G ${GTF_FILE} -f 0.01 data/${SPECIES}/${ERR}.${GENE_ID}.bam >> logs/${ERR}.mapping.log 2>&1
		  
		  # If stringtie fails, the usual reason is failure to parse the gtf file
		  # - most common reason is rows with an empty transcript_id field. Create a
		  # new GTF if needed and try again.
		  if [ $? -ne 0 ]; then
		    if [ ! -e  ${GTF_FILE}.no.gene.gtf ]; then
		      awk '$3 != "gene" ' ${GTF_FILE} > ${GTF_FILE}.no.gene.gtf
		    fi
		    stringtie -o data/${SPECIES}/${ERR}.${GENE_ID}.gtf -p 1 -l ${SPECIES} -G ${GTF_FILE}.no.gene.gtf -f 0.01 data/${SPECIES}/${ERR}.${GENE_ID}.bam >> logs/${ERR}.mapping.log 2>&1
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


# Ensure samples selected and genome indexes available.
# We may have parallel scripts running, so ensure only one runs this step
if [ ! -e data/preMapping.lck ]; then
	touch data/preMapping.lck
	echo "`date '+%Y-%m-%d %X'` Running genome indexing and sample selection" > logs/preMapping.log 2>&1
	# Ensure all genome and annotations are present
	bash src/makeIndexedGenomes.sh >> logs/preMapping.log 2>&1
	if [ $? -ne 0 ]; then
		echo "`date '+%Y-%m-%d %X'` Error making genome indexes, exiting" >> logs/preMapping.log 2>&1
		rm data/preMapping.lck
		exit 1
	fi
	# Select samples to map from metadata by loading all functions and data
	Rscript src/functions.R >> logs/preMapping.log 2>&1
	if [ $? -ne 0 ]; then
		echo "`date '+%Y-%m-%d %X'` Error running sample selection, exiting" >> logs/preMapping.log 2>&1
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

# Remap samples and generate bam files for splice junction detection
# Select the Run, Library type, species and genome columns
SE_SAMPLES=$(cat data/mapping.samples.csv | cut -f 1,2,3,4,5,6,7 -d , | grep -e '[S|E|D]RR' | grep -e 'SINGLE')
PE_SAMPLES=$(cat data/mapping.samples.csv | cut -f 1,2,3,4,5,6,7 -d , | grep -e '[S|E|D]RR' | grep -e 'PAIRED')
LR_SAMPLES=$(cat data/mapping.samples.csv | cut -f 1,2,3,4,5,6,7 -d , | grep -e '[S|E|D]RR' | grep -e 'LONG_READ')

echo "`date '+%Y-%m-%d %X'` Processing single end samples" >> logs/mapping.log 2>&1
for LINE in ${SE_SAMPLES}; do
	ERR=$(echo ${LINE} | cut -f 1 -d , )
	SPECIES=$(echo ${LINE} | cut -f 3 -d , )
	GENOME=$(echo ${LINE} | cut -f 4 -d , )
	GTF_FILE=$(echo ${LINE} | cut -f 5 -d , )
	GENE_ID=$(echo ${LINE} | cut -f 6 -d , )
	LOCATION=$(echo ${LINE} | cut -f 7 -d , )
	map_se_sample ${SPECIES} ${ERR} ${GENOME} ${GTF_FILE} ${GENE_ID} ${LOCATION}
done

echo "`date '+%Y-%m-%d %X'` Processing paired end samples" >> logs/mapping.log 2>&1
for LINE in ${PE_SAMPLES}; do
	ERR=$(echo ${LINE} | cut -f 1 -d , )
	SPECIES=$(echo ${LINE} | cut -f 3 -d , )
	GENOME=$(echo ${LINE} | cut -f 4 -d , )
	GTF_FILE=$(echo ${LINE} | cut -f 5 -d , )
	GENE_ID=$(echo ${LINE} | cut -f 6 -d , )
	LOCATION=$(echo ${LINE} | cut -f 7 -d , )
	map_pe_sample ${SPECIES} ${ERR} ${GENOME} ${GTF_FILE} ${GENE_ID} ${LOCATION}
done

echo "`date '+%Y-%m-%d %X'` Processing long read samples" >> logs/mapping.log 2>&1
for LINE in ${LR_SAMPLES}; do
	ERR=$(echo ${LINE} | cut -f 1 -d , )
	SPECIES=$(echo ${LINE} | cut -f 3 -d , )
	GENOME=$(echo ${LINE} | cut -f 4 -d , )
	GTF_FILE=$(echo ${LINE} | cut -f 5 -d , )
	GENE_ID=$(echo ${LINE} | cut -f 6 -d , )
	LOCATION=$(echo ${LINE} | cut -f 7 -d , )
	map_long_read_sample ${SPECIES} ${ERR} ${GENOME} ${GTF_FILE} ${GENE_ID} ${LOCATION}
done
