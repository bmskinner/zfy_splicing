#!/bin/bash

# Create HISAT indexed genomes
# Invoke from project base directory

# Genome indexes created for hisat2. Note we can't use gz fa file in hisat2-build

# Create a genome index
# 1 - the name of the genome e.g. GRCg7b
# 2 - the URL of the FASTA sequences
# 3 - the URL of the GTF annotation
# 4 - the desired output file name for the FASTA sequences
# 5 - the desired output file name for the GTF
build_genome_index () {
	GENOME=$(echo $1 | tr -d '"')
	FASTA_URL=$(echo $2 | tr -d '"')
	GTF_URL=$(echo $3 | tr -d '"')
	FASTA_FILE=$(echo $4 | tr -d '"') # should be .gz
	FASTA_FILE=$(basename $FASTA_FILE)
	GTF_FILE=$(echo $5 | tr -d '"') # should not be .gz
	GTF_FILE=$(basename $GTF_FILE)
	
	# The name of the unzipped files - may not be the desired final name
	# e.g. genes.gtf
	FASTA_GZ_FILE=$(basename ${FASTA_URL})
	FASTA_RAW_FILE=$(echo $FASTA_GZ_FILE | sed -e 's/.gz//')
	
	GTF_GZ_FILE=$(basename ${GTF_URL})
	GTF_RAW_FILE=$(echo $GTF_GZ_FILE | sed -e 's/.gz//')
	
	# Ensure we have genome GTF
	if [ ! -e ${GTF_FILE} ]; then
	  # download and unzip the GTF
	  echo "${GENOME}: Downloading GTF from ${GTF_URL}"
		wget --quiet ${GTF_URL}
		gunzip $GTF_GZ_FILE
		# Downloaded filename may not match desired name - move. Noop if name matches
		echo "${GENOME}: Moving GTF file from ${GTF_RAW_FILE} to ${GTF_FILE}"
		mv ${GTF_RAW_FILE} ${GTF_FILE}
	fi
	
	# Ensure we have a gzipped genome FASTA.
  if [ ! -e ${FASTA_FILE} ]; then
	  echo "${GENOME}: Downloading FASTA from ${FASTA_URL}"

	  wget --quiet ${FASTA_URL}
	  # Ensure gz file is correctly named
	  mv ${FASTA_GZ_FILE} ${FASTA_FILE}
  fi
	 
	# Create index for histat2
	if [ ! -e ${GENOME}.1.ht2 ] && [ ! -e ${GENOME}.1.ht2l ]; then # could be .ht2 or .ht2l for large genomes
		echo "${GENOME}: Creating hisat genome index"
		
		# -c keeps original file
		gunzip -c ${FASTA_FILE} > ${FASTA_RAW_FILE}

		# Create splice site and exon names
		SSFILE=$(echo ${GTF_FILE} | sed -e 's/gtf/ss/')
		EXONFILE=$(echo ${GTF_FILE} | sed -e 's/gtf/exons/')

		# Extract splice and exon coordinates from GTF annotations
		hisat2_extract_splice_sites.py ${GTF_FILE} > ${SSFILE}
		hisat2_extract_exons.py ${GTF_FILE} > ${EXONFILE}

		# Make the genome index - note we can't use fa.gz file in hisat2-build
		# so rezip once complete and move to final name
		hisat2-build --ss ${SSFILE} --exon ${EXONFILE} ${FASTA_RAW_FILE} ${GENOME}
	  rm ${FASTA_RAW_FILE}
	fi
	
	# Create index for minimap2 long read mapping
	if [ ! -e ${GENOME}.mmi ]; then
	  echo "${GENOME}: Creating minimap2 genome index"
	  gunzip -c ${FASTA_FILE} > ${FASTA_RAW_FILE}
		minimap2 -d ${GENOME}.mmi ${FASTA_RAW_FILE}
		if [ $? -ne 0 ]; then
		  echo "${GENOME}: Could not run genome indexing, exiting"
			exit 1
		fi
		rm  ${FASTA_RAW_FILE}
	fi
}

# Switch working directory for downloads
mkdir -p genomes
cd genomes
	
 # Check metadata exists
if [ ! -e "../metadata/genomes.csv" ]; then
	echo "Could not find relative file ../metadata/genomes.csv"
	echo "Running from: `pwd`"
	exit 1
fi

# Read the genomes metadata file and create indexes if missing
# Expecting csv
while read LINE; do
	GENOME=$(echo ${LINE} | cut -f 1 -d , | tr -d '"')
	FASTA_URL=$(echo ${LINE} | cut -f 4 -d , )
	GTF_URL=$(echo ${LINE} | cut -f 5 -d , )
	FASTA_FILE=$(echo ${LINE} | cut -f 6 -d , )
	GTF_FILE=$(echo ${LINE} | cut -f 7 -d , )
	if [ ${GENOME} == "Genome" ]; then
		continue # skip header line
	fi
	build_genome_index "${GENOME}" "${FASTA_URL}" "${GTF_URL}" "${FASTA_FILE}" "${GTF_FILE}"
done < "../metadata/genomes.csv"
