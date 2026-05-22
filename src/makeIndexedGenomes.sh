#!/bin/bash

# Create HISAT indexed genomes
# Invoke from project base directory

# Genome indexes created for hisat2. Note we can't use gz fa file in hisat2-build

# Create a genome index
# 1 - the name of the genome e.g. GRCg7b
# 2 - the URL of the FASTA sequences
# 3 - the URL of the GTF annotation
build_genome_index () {
	GENOME=$(echo $1 | tr -d '"')
	FASTAURL=$(echo $2 | tr -d '"')
	GTFURL=$(echo $3 | tr -d '"')

	if [ ! -e ${GENOME}.1.ht2 ] && [ ! -e ${GENOME}.1.ht2l ]; then # could be .ht2 or .ht2l for large genomes
	
		FASTAGZFILE=$(basename ${FASTAURL})
		GTFGZFILE=$(basename ${GTFURL})

		GTFFILE=$(echo $GTFGZFILE | sed -e 's/.gz//')
		SSFILE=$(echo $GTFFILE | sed -e 's/gtf/ss/')
		EXONFILE=$(echo $GTFFILE | sed -e 's/gtf/exons/')
		FASTAFILE=$(echo $FASTAGZFILE | sed -e 's/.gz//')
  	
		# Get the annotations and sequence
		wget ${FASTAURL}
		wget ${GTFURL}

		gunzip $FASTAGZFILE
		gunzip $GTFGZFILE

		# Extract splice and exon coordinates from annotations
		hisat2_extract_splice_sites.py $GTFFILE >  $SSFILE
		hisat2_extract_exons.py $GTFFILE > $EXONFILE

		# Make the genome - note we can't use fa.gz file in hisat2-build
		# so rezip once complete
		hisat2-build --ss ${SSFILE} --exon ${EXONFILE} ${FASTAFILE} ${GENOME}
		gzip $FASTAFILE
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
	GENOME=$(echo ${LINE} | cut -f 1 -d , )
	FASTA_URL=$(echo ${LINE} | cut -f 4 -d , )
	GTF_URL=$(echo ${LINE} | cut -f 5 -d , )
	build_genome_index "${GENOME}" "${FASTA_URL}" "${GTF_URL}"
done < "../metadata/genomes.csv"
