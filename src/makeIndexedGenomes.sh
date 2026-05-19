#!/bin/bash

# Create HISAT indexed genomes
# Invoke from project base directory

# Genome indexes created for hisat2. Note we can't use gz fa file in hisat2-build

mkdir -p genomes
cd genomes

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
	  hisat2-build --ss ${SSFILE} --exon ${EXONFILE} ${FASTAFILE} ${GENOME}
    gzip $FASTAFILE
	fi
	
}	

# Chicken GRCg7b
build_genome_index GRCg7b https://ftp.ensembl.org/pub/release-112/fasta/gallus_gallus/dna/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa.gz https://ftp.ensembl.org/pub/release-112/gtf/gallus_gallus/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.112.gtf.gz

# Opossum ASM229v1
build_genome_index ASM229v1 https://ftp.ensembl.org/pub/release-112/fasta/monodelphis_domestica/dna/Monodelphis_domestica.ASM229v1.dna.toplevel.fa.gz https://ftp.ensembl.org/pub/release-112/gtf/monodelphis_domestica/Monodelphis_domestica.ASM229v1.112.gtf.gz

# Mouse GRCm39
build_genome_index GRCm39 https://ftp.ensembl.org/pub/release-112/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.toplevel.fa.gz https://ftp.ensembl.org/pub/release-112/gtf/mus_musculus/Mus_musculus.GRCm39.112.gtf.gz

# Human GRCh38
build_genome_index GRCh38 https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz

# Rat mRatBN7.2
build_genome_index mRatBN7.2 https://ftp.ensembl.org/pub/release-112/fasta/rattus_norvegicus/dna/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.gz https://ftp.ensembl.org/pub/release-112/gtf/rattus_norvegicus/Rattus_norvegicus.mRatBN7.2.112.gtf.gz

# Macaque Mmul_10
build_genome_index Mmul_10 https://ftp.ensembl.org/pub/release-112/fasta/macaca_mulatta/dna/Macaca_mulatta.Mmul_10.dna.toplevel.fa.gz https://ftp.ensembl.org/pub/release-112/gtf/macaca_mulatta/Macaca_mulatta.Mmul_10.112.gtf.gz

# Playtpus mOrnAna1.p.v1
build_genome_index mOrnAna1.p.v1 https://ftp.ensembl.org/pub/release-112/fasta/ornithorhynchus_anatinus/dna/Ornithorhynchus_anatinus.mOrnAna1.p.v1.dna.toplevel.fa.gz  https://ftp.ensembl.org/pub/release-112/gtf/ornithorhynchus_anatinus/Ornithorhynchus_anatinus.mOrnAna1.p.v1.112.gtf.gz

# Zebra finch bTaeGut1_v1.p
build_genome_index bTaeGut1_v1.p https://ftp.ensembl.org/pub/release-112/fasta/taeniopygia_guttata/dna/Taeniopygia_guttata.bTaeGut1_v1.p.dna.toplevel.fa.gz  https://ftp.ensembl.org/pub/release-112/gtf/taeniopygia_guttata/Taeniopygia_guttata.bTaeGut1_v1.p.112.gtf.gz

# Xenopus bTaeGut1_v1.p
build_genome_index UCB_Xtro_10.0 https://ftp.ensembl.org/pub/release-112/fasta/xenopus_tropicalis/dna/Xenopus_tropicalis.UCB_Xtro_10.0.dna.toplevel.fa.gz  https://ftp.ensembl.org/pub/release-112/gtf/xenopus_tropicalis/Xenopus_tropicalis.UCB_Xtro_10.0.112.gtf.gz

# Green Anole AnoCar2.0v2
build_genome_index AnoCar2.0v2 https://ftp.ensembl.org/pub/release-112/fasta/anolis_carolinensis/dna/Anolis_carolinensis.AnoCar2.0v2.dna.toplevel.fa.gz https://ftp.ensembl.org/pub/release-112/gtf/anolis_carolinensis/Anolis_carolinensis.AnoCar2.0v2.112.gtf.gz

# Zebra fish GRCz11
build_genome_index GRCz11 https://ftp.ensembl.org/pub/release-112/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.toplevel.fa.gz https://ftp.ensembl.org/pub/release-112/gtf/danio_rerio/Danio_rerio.GRCz11.112.gtf.gz

# Pig Sscrofa11.1
build_genome_index Sscrofa11.1 https://ftp.ensembl.org/pub/release-112/fasta/sus_scrofa/dna/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz https://ftp.ensembl.org/pub/release-112/gtf/sus_scrofa/Sus_scrofa.Sscrofa11.1.112.gtf.gz

# Koala phaCin_unsw_v4.1
build_genome_index phaCin_unsw_v4.1 https://ftp.ensembl.org/pub/release-115/fasta/phascolarctos_cinereus/dna/Phascolarctos_cinereus.phaCin_unsw_v4.1.dna.toplevel.fa.gz https://ftp.ensembl.org/pub/release-115/gtf/phascolarctos_cinereus/Phascolarctos_cinereus.phaCin_unsw_v4.1.115.gtf.gz

# Echidna mTacAcu1
build_genome_index mTacAcu1 https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/015/852/505/GCF_015852505.1_mTacAcu1.pri/GCF_015852505.1_mTacAcu1.pri_genomic.fna.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/015/852/505/GCF_015852505.1_mTacAcu1.pri/GCF_015852505.1_mTacAcu1.pri_genomic.gtf.gz