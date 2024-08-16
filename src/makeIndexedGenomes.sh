#!/bin/bash

# Create HISAT indexed genomes
# Invoke from project base directory

# Genome indexes created for hisat2. Note we can't use gz fa file in hisat2-build

mkdir -p genomes
cd genomes

# Chicken GRCg7b
if [ ! -e GRCg7b.1.ht2 ]; then
	# Get the annotations and sequence
	wget https://ftp.ensembl.org/pub/release-112/fasta/gallus_gallus/dna/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa.gz
	wget https://ftp.ensembl.org/pub/release-112/gtf/gallus_gallus/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.112.gtf.gz
	gunzip Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.112.gtf.gz
	gunzip Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa.gz

	# Extract splice and exon coordinates from annotations
	hisat2_extract_splice_sites.py Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.112.gtf >  Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.112.ss
	hisat2_extract_exons.py Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.112.gtf >  Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.112.exons

	# Make the genome - note we can't use gz file in hisat2-build

	hisat2-build --ss Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.112.ss --exon Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.112.exons Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa GRCg7b
	gzip Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa
fi

# Opossum ASM229v1
if [ ! -e ASM229v1.1.ht2 ]; then
	# Get the annotations and sequence
	wget https://ftp.ensembl.org/pub/release-112/fasta/monodelphis_domestica/dna/Monodelphis_domestica.ASM229v1.dna.toplevel.fa.gz
	wget https://ftp.ensembl.org/pub/release-112/gtf/monodelphis_domestica/Monodelphis_domestica.ASM229v1.112.gtf.gz
	gunzip Monodelphis_domestica.ASM229v1.112.gtf.gz
	gunzip Monodelphis_domestica.ASM229v1.dna.toplevel.fa.gz

	# Extract splice and exon coordinates from annotations
	hisat2_extract_splice_sites.py Monodelphis_domestica.ASM229v1.112.gtf >  Monodelphis_domestica.ASM229v1.112.ss
	hisat2_extract_exons.py Monodelphis_domestica.ASM229v1.112.gtf >  Monodelphis_domestica.ASM229v1.112.exons

	# Make the genome
	hisat2-build --ss Monodelphis_domestica.ASM229v1.112.ss --exon Monodelphis_domestica.ASM229v1.112.exons Monodelphis_domestica.ASM229v1.dna.toplevel.fa ASM229v1
	gzip Monodelphis_domestica.ASM229v1.dna.toplevel.fa
fi


# Mouse GRCm39
if [ ! -e GRCm39.1.ht2 ]; then
	# Get the annotations and sequence
	wget https://ftp.ensembl.org/pub/release-112/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.toplevel.fa.gz
	wget https://ftp.ensembl.org/pub/release-112/gtf/mus_musculus/Mus_musculus.GRCm39.112.gtf.gz
	gunzip Mus_musculus.GRCm39.112.gtf.gz
	gunzip Mus_musculus.GRCm39.dna.toplevel.fa.gz

	# Extract splice and exon coordinates from annotations
	hisat2_extract_splice_sites.py Mus_musculus.GRCm39.112.gtf >  Mus_musculus.GRCm39.112.ss
	hisat2_extract_exons.py Mus_musculus.GRCm39.112.gtf >  Mus_musculus.GRCm39.112.exons

	# Make the genome
	hisat2-build --ss Mus_musculus.GRCm39.112.ss --exon Mus_musculus.GRCm39.112.exons Mus_musculus.GRCm39.dna.toplevel.fa GRCm39
	gzip Mus_musculus.GRCm39.dna.toplevel.fa
fi

# Human GRCh38
if [ ! -e GRCh38.1.ht2 ]; then
	# Get the annotations and sequence
	wget https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
	wget https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz
	gunzip Homo_sapiens.GRCh38.112.gtf.gz
	gunzip Homo_sapiens.GRCh38.dna.toplevel.fa.gz

	# Extract splice and exon coordinates from annotations
	hisat2_extract_splice_sites.py Homo_sapiens.GRCh38.112.gtf >  Homo_sapiens.GRCh38.112.ss
	hisat2_extract_exons.py Homo_sapiens.GRCh38.112.gtf >  Homo_sapiens.GRCh38.112.exons

	# Make the genome
	hisat2-build --ss Homo_sapiens.GRCh38.112.ss --exon Homo_sapiens.GRCh38.112.exons Homo_sapiens.GRCh38.dna.toplevel.fa GRCh38
	gzip Homo_sapiens.GRCh38.dna.toplevel.fa
fi

# Rabbit OryCun2.0
if [ ! -e OryCun2.0.1.ht2 ]; then
	# Get the annotations and sequence
	wget https://ftp.ensembl.org/pub/release-112/fasta/oryctolagus_cuniculus/dna/Oryctolagus_cuniculus.OryCun2.0.dna.toplevel.fa.gz
	wget https://ftp.ensembl.org/pub/release-112/gtf/oryctolagus_cuniculus/Oryctolagus_cuniculus.OryCun2.0.112.gtf.gz
	gunzip Oryctolagus_cuniculus.OryCun2.0.112.gtf.gz
	gunzip Oryctolagus_cuniculus.OryCun2.0.dna.toplevel.fa.gz

	# Extract splice and exon coordinates from annotations
	hisat2_extract_splice_sites.py Oryctolagus_cuniculus.OryCun2.0.112.gtf >  Oryctolagus_cuniculus.OryCun2.0.112.ss
	hisat2_extract_exons.py Oryctolagus_cuniculus.OryCun2.0.112.gtf >  Oryctolagus_cuniculus.OryCun2.0.112.exons

	# Make the genome
	hisat2-build --ss Oryctolagus_cuniculus.OryCun2.0.112.ss --exon Oryctolagus_cuniculus.OryCun2.0.112.exons Oryctolagus_cuniculus.OryCun2.0.dna.toplevel.fa OryCun2.0
	gzip Oryctolagus_cuniculus.OryCun2.0.dna.toplevel.fa
fi

# Rat mRatBN7.2
if [ ! -e mRatBN7.2.1.ht2 ]; then
	# Get the annotations and sequence
	wget https://ftp.ensembl.org/pub/release-112/fasta/rattus_norvegicus/dna/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.gz
	wget https://ftp.ensembl.org/pub/release-112/gtf/rattus_norvegicus/Rattus_norvegicus.mRatBN7.2.112.gtf.gz
	gunzip Rattus_norvegicus.mRatBN7.2.112.gtf.gz
	gunzip Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.gz

	# Extract splice and exon coordinates from annotations
	hisat2_extract_splice_sites.py Rattus_norvegicus.mRatBN7.2.112.gtf >  Rattus_norvegicus.mRatBN7.2.112.ss
	hisat2_extract_exons.py Rattus_norvegicus.mRatBN7.2.112.gtf >  Rattus_norvegicus.mRatBN7.2.112.exons

	# Make the genome
	hisat2-build --ss Rattus_norvegicus.mRatBN7.2.112.ss --exon Rattus_norvegicus.mRatBN7.2.112.exons Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa mRatBN7.2
	gzip Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa
fi

# Macaque Mmul_10
if [ ! -e Mmul_10.1.ht2 ]; then
	# Get the annotations and sequence
	wget https://ftp.ensembl.org/pub/release-112/fasta/macaca_mulatta/dna/Macaca_mulatta.Mmul_10.dna.toplevel.fa.gz
	wget https://ftp.ensembl.org/pub/release-112/gtf/macaca_mulatta/Macaca_mulatta.Mmul_10.112.gtf.gz
	gunzip Macaca_mulatta.Mmul_10.112.gtf.gz
	gunzip Macaca_mulatta.Mmul_10.dna.toplevel.fa.gz

	# Extract splice and exon coordinates from annotations
	hisat2_extract_splice_sites.py Macaca_mulatta.Mmul_10.112.gtf > Macaca_mulatta.Mmul_10.112.ss
	hisat2_extract_exons.py Macaca_mulatta.Mmul_10.112.gtf > Macaca_mulatta.Mmul_10.112.exons

	# Make the genome
	hisat2-build --ss Macaca_mulatta.Mmul_10.112.ss --exon Macaca_mulatta.Mmul_10.112.exons Macaca_mulatta.Mmul_10.dna.toplevel.fa Mmul_10
	gzip Macaca_mulatta.Mmul_10.dna.toplevel.fa
fi
