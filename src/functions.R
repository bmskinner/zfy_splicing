library(parallel)
library(xlsx)

TIME.ORDER <- factor(c("birth", "mid-meiosis", "adult"), levels = c("birth", "mid-meiosis", "adult"))

# Write the given data frame to an Excel file
create.xlsx = function(data, file.name){
  
  data <- as.data.frame(data) # ensure not a tibble
  
  oldOpt = options()
  options(xlsx.date.format="yyyy-mm-dd") # change date format
  wb = xlsx::createWorkbook(type = "xlsx")
  sh = xlsx::createSheet(wb)
  xlsx::addDataFrame(data, sh, row.names = F)
  xlsx::createFreezePane(sh, 2, 2, 2, 2) # freeze top row and first column
  xlsx::autoSizeColumn(sh, 1:ncol(data))
  xlsx::saveWorkbook(wb, file=file.name)
  options(oldOpt)
}

# Get the names of GTF files for a genome
get.genome.data <- function(){
  genomes <- matrix(c("chicken", "Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.112.gtf",
                      "opossum", "Monodelphis_domestica.ASM229v1.112.gtf",
                      "mouse", "Mus_musculus.GRCm39.112.gtf",
                      "human", "Homo_sapiens.GRCh38.112.gtf",
                      "rabbit", "Oryctolagus_cuniculus.OryCun2.0.112.gtf",
                      "rat", "Rattus_norvegicus.mRatBN7.2.112.gtf",
                      "macaque", "Macaca_mulatta.Mmul_10.112.gtf"),
                    byrow = TRUE, ncol = 2 )
  colnames(genomes) <- c("CommonName", "GTF")
  genomes
}

# Get the identifiers for genes of interest
get.gene.locations <- function(){
  # Manual locations from Ensembl
  zfx.y.locations <- matrix(c("chicken", "ZFX",  "ENSGALG00010003052", "ENSGALT00010007119",
                              "opossum", "ZFX",  "ENSMODG00000007512", "ENSMODT00000009508",
                              "mouse",   "Zfx",  "ENSMUSG00000079509", "ENSMUST00000088102",
                              "mouse",   "Zfy1", "ENSMUSG00000053211", "ENSMUST00000189888",
                              "mouse",   "Zfy2", "ENSMUSG00000000103", "ENSMUST00000187148",
                              "mouse",   "Zfa", "ENSMUSG00000121690",  "ENSMUST00000180673",
                              "human",   "ZFX",  "ENSG00000005889",    "ENST00000304543",
                              "human",   "ZFY",  "ENSG00000067646",    "ENST00000155093",
                              "macaque", "ZFX",  "ENSMMUG00000009801", "ENSMMUT00000013690",
                              "macaque", "ZFY",  "ENSMMUG00000046378", "ENSMMUT00000057467",
                              "rabbit",  "ZFX",  "ENSOCUG00000003815", "ENSOCUT00000003815",
                              "rat",     "Zfx",  "ENSRNOG00000005624", "ENSRNOT00000076613",
                              "rat",     "Zfy2", "ENSRNOG00000053042", "ENSRNOT00000077708"), 
                            byrow = TRUE, ncol = 4)
  colnames(zfx.y.locations) <- c("CommonName", "Gene","EnsemblId", "CanonicalTranscript")
  zfx.y.locations
}

# Global data frame with gene ids for all species
GENE.LOCATIONS <- merge(get.gene.locations(), get.genome.data(), by="CommonName")

# Read all GTF files. Parallel on Unix.
read.gtf.data <- function(){
  genome.data <- get.genome.data()
  cat("Reading GTF files\n")
  gtf.data <- mclapply(genome.data[,2], \(f){ 
    cat("Reading GTF file", f, "\n")
    rtracklayer::import( paste0("genomes/", f))
  }, mc.cores = ifelse(installr::is.windows(), 1, 4))
  names(gtf.data) <- genome.data[,1]
  cat("Read GTF files\n")
  gtf.data
}

# Read the metadata to find samples. Aggregate to groups based on tissue type
# and note which samples still need processing
make.sample.groups <- function(){
  filtered.samples <- do.call(rbind, lapply(list.files(path="metadata", pattern = "*.filt.csv", full.names = TRUE), 
                                            \(f) read.csv(f) %>% dplyr::mutate(CommonName = str_replace(str_replace(f, "metadata/", ""), ".filt.csv", ""))))
  
  # Create command to merge bams in groups
  filtered.samples %>% 
    dplyr::group_by(Organism, Organism_part, Timepoint, CommonName) %>% # not by sex - no difference seen in first pass
    dplyr::mutate(bam.file = paste0("data/", CommonName, "/", Run, ".bam")) %>%
    dplyr::summarise(bams = paste(bam.file, collapse = " "),
                     all.bams.present = all(file.exists(bam.file)),
                     count = n()) %>%
    dplyr::mutate(merged.bam = paste0("data/merged/", CommonName, ".", Organism_part, ".", Timepoint, ".bam"),
                  samtools.merge.arguments = paste("merge -@ 7 -r -o", merged.bam,  bams))  %>%
    dplyr::mutate(merged.bam.exists = file.exists(merged.bam),
                  index.exists = file.exists(paste0(merged.bam, ".csi"))) 
}
