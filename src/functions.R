library(parallel)
library(xlsx)
library(tidyverse)
library(GenomicRanges)

TIME.ORDER <- factor(c("birth", "mid-meiosis", "adult", "Day_00-06",  "Day_07-13", "Day_14-20", "Day_21-27"), 
                     levels = c("birth", "mid-meiosis", "adult","Day_00-06",  "Day_07-13", "Day_14-20", "Day_21-27"))

#### Common functions ####

save.double.width <- function(filename, plot, width=170, height=170){
  ggsave(filename, plot, units = "mm", height = height, width = width, dpi = 300)
}

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
  genomes <- matrix(c(
    "human",      "Homo_sapiens.GRCh38.112.gtf",                       "GRCh38",        "Homo sapiens",
    "mouse",      "Mus_musculus.GRCm39.112.gtf",                       "GRCm39",        "Mus musculus",
    "chicken",    "Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.112.gtf", "GRCg7b",        "Gallus gallus",
    "rat",        "Rattus_norvegicus.mRatBN7.2.112.gtf",               "mRatBN7.2",     "Rattus norvegicus",
    "zebrafinch", "Taeniopygia_guttata.bTaeGut1_v1.p.112.gtf",         "bTaeGut1_v1.p", "Taeniopygia guttata",
    "xenopus",    "Xenopus_tropicalis.UCB_Xtro_10.0.112.gtf",          "UCB_Xtro_10.0", "Xenopus tropicalis",
    "anole",      "Anolis_carolinensis.AnoCar2.0v2.112.gtf",           "AnoCar2.0v2",   "Anolis carolinensis",
    "zebrafish",  "Danio_rerio.GRCz11.112.gtf",                        "GRCz11",        "Danio rerio",
    "opossum",    "Monodelphis_domestica.ASM229v1.112.gtf",            "ASM229v1",      "Monodelphis domestica",
    "platypus",   "Ornithorhynchus_anatinus.mOrnAna1.p.v1.112.gtf",    "mOrnAna1.p.v1", "Ornithorhynchus anatinus",
    "macaque",    "Macaca_mulatta.Mmul_10.112.gtf",                    "Mmul_10",       "Macaca mulatta"),
                    byrow = TRUE, ncol = 4 )
  colnames(genomes) <- c("CommonName", "GTF", "Genome", "Species")
  genomes
}

# Identify the coordinates of a given gene id from GTF. Expand by size on each flank if desired
get.gene.coordinates <- function(common.name, gene.id, gtf.data, size=1000){
  cat("Finding gene coordinates for", common.name, gene.id, "\n")
  # print(str(gtf.data))
  gtf <- gtf.data[[common.name]]
  # cat("Getting gtf coordinates for", gene.id, "\n")
  data <- gtf[gtf$gene_id==gene.id]
  # Expand by size on each flank
  if(size>0) data <- GenomicRanges::resize(data, size*2, fix = "center")
  paste0(unique(GenomicRanges::seqnames(data)), ":", min(GenomicRanges::start(data)), "-", max(GenomicRanges::end(data)))
}

# Get the identifiers for genes of interest
get.gene.locations <- function(gtf.data){
  cat("Finding gene ids\n")
  # Key gene ids from Ensembl
  zfx.y.locations <- matrix(c("chicken", "ZFX",  "ENSGALG00010003052", "ENSGALT00010007119",
                              "zebrafinch", "ZFX", "ENSTGUG00000007219" ,"ENSTGUT00000021043",
                              "opossum", "ZFX",  "ENSMODG00000007512", "ENSMODT00000009508",
                              "xenopus", "ZFX",  "ENSXETG00000007785", "ENSXETT00000017004",
                              "anole",   "ZFX",  "ENSACAG00000007227", "ENSACAT00000007264",
                              "zebrafish", "ZFX", "ENSDARG00000074453", "ENSDART00000110652",
                              "platypus", "ZFX", "ENSOANG00000046710", "ENSOANT00000073933",
                              "mouse",   "Zfx",  "ENSMUSG00000079509", "ENSMUST00000088102",
                              "mouse",   "Zfy1", "ENSMUSG00000053211", "ENSMUST00000189888",
                              "mouse",   "Zfy2", "ENSMUSG00000000103", "ENSMUST00000187148",

                              "human",   "ZFX",  "ENSG00000005889",    "ENST00000304543",
                              "human",   "ZFY",  "ENSG00000067646",    "ENST00000155093",
                              "macaque", "ZFX",  "ENSMMUG00000009801", "ENSMMUT00000013690",
                              "macaque", "ZFY",  "ENSMMUG00000046378", "ENSMMUT00000057467",

                              "rat",     "Zfx",  "ENSRNOG00000005624", "ENSRNOT00000076613",
                              "rat",     "Zfy2", "ENSRNOG00000053042", "ENSRNOT00000077708"), 
                            byrow = TRUE, ncol = 4)
  colnames(zfx.y.locations) <- c("CommonName", "Gene","EnsemblId", "CanonicalTranscript")
  
  zfx.y.locations <- as.data.frame(zfx.y.locations)
  
  cat("Merging gene coordinates\n")
  # Add gene locations from the GTF files
  zfx.y.locations$FlankedLocations <- mapply(get.gene.coordinates, 
                                             common.name = zfx.y.locations$CommonName, 
                                             gene.id     = zfx.y.locations$EnsemblId, 
                                             MoreArgs    = list(gtf.data=gtf.data,
                                                                size=1000), # ensure flanking lncRNAs will be detected
                                             SIMPLIFY = TRUE)
  
  
  zfx.y.locations
}


# Annotate which exons contain interesting features for labelling plots
get.annotated.exons <- function(){
  features <- matrix(c(
                       # Mouse Zfy1
                       "ENSMUSE00001036800", "1", "",
                       "ENSMUSE00001015775", "2", "SP",
                       "ENSMUSE00000993470", "3", "",
                       "ENSMUSE00000568758", "4", "",
                       "ENSMUSE00001050190", "5", "",
                       "ENSMUSE00001068922", "6", "",
                       "ENSMUSE00001324629", "7", "DBD",
                       
                       # Mouse Zfy2
                       "ENSMUSE00000992753", "1", "",
                       "ENSMUSE00000984092", "2", "SP",
                       "ENSMUSE00001047922", "3", "",
                       "ENSMUSE00000704557", "4", "",
                       "ENSMUSE00001005853", "5", "",
                       "ENSMUSE00001089521", "6", "",
                       "ENSMUSE00001334440", "7", "DBD",
                       
                       # Mouse Zfx
                       "ENSMUSE00001269295", "1", "",
                       "ENSMUSE00000149474", "2", "",
                       "ENSMUSE00000149476", "3", "",
                       "ENSMUSE00000149475", "4", "",
                       "ENSMUSE00000477715", "5", "",
                       "ENSMUSE00000149469", "6", "",
                       "ENSMUSE00000744375", "7", "DBD",

                       # Chicken ZFX
                       "ENSGALE00010028768", "1", "",
                       "ENSGALE00010028778", "2", "SP",
                       "ENSGALE00010028781", "3", "",
                       "ENSGALE00010028783", "4", "",
                       "ENSGALE00010028784", "5", "",
                       "ENSGALE00010028785", "6", "",
                       "ENSGALE00010028786", "7", "DBD",
                       
                       # Opossum ZFX
                       "ENSMODE00000078848", "1", "",
                       "ENSMODE00000078849", "2", "SP",
                       "ENSMODE00000078850", "3", "",
                       "ENSMODE00000078860", "4", "",
                       "ENSMODE00000078869", "5", "",
                       "ENSMODE00000307185", "6", "",
                       "ENSMODE00000373268", "7", "DBD",
                       
                       # Human ZFX
                       "ENSE00002688681", "1", "",
                       "ENSE00001176296", "2", "SP",
                       "ENSE00003592125", "3", "",
                       "ENSE00003471331", "4", "",
                       "ENSE00001598623", "5", "",
                       "ENSE00002732075", "6", "",
                       "ENSE00001708883", "7", "DBD",
                       
                       # Human ZFY
                       "ENSE00003889480", "1", "",
                       "ENSE00003895848", "2", "SP",
                       "ENSE00003764421", "3", "",
                       "ENSE00003768468", "4", "",
                       "ENSE00003889859", "5", "",
                       "ENSE00003891660", "6", "",
                       "ENSE00003895708", "7", "DBD",
                       
                       # Rat Zfx
                       "ENSRNOE00000647045", "1", "",
                       "ENSRNOE00000599588", "2", "",
                       "ENSRNOE00000054034", "3", "SP",
                       "ENSRNOE00000613060", "4", "",
                       "ENSRNOE00000053082", "5", "",
                       "ENSRNOE00000053178", "6", "",
                       "ENSRNOE00000296274", "7", "",
                       "ENSRNOE00000514747", "8", "DBD",
                       
                       # Rat Zfy2
                       "ENSRNOE00000541567", "1", "",
                       "ENSRNOE00000544130", "2", "SP",
                       "ENSRNOE00000537273", "3", "",
                       "ENSRNOE00000546137", "4", "",
                       "ENSRNOE00000568193", "5", "",
                       "ENSRNOE00000517140", "6", "",
                       "ENSRNOE00000539578", "7", "DBD",
                       
                       # Macaque ZFX
                       "ENSMMUE00000095158", "1", "",
                       "ENSMMUE00000413272", "2", "",
                       "ENSMMUE00000388194", "3", "SP",
                       "ENSMMUE00000095132", "4", "",
                       "ENSMMUE00000095135", "5", "",
                       "ENSMMUE00000095138", "6", "",
                       "ENSMMUE00000095142", "7", "",
                       "ENSMMUE00000095147", "8", "DBD",
                       
                       # Macaque ZFY
                       "ENSMMUE00000393097", "1", "",
                       "ENSMMUE00000415995", "2", "SP",
                       "ENSMMUE00000337864", "3", "",
                       "ENSMMUE00000095156", "4", "",
                       "ENSMMUE00000355682", "5", "",
                       "ENSMMUE00000095161", "6", "",
                       "ENSMMUE00000407322", "7", "DBD",
                       
                       # Platypus ZFX
                       "ENSOANE00000261828", "1", "",
                       "ENSOANE00000124688", "2", "SP",
                       "ENSOANE00000124690", "3", "",
                       "ENSOANE00000124691", "4", "",
                       "ENSOANE00000124692", "5", "",
                       "ENSOANE00000124693", "6", "",
                       "ENSOANE00000249399", "7", "DBD",
                       
                       # Zebra finch ZFX
                       "ENSTGUE00000074604", "1", "",
                       "ENSTGUE00000074612", "2", "SP",
                       "ENSTGUE00000074659", "3", "",
                       "ENSTGUEE00000216832", "4", "",
                       "ENSTGUE00000074731", "5", "",
                       "ENSTGUE00000074763", "6", "",
                       "ENSTGUEE00000199619", "7", "DBD",
                       
                       # Xenopus ZFX
                       "ENSXETE00000095711", "1", "",
                       "ENSXETE00000095712", "2", "SP",
                       "ENSXETE00000095713", "3", "",
                       "ENSXETE00000095714", "4", "",
                       "ENSXETE00000095715", "5", "",
                       "ENSXETE00000095716", "6", "",
                       "ENSXETE00000570894", "7", "DBD",
                       
                       # Anole ZFX
                       "ENSACAE00000067344", "1", "",
                       "ENSACAE00000067407", "2", "SP",
                       "ENSACAE00000067465", "3", "",
                       "ENSACAE00000067528", "4", "",
                       "ENSACAE00000067604", "5", "",
                       "ENSACAE00000067676", "6", "",
                       "ENSACAE00000067773", "7", "DBD",
                       
                       # Zebrafish ZFX
                       "ENSDARE00001112692", "1", "",
                       "ENSDARE00000798648", "2", "SP",
                       "ENSDARE00000808537", "3", "",
                       "ENSDARE00000846717", "4", "",
                       "ENSDARE00000822128", "5", "",
                       "ENSDARE00000794351", "6", "",
                       "ENSDARE00001272363", "7", "DBD"
                       
                       ),
                     
                     byrow=TRUE, ncol = 3)
  colnames(features) <- c("ExonId", "CodingExonNumber", "Feature")
  features
}

# Read GTF files. Parallel on Unix.
# Reads all GTF files from ./genomes
read.gtf.data <- function(){
  genome.data <- get.genome.data()
  cat("Reading GTF files\n")
  gtf.data <- mclapply(genome.data[,2], \(f){ 
    cat("Reading GTF file", f, "\n")
    rtracklayer::import( paste0("genomes/", f))
  }, mc.cores = ifelse(installr::is.windows(), 1, 5))
  names(gtf.data) <- genome.data[,1]
  cat("Read GTF files\n")
  gtf.data
}

# Read all selected samples from ./metadata 
# i.e. all files with .filt. in the name
read.selected.samples <- function(){
  cat("Reading selected samples\n")
  # Make a factor of times to allow ordering of plots
  do.call(rbind, lapply(list.files(path="metadata", pattern = "*.filt.csv", full.names = TRUE), read.csv)) %>%
    dplyr::mutate(Timepoint = factor(Timepoint, levels = TIME.ORDER)) %>%
    dplyr::arrange(CommonName, Run)
  
}


# Read the metadata to find samples. Aggregate to groups based on tissue type
# and note which samples still need processing
make.sample.groups <- function(){
  # Create command to merge bams in groups
  SELECTED.SAMPLES %>% 
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

#### Global variables ####

cat("Making global variables\n")

# Read all GTF files
GTF.DATA <- read.gtf.data()

# Global data frame with gene ids for all species
GENE.LOCATIONS <- merge(get.gene.locations(GTF.DATA), get.genome.data(), by="CommonName")

# Read the filtered samples, match folder names
SELECTED.SAMPLES <- read.selected.samples()

cat("Common functions and global variables loaded\n")