# Assemble transcripts usign StringTie
library(fs)
source("src/functions.R")


fs::dir_create("data/stringtie")
fs::file_delete("data/stringtie/ratios.txt")

# Get the distinct groups
sample.groups <- merge(make.sample.groups(), GENE.LOCATIONS, by="CommonName") 

assemble.transcript <- function(common.name, tissue, gene.id, coordinates, merged.bam.file, gtf.file){
  
  cat("Assembling transcripts from", common.name, tissue, gene.id, "\n")
  
  bam.file <- paste0("data/stringtie/", common.name, ".", tissue, ".", gene.id, ".gtf")
  out.file <- paste0("data/stringtie/", common.name, ".", tissue, ".", gene.id, ".gtf")
  
  system2("samtools", paste0("view -o ", bam.file, " ", merged.bam.file, " '", coordinates, "'"))
  
  system2("~/bin/stringtie-2.2.3.Linux_x86_64/stringtie", paste("-o ", out.file, 
                                                                "-p 1 -l", common.name, 
                                                                "-G", gtf.file,
                                                                "-f 0.01",
                                                                bam.file))
  
  # Count the reads on sense and antisense strands
  system2("bash", paste("src/countStrandRatio.sh", bam.file))

}

mapply(assemble.transcript, sample.groups$CommonName, sample.groups$Organism_part, 
       sample.groups$EnsemblId, sample.groups$FlankedLocations, 
       sample.groups$merged.bam, sample.groups$GTF)

