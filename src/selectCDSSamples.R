# Select CDS mapping samples

#These are samples for which there is not a fully annotated genome assembly. In
#these cases, we map to just the ZFX/ZFY CDS or mRNA. The FASTA sequence is in 
# cds/fasta.
library(tidyverse)

#### Functions ####
CDS.MAPPING.FILE <- "cds/cds.samples.csv"

if(file.exists(CDS.MAPPING.FILE)) file.remove(CDS.MAPPING.FILE)

.append.cds.data <- function(data){
  data <- data |>
    
    # Adjust the library field to note long read sequencing distinct from short
    # read paired
    dplyr::mutate(,
                  LibraryLayout = ifelse(Platform %in% c("PACBIO_SMRT"), "LONG_READ", LibraryLayout)) |>
    dplyr::select(Run, CommonName, Genome, LibraryLayout,
                  BioProject, DevStage, Sex, Timepoint, Tissue, Organism, LibrarySelection, 
                  LibrarySource, Bases
    )
  
  if(file.exists(CDS.MAPPING.FILE)){
    write.table(data, file = CDS.MAPPING.FILE, row.names = FALSE, quote = TRUE, 
                append = TRUE, sep = ",", col.names = FALSE)
  } else {
    write.table(data, file = CDS.MAPPING.FILE, row.names = FALSE, quote = TRUE, 
                append = FALSE, sep = ",", col.names = TRUE)
  }
  data
}
#### Select samples ####

##### Tokudaia #####
read.csv("metadata/Tokudaia_osimensis.csv") |>
  dplyr::filter(tissue!="Tail") |>
  dplyr::mutate(
    Sex = str_replace_all(sex, " ", "-"),
    Tissue = tissue,
    Timepoint = AGE,
    DevStage = dev_stage,
    CommonName = "tokudaia_osimensis",
    Genome = "Tokudaia_osimensis_CDS"
  )|>
  .append.cds.data()

##### Cat #####
read.csv("metadata/cat.testis.csv")|>
  dplyr::filter(BioProject=="PRJEB33381") |>
  dplyr::mutate(
    Sex = sex,
    Tissue = tissue,
    Timepoint = Developmental_Stage,
    DevStage = AGE,
    CommonName = "cat",
    Genome = "Felis_catus_mRNA"
  )|>
  .append.cds.data()

##### Dog #####

read.csv("metadata/dog.testis.csv")|>
  dplyr::filter(BioProject=="PRJEB33381") |>
  dplyr::mutate(
    Sex = sex,
    Tissue = tissue,
    Timepoint = Developmental_Stage,
    DevStage = AGE,
    CommonName = "dog",
    Genome = "Canis_lupus_mRNA"
  )|>
  .append.cds.data()


##### Peromyscus #####
read.csv("metadata/PRJNA1068126.csv")|>
  dplyr::filter(tissue=="testis") |>
  dplyr::mutate(
    Sex = sex,
    Tissue = tissue,
    Timepoint = dev_stage,
    DevStage = dev_stage,
    CommonName = "peromyscus",
    Genome = "Peromyscus_maniculatus_CDS"
  )|>
  .append.cds.data()
