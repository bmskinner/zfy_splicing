library(gridExtra)
library(png)
library(grid)
source("src/functions.R")

png.files <- data.frame(path = list.files(path = "report/raw_sashimi", pattern = "[SDE]RR.*.png$", 
                                          full.names = TRUE, recursive=TRUE)) |>
  dplyr::mutate(file = basename(path)) |>
  tidyr::separate_wider_delim(file,
                              delim = ".", names = c("Run", "GeneId", "condensed", "ext")
  ) |>
  merge(GENE.LOCATIONS, by = "GeneId") |>
  merge(GENOME.DATA, by=c("CommonName", "GTF_FILE")) |>
  merge(SELECTED.SAMPLES, by = c("CommonName", "Run",  "Genome", "GTF_FILE"))

unique.species <- sort(unique(png.files$CommonName))
tissues <- sort(unique(png.files$Tissue))
timepoints <- sort(unique(png.files$Timepoint))
genes <- sort(unique(png.files$Gene))
sexes <- sort(unique(png.files$Sex))

for(species in unique.species){
  for (tissue in tissues) {
    for (timepoint in timepoints) {
      for (gene in genes) {
        for(sex in sexes){
          
          out.file <- paste0("report/tissues/",species, ".", sex, ".", tissue, ".", timepoint, ".", gene , ".pdf")
          if(file.exists(out.file)) next
          
          srrs <- png.files |> 
            dplyr::filter(CommonName==species, Tissue==tissue,
                          Timepoint==timepoint, Gene==gene, Sex==sex) |>
            dplyr::select(path) |>
            dplyr::pull()
          
          srrs <- srrs[file.size(srrs)>0] # there are some placeholder images of zero size when no reads were in a sample
          
          if(length(srrs)==0) next

          plots <- lapply(srrs, \(x){
            img <- as.raster(readPNG(x))
            rasterGrob(img, interpolate = FALSE, default.units = "mm", width = 170, height = 50)
          })
          
          ggsave(
            out.file,
            arrangeGrob(grobs = plots, nrow = length(plots), ncol = 1),
            dpi = 300, units = "mm", width = 170, height = min(1000, 50 * length(plots))
          )
        }
      }
    }
  }
}
