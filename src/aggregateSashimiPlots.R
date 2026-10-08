source("src/functions.R")

#### Individual runs ####
png.files <- data.frame(path = list.files(path = "report/raw_sashimi", pattern = "[SDE]RR.*.png$", 
                                          full.names = TRUE, recursive=TRUE)) |>
  dplyr::mutate(file = basename(path)) |>
  tidyr::separate_wider_delim(file,
                              delim = ".", names = c("Run", "GeneId", "condensed", "ext")
  ) |>
  merge(SELECTED.SAMPLES, by = c("Run", "GeneId"))



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
          
          out.file <- paste0("report/tissues/",species, ".", sex, ".", tissue, ".", timepoint, ".", gene , ".png")
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

#### Merged ####

png.files <- data.frame(path = list.files(path = "report/merged_sashimi", pattern = ".*.png$", 
                                          full.names = TRUE, recursive=TRUE)) |>
  dplyr::mutate(file = basename(path)) |>
  tidyr::separate_wider_delim(file,
                              delim = ".", names = c("CommonName", "Tissue", "Timepoint", "Sex", "GeneId", "Gene","condensed", "ext")
  )

unique.species <- sort(unique(png.files$CommonName))
tissues <- sort(unique(png.files$Tissue))
timepoints <- sort(unique(png.files$Timepoint))
genes <- sort(unique(png.files$Gene))
sexes <- sort(unique(png.files$Sex))

for(species in unique.species){
  
  out.file <- paste0("report/species/",species, ".png")
  if(file.exists(out.file)) next
  
  pngs <- png.files |> 
    dplyr::filter(CommonName==species) |>
    dplyr::select(path) |>
    dplyr::pull()
  
  pngs <- pngs[file.size(pngs)>0] # there are some placeholder images of zero size when no reads were in a sample
  
  if(length(pngs)==0) next
  
  plots <- lapply(pngs, \(x){
    img <- as.raster(readPNG(x))
    rasterGrob(img, interpolate = FALSE, default.units = "mm", width = 170, height = 50)
  })
  
  ggsave(
    out.file,
    arrangeGrob(grobs = plots, nrow = length(plots), ncol = 1),
    dpi = 300, units = "mm", width = 170, height = min(1200, 50 * length(plots))
  )
}

for (tissue in tissues) {
  out.file <- paste0("report/tissues/", tissue, ".png")
  if(file.exists(out.file)) next
  
  pngs <- png.files |> 
    dplyr::filter(Tissue==tissue) |>
    dplyr::select(path) |>
    dplyr::pull()
  
  pngs <- pngs[file.size(pngs)>0] # there are some placeholder images of zero size when no reads were in a sample
  
  if(length(pngs)==0) next
  
  plots <- lapply(pngs, \(x){
    img <- as.raster(readPNG(x))
    rasterGrob(img, interpolate = FALSE, default.units = "mm", width = 170, height = 50)
  })
  
  ggsave(
    out.file,
    arrangeGrob(grobs = plots, nrow = length(plots), ncol = 1),
    dpi = 300, units = "mm", width = 170, height = min(1200, 50 * length(plots))
  )
}
