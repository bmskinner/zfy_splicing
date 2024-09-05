# Make sashimi plot based on ggsashimi output.
# ggsashimi.py modified to write data objects to Rds when run
# This custom sashimi plot ensures the transcripts are always left to right 
# irrespective of strand
library(tidyverse)
library(data.table)
library(patchwork)
library(grid)
library(scales)
library(fs)
source("src/functions.R")

# Ensure output dirs exist
fs::dir_create(c("report/species", "report/timepoints", "report/tissues"))
gtf.data <- read.gtf.data()

#### Main functions #### 

is.reverse.strand <- function(transcript.id, gtf.data){
  exons <- gtf.data[gtf.data$type=="exon" & gtf.data$transcript_id==transcript.id]
  any(GenomicRanges::strand(exons)=="-")
}

# Given a canonical transcript id, find the splice junctions
get.canonical.junctions <- function(canonical.transcript.id, gtf.data){
  exons <- gtf.data[gtf.data$type=="exon" & gtf.data$transcript_id==canonical.transcript.id]

  if( is.reverse.strand(canonical.transcript.id, gtf.data) ) {
    return( na.omit(data.frame(
      j1 = lead(GenomicRanges::end(exons))+1, # end of intron is start of next exon (end coordinate rev strand)
      j2 = GenomicRanges::start(exons)        # first base of intron is start of exon coordinate (rev strand)
    )))  
  }
  
  data.frame(
    j1 = GenomicRanges::end(exons)+1,         # start of the intron is the end of the exon
    j2 = lead(GenomicRanges::start(exons))) # end of the intron is the start of the next exon
}

# given a splice junction, determine if this is canonical. i.e. is the junction
# found in the canonical transcript?
# is.junction.canonical("ENSGALT00010007119",chicken.gtf, 118299992 ,118301364)
is.junction.canonical <- function(canonical.transcript.id, gtf.data, start, end){
  junctions <- get.canonical.junctions(canonical.transcript.id,gtf.data)
  any(junctions$j1==start & junctions$j2==end)
}

# Create the panel of transcripts
make.gene.track <- function(sashimi.data){
  # cat("Making gene track\n")
  data <- sashimi.data$density_list
  junctions <- sashimi.data$junction_list
  anns <- sashimi.data$ann_list
  
  # Set coordinates for the x axis, reversing if on reverse strand
  xmin <- ifelse(sashimi.data$is.reverse.strand, max(data$x), min(data$x))
  xmax <- ifelse(sashimi.data$is.reverse.strand, min(data$x), max(data$x))
  
  # Create the arrows to draw on introns
  make.tx.arrows <- function(){
    # Create data table for strand arrows
    txarrows = data.table()
    introns = sashimi.data$ann_list[['introns']]
    # Add right-pointing arrows for plus strand
    if ("+" %in% introns$strand && nrow(introns[strand=="+" & end-start>5, ]) > 0) {
      txarrows = rbind(
        txarrows,
        introns[strand=="+" & end-start>5, list(
          seq(start+4,end,by=453)-1,
          seq(start+4,end,by=453)
        ), by=.(tx,start,end)
        ]
      )
    }
    # Add left-pointing arrows for minus strand
    if ("-" %in% introns$strand && nrow(introns[strand=="-" & end-start>5, ]) > 0) {
      txarrows = rbind(
        txarrows,
        introns[strand=="-" & end-start>5, list(
          seq(start,max(start+1, end-4), by=453),
          seq(start,max(start+1, end-4), by=453)-1
        ), by=.(tx,start,end)
        ]
      )
    }
    txarrows
  }
  
  # Make the gene track
  ggplot()+
    # Introns and arrows
    geom_segment(data=sashimi.data$ann_list$introns, aes(x=start, xend=end, y=tx, yend=tx), linewidth=0.3)+
    # geom_segment(data=make.tx.arrows(), aes(x=V1,xend=V2,y=tx,yend=tx), arrow=arrow(length=unit(0.02,"npc")))+
    
    # Exons
    geom_segment(data=sashimi.data$ann_list$exons, aes(x=start, xend=end, y=tx, yend=tx, col=tx==sashimi.data$canonical.transcript.id), size=5, alpha=1)+
    scale_color_manual(values = c(`TRUE`="blue", `FALSE`="grey"))+
    scale_y_discrete(expand=c(0,0.5))+
    coord_cartesian(xlim = c(xmin,xmax))+
    scale_x_continuous(expand=c(0,0.25), labels=comma)+
    theme_minimal()+
    theme(axis.line.y = element_blank(),
          axis.line.x = element_line(),
          axis.ticks.x = element_line(),
          axis.title = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid = element_blank(),
          legend.position = "none")
}

# Create a sashimi panel for the given ggsashimi data
# min.spanning.reads - the minimum number of junction-spanning reads needed to display on the chart
# y.label - the data value to be used for the y axis label
make.sashimi.panel <- function(sashimi.data, min.spanning.reads=5, label="tissue"){
  # cat("Making sashimi panel\n")
  data <- sashimi.data$density_list
  junctions <- sashimi.data$junction_list
  anns <- sashimi.data$ann_list

  # Set coordinates for the x axis
  is.minus.strand <- sashimi.data$is.reverse.strand
  xmin <- min(data$x)
  xmax <- max(data$x)
  
  xtmp <- xmin
  xmin <- ifelse(is.minus.strand, xmax, xmin)
  xmax <- ifelse(is.minus.strand, xtmp, xmax)

  ymax <- max(data$y)
  
  # Create the coverage plot
  splot <- ggplot() + 
    geom_bar(data = data, aes(x, y), width=1, position='identity', stat='identity', alpha=0.5)+
    coord_cartesian(ylim = c(0-(ymax), ymax*1.5), 
                    xlim = c(xmin, xmax))+
    scale_x_continuous(expand=c(0,0.25))+
    labs(y = label)
  
  add.junction <- function(splot, xmin, xmax, ymin, ymax, is.even, count, is.canonical){
    
    # Define the spline shapes that make the junction lines
    spline.color <- ifelse(is.canonical, "blue", "black")
    spline.size  <- ifelse(is.canonical, 1, 2)
    
    l.spline.btm <- grid::xsplineGrob(x=c(0, 0, 1, 1), y=c(1, 0, 0, 0), shape=1, gp=gpar(lwd=spline.size, col=spline.color))
    l.spline.top <- grid::xsplineGrob(x=c(0, 0, 1, 1), y=c(0, 1, 1, 1), shape=1, gp=gpar(lwd=spline.size, col=spline.color))
    r.spline.btm <- grid::xsplineGrob(x=c(1, 1, 0, 0), y=c(1, 0, 0, 0), shape=1, gp=gpar(lwd=spline.size, col=spline.color))
    r.spline.top <- grid::xsplineGrob(x=c(1, 1, 0, 0), y=c(0, 1, 1, 1), shape=1, gp=gpar(lwd=spline.size, col=spline.color))
    
    # Determine which splines to use for minus strand versus plus strand transcripts
    l.grob.btm <- l.spline.btm 
    if(is.minus.strand) l.grob.btm <- r.spline.btm
    
    l.grob.top <- l.spline.top
    if(is.minus.strand) l.grob.top <- r.spline.top
    
    r.grob.btm <- r.spline.btm 
    if(is.minus.strand) r.grob.btm <- l.spline.btm
    
    r.grob.top <- r.spline.top
    if(is.minus.strand) r.grob.top <- l.spline.top
    
    xmid <- (xmin + xmax) / 2
    ymid <- ((ymin + ymax) / 2) * 1.25

    # Left arc
    if(is.even){ # Junctions below zero
      splot <- splot+annotation_custom(grob = l.grob.btm, 
                                       xmin = xmin, 
                                       xmax = xmid, 
                                       ymin = -ymid*0.5, 
                                       ymax = 0)
    } else {
      splot <- splot+annotation_custom(grob = l.grob.top, 
                                       xmin = xmin, 
                                       xmax = xmid, 
                                       ymax = ymid, 
                                       ymin = 0)
    }
    
    # Right arc    
    if(is.even){ # Junctions below zero
      splot <- splot+annotation_custom(grob = r.grob.btm, 
                                       xmin = xmid, 
                                       xmax = xmax, 
                                       ymin = -ymid*0.5, 
                                       ymax = 0)
    } else {
      splot <- splot+annotation_custom(grob = r.grob.top, 
                                       xmin = xmid, 
                                       xmax = xmax, 
                                       ymax = ymid, 
                                       ymin = 0)
    }
    
    splot <- splot + annotate("label", x = xmid, 
                              y = ifelse(is.even, -ymid*1.1, ymid*1.2), 
                              label = as.character(count),
                              size=2, col = spline.color, fill=NA, label.size=NA)
    
    return(splot)
  }
  
  if(nrow(junctions)>0){
    junctions <- junctions %>% 
      dplyr::filter(count>=min.spanning.reads)
  }
  
  # Add the junctions, adjusting for plus vs minus strand
  if(nrow(junctions)>0){

    # Calculate charting coordinates
    junctions$total.count <- sum(junctions$count)
    junctions$isEven <- sapply(1:nrow(junctions), \(x) x%%2==0)
    
    for(i in 1:nrow(junctions)){
      splot <- add.junction(splot, junctions[i,]$x, junctions[i,]$xend, junctions[i,]$y, 
                            junctions[i,]$yend, junctions[i,]$isEven, junctions[i,]$count,
                            junctions[i,]$is.canonical)
    }
  }
  
  # Format the final plot
  splot <- splot +
    
    theme_minimal()+
    theme(axis.line = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(angle = 0, hjust = 1, vjust = 0.5),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  
  list(plot=splot, junctions=junctions)
}

# Read a junction file and note if each junction is in the canonical transcript
read.rds.file <- function(rds.file){

  tryCatch({
    rds.data <- readRDS(rds.file)
    
    file.name.parts <- str_split_1(basename(rds.file), "\\.")
    rds.data$species <- file.name.parts[1]
    rds.data$tissue <- file.name.parts[2]
    rds.data$timepoint <- file.name.parts[3]
    rds.data$gene.id <- file.name.parts[4]
    rds.data$gene.name <- GENE.LOCATIONS[GENE.LOCATIONS$EnsemblId==rds.data$gene.id, "Gene"]
    rds.data$filename <- basename(rds.file)
    
    rds.data$density_list <- rds.data$density_list[[1]]
    rds.data$junction_list <- rds.data$junction_list[[1]]
    
    rds.data$canonical.transcript.id <- GENE.LOCATIONS[GENE.LOCATIONS$EnsemblId==rds.data$gene.id, "CanonicalTranscript"]
    rds.data$gtf.data <- gtf.data[[rds.data$species]]
    
    # Is the gene on the forward or reverse strand? Note - this is the gene, not the junctions or reads
    rds.data$is.reverse.strand <- is.reverse.strand(rds.data$canonical.transcript.id, rds.data$gtf.data)
    rds.data$strand <- ifelse(rds.data$is.reverse.strand, "-", "+")
    
    rds.data$junction_list$is.canonical <- mapply(is.junction.canonical, 
                                                  start = rds.data$junction_list$x, 
                                                  end   = rds.data$junction_list$xend,
                                                  MoreArgs = list(canonical.transcript.id = rds.data$canonical.transcript.id,
                                                                  gtf.data      = rds.data$gtf.data))
    
    
    return(rds.data)
  }, error=\(e) {cat("Error making data from", rds.file,"\n", paste(e)); e})
}

# Create a sashimi panel plot for all tissues of the given species and timepoint
make.species.panel <- function(species, timepoint, gene.id){
  data.files <- list.files(path = "data/merged", pattern = paste0(species, ".*\\.", timepoint, "\\.", gene.id, "\\.Rds_*"), full.names = TRUE)
  if(length(data.files)==0) return()
  data <- lapply(data.files, read.rds.file)
  
  gene.name <- GENE.LOCATIONS[GENE.LOCATIONS$EnsemblId==gene.id,]$Gene
  
  out.png.file <- paste0("report/species/", species, ".", timepoint, ".", gene.id, ".", gene.name, ".png")
  
  # if(!file.exists(out.png.file)){
    plots <- lapply(data, \(x)  make.sashimi.panel(x,label=paste0(x$tissue))$plot)
    track <- make.gene.track(data[[1]]) # only one gene, only need one track
    plots[[length(plots)+1]] <- track
    
    patchwork::wrap_plots(plots, nrow = length(plots))
    ggsave(plot = last_plot(), filename =out.png.file, dpi = 300, units = "mm", width = 170, height = 240)
  # }
}

# Create a sashimi panel plot for all species of the given tissue and timepoint
make.tissue.panel <- function(tissue, timepoint){
  cat("Making", tissue, "at", timepoint, "\n")
  data.files <- list.files(path = "data/merged", pattern = paste0(".*\\.", tissue, "\\.", timepoint, "\\..*.Rds_*"), full.names = TRUE)
  if(length(data.files)==0) return()
  data <- lapply(data.files, read.rds.file)
  
  out.png.file <- paste0("report/tissues/", tissue, ".", timepoint, ".png")
  
  plots <- lapply(data, \(x) make.sashimi.panel(x, label= paste0(x$species, "\n", x$gene.name))$plot)
  tracks <- lapply(data, make.gene.track)
  plots <- c(rbind(plots, tracks))
  
  patchwork::wrap_plots(plots, nrow = length(plots))
  ggsave(plot = last_plot(), filename = out.png.file, dpi = 300, units = "mm", width = 170, height = 240)
}

# Create a sashimi panel plot for all timepoint of the given tissue and species
make.timepoint.panel <- function(species, tissue, gene.id){
  data.files <- list.files(path = "data/merged", pattern = paste0(species, ".*\\.", tissue, "\\..*", gene.id, ".*.Rds_*"), full.names = TRUE)
  if(length(data.files)==0) return()
  
  # Ensure files are plotted in time order
  ordered.files <- list()
  for(t in TIME.ORDER){
    for(f in data.files){
      if(str_detect(f, t)){
        ordered.files <- c(ordered.files, f)
      }
    }
  }
  
  data <- lapply(ordered.files, read.rds.file)
  
  gene.name <- GENE.LOCATIONS[GENE.LOCATIONS$EnsemblId==gene.id,]$Gene
  
  out.png.file <-  paste0("report/timepoints/", species, ".", tissue, ".", gene.id, ".", gene.name, ".png")

  plots <- lapply(data, \(x) make.sashimi.panel(x, label=paste0(x$species, " ", gene.name, "\n",  x$timepoint))$plot)
  track <- make.gene.track(data[[1]])
  plots[[length(plots)+1]] <- track
  
  patchwork::wrap_plots(plots, nrow = length(plots))
  ggsave(plot = last_plot(), filename =out.png.file, dpi = 300, units = "mm", width = 170, height = 240)
}

#### Select groups for plotting ####

# find all species combinations for plotting
all.samples <- merge(make.sample.groups(), GENE.LOCATIONS, by="CommonName") 

species.groups <- all.samples %>% 
  dplyr::group_by(CommonName, Timepoint, EnsemblId) %>%
  dplyr::summarise(SampleCount = n())

tissue.groups <- all.samples %>% 
  dplyr::group_by(Organism_part, Timepoint) %>%
  dplyr::summarise(SampleCount = n())


timepoint.groups <- all.samples %>% 
  dplyr::group_by(CommonName, Organism_part, EnsemblId) %>%
  dplyr::summarise(SampleCount = n())

#### Generate results figures ####

# Make species plots showing variation over tissues as a specific timepoint
mapply(make.species.panel, species.groups$CommonName, species.groups$Timepoint,species.groups$EnsemblId)

# Make tissue plots showing variation over species as a specific timepoint
mapply(make.tissue.panel, tissue.groups$Organism_part, tissue.groups$Timepoint)

# Make timepoint plots showing variation over times in a specific tissue
mapply(make.timepoint.panel, timepoint.groups$CommonName, timepoint.groups$Organism_part, timepoint.groups$EnsemblId)
# make.timepoint.panel("chicken", "testis", "ENSGALG00010003052")

#### Condense introns for neater plotting #### 

# Plot junctions directly on exon track
make.gene.track.sashimi.panel <- function(sashimi.data, min.spanning.reads=10, label="tissue", 
                                          show.x.axis=TRUE, is.collapse.introns=FALSE){
  
  # Annotatable exon features
  feature.data <- sashimi.data$gtf.data %>%
    as.data.frame %>%
    dplyr::filter(gene_id==sashimi.data$gene.id & !is.na(exon_id)) %>%
    merge(., get.annotated.exons(), by.x = "exon_id", by.y = "ExonId", all.x=TRUE) %>%
    dplyr::mutate(xmid = (start+end)/2)

  collapse.introns <- function(){
    
    # Calculate offsets to make all introns at most 500bp
    calculate.coordinate.conversion <- function(){
      
      # Ensure we only look at transcripts that are part of the gene of interest
      transcript.ids <- sashimi.data$gtf.data %>%
        as.data.frame %>%
        dplyr::filter(gene_id==sashimi.data$gene.id & !is.na(transcript_id)) %>%
        dplyr::select(transcript_id) %>%
        dplyr::distinct()
      
      # Get reduced ranges for the exons
      exons <- sashimi.data$ann_list$exons %>%
        dplyr::filter(strand==sashimi.data$strand & tx %in% transcript.ids$transcript_id)
      
      introns <- sashimi.data$ann_list$introns %>%
        dplyr::filter(strand==sashimi.data$strand& tx %in% transcript.ids$transcript_id)
      
      # print(introns)
      
      exon.ranges <- GenomicRanges::GRanges(seqnames=rep("test", nrow(exons)), 
                                            ranges = IRanges::IRanges(start=exons$start, 
                                                                      end = exons$end),
                                            strand = exons$strand)
      
      exon.ranges <- GenomicRanges::reduce(exon.ranges)
    
      intron.ranges <- GenomicRanges::GRanges(seqnames=rep("test", nrow(introns)), 
                                            ranges = IRanges::IRanges(start=introns$start, 
                                                                      end = introns$end),
                                            strand = introns$strand)
      
      # print(intron.ranges)
      
      # Break introns apart, since they can be part of an exon for a different transcript
      intron.ranges <- GenomicRanges::disjoin(intron.ranges)
      intron.ranges$overlappingExons <- GenomicRanges::countOverlaps(intron.ranges,exon.ranges, minoverlap = 10)
      
      # print(intron.ranges)
      
      # Remove introns that overlap an exon
      intron.ranges <- intron.ranges[intron.ranges$overlappingExons==0,]
      
      range.start <- min(exons$start)
      
      # May be some introns missing. Fill in gaps that are not covered by intron or exons
      missing.introns <- GenomicRanges::gaps(GenomicRanges::reduce(c(intron.ranges, exon.ranges)),
                                             start = range.start)#,
                                             # end = max(end(exon.ranges)))
      
      
      # print(missing.introns)
      intron.ranges <- c(intron.ranges, missing.introns)
      
      intron.ranges <- intron.ranges %>%
        as.data.frame %>%
        dplyr::select(start, end, strand) %>%
        dplyr::mutate(Type = "intron")
      
      exon.ranges <- exon.ranges %>%
        as.data.frame %>%
        dplyr::select(start, end, strand) %>%
        dplyr::mutate(Type = "exon")
      
      # Calculate how to convert coordinate ranges for introns and exons
      conversion.coords <- rbind(intron.ranges, exon.ranges)
      

      
      conversion.coords <- conversion.coords %>%
        dplyr::arrange(start, end) %>%
        dplyr::distinct() %>%
        dplyr::mutate(length = abs(start-end),
                      # offset = ifelse(length>500 & Type=="intron", 500, 0),
                      new.length = ifelse(length>500 & Type=="intron", 500, length),
                      new.end = min(start)+cumsum(new.length),
                      new.start = new.end - new.length) # how much offset to apply
      # print(conversion.coords)
      conversion.coords
    }
    
    # Determine the offset conversions
    coord.conversion.table <- calculate.coordinate.conversion()

    # Convert a coordinate to collapsed space using the conversion table
    convert.coordinates <- function(coordinate){

      # Take the first matching row
      feature <- coord.conversion.table %>% 
        dplyr::filter(start<=coordinate & end>=coordinate) %>%
        dplyr::slice_head(n=1)
      
      # How far along the feature are we?
      f.feature <- (coordinate-feature$start)/feature$length
      
      # Nearest integer to the same fraction of the new coordinate space
      
      result <- round(f.feature*feature$new.length + feature$new.start)
      unique(result)
    }
    
    # Copy the existing coordinates
    junctions <- sashimi.data$junction_list
    anns <- sashimi.data$ann_list
    
    # Apply offsets to coordinates
    anns$exons <- as.data.frame(anns$exons)
    anns$introns <- as.data.frame(anns$introns)

    anns$exons$old.start <- anns$exons$start
    anns$exons$old.end <- anns$exons$end
    anns$exons$start <- sapply(anns$exons$old.start, convert.coordinates)
    anns$exons$end <- sapply(anns$exons$old.end, convert.coordinates)
    
    anns$introns$old.start <- anns$introns$start
    anns$introns$old.end <- anns$introns$end
    anns$introns$start <- sapply(anns$introns$old.start, convert.coordinates)
    anns$introns$end <- sapply(anns$introns$old.end, convert.coordinates)
    
    # print(anns$introns)
    
    sashimi.data$ann_list <<- anns
    
    junctions$old.x <- junctions$x
    junctions$old.xend <- junctions$xend
    junctions$x <- sapply(junctions$old.x, convert.coordinates)
    junctions$xend <- sapply(junctions$old.xend, convert.coordinates)
    
    sashimi.data$junction_list <<- junctions
    
    feature.data$old.xmid <<- feature.data$xmid
    feature.data$xmid <<- sapply(feature.data$old.xmid, convert.coordinates)
  }
 
  if(is.collapse.introns){
    collapse.introns()
  }

  junctions <- sashimi.data$junction_list
  anns <- sashimi.data$ann_list

  # Set coordinates for the x axis, reversing if on reverse strand
  xmin <- ifelse(sashimi.data$is.reverse.strand, max(anns$exons$end)+500, min(anns$exons$start)-500)
  xmax <- ifelse(sashimi.data$is.reverse.strand, min(anns$exons$start)-500, max(anns$exons$end)+500)
  
  # Canonical exons
  canonical.exons <- sashimi.data$ann_list$exons %>% dplyr::filter(tx==sashimi.data$canonical.transcript.id)
  non.canonical.exons <- sashimi.data$ann_list$exons %>% dplyr::filter(tx!=sashimi.data$canonical.transcript.id)

  # Make the gene track
  splot <- ggplot()+
    # Introns
    geom_segment(data=sashimi.data$ann_list$introns, aes(x=start, xend=end, y=0, yend=0), linewidth=0.3)+
    
    # Exons
    geom_rect(data=non.canonical.exons, aes(xmin=start, xmax=end, ymin=-0.5, ymax=0.5), fill="grey",  alpha=1)+
    
    # Canonical exons
    geom_rect(data=canonical.exons, aes(xmin=start, xmax=end, ymin=-0.5, ymax=0.5), fill="blue", alpha=1)+
    
    scale_y_discrete(expand=c(0.5,0.5))+
    coord_cartesian(xlim = c(xmin,xmax))+
    scale_x_continuous(expand=c(0,0.25), labels=comma)+
    theme_minimal()+
    labs(y = label)+
    theme(axis.line.y = element_blank(),
          axis.line.x = element_line(),
          axis.ticks.x = element_line(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(angle = 0, hjust = 1, vjust = 0.5),
          axis.ticks.y = element_blank(),
          panel.grid = element_blank(),
          legend.position = "none")
  
  # Hide x-axis if needed
  if(!show.x.axis){
    splot <- splot+
      theme(axis.line.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank())
  }
  
  # Add feature annotations of interest
  add.features <- function(splot){
    
    feature.label.y.offset <- 0.75

    splot <- splot +
      # Number the coding exons
      geom_label(data=feature.data, aes(x = xmid, y = 0, label=CodingExonNumber), 
                 size=2,fill=NA, label.size=NA, col="white" )
      # geom_label(data=feature.data, aes(x = xmid, y = feature.label.y.offset, label=Feature), 
      #            size=2,fill=NA, label.size=NA, col="black" )
      
    splot
    
  }
  
  splot <- add.features(splot)
  
  # Add the sashimi splines
  add.junction <- function(splot, xmin, xmax, ymin, ymax, is.even, count, is.canonical, f.length){
    
    # Define the spline shapes that make the junction lines
    spline.color <- ifelse(is.canonical, "blue", "black")
    spline.size  <- ifelse(is.canonical, 1, 2)
    spline.alpha <- ifelse(is.canonical, 0.5, 1)
    
    l.spline.btm <- grid::xsplineGrob(x=c(0, 0, 1, 1), y=c(1, 0, 0, 0), shape=1, gp=gpar(lwd=spline.size, col=spline.color, alpha=spline.alpha))
    l.spline.top <- grid::xsplineGrob(x=c(0, 0, 1, 1), y=c(0, 1, 1, 1), shape=1, gp=gpar(lwd=spline.size, col=spline.color, alpha=spline.alpha))
    r.spline.btm <- grid::xsplineGrob(x=c(1, 1, 0, 0), y=c(1, 0, 0, 0), shape=1, gp=gpar(lwd=spline.size, col=spline.color, alpha=spline.alpha))
    r.spline.top <- grid::xsplineGrob(x=c(1, 1, 0, 0), y=c(0, 1, 1, 1), shape=1, gp=gpar(lwd=spline.size, col=spline.color, alpha=spline.alpha))
    
    spline.y.offset = 0.25 # separation between exon and spline
    label.y.offset = 1.5 # separation between spline and label
    
    # Determine which splines to use for minus strand versus plus strand transcripts
    l.grob.btm <- l.spline.btm 
    if(sashimi.data$is.reverse.strand) l.grob.btm <- r.spline.btm
    
    l.grob.top <- l.spline.top
    if(sashimi.data$is.reverse.strand) l.grob.top <- r.spline.top
    
    r.grob.btm <- r.spline.btm 
    if(sashimi.data$is.reverse.strand) r.grob.btm <- l.spline.btm
    
    r.grob.top <- r.spline.top
    if(sashimi.data$is.reverse.strand) r.grob.top <- l.spline.top
    
    xmid <- (xmin + xmax) / 2
    
    # Left arc
    if(is.even){ # Junctions below zero
      splot <- splot+annotation_custom(grob = l.grob.btm, 
                                       xmin = xmin, 
                                       xmax = xmid, 
                                       ymin = -(1 + f.length), 
                                       ymax = -spline.y.offset)
    } else {
      splot <- splot+annotation_custom(grob = l.grob.top, 
                                       xmin = xmin, 
                                       xmax = xmid, 
                                       ymax = 1 + f.length, 
                                       ymin = spline.y.offset)
    }
    
    # Right arc    
    if(is.even){ # Junctions below zero
      splot <- splot+annotation_custom(grob = r.grob.btm, 
                                       xmin = xmid, 
                                       xmax = xmax, 
                                       ymin = -(1 + f.length), 
                                       ymax = -spline.y.offset)
    } else {
      splot <- splot+annotation_custom(grob = r.grob.top, 
                                       xmin = xmid, 
                                       xmax = xmax, 
                                       ymax = 1 + f.length, 
                                       ymin = spline.y.offset)
    }
    
    splot <- splot + annotate("label", x = xmid, 
                              y = ifelse(is.even, -label.y.offset - f.length, label.y.offset + f.length), 
                              label = as.character(count),
                              size=2, col = spline.color, fill=NA, label.size=NA)
    
    return(splot)
  }
  
  # Check junctions meet plot criteria
  if(nrow(junctions)>0){
    
    # print(junctions)
    # str(junctions)
    junctions <- junctions %>% 
      dplyr::arrange(x, xend) %>%
      dplyr::mutate(length = abs(xend-x))%>%
      dplyr::filter(count>=min.spanning.reads) 
  }
  
  # Add the junctions, adjusting for plus vs minus strand
  if(nrow(junctions)>0){
    
    # Calculate charting coordinates
    junctions$isEven <- sapply(1:nrow(junctions), \(x) x%%2==0)
    junctions$y.offset <- rep(seq(0, 1, 0.25), length.out = nrow(junctions)) # give each junction a separate y offset
    
    for(i in 1:nrow(junctions)){
      splot <- add.junction(splot, junctions[i,]$x, junctions[i,]$xend, junctions[i,]$y, 
                            junctions[i,]$yend, junctions[i,]$isEven, junctions[i,]$count,
                            junctions[i,]$is.canonical, junctions[i,]$y.offset)
    }
  }

  list(plot=splot, junctions=junctions)
}

# Make species plot using combined panels
make.condensed.species.panels <- function(species, timepoint, gene.id){
  data.files <- list.files(path = "data/merged", pattern = paste0(species, ".*\\.", timepoint, "\\.", gene.id, "\\sense.Rds_*"), full.names = TRUE)
  if(length(data.files)==0) return()
  data <- lapply(data.files, read.rds.file)
  
  gene.name <- GENE.LOCATIONS[GENE.LOCATIONS$EnsemblId==gene.id,]$Gene
  
  out.png.file <- paste0("report/species/", species, ".", timepoint, ".", gene.id, ".", gene.name, ".condensed.png")
  cat("Making", out.png.file, "\n")
  
  plots <- lapply(data, \(x)  make.gene.track.sashimi.panel(x,label=paste0(species, " ", gene.name, "\n", x$tissue), show.x.axis = FALSE, is.collapse.introns = TRUE)$plot)
  
  patchwork::wrap_plots(plots, nrow = length(plots))
  ggsave(plot = last_plot(), filename =out.png.file, dpi = 300, units = "mm", width = 170, height = 240)

  data
}

# Make tissue plot using condensed panels
make.condensed.tissue.panels <- function(tissue, timepoint){
  
  cat("Making", tissue, "at", timepoint, "\n")
  data.files <- list.files(path = "data/merged", pattern = paste0(".*\\.", tissue, "\\.", timepoint, "\\..*.Rds_*"), full.names = TRUE)
  if(length(data.files)==0) return()
  data <- lapply(data.files, read.rds.file)
  
  out.png.file <- paste0("report/tissues/", tissue, ".", timepoint, ".condensed.png")
  
  plots <- lapply(data, \(x) make.gene.track.sashimi.panel(x, label= paste0(x$species, "\n", x$gene.name), 
                                                           show.x.axis = FALSE, is.collapse.introns = TRUE)$plot)

  patchwork::wrap_plots(plots, nrow = length(plots))
  ggsave(plot = last_plot(), filename = out.png.file, dpi = 300, units = "mm", width = 170, height = 240)
}

# Make species plot using combined panels
make.condensed.timepoint.panels <- function(species, tissue, gene.id){
  cat("Making", species, tissue, "for", gene.id, "\n")
  data.files <- list.files(path = "data/merged", pattern = paste0(species, ".*\\.", tissue, "\\..*", gene.id, ".*.Rds_*"), full.names = TRUE)
  if(length(data.files)==0) return()
  
  # Ensure files are plotted in time order
  ordered.files <- list()
  for(t in TIME.ORDER){
    for(f in data.files){
      if(str_detect(f, t)){
        ordered.files <- c(ordered.files, f)
      }
    }
  }
  
  data <- lapply(ordered.files, read.rds.file)
  
  gene.name <- GENE.LOCATIONS[GENE.LOCATIONS$EnsemblId==gene.id,]$Gene
  
  out.png.file <-  paste0("report/timepoints/", species, ".", tissue, ".", gene.id, ".", gene.name, "condensed.png")
  
  plots <- lapply(data, \(x) make.gene.track.sashimi.panel(x, label=paste0(x$species, " ", gene.name, "\n",  x$timepoint),
                                                           show.x.axis = FALSE, is.collapse.introns = TRUE)$plot)

  patchwork::wrap_plots(plots, nrow = length(plots))
  ggsave(plot = last_plot(), filename =out.png.file, dpi = 300, units = "mm", width = 170, height = 240)
}

# test.data <- make.combined.panels("macaque", "adult", "ENSMMUG00000046378")
mapply(make.condensed.species.panels, species.groups$CommonName, species.groups$Timepoint,species.groups$EnsemblId)

# make.condensed.tissue.panels("forebrain", "adult")
mapply(make.condensed.tissue.panels, tissue.groups$Organism_part, tissue.groups$Timepoint)

# make.condensed.timepoint.panels("mouse", "testis", "ENSMUSG00000053211")
mapply(make.condensed.timepoint.panels, timepoint.groups$CommonName, timepoint.groups$Organism_part, timepoint.groups$EnsemblId)

