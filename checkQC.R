# Check sample QC

library(tidyverse)
library(ggbeeswarm)
library(patchwork)

filtered.samples <- read.filtered.samples()

extract.pct <- function(x){
  x <- stringr::str_extract(x, "\\(.*\\)")
  x <- stringr::str_replace(x, "\\(", "")
  x <- stringr::str_replace(x, "%\\)", "")
  as.numeric(x)
}

extract.val <- function(x){
  as.numeric(stringr::str_replace(x, " \\(.*\\)" , ""))
}

# Read the mapping efficiencies

# Extract the mapping summary from stdout files
# cat bash.o* | grep -w -e 'mapping' -e 'Aligned' -e 'rate' | tr -d '\t' > report/mapping.txt

map.data <- read.table("report/mapping.txt", sep="$")

# Spread to 5 columns
map.values <- as.data.frame(do.call(cbind, lapply(1:5, \(x)  map.data[seq(x, nrow(map.data), 5),1])))
colnames(map.values) <- c("Run", "Zero", "One", "Multiple", "Overall")

map.values <- map.values %>% 
  dplyr::mutate(Run = stringr::str_remove_all(Run, ": mapping"),
                Zero = stringr::str_remove_all(Zero, "Aligned 0 time: "),
                One = stringr::str_remove_all(One, "Aligned 1 time: "),
                Multiple = stringr::str_remove_all(Multiple, "Aligned >1 times: "),
                Overall = stringr::str_remove_all(Overall, "Overall alignment rate: "),
                ZeroPct = extract.pct(Zero),
                OnePct = extract.pct(One),
                MultiplePct = extract.pct(Multiple),
                OverallPct = as.numeric(stringr::str_replace(Overall, "%", "")),
                ZeroVal  = extract.val(Zero),
                OneVal  = extract.val(One),
                MultipleVal  = extract.val(Multiple),
                OverallVal  = ZeroVal+OneVal+MultipleVal
                
                ) %>%
  # Merge in the sample info,
  merge(., filtered.samples, by="Run")

# Plot the mapping efficiencies

p2 <- ggplot(map.values, aes(x = CommonName, y = ZeroPct, col=CommonName))+
  geom_beeswarm()+
  labs(y="Unmapped reads (%)", col="Species")+
  facet_wrap(~Organism_part)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.title.x = element_blank())

p3 <- ggplot(map.values, aes(x = CommonName, y = OnePct, col=CommonName))+
  geom_beeswarm()+
  labs(y="Single mapping reads (%)", col="Species")+
  facet_wrap(~Organism_part)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.title.x = element_blank())

p4 <- ggplot(map.values, aes(x = CommonName, y = MultiplePct, col=CommonName))+
  geom_beeswarm()+
  labs(y="Multiple mapping reads (%)", col="Species")+
  facet_wrap(~Organism_part)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.title.x = element_blank())

p5 <- ggplot(map.values, aes(x = CommonName, y = OverallVal, col=CommonName))+
  geom_beeswarm()+
  labs(y="Number of reads", col="Species")+
  facet_wrap(~Organism_part)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.title.x = element_blank())


p1 <- ggplot(map.values, aes(x = CommonName, y = OverallPct, col=CommonName))+
  geom_beeswarm()+
  labs(y="Overall mapping (%)", col="Species")+
  theme_bw()+
  facet_wrap(~Organism_part)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.title.x = element_blank())


p2 + p3 +p4 + patchwork::plot_layout(guides = "collect") & theme(legend.position = "none")

ggsave(plot=last_plot(), filename = "report/mapping.qc.png", dpi=300, units="mm",
       width=200, height = 170)


p1 +p5 + patchwork::plot_layout(guides = "collect")
ggsave(plot=last_plot(), filename = "report/mapping.qc.a.png", dpi=300, units="mm",
       width=200, height = 170)


