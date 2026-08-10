library(tidyverse)

# Read files
bed <- read.table("cactus4_12_merged.bed", header = FALSE)
colnames(bed) <- c("chr", "start", "end")

genome <- read.table("genome.txt", header = FALSE)
colnames(genome) <- c("chr", "length")

# Order chromosomes
genome$chr <- factor(genome$chr, levels = rev(genome$chr))
bed$chr <- factor(bed$chr, levels = levels(genome$chr))

# Plot
ggplot() +
  geom_segment(data = genome,
               aes(x = 0, xend = length,
                   y = chr, yend = chr),
               linewidth = 2,
               color = "grey70") +
  geom_segment(data = bed,
               aes(x = start, xend = end,
                   y = chr, yend = chr),
               linewidth = 3,
               color = "red") +
  theme_bw() +
  labs(x = "Position (bp)",
       y = "Chromosome",
       title = "Distribution of cactus regions across the genome")

ggsave("Distribution_cactus_regions.svg", plot = get_last_plot(), dpi = 300, height = 10, width = 15)

