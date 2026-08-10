# Load necessary libraries
library(ggplot2)
library(reshape2)
library(scales)
library(gridExtra)
library(RColorBrewer)
library(grid)

# Define a better color palette
color_palette <- brewer.pal(5, "Set2")

# List to store plots
plots <- list()

# Loop through K values
for (i in 1:5) {
  filename <- paste0("Lcas_NGSadmix_k", i, ".1.qopt")
  pop <- read.table("Lcas_bam.filelist")
  df <- read.table(filename)
  df$Individual <- pop$V1
  mdfr <- melt(df, id.vars = "Individual")
  
  # Create the plot
  p <- ggplot(mdfr, aes(Individual, value, fill = variable)) +
    geom_bar(position = "fill", stat = "identity") +
    scale_y_continuous(labels = percent) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      legend.position = "none"
    ) +
    scale_fill_manual(values = color_palette) +
    ggtitle(paste0("Admixture proportions for K=", i))
  
  # Store plot in the list
  plots[[i]] <- p
}

# Combine all plots in a grid layout
combined_plot <- arrangeGrob(grobs = plots, ncol = 2)  # Adjust ncol if needed

# Save the combined plot
ggsave("Lcas_Admix_paper.svg", plot = combined_plot, width = 15, height = 10, dpi = 300)
