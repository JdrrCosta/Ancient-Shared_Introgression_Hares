# Load saguaro output
#change this to how many cacti to plot
for (i in 0:12) {
  cactus <- paste0("cactus", i)
  matrix_data <- read.table(paste0(cactus,".txt"), header = TRUE, row.names = 1)
  
  # Transform into distance matrix
  dist_matrix <- as.dist(as.matrix(matrix_data))
  
  # Check matrix
  print(dist_matrix)
  
  # Install and load necessary libraries
  library(ggplot2)
  library(ggtree)
  library(ape)
  library(dplyr)
  
  
  # Create the tree using neighbor-joining method
  tree <- nj(dist_matrix)
  # Plot the tree with ggtree
  p_unrooted <- ggtree(tree, layout = "unrooted") +
    geom_tiplab2(size = 3, hjust = 0, aes(color = case_when(
      grepl("Lcor", label) ~ "Lcor",   # Lcor labels will be red
      grepl("Leur", label) ~ "Leur",
      grepl("Ltim", label) ~ "Ltim",
      grepl("Lcas", label) ~ "Lcas",
      grepl("Lame", label) ~ "Lame",
      grepl("Lgra", label) ~ "Lgra"))) +  # Tip labels, adjust for left alignment
    scale_color_manual(values = c("Lcor" = "#e6ab02", #chose colors
                                  "Leur" = "#7570b3", 
                                  "Ltim" = "#66a61e", 
                                  "Lcas" = "#d95f02", 
                                  "Lame" = "#e7298a", 
                                  "Lgra" = "#1b9e77")) +
    theme_tree2() +          # Enhanced theme
    theme(
      axis.title = element_blank(),   # Remove axis titles
      axis.text.x = element_blank(),  # Remove x-axis text
      axis.ticks.x = element_blank(), # Remove x-axis ticks
      axis.line = element_blank(),    # Remove axis lines
      plot.margin = margin(0.05, 0.05, 0.05, 0.1, "cm"),  # Expand margins further
      plot.title.position = "plot",    # Position title at the top of the plot area
      plot.title = element_text(hjust = 0)  # Left-align the title
    ) +
    labs(
      title = "Saguaro Iter 12", 
      subtitle = paste0(cactus, " unrooted")
    ) + 
    theme(
      plot.title = element_text(hjust = 0),  # Align title to the left
      plot.margin = margin(0.1, 0.1, 0.1, 0.15, "cm"), # Extra margin to ensure space for labels
      plot.background = element_rect(fill = "white")  # Ensure plot background doesn't interfere
    ) + 
    coord_cartesian(clip = "off")  # Prevent clipping of tip labels
  
  # Remove the scale bar
  p_unrooted <- p_unrooted + theme(legend.position = "none")
  
  # Print the plot
  print(p_unrooted)
  
  # Save the plot with large enough dimensions
  ggsave(paste0(cactus, "_unrooted_tree_color.png"), plot = p_unrooted, width = 15, height = 10, dpi = 400)
}
