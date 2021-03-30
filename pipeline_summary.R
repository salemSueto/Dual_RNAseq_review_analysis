#!/pipeline_summary.R

#The script is used to plot the results of the mapping and the quantification for each step of the dual RNA-seq pipeline. 
#It is possibly to use the number of reads or the percentage of the reads present in that particular step.

#The inputs are the following:
#1)file_path -> file's path where all the information are saved.
#2)sample_group_color -> list of colors used to select the different sample groups.

#The output is a PDF file where the plots are saved.


library(tidyverse)
library(gridExtra)

rm(list = ls()) #Remove all saved objects

##### USER's INPUT ######
input_path <- "/Users/salem_sueto/Desktop/Dual_RNAseq_review_methods/Zimmermann/metadata/Zimmermann_summary.txt"
sample_group_color <- c("coral", "brown1", "aquamarine", "dodgerblue", "gold", "blue", "grey", "red", "green", "magenta")
#########################

input_df <- read.table(input_path, header=TRUE, fill=TRUE, sep="\t")
input_no_group_df <- input_df %>% dplyr::select(-(Group))

plot_list <- list()

for (sp in unique(input_no_group_df$Species)) 
{
  #Melt the table based on the species
  sp_melt <- input_no_group_df %>% filter(Species==sp) %>% dplyr::select(-Species)
  sp_melt <- reshape2::melt(sp_melt,id.vars=c("Sample"))
  
  #Plot 
  plot_list[[sp]] <- sp_melt %>%
    ggplot (aes(x=variable, y=value, shape=Sample, colour=Sample, group=Sample)) + 
    geom_line () +
    geom_point (size = 2) +
    scale_shape_manual(values = c(1:length(unique(sp_melt$Sample)))) +
    xlab("") + ylab("") + ggtitle(sprintf("%s - Percentage Result", sp)) +
    scale_y_continuous(breaks = seq(0, 100, by = 10))
}

#Sample by Sample
melt_list <- list()

for (sp in unique(input_no_group_df$Species))
{
  melt_list[[sp]] <- input_no_group_df %>% 
    filter(Species==sp) %>% 
    dplyr::select(-c(Species)) %>% 
    reshape2::melt(id.vars=c("Sample")) %>%
    add_column(Species = sp)
}

input_melt_df <- plyr::ldply(melt_list, data.frame) %>% dplyr::select(-c(.id))

#Facet_wrap for each group of sample
groups <- unique(input_df$Group)

for (iG in 1:length(groups))
{
  sample_list <- input_df %>% filter(Group==groups[iG]) %>% pull(Sample)
  
  plot_list[[groups[iG]]] <- input_melt_df %>%
    filter(Sample %in% sample_list) %>%
    ggplot (aes(x=variable, y=value, group=Species, color=Species, shape=Species)) + 
    geom_line (alpha=0.3) +
    geom_point (size = 2) +
    facet_wrap(~ Sample, ncol = 1) +
    xlab("") + ylab("") +
    scale_y_continuous(breaks = seq(0, 100, by = 20)) +
    theme(strip.background = element_rect(fill=sample_group_color[iG])) +
    labs(title=groups[iG])
  
}

#All samples together
plot_list[["All_samples"]] <- input_melt_df %>%
  ggplot (aes(x=variable, y=value, group=Species, color=Species, shape=Species)) + 
  geom_line (alpha=0.3) +
  geom_point (size = 2) +
  facet_wrap(~ Sample, ncol = 4) +
  xlab("") + ylab("") #+
  scale_y_continuous(breaks = seq(0, 100, by = 25))

##### Save Heatmap plots #####
save(plot_list, file = sprintf("%s/summary_plots.RData", dirname(input_path)))

pdf(sprintf("%s/summary_plots.pdf", dirname(input_path)), onefile=TRUE)
for (i in seq(length(plot_list))) {grid.arrange (plot_list[[i]])}
dev.off()
