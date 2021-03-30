#!/scotty_plot.R

#The script takes Scotty results and create heatmap plots for the Power Plot and Power Bias analysis.

#The inputs are:
#1)scotty_dir -> directory's path where the Scotty results are saved.
#Due to the high number of outputs from one Scotty analysis, the files are arranged in the following way:
#Create one main Directory called for example "Scotty". And create sub-directories based on the following scheme:
#Main Directory -> Tool -> Species -> ReadType -> Comparison -> Scotty's result
#The "Comparison" directory name has the following scheme: Group1_Group2.
#Example: Scotty -> STAR-HTSeq -> Hsapiens -> Gene -> wildType_control -> Scotty's result files

#2)scotty_info -> path to the count table information file which has the following columns:
#Comparison -> it shows the group comparison (Group1_Group2) as saved in the Scotty directories.
#minReplicate -> it sets the minimum number of the replicate in the group comparison.
#Host -> column name is the species name and it shows the minimum sequencing depth.
#Pathogen -> column name is the species name and it shows the minimum sequencing depth.

#3)outDir -> directory path where to save the output files. 

#The output is one pdf file called "scotty_plot.pdf" 
#and a text file called "scotty_data.txt" 
#saved in the specified ouput-directory 


#Load Library
library(tidyverse) #Collection of R packages designed for data science: ggplot2, dplyr, tidyr, readr, purrr, tibble, stringr, forcats.
library(ggforce) #Aid in visual data investigations.

rm(list=ls(all=TRUE))

####### USER's Input ######
scotty_dir <- "/Users/salem_sueto/Desktop/Dual_RNAseq_review_methods/Zimmermann/Scotty_Group_Comparisons"
scotty_info <- "/Users/salem_sueto/Desktop/Dual_RNAseq_review_methods/Zimmermann/metadata/scotty_minReplicate_SeqDepth.txt" #Comparison #MinReplicate #SeqDepth
outDir <- "/Users/salem_sueto/Desktop/Dual_RNAseq_review_methods/Zimmermann/Scotty_Group_Comparisons"
###########################

###Global Objects
scotty_info_df <- read.table(scotty_info, header=TRUE, sep = "\t")
plot_list <- list()  #Save all the plots

###Get all files in the Scotty directory and in all its sub-directories
powerPlot_list <- list.files (path=scotty_dir, pattern = "powerOptimization.csv", recursive = TRUE)
powerBias_list <- list.files (path=scotty_dir, pattern = "poissonNoise.csv", recursive = TRUE)
nReps_list <- list.files (path=scotty_dir, pattern = "nReps.csv", recursive = TRUE)
seqDepth_list <- list.files (path=scotty_dir, pattern = "seqDepth.csv", recursive = TRUE)

###Get all information based on the working condition
nInfos <- str_count(powerPlot_list[1], "/")   #Number of Information based on the count of directories present by "/"
powerPoisson_df <- as.data.frame(matrix(data = NA, nrow = 0, ncol = nInfos + 5))

#Check that all the 4 files are present in each directory
if (var(c (length(powerPlot_list), length(powerBias_list), length(nReps_list), length(seqDepth_list))) == 0)
{
  for (iPP in 1:length(powerPlot_list)) 
  {
    #Get Power Plot & Power Bias Tables
    now_PowerPlot <- read.table(sprintf("%s/%s", scotty_dir, powerPlot_list[iPP]), header=FALSE, sep = ",") #Power Plot table
    now_PowerBias <- read.table(sprintf("%s/%s", scotty_dir, powerBias_list[iPP]), header=FALSE, sep = ",") #Power Bias table
    now_nReps <- read.table(sprintf("%s/%s", scotty_dir, nReps_list[iPP]), header=FALSE, sep = ",")         #Replicate table
    now_seqDepth <- read.table(sprintf("%s/%s", scotty_dir, seqDepth_list[iPP]), header=FALSE, sep = ",")   #Sequencing Depth table
    
    #Get the working conditions
    nowTableName <- str_replace_all(powerPlot_list[iPP], pattern = "/", "564752xuqyto4e488688098")
    nowTableName <- head (str_split (nowTableName, pattern = "564752xuqyto4e488688098")[[1]], -1) #Tool #Species #Gene/Isoform #Comparison
    
    ###ScottyState -> "WorkingCondition"###
    nowComparison_minRep <- scotty_info_df %>% dplyr::filter(Comparison==nowTableName[4]) %>% dplyr::pull(minReplicate) #Minimum Replicate
    nowComparison_minDepth <- scotty_info_df %>% dplyr::filter(Comparison==nowTableName[4]) %>% dplyr::pull(nowTableName[2]) #Minimum Sequencing Depth
    
    getMinColumn <- round(nowComparison_minDepth, digits = -1)/10
    if(getMinColumn==0){getMinColumn<- 1}
    
    powerOptimization <- now_PowerPlot[nowComparison_minRep-1, getMinColumn] #PowerOptimization
    poissonNoise <- now_PowerBias[nowComparison_minRep-1, getMinColumn] #PoissonNoise
    
    powerPoisson_df <- rbind(powerPoisson_df, c(nowTableName, nowComparison_minRep, nowComparison_minDepth, powerOptimization, poissonNoise, "WorkingCondition"))
    
    ###ScottyState -> "MaxReplicate"###
    powerOptimization <- now_PowerPlot[max(now_nReps)-1, getMinColumn] #PowerOptimization
    poissonNoise <- now_PowerBias[max(now_nReps)-1, getMinColumn] #PoissonNoise
    powerPoisson_df <- rbind(powerPoisson_df, c(nowTableName, max(now_nReps), nowComparison_minDepth, powerOptimization, poissonNoise, "MaxReplicate"))

    #ScottyState -> "MaxSeqDepth"
    powerOptimization <- now_PowerPlot[nowComparison_minRep-1, max(now_seqDepth)/10] #PowerOptimization
    poissonNoise <- now_PowerBias[nowComparison_minRep-1, max(now_seqDepth)/10] #PoissonNoise
    powerPoisson_df <- rbind(powerPoisson_df, c(nowTableName, nowComparison_minRep, max(now_seqDepth), powerOptimization, poissonNoise, "MaxSeqDepth"))
    
    #ScottyState -> "MaxReplicateDepth"
    powerOptimization <- now_PowerPlot[max(now_nReps)-1, max(now_seqDepth)/10] #PowerOptimization
    poissonNoise <- now_PowerBias[max(now_nReps)-1, max(now_seqDepth)/10] #PoissonNoise
    powerPoisson_df <- rbind(powerPoisson_df, c(nowTableName, max(now_nReps), max(now_seqDepth), powerOptimization, poissonNoise, "MaxRep_MaxSeqDepth"))
  }
} else(print("The four files (powerPlot.csv, powerBias.csv, nReps.csv, seqDepth.csv) are not present together in all directories."))

#Set the order of the Pipelines
colnames(powerPoisson_df) <- c("Tool", "Species", "Type", "Comparison", "Replicate", "SeqDepth", "PowerOptimization", "PoissonNoise", "ScottyState")
powerPoisson_df$Tool <- factor(powerPoisson_df$Tool, levels = c("STAR-HTSeq", "STAR-RSEM", "Kallisto", "Salmon", "Paper"))
powerPoisson_df$Comparison <- factor(powerPoisson_df$Comparison, levels = unique(scotty_info_df$Comparison))

#Convert to numeric: Replicate - FalsePositive - TruePositiveRate - TotalDEs
iNum <- c(5, 6, 7, 8)
powerPoisson_df [ , iNum] <- apply(powerPoisson_df [ , iNum], 2, function(x) as.numeric(as.character(x)))


###Create Heatmap plot -> Tools used against the group comparisons
for (iSpecies in unique(powerPoisson_df$Species)) 
{
  for (iState in unique(powerPoisson_df$ScottyState)) 
  {
    plot_list [[sprintf("%s - %s - %s", iSpecies, iState, "PowerOptimization")]] <- powerPoisson_df %>% 
      filter (Species==iSpecies & ScottyState==iState) %>%
      ggplot(aes(x=Tool, y=Comparison, fill=PowerOptimization)) + 
      geom_tile() + 
      geom_text(aes(label = format(round(PowerOptimization, 2), nsmall = 2))) +
      scale_fill_gradient(low="blue", high="red") +
      coord_equal() +
      labs(x="", y="", title = sprintf("%s - %s - PowerOptimization", iSpecies, iState)) +
      theme(axis.text.x = element_text(size=9, face="bold", colour = "black"), axis.text.y = element_text(size=11, face="bold", colour = "black"))
    
    plot_list [[sprintf("%s - %s - %s", iSpecies, iState, "PoissonNoise")]] <- powerPoisson_df %>% 
      filter (Species==iSpecies & ScottyState==iState) %>%
      ggplot(aes(x=Tool, y=Comparison, fill=PoissonNoise)) + 
      geom_tile() + 
      geom_text(aes(label = format(round(PoissonNoise, 2), nsmall = 2))) +
      scale_fill_gradient(low="blue", high="red") +
      coord_equal() +
      labs(x="", y="", title = sprintf("%s - %s - PoissonNoise", iSpecies, iState)) +
      theme(axis.text.x = element_text(size=9, face="bold", colour = "black"), axis.text.y = element_text(size=11, face="bold", colour = "black"))
  }
}


###Plot -> Power-Optimization Vs. Poisson-Noise
plot_list [["Power-Optimization - Poisson-Noise"]] <- powerPoisson_df %>%
  ggplot(aes(x=PowerOptimization, y=PoissonNoise, shape=Tool, color=ScottyState)) +
  geom_point(size=3) +
  facet_grid(cols = vars(Comparison), rows = vars(Species)) +
  scale_shape_manual (values=c(0:19)) +
  scale_color_manual(values=c("red", "blue", "purple", "green"))


###Save plots in one file
write.table(powerPoisson_df, file= sprintf("%s/scotty_data.txt", outDir), row.names=FALSE, sep="\t", quote = FALSE)
save(plot_list, powerPoisson_df, file = sprintf("%s/scotty_data_plots.RData", outDir))

pdf(sprintf("%s/scotty_plot.pdf", outDir), onefile=TRUE)
plot_list
dev.off()
