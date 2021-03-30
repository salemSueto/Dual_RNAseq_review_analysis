#!/htseq_merge_table.R

#The script takes the gene counts obtained from STAR in-built HTSeq-count, during the mapping step with the flag "--quantMode GeneCounts", 
#and merges them together to form one file for each of the three sequencing type.

#The inputs are:
#1)pathFiles -> list of directories where the files reside
#2)patternFiles -> one pattern used to select the right files in the directory

#The output are:
#1)Count tables for each strand type (unstranded, strand, and reverse strand).
#2)summary\_count.txt -> text file which has the summary info for all three strand types:
#number of unmapped, number of multimapping, number of no features, number of ambiguous, and number of mapped reads.
#3)One pdf file of the summary plots information.

#Load Library
library (tidyverse) #Collection of R packages designed for data science: ggplot2, dplyr, tidyr, readr, purrr, tibble, stringr, forcats.
library(gridExtra) #Provides the ability to arrange multiple grid-based plots on a page, and draw tables.

rm(list = ls()) #Remove all saved objects

####### USER's Input ######
pathFiles <- c("Hsapiens_Spneumoniae_SingleStrand_SE/3_star_spneumoniae/htseq") #Get list of all files in directory
patternFiles <- "sample" #Pattern used to select the right files
###########################

#Loop through each directory
for (iDir in pathFiles) 
{
  #Get all files which has the patternFiles character
  files <- list.files(path=iDir, pattern=patternFiles)
  list_of_files <- list()
  
  if (length(files) > 0)
  {
    # Save all files
    for (i in files) 
    {list_of_files[[i]] <- read.table(paste(iDir, "/", i, sep = ""), header=FALSE)}
    
    #Summary Info Dataframe
    summary_plot_df <- as.data.frame (matrix(nrow = 0, ncol = 4))
    colnames(summary_plot_df) <- c("ID", "variable", "value", "MapType")
    
    #Extract count for each sequencing type
    selectColID <- c("unstranded", "strand_yes", "strand_rev")
    
    for (selectCol in c(2,3,4)) #2: unstranded - #3: 1st read strand aligned with RNA (htseq-count -s yes) - #4: 1st read strand aligned with RNA (htseq-count -s reverse)
    {
      # Select the column for each files
      merge_htseq_df <- as.data.frame (matrix(data = NA, nrow = nrow(list_of_files[[1]])-4, ncol = length(list_of_files)))
      summary_htseq_df <- as.data.frame (matrix(data = NA, nrow = 4, ncol = length(list_of_files)))
      
      colnames(merge_htseq_df) <- names(list_of_files)
      colnames(summary_htseq_df) <- names(list_of_files)
      
      rownames(merge_htseq_df) <- list_of_files[[1]]$V1 [-(1:4)]
      rownames(summary_htseq_df) <- list_of_files[[1]]$V1 [1:4]
      
      for (i in 1:length(list_of_files)) 
      {
        sFile <- list_of_files[[i]] %>% dplyr::select (contains(as.character(selectCol)))
        summary_htseq_df[,i] <- sFile[1:4, 1]
        merge_htseq_df [,i] <- tail(sFile, -4)
      }
      
      merge_htseq_colsums <- colSums(merge_htseq_df)
      merge_htseq_df <- merge_htseq_df %>% rownames_to_column(var = "ID")
      colnames(merge_htseq_df) <- sapply (str_split(colnames(merge_htseq_df), "\\."), `[`, 1)
      
      summary_htseq_df <- rbind (summary_htseq_df, N_mapped_Quantified = merge_htseq_colsums) %>% rownames_to_column(var = "ID")
      summary_htseq_df$MapType <- selectColID[selectCol-1]
      colnames(summary_htseq_df) <- sapply (str_split(colnames(summary_htseq_df), "\\."), `[`, 1)
      
      temp_summary <- reshape2::melt(summary_htseq_df,id.vars=c("ID"))
      temp_summary$MapType <- selectColID[selectCol-1]
      summary_plot_df <- rbind(summary_plot_df, temp_summary)
      
      fileName <- paste("htseq-count", selectColID[selectCol-1], sep = "_")
      write.table(merge_htseq_df, paste(iDir, "/", fileName, "_count.txt", sep = ""), row.names=F, sep="\t", quote=F)
    }
    
    ####PLOT SUMMARY INFO #####
    summary_plot_df <- summary_plot_df %>% filter(value!=MapType)
    write.table(summary_plot_df, paste(iDir, "/summary_count.txt", sep = ""), row.names=F, sep="\t", quote=F)
    
    plots_list <- list()
    for (i in unique(summary_plot_df$ID))
    {
      plots_list[[paste("Variable:", i)]] <- summary_plot_df %>%
        filter(ID==i) %>%
        ggplot(aes(x=variable, y=value, fill=MapType)) + 
        geom_bar(position = "dodge", stat="identity") +
        labs(title=paste("Variable:", i)) + 
        theme(axis.text.x = element_text(angle=65, hjust=1.0))
    }
    
    pdf(paste(iDir, "/htseq_summary_plots.pdf", sep = ""), onefile=TRUE)
    for (i in seq(length(plots_list))) {grid.arrange (plots_list[[i]])}
    dev.off()
  }
}
