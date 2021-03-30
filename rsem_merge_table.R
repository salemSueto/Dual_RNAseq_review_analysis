#!/rsem_merge_table.R

#The script takes RSEM outputs from the transcriptome alignment of STAR and merges them together to form one file at the gene and transcript level.

#The user's input is pathFiles which is the list of directories where the RSEM output files reside.

#The outputs are two txt files for each directory:
#1)"rsem_gene_count.txt" -> merges all the gene count
#1)"rsem_isoform_count.txt" -> merges all the transcripts count

#Load Library
library (tidyverse) #Collection of R packages designed for data science: ggplot2, dplyr, tidyr, readr, purrr, tibble, stringr, forcats.

#Remove all saved objects
rm(list = ls())

################
#USER's INPUT
pathFiles <- c("Hsapiens_Spneumoniae_SingleStrand_SE/2_star_hsapiens/rsem") #Get list of all files in directory
################

#Get the samples
for (iDir in pathFiles) #Loop through the directories
{
  for (iPat in c("gene", "isoform")) #Loop through the patterns
  {
    # Save all files
    files <- list.files(path=iDir, pattern=iPat)
    list_of_files <- list()
    
    if (length(files)>0)
    {
      # Save all files
      for (i in files) 
      {list_of_files[[i]] <- read.table(paste (iDir, "/", i, sep = ""), header=TRUE)}
      
      # Select the "expected_count" column for each files
      merge_rsem_df <- as.data.frame (matrix(data = NA, nrow = nrow(list_of_files[[1]]), ncol = length(list_of_files)+4))
      merge_rsem_df[,c(1,2,3,4)] <- list_of_files[[1]][,1:4]
      for (i in 1:length(list_of_files)) {merge_rsem_df [,i+4] <- list_of_files[[i]]$expected_count}
      merge_rsem_df <- merge_rsem_df %>% dplyr::select(-c(2,3,4))
      colnames(merge_rsem_df) <- c("ID", sapply(str_split(names(list_of_files), "\\."), `[`, 1))
      write.table (merge_rsem_df, file = sprintf("%s/rsem_%s_count.txt", iDir, iPat), sep="\t", row.names=F, quote=F)
    }
  }
}
