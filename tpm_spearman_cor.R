#!/tpm_spearman_cor.R

#The script takes a series of count tables from a directory.
#Calculate for each gene/tx its TPM value.
#Calculate the Spearman correlation among each count table duo combination.
#Plot a heatmap of the Spearman correlation.

#The inputs are:
#1)pathFiles -> list of directories where the count tables reside.
#The count tables, from the same directory, need to have the same headers and in the same order.
#The count tables, from the same directory, need to have the first column as their Gene/Trasnscript ID
#2)patternFiles -> select the right count tables, in the directory, based on the pattern.
#3)gene_tx_Length -> file path to the file where the gene/transcript and its length are saved (#ID - #Length).
#4)mean_read_seq -> mean of the fragment length sequenced.

#There are two output for each directory-pattern combination in the analysis and they are saved in the same directory as the input files:
#1)pattern_tpm_spearman_df.txt -> text file that shows the TPM Spearman correlation between each count table.
#2)pattern_tpm_spearman_plot.pdf -> heatmap plot of the TPM correlation.

#Load Library
library (tidyverse)   #Collection of R packages designed for data science: ggplot2, dplyr, tidyr, readr, purrr, tibble, stringr, forcats.
library (forestmangr) #Processing forest inventory data with methods such as simple random sampling, stratified random sampling and systematic sampling.

rm(list = ls()) #Remove all saved objects

########### USER's INPUT ###########
pathFiles <- c("/Users/salem_sueto/Desktop/Polyester/Hsapiens_Spneumoniae_SingleStrand_SE/6_count_table/spneumoniae")  #Directories where the count tables reside.
patternFiles <- c("gene_count", "isoform_count")           #Select the right count tables based on a specific pattern
ID_Length <- "/Users/salem_sueto/Desktop/Polyester/Hsapiens_Spneumoniae_SingleStrand_PE/reference/spneumoniae/txome_annotation/tx_gene_length.txt" #Filepath with the gene/tx and its length (#1 ID - #2 Length)
mean_read_seq <- 75                               #Mean of the Fragment Length sequenced
####################################

ID_Length_df <- read.table(ID_Length, header = TRUE) #Upload the ID-Length file
tpm_spearman_plots <- list() #tpm_cor

######Loop through each directory
for (iDir in pathFiles) 
{
  #Loop through each pattern
  for (iPat in patternFiles) 
  {
    #Get the file names
    files <- list.files(path=iDir, pattern=iPat)
    list_of_files <- list()
    uniqueID <- c()
    
    #Save all count tables inside the list
    for (i in files) 
    {
      nowCT <- read.table(paste(iDir, "/", i, sep = ""), header=TRUE)
      nowCT[,1] <- str_replace(nowCT[,1],"\\..*","") #Remove the eventual presence of the Gene/Transcript version after the dot
      list_of_files[[i]] <- nowCT
      uniqueID <- append(uniqueID, nowCT[,1]) #Get the Gene/Transcript IDs from each count table
    }
    uniqueID <- unique(uniqueID) #Get the unique ID from all count tables
    
    #Add rows from the missing IDs
    for (k in 1:length (list_of_files))
    {
      #Get the missing IDs from the specific table based on all unique IDs
      indexID <- list_of_files[[k]] [,1]
      missingIDs <- setdiff(uniqueID, indexID)
      
      #Add the missing IDs to the table with 0 reads in each samples
      if (length(missingIDs) > 0) 
      {
        addRowsDF <- as.data.frame(matrix(data = 0, nrow = length(missingIDs), ncol = ncol(list_of_files[[k]])))
        addRowsDF$V1 <- missingIDs
        colnames(addRowsDF) <- colnames(list_of_files[[k]])
        
        list_of_files[[k]] <- bind_rows(list_of_files[[k]], addRowsDF) %>%  #Row-bind the old table with the new dataframe for the missingIDs
          arrange_at(1) %>%                                                 #Sort the first column (ID)
          column_to_rownames(var = colnames(list_of_files[[k]]) [1])        #Set the first column (ID) as rownames
      } else {list_of_files[[k]] <- list_of_files[[k]] %>% arrange_at(1)  %>% column_to_rownames(var = colnames(list_of_files[[k]]) [1])}
    }
    
    ########### Calculate TPM values ###########
    
    #Filter for IDs present in uniqueIDs
    ID_Length_df_filter <- ID_Length_df %>% filter(ID %in% uniqueID) %>% distinct() %>% arrange_at(1) #Filter for values in uniqueID & Order based on the ID
    
    #TPM calculation function
    counts_to_tpm <- function(counts, featureLength, meanFragmentLength) {
      
      # Ensure valid arguments.
      stopifnot(length(featureLength) == nrow(counts))
      stopifnot(length(meanFragmentLength) == ncol(counts))
      
      # Compute effective lengths of features in each library.
      effLen <- do.call(cbind, lapply(1:ncol(counts), function(i) {
        featureLength - meanFragmentLength[i] + 1
      }))
      
      # Exclude genes with length less than the mean fragment length (NOT IN THIS CASE)
      idx <- apply(effLen, 1, function(x) min(x) > 1) 
      
      idx_true <- ID_Length_df_filter$ID [which(idx == "TRUE")]   #IDs where it is true
      idx_false <- ID_Length_df_filter$ID [which(idx == "FALSE")] #IDs where it is false
      
      counts <- counts[idx,]
      effLen <- effLen[idx,]
      featureLength <- featureLength[idx]
      
      # Process one column at a time.
      tpm <- do.call(cbind, lapply(1:ncol(counts), function(i) {
        rate = log(counts[,i]) - log(effLen[,i])
        denom = log(sum(exp(rate)))
        exp(rate - denom + log(1e6))
      }))
      
      #Prepare TPM output
      colnames(tpm) <- colnames(counts) # Copy the column names from the original matrix.
      rownames(tpm) <- idx_true
      
      idx_false_df <- matrix(data = 0, nrow = length(idx_false), ncol = ncol(counts)) 
      colnames(idx_false_df) <- colnames(counts)
      rownames(idx_false_df) <- idx_false
      
      final_tpm <- as.data.frame(rbind(tpm, idx_false_df)) %>% rownames_to_column(var = "ID") %>% arrange_at(1) %>% column_to_rownames("ID")
      return(final_tpm)
    }
    
    #Calculate TPM 
    list_count_table_tpm <- list()
    
    for (k in 1:length (list_of_files))
    {
      #Calculate TPM from the count table
      meanFragmentLength <- rep(mean_read_seq, ncol(list_of_files[[k]]))
      list_count_table_tpm[[files[k]]] <- counts_to_tpm (list_of_files[[k]], ID_Length_df_filter$Length, meanFragmentLength)
    }
    
    
    ##### Get all the different combinations of count table comparisons
    args <- list(group1 = names(list_of_files), group2 = names(list_of_files))
    args <- args %>% cross_df() %>% filter(group1 != group2) #Get the combinations but remove the ones where the same table name is present
    removeRow <- c()
    
    for (i in 1:nrow(args)) 
    {
      indexRow <- which (args$group2 == args$group1[i] & args$group1 == args$group2[i])
      if (i < indexRow) removeRow <- append(removeRow, indexRow)
    }
    
    args <- args[-removeRow, ] 
    
    
    #Get table with all the Spearman correlation values for each sample and table combinations
    tpm_model_df <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 5))
    colnames(tpm_model_df) <- c("Sample", "Group1", "Group2", "Type", "Spearman")
    
    for (i in colnames (list_count_table_tpm[[1]]))
    {
      for (j in 1:nrow(args))
      {
        #First Table
        tableIndex <- which (names(list_of_files) == args$group1[j])
        singleCol <- list_count_table_tpm[[tableIndex]] %>% dplyr::select (i) #Get the single column (sample) for that count table
        colnames(singleCol) <- paste(str_split(args$group1[j], pattern = "_") [[1]][1], i, sep = "_") #Change the column name (e.i. Tool_sample)
        
        #Second Table
        tableIndex <- which (names(list_of_files) == args$group2[j])
        singleCol2 <- list_count_table_tpm[[tableIndex]] %>% dplyr::select (i) #Get the single column (sample) for that count table
        colnames(singleCol2) <- paste(str_split(args$group2[j], pattern = "_") [[1]][1], i, sep = "_") #Change the column name (e.i. Tool_sample)
        
        #Full TPM table & Calculate Spearman correlation
        tpm_table <- cbind(singleCol, singleCol2)
        spearmanCor <- round(cor(log2(tpm_table[,1]), log2(tpm_table[,2]), method="spearman"), digits = 5)
        
        #Save inside the "tpm_model_df" table
        g1 <- str_split(args$group1[j], pattern = "_")[[1]][1]
        g2 <- str_split(args$group2[j], pattern = "_")[[1]][1]
        type <- paste(g1, g2, sep = "_")
        tpm_model_df [nrow(tpm_model_df) + 1,] <- c(i, g1, g2, type, spearmanCor)
      }
    }
    tpm_model_df$Spearman <- as.numeric(as.character(tpm_model_df$Spearman))
    
    #Save results inside the directory: TPM_Spearman_Result
    dir.create(file.path(iDir, "TPM_Spearman_Result"), showWarnings = FALSE)
    write.table (tpm_model_df, file = sprintf("%s/TPM_Spearman_Result/%s_tpm_spearman_df.txt", iDir, iPat), sep="\t", row.names = FALSE, quote=F)
    
    
    #Heatmap Plot
    tpm_spearman_plots[[iPat]] <- tpm_model_df %>% 
      ggplot(aes(x=Sample, y=Type, fill=Spearman)) + 
      geom_tile(colour = "white") + 
      geom_text(aes(label = round(Spearman, 4))) +
      #facet_grid(rows = vars(Type)) +
      scale_fill_gradient(low="blue", high="red") +
      coord_equal() +
      labs(x="", y="", title = sprintf("Count Table TPM Analysis: %s", iPat), fill="Spearman") +
      theme(axis.text.x = element_text(size=9, face="bold", colour = "black"), axis.text.y = element_text(size=11, face="bold", colour = "black"))
  }
  
  #Print plots
  pdf(sprintf("%s/TPM_Spearman_Result/tpm_spearman_plot.pdf", iDir), onefile=TRUE)
  tpm_spearman_plots
  dev.off()
}
