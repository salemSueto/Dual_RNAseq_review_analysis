#!/polyester_sim.R

#The script takes as input the complete transcriptome of two, or more, species (Host and Pathogen).
#Then, it selects all the transcripts isoforms from the list of genes given by the user.
#Finally, it creates simulated RNA-seq samples from each transcriptome fasta.

#The script has the following inputs:
#1)txome_list -> list of path to the the transcriptome fasta files.
#2)gene_selection -> path to the file with the list of genes used to filter the transcriptome.
#The file contains one column for each transcriptome with header having the same name of the transcriptome.
#Each column has the list of genes which is used for the filtering step.
#3)polyester_settings -> path the file which has all the settings for Polyester to simulate data.
#The file has the following columns:
#DE_TXs (Number of differentially expressed transcripts)
#Coverage (average coverage for all transcript)
#Min_FC (minimum value of fold change for differentially expressed transcripts)
#Max_FC (maximum value of fold change for differentially expressed transcripts)
#Read_Length (length of the simulated reads)
#Paired_End (TRUE for paired-end; FALSE for single-end)
#Strand_Specific (TRUE for strand specific, FALSE for unstranded)
#Replicate_Group1 (number of replicates for group1)
#Replicate_Group2 (number of replicates for group2).
#The rownames of the file has the same name of the transcriptome fasta files.
#4)outputDir -> name of the directory where to save all the output files of the script.

#The script has different outputs and saved inside the output-directory (outputDir):

#1)One output directory for each transcriptome file where:
#1A)List of simulated compressed fasta file equal to the sum of replicates for both groups.
#1B)id_length.txt -> shows the gene/transcript ID and its length.
#The transcript's length is obtained from the transcriptome fasta file. And the gene's length is the sum of all the gene's isoforms length.
#1C)sim_rep_info.txt -> shows the sample's group membership.
#1D)sim_tx_info.txt -> shows the fold change and its DE status for each transcript simulated.
#1E)sim_counts_matrix.rda -> shows the correct count table for all the simulated transcripts in rda format.
#1F)sim_isoform_count.txt -> shows the correct transcript count for the simulated samples.
#1G)sim_gene_count.txt -> shows the correct gene count for the simulated samples.

#2)The Polyester merged simulation (polyester_merged_sim) directory has the information from each simulated data is merged together.
#2A)id_length.txt -> merges all the gene/transcript's length (1B) files together.
#2B)sim_count_gene.txt -> merges the gene count for each simulated transcriptome.
#2C)sim_count_isoform.txt -> merges the transcript count for each simulated transcriptome.
#2D)sim_count_percentage.txt -> shows the percentage of each transcriptome based on the total number of reads.
#2E)sim_count_summary.txt -> shows the total number of reads for each simulated transcriptome.
#2F)sim_de_gene.txt -> shows the differential expressed genes for each simulated transcriptome.
#2G)sim_de_tx.txt -> shows the differential expressed transcripts for each simulated transcriptome.
#2H)sim_rep_info.txt -> shows the sample's group membership for each simulated transcriptome.
#2I)sim_tx_info.txt -> shows the fold change and its DE status for each simulated transcriptome.

#Load Library
library(polyester)  #Simulate RNA-seq reads from differential expression experiments with replicates.
library(Biostrings) #Fast manipulation of large biological sequences or sets of sequences.
library(tidyverse)  ##Collection of R packages designed for data science: ggplot2, dplyr, tidyr, readr, purrr, tibble, stringr, forcats.

rm(list = ls()) #Remove all saved objects

##### ##### ##### USER's INPUT ##### ##### ##### 
txome_list <- c("Hsapiens_Spneumoniae_SingleStrand_SE/reference/hsapiens/Homo_sapiens.GRCh38.cdna.all.fa", "Hsapiens_Spneumoniae_SingleStrand_SE/reference/spneumoniae/Streptococcus_pneumoniae_d39_gca_000014365.ASM1436v2.cdna.all.fa") #Transcriptome List
gene_selection <- "Aprianto_gene_selection.txt" #Gene selection file
polyester_settings <- "Aprianto_polyester_settings.txt" #Polyester parameters setting for each species in the analysis
outputDir <- "Polyester_dual_RNAseq_Hs_Sp" #Path to the directory where the outputs will be saved (Do not add the last "/" in the path)

######General Objects
gene_selection_df <- read.table(gene_selection, header=TRUE, fill = TRUE) #Dataframe with the fasta filename where the genes/transcripts to select for a listed for each fasta file
polyester_settings_df <- read.table(polyester_settings, header=TRUE, fill = TRUE) #Polyester settings for each fasta file
txome_selectedTx_list <- list() #The list saves the filtered version of the transcriptome files
out_dirs <- c() #Saves the path destination of the simulated data
polyester_outputs_list <- list() #Saves all the polyester output dataframes
polyester_de_gene_list <- list() #Saves the list of differential expressed genes from the polyester simulated data
polyester_de_tx_list <- list() #Saves the list of differential expressed transcripts from the polyester simulated data


###### Create a filtered fasta version of the species' transcriptome based on the gene list saved in the gene_selection_df dataframe
dir.create(file.path(getwd(), outputDir), showWarnings = FALSE) #Create output directory where to save all the output files

for (iTx in txome_list) 
{
  iTx_fa <- readDNAStringSet(iTx) #Read the fasta file

  if(length(iTx_fa) > 0)
  {
    iCol <- grep(basename(iTx), colnames(gene_selection_df))
    
    if(length(iCol)==1)
    {
      gene_sel <- gene_selection_df[,iCol]
      gene_sel <- gene_sel[gene_sel != ""]
      gene_sel <- sapply(strsplit(gene_sel, "\\."), `[`, 1)
      
      #Save the selected transcripts in a new Fasta file
      all_search <- paste(gene_sel, collapse="|")
      
      txID <- sapply(strsplit(names(iTx_fa), " "), `[`, 1)
      txID <- sapply(strsplit(txID, "\\."), `[`, 1)
      
      geneID <- sapply(str_split(names(iTx_fa), "gene:"), `[`, 2)
      geneID <- sapply(str_split(geneID, " "), `[`, 1)
      geneID <- sapply(str_split(geneID, "\\."), `[`, 1)
      
      iTx_fa_df <- as.data.frame(cbind(txID, geneID))
      selIndex <- unique (c(which(iTx_fa_df$txID %in% gene_sel), which(iTx_fa_df$geneID %in% gene_sel)))
      txome_selectedTx_list [[basename(iTx)]] <- iTx_fa[selIndex]
      writeXStringSet (iTx_fa[selIndex], sprintf("%s/%s/%s.selectedTXs.fa", getwd(), outputDir, basename(iTx)), append=F, compress=F, compression_level=NA, format="fasta")
    }else{print(paste("The ", iTx, " fasta is not saved in the gene selection file.", sep = ""))}
  }else{print(paste("The ", iTx, " fasta file has no sequence saved.", sep = ""))}
}


###### Simulate RNA-seq reads based on the input from the polyester_species_settings_df dataframe
for (iTx in txome_list) 
{
  iTx_basename <- basename(iTx)
  iRow <- grep(iTx_basename, rownames (polyester_settings_df))
  cov_reads <- round(polyester_settings_df$Coverage[iRow] * width(txome_selectedTx_list[[iTx_basename]]) / 100) #Coverage Reads per Transcript
  
  #Fold Change
  fc <- matrix (data = 1, nrow = length(txome_selectedTx_list[[iTx_basename]]), ncol = 2) 
  
  for (i in 1:polyester_settings_df$DE_TXs[iRow]) 
  {
    colSelect <- floor(runif(1, min=1, max=3)) #Randomly select the columnn between 1 & 2 (max=3 because it will never actually equal 3)
    fc [i, colSelect] <- runif(1, min=polyester_settings_df$Min_FC[iRow], max=polyester_settings_df$Max_FC[iRow]) #Choose random number between the min and max value of the FC
  }
  fc <- fc[sample(1:nrow(fc)), ] #Randomise the row order
  
  #Simulations Reads
  out_name <- sprintf("%s.selectedTXs_cov_%d_rep_%d_%d", iTx_basename, polyester_settings_df$Coverage[iRow], polyester_settings_df$Replicate_Group1[iRow], polyester_settings_df$Replicate_Group2[iRow])
  out_dirs <- append(out_dirs, sprintf("%s/%s", outputDir, out_name))

  simulate_experiment (sprintf("%s/%s.selectedTXs.fa", outputDir, iTx_basename), 
                       reads_per_transcript = cov_reads, 
                       num_reps = c(polyester_settings_df$Replicate_Group1[iRow], polyester_settings_df$Replicate_Group2[iRow]),
                       fold_changes = fc, 
                       outdir = sprintf("%s/%s", outputDir, out_name), 
                       readlen = polyester_settings_df$Read_Length[iRow],
                       paired = polyester_settings_df$Paired_End[iRow],
                       strand_specific = polyester_settings_df$Strand_Specific[iRow], 
                       seed = 200,
                       gzip=TRUE)
  
  #Get the Length of Gene & Transcript 
  txome_tx_len <- width(txome_selectedTx_list [[grep(iTx_basename, names(txome_selectedTx_list))]])
  txome_name <- names(txome_selectedTx_list [[grep(iTx_basename, names(txome_selectedTx_list))]])
  
  txome_name_tx <- sapply(strsplit(txome_name, " "), `[`, 1)
  txome_name_tx <- sapply(strsplit(txome_name_tx, "\\."), `[`, 1)
  
  txome_name_gene <- sapply(strsplit(txome_name, " "), `[`, 4) %>% str_replace("gene:", "")
  txome_name_gene <- sapply(strsplit(txome_name_gene, "\\."), `[`, 1)
  
  tx_gene_len <- as.data.frame(cbind(txome_name_tx, txome_name_gene, txome_tx_len))
  
  gene_len <- tx_gene_len %>% 
    mutate_at(3, as.integer) %>%
    group_by(txome_name_gene) %>% 
    summarise(txome_tx_len = sum(txome_tx_len)) %>%
    rename(ID=txome_name_gene, Length=txome_tx_len)
  
  tx_len <- tx_gene_len[,c(1,3)] %>% rename(ID=txome_name_tx, Length=txome_tx_len)
  
  id_length_df <- rbind(gene_len, tx_len)
  write.table(id_length_df, file= sprintf("%s/id_length.txt", sprintf("%s/%s", outputDir, out_name)), row.names=FALSE, sep="\t", quote = FALSE)
  
  #Transcript Information
  nowInfo <- read.csv2 (paste(sprintf("%s/%s", outputDir, out_name), "/", "sim_tx_info.txt", sep = ""), header=TRUE, sep = "\t")
  nowInfo$Species <- str_replace(iTx_basename, pattern = ".cdna.all.fa", "")
  nowInfo$DEstatus <- "FALSE"
  nowInfo$DEstatus [c(which(nowInfo$DEstatus.1 == "TRUE"), which(nowInfo$DEstatus.2 == "TRUE"))] <- "TRUE"
  write.table(nowInfo, file= sprintf("%s/sim_tx_info.txt", sprintf("%s/%s", outputDir, out_name)), row.names=FALSE, sep="\t", quote = FALSE)
  
  #Count Transcript Table
  load (paste (sprintf("%s/%s", outputDir, out_name), "/sim_counts_matrix.rda", sep = ""))
  sim_count_matrix <- as.data.frame(counts_matrix) %>% rownames_to_column(var = "ID")
  
  sim_count_matrix$TranscriptID <- sapply(strsplit(sim_count_matrix$ID, " "), `[`, 1)
  sim_count_matrix$TranscriptID <- sapply(strsplit(sim_count_matrix$TranscriptID, "\\."), `[`, 1)
  
  sim_count_matrix$GeneID <- sapply(str_split(sim_count_matrix$ID, "gene:"), `[`, 2)
  sim_count_matrix$GeneID <- sapply(str_split(sim_count_matrix$GeneID, " "), `[`, 1)
  sim_count_matrix$GeneID <- sapply(str_split(sim_count_matrix$GeneID, "\\."), `[`, 1)
  
  sim_count_matrix <- subset(sim_count_matrix, select = -c(ID))
  sim_count_matrix <- sim_count_matrix %>% dplyr::select(TranscriptID, GeneID, everything())
  write.table(sim_count_matrix[,c(-2)], file= sprintf("%s/sim_isoform_count.txt", sprintf("%s/%s", outputDir, out_name)), row.names=FALSE, sep="\t", quote = FALSE)
  
  #Count Gene Table
  sim_gene_count <- sim_count_matrix %>% dplyr::select(c(-1)) %>% group_by(GeneID) %>% summarise(across(everything(), sum))
  write.table(sim_gene_count, file= sprintf("%s/sim_gene_count.txt", sprintf("%s/%s", outputDir, out_name)), row.names=FALSE, sep="\t", quote = FALSE)
}


###### Merge all the Polyester output from each simulated RNA-seq data

#Check that all organisms have the same replicates in each group
if (setequal(polyester_settings_df$Replicate_Group1, polyester_settings_df$Replicate_Group2) == TRUE) 
{
  #Create the merged directory
  merged_info_dir <- sprintf("%s/polyester_merged_sim", outputDir)
  ifelse(!dir.exists(merged_info_dir), dir.create(merged_info_dir), FALSE)
  
  #Save the transcript info from each organism
  all_tx_info <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 9))
  de_gene_list <- list()
  de_tx_list <- list()
  all_tx_count <- as.data.frame(matrix(data = NA, nrow = 0, ncol = polyester_settings_df$Replicate_Group1[1] + polyester_settings_df$Replicate_Group1[2] + 1))
  all_tx_summary <- as.data.frame(matrix(data = NA, nrow = 0, ncol = polyester_settings_df$Replicate_Group1[1] + polyester_settings_df$Replicate_Group1[2] + 1))
  all_gene_count <- as.data.frame(matrix(data = NA, nrow = 0, ncol = polyester_settings_df$Replicate_Group1[1] + polyester_settings_df$Replicate_Group1[2] + 1))
  all_id_length <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 2))
  
  for (iDir in out_dirs)
  {
    #Save Group Information
    sim_rep_info <- read.csv2 (paste(iDir, "/sim_rep_info.txt", sep = ""), header=TRUE, sep = "\t")
    if (file.exists(sprintf("%s/sim_rep_info.txt", merged_info_dir))==FALSE)
    {write.table(sim_rep_info, file= sprintf("%s/sim_rep_info.txt", merged_info_dir), row.names=FALSE, sep="\t", quote = FALSE)}
    
    #Save Transcript Info
    nowInfo <- read.csv2 (paste(iDir, "/", "sim_tx_info.txt", sep = ""), header=TRUE, sep = "\t")
    nowInfo$TxID <- sapply(strsplit(nowInfo$transcriptid, " "), `[`, 1)
    nowInfo$TxID <- sapply(strsplit(nowInfo$TxID, "\\."), `[`, 1)
    
    nowInfo$GeneID <- sapply(strsplit(nowInfo$transcriptid, " "), `[`, 4) %>% str_replace("gene:", "")
    nowInfo$GeneID <- sapply(strsplit(nowInfo$GeneID, "\\."), `[`, 1)
    all_tx_info <- rbind(all_tx_info, nowInfo)
    
    #Save DE Genes & DE Transcripts
    nowInfo_de <- nowInfo %>% filter(DEstatus.1==TRUE | DEstatus.2==TRUE)
    
    de_gene_list[[tail (str_split(iDir, "/")[[1]], 1)]] <-  sapply(str_split(unique(nowInfo_de$GeneID), pattern = "\\."), `[`, 1)
    de_tx_list[[tail (str_split(iDir, "/")[[1]], 1)]] <- sapply(str_split(unique(nowInfo_de$TxID), pattern = "\\."), `[`, 1)
    
    #Save Count tables
    sim_isoform_count <- read.csv2 (paste(iDir, "/sim_isoform_count.txt", sep = ""), header=TRUE, sep = "\t") #Transcript count data
    all_tx_count <- rbind(all_tx_count, sim_isoform_count)
    all_tx_summary <- rbind(all_tx_summary, c(tail (str_split(iDir, "/")[[1]], 1), colSums(Filter(is.numeric, sim_isoform_count))))
    
    sim_gene_count <- read.csv2 (paste(iDir, "/sim_gene_count.txt", sep = ""), header=TRUE, sep = "\t") #Gene count data
    all_gene_count <- rbind(all_gene_count, sim_gene_count)

    #Save the id_length file
    id_length_now <- read.table (paste(iDir, "/", "id_length.txt", sep = ""), header=TRUE, sep = "\t")
    all_id_length <- rbind(all_id_length, id_length_now)
  }
  
  #Summary Table
  colnames(all_tx_summary) <- c("ID", colnames(sim_isoform_count)[-1])
  iNum <- 2:ncol(all_tx_summary) #All columns without the Column "ID"
  all_tx_summary [ , iNum] <- apply(all_tx_summary [ , iNum], 2, function(x) as.numeric(as.character(x)))
  tot_sample <- colSums (Filter(is.numeric, sim_isoform_count))
  
  species_sample_per <- as.data.frame(matrix(data=NA, nrow = nrow(all_tx_summary), ncol = ncol(all_tx_summary)))
  colnames(species_sample_per) <- colnames(all_tx_summary)
  species_sample_per[,1] <- all_tx_summary[,1]
  
  for (iCol in 2:ncol(all_tx_summary)) 
  {
    sampleSum <- sum(all_tx_summary[,iCol])
    for (iRow in 1:nrow(all_tx_summary)) {species_sample_per[iRow, iCol] <- all_tx_summary[iRow, iCol] / sampleSum * 100}
  }
  
  #Save all the Polyester results table
  write.table(all_tx_info, file= sprintf("%s/sim_tx_info.txt", merged_info_dir), row.names=FALSE, sep="\t", quote = FALSE)
  write.table(data.frame(lapply(de_gene_list, "length<-", max(lengths(de_gene_list)))), file= sprintf("%s/sim_de_gene.txt", merged_info_dir), row.names=FALSE, sep="\t", quote = FALSE)
  write.table(data.frame(lapply(de_tx_list, "length<-", max(lengths(de_tx_list)))), file= sprintf("%s/sim_de_tx.txt", merged_info_dir), row.names=FALSE, sep="\t", quote = FALSE)
  write.table(all_tx_count, file= sprintf("%s/sim_count_isoform.txt", merged_info_dir), row.names=FALSE, sep="\t", quote = FALSE)
  write.table(all_tx_summary, file= sprintf("%s/sim_count_summary.txt", merged_info_dir), row.names=FALSE, sep="\t", quote = FALSE) ###########
  write.table(species_sample_per, file= sprintf("%s/sim_count_percentage.txt", merged_info_dir), row.names=FALSE, sep="\t", quote = FALSE)
  write.table(all_gene_count, file= sprintf("%s/sim_count_gene.txt", merged_info_dir), row.names=FALSE, sep="\t", quote = FALSE)
  write.table(all_id_length, file= sprintf("%s/id_length.txt", merged_info_dir), row.names=FALSE, sep="\t", quote = FALSE)
}
