#!/salmon_merge_table.R

#The script merges the Salmon reference-free alignments outputs into one file for both gene and transcript level.

#The script has 2 inputs:
#1)fileDirectory -> list of directories where the salmon's outputs are saved.
#2)gene_tx -> path to a file which has atleast two columns called "GeneID" and "TranscriptID".

#The script outputs two txt files for each directory
#1)"salmon_gene_count.txt" -> merges all the gene count
#2)"salmon_isoform_count.txt" -> merges all the transcripts count

#Load Library
library(tidyverse) #Collection of R packages designed for data science: ggplot2, dplyr, tidyr, readr, purrr, tibble, stringr, forcats.
library(tximport) #Import and summarize transcript-level estimates for transcript- and gene-level analysis

rm(list=ls(all=TRUE))

##### ##### ##### User's Input ##### ##### ##### 
fileDirectory <- c("Hsapiens_Spneumoniae_SingleStrand_SE/4_salmon")  #Directory where the Salmon results are saved (DO NOT add the last "/")
gene_tx <- "Hsapiens_Spneumoniae_SingleStrand_SE/metadata/Hsapiens_Spneumoniae_tx_gene_rel.txt"  #TranscriptID - GeneID - GTF - CDNA file's path
################################################## 

##### Get Polyester Output Data
gene_tx_df <- read.csv2(gene_tx, sep="\t", header = TRUE)
tx2tx <- gene_tx_df[,c(1,1)]               #Get count for transcripts
tx2gene <- gene_tx_df[,c(1,2)]             #Get count for Genes
colnames(tx2tx) <- c("TXNAME", "GENEID")
colnames(tx2gene) <- c("TXNAME", "GENEID")

#Loop through each directory
for (iDir in fileDirectory) 
{
  #Get all the sample names 
  files <- list.dirs(iDir, recursive=FALSE)
  files <- paste0(files, "/quant.sf")
  dirNames <- str_split(list.dirs(iDir, recursive=FALSE), pattern = "/")
  dirNames <- unlist (lapply(dirNames, function(x) x[length(x)]))
  names(files) <- dirNames
  
  #Remove the transcript version from the transcript ID column from the abundance.tsv files
  for (iTxt in files) 
  {
    iTxT_df <- read.csv2(iTxt, sep="\t", header = TRUE)
    iTxT_df$Name <- sapply (str_split(iTxT_df$Name, "\\."), `[`, 1)
    write.table(iTxT_df, file = unname(iTxt), row.names=F, sep="\t", quote=F)
  }
  
  #Get Gene Count Table
  txi_salmon <- tximport(files, type = "salmon", tx2gene = tx2gene)
  count <- as.data.frame(txi_salmon$counts) %>% rownames_to_column("ID")
  write.table(count, file = sprintf("%s/salmon_gene_count.txt", iDir), row.names=F, sep="\t", quote=F)
  
  ##### Get Transcript Count Table
  txi_salmon <- tximport(files, type = "salmon", tx2gene = tx2tx)
  count <- as.data.frame(txi_salmon$counts) %>% rownames_to_column("ID")
  count$ID <- str_replace (count$ID, "\\..*", "")
  write.table(count, file = sprintf("%s/salmon_isoform_count.txt", iDir), row.names=F, sep="\t", quote=F)
}
