#!/de_countTables_comparison.R

#The script takes a group of count tables and identifies its differential expressed (DE) elements with DESeq2, edgeR, and NOISeq at different values of cut-off, replicates, and adjusted p-values.
#Then, it compares the identified DE elements to a list furnished by the user.
#Then, it plots different figures.

#The inputs are:
#1)count_table_dir -> Path to the directory where the count tables are saved. DO NOT add the last "/" on the path.
#2)patternFile -> Pattern used to select the right count tables.
#3)compare_table -> Path to the file the genes/isoforms' and their LFC are saved for the comparison. The column names are ID and LFC.
#4)sample_group -> Path to the file where the sample's group are saved as two columns called "Sample" and "Group".
#5)output_dir -> output directory name where the outputs will be save.
#6)groupComparison -> List of group comparison written as "Group1_Group2".
#7)selectReplicate -> List of values for the number of replicates to use per group during the differential expressed analysis.
#8)cutOff -> List of cut-off values for the filtering row (genes/isoforms) step which has a row-sum higher than the cut-off value.
#9)padj_FDR_value -> List of adjusted p-values used to filter for differential expressed genes/isoforms.
#10)log2FC_value -> Single value for the Log2 Fold Change used to select differential expressed genes/isoforms.

#The outputs are saved in a new created directory called from the value of "output_dir" saved in the "count_table_dir" directory.
#1)comparison_tpr_fp.txt -> table with the comparison results with the columns:
#Table (table's name) #Comparison (comparison name) #Tool (tool used to find DE) #Replicate (number of replicates) #CutOff (cut-off value) #AdjPvalue (adjusted p-value)
#TruthDE (number of DE from the user's input) 
#TruePositive (number of DE elements in common between the user's input and the identified DE)
#FalsePositive (number of DE elements found in the identified DE but absent in the user's input)
#FalseNegative (number of DE elements found in the user's input but absent in the identified DE)
#TruePositiveRate (number of True Positive divided by the user's input DE)
#2)comparison_truthDE_lfc.txt -> table shows the Spearman correlation between the real LFC and the value from from the analysis with the columns:
#Table (table's name) #Comparison (comparison name) #Tool (tool used to find DE) #Replicate (number of replicates) #CutOff (cut-off value) #AdjPvalue (adjusted p-value)
#Correlation (Spearman's correlation between the truth LFC value and the LFC found from our analysis)
#3)plots.pdf -> one PDF plot file for all plots.
#4)plots.rds -> one RData where the plots are saved.

#Load Library
library(tidyverse)  #Collection of R packages designed for data science: ggplot2, dplyr, tidyr, readr, purrr, tibble, stringr, forcats.
library(DESeq2)     #Differential gene expression analysis based on the negative binomial distribution
library(edgeR)      #Empirical Analysis of Digital Gene Expression Data in R      
library(NOISeq)     #Exploratory analysis and differential expression for RNA-seq data
library(gridExtra)  #Provides functions to work with "grid" graphics, notably to arrange multiple grid-based plots on a page, and draw tables.
library(corrplot)

rm(list = ls()) #Remove all saved objects

####### USER's Input ######
#Files selection
count_table_dir <- "Hsapiens_Spneumoniae_SingleStrand_SE/6_count_table/hsapiens"                   #Path to the directory where the count tables are saved.
patternFile <- "isoform"                                  #Pattern used to select the right count tables from the directory
compare_table <- "Hsapiens_Spneumoniae_SingleStrand_SE/metadata/truth_hsapiens_2_1_de_isoform.txt"  #Path to the directory where the compare tables are present
sample_group <- "Hsapiens_Spneumoniae_SingleStrand_SE/metadata/sim_rep_info.txt"     #Path to the file where the sample's group are saved as two columns called "Sample" and "Group".
output_dir <- "Comparison_hsapiens_isoform" #Create Directory in the "count_table_dir" directory where to save all results

#DE Information
groupComparison <- c("2_1")                          #List of group comparison written as "Group1_Group2".
selectReplicate <- c(2)                                #The highest number of replicates to use per group during the differential expressed analysis
cutOff <- c(0)                                     #Filter for genes/isoforms which has a row-sum higher than the cut-off value.
padj_FDR_value <- c(0.05)          #Adjusted p-values used to filter for differential expressed genes/isoforms.
log2FC_value <- 1                                             #Single value for the Log2 Fold Change used to select differential expressed genes/isoforms.
###########################

###Global Objects
list_count_table <- list() #Save all the count tables

de_list <- list() #Save the list of DEGs found for each count table at different settings of cutOff, replicates, etc..
de_allInfo_list <- list() #Save all the information at the end of the DEGs identification step.

ggplot_list <- list() #Save the plots

#Get the compare file with column names as "ID" and "LFC"
compare_table_df <- read.table(compare_table, header=TRUE, sep = "\t") %>% 
                    drop_na() %>%                                       #Get the compare DE list file
                    mutate(LFC = as.numeric(LFC)) %>%                   #Change the LFC column to numeric type
                    dplyr::group_by(ID) %>%                             #Group by the ID
                    dplyr::summarize(LFC = mean(LFC, na.rm=TRUE))       #In case of duplicated IDs get the mean LFC value.


sample_group_df <- read.table(sample_group, header=TRUE, sep = "\t")  #Get the sample and its group membership file
comparison_de_df <- as.data.frame(matrix(data=NA, ncol=12, nrow=0))   #Comparison Table Result

#Get the count tables
files <- list.files(path=count_table_dir, pattern=patternFile)
for (i in files) {list_count_table[[i]] <- read.table(sprintf("%s/%s", count_table_dir, i), header=TRUE, sep = "\t") %>% drop_na()}

### Differential expressed analysis ###
iGC <- 1
iCT <- 1
iRep <- selectReplicate[1]
iCO <- cutOff[1]
pv <- padj_FDR_value[1]

for (iGC in 1:length(groupComparison))   #Loop through the group comparisons
{
  #Get the list of DE to compare to
  de_compare_list <- unique (compare_table_df$ID)
  
  #Dataframe for LFC analysis
  LFC_analysis_df <- compare_table_df %>% 
    dplyr::select(ID, LFC) %>% 
    dplyr::rename(Truth=LFC) %>% 
    dplyr::distinct(ID, .keep_all = TRUE) %>%
    dplyr::arrange(ID)
  
  for (iCT in 1:length(list_count_table)) #Loop through the count tables
  {
    #Tidy the Count Table
    nowTable <- list_count_table[[iCT]]
    rowID <- nowTable[,1]
    nowTable <- nowTable %>% dplyr::select (-c(colnames(nowTable)[1])) %>% mutate_if (is.character,as.numeric)
    rownames(nowTable) <- rowID
    
    for (iRep in selectReplicate)  #Loop through the number of replicates from 2 to the maximum number of replicate
    {
      for (iCO in cutOff) #Loop through the cut-off values
      {
        print(sprintf("Table: %s * Group Comparison: %s * Replicate: %s * CutOff: %s", files[iCT], iGC, iRep, iCO))
        
        #Get the info for the two groups
        groupInfo <- str_split(groupComparison, pattern="_")[[1]]   
        
        if(length(groupInfo)==2)
        {
          #Select the right samples for Group1 based on the number of replicates
          group1Samples <- sample_group_df %>% filter(Group==groupInfo[1])
          group1Samples <- sample (group1Samples$Sample, iRep)  #Randomly select the samples in the group
          
          #Select the right samples for Group2 based on the number of replicates
          group2Samples <- sample_group_df %>% filter(Group==groupInfo[2])
          group2Samples <- sample (group2Samples$Sample, iRep)  #Randomly select the samples in the group
          
          #Select the samples selecteed and tidy the table
          compaDF <- nowTable %>% 
            dplyr::select(all_of(c(group1Samples, group2Samples))) %>% 
            rownames_to_column() %>%
            mutate(Sum=rowSums(dplyr::select(., -starts_with("rowname")))) %>% 
            filter(Sum >= iCO) %>% 
            dplyr::select(-Sum) %>%
            column_to_rownames()
          
          # Group Factors
          sampleGroups <- append (rep(1,iRep), rep(2,iRep))
          myfactors_D_E <- factor(sampleGroups)             #Factors (DESeq2 & edgeR)
          myfactors_N <- data.frame(group = sampleGroups)    #Factors for NOISeq
          
          #------- DESeq2 ------#
          matrixCountTable <- round (as.matrix(compaDF))
          coldata <- data.frame(row.names=colnames(matrixCountTable), myfactors_D_E)
          dds <- DESeqDataSetFromMatrix(countData=matrixCountTable, colData=coldata, design=~myfactors_D_E)
          dds <- DESeq(dds)
          res <- results(dds)
          res <- as.data.frame(res)
          
          #------- edgeR ------#
          y <- DGEList(counts = compaDF, group = myfactors_D_E)
          y <- calcNormFactors(y) #TMM normalisation
          y <- estimateDisp(y)
          et <- exactTest(y)
          edgeR_res <- topTags(et, n = nrow(compaDF), sort.by = "none")
          
          #------- NOISeq (biological replicates ------#
          noiseq_data <- readData(data = compaDF, factors = myfactors_N) #Converting data into a NOISeq object
          gc() #Running Garbage Collection
          mynoiseqbio <- noiseqbio(noiseq_data, k = 0.5, norm = "tmm", factor = "group", r = 20, adj = 1.5, plot = FALSE, a0per = 0.9, random.seed = 12345, filter = 0)
          
          #Adjusted P-values
          for (pv in padj_FDR_value) 
          {
            ###DESeq2
            deseq2_deg <- res %>% filter(padj < pv, abs(log2FoldChange) >= log2FC_value)
            de_list[[sprintf("%s*%s*DESeq2*%s*%s*%s", files[iCT], groupComparison[iGC], iRep, iCO, pv)]] <- unique(rownames(deseq2_deg))
            
            #Add DESeq2 results to the Comparison result dataframe
            tp <- intersect(de_compare_list, rownames(deseq2_deg))
            fp <- setdiff(rownames(deseq2_deg), de_compare_list)
            fn <- setdiff(de_compare_list, rownames(deseq2_deg))
            tpr <- length(tp)/length(de_compare_list)
            comparison_de_df <- rbind(comparison_de_df, c(files[iCT], groupComparison[iGC], "DESeq2", iRep, iCO, pv, length(de_compare_list), length(tp), length(fp), length(fn), tpr, length(unique(rownames(deseq2_deg)))))
            
            #Add the LFC
            deseq2_lfc <- deseq2_deg %>% rownames_to_column(var="ID") %>% filter(ID %in% LFC_analysis_df$ID) %>% dplyr::select(ID, log2FoldChange) %>% dplyr::rename(LFC=log2FoldChange)
            notFoundID <- data.frame(ID=setdiff(LFC_analysis_df$ID, rownames(deseq2_deg)))
            notFoundID$LFC <- 0
            deseq2_lfc_toAdd <- rbind(deseq2_lfc, notFoundID) %>% arrange(ID) %>% dplyr::select(LFC)
            LFC_analysis_df$deseq2_lfc_toAdd <- deseq2_lfc_toAdd
            colnames(LFC_analysis_df)[which(colnames(LFC_analysis_df)=="deseq2_lfc_toAdd")] <- sprintf("%s*%s*DESeq2*%s*%s*%s", files[iCT], groupComparison[iGC], iRep, iCO, pv)
            
            ###edgeR
            edgeR_deg <- tibble::rownames_to_column(edgeR_res$table, "Gene") %>% filter (FDR < pv, abs(logFC) >= log2FC_value)
            de_list[[sprintf("%s*%s*edgeR*%s*%s*%s", files[iCT], groupComparison[iGC], iRep, iCO, pv)]] <- unique(edgeR_deg$Gene)
            
            #Add edgeR results to the Comparison result dataframe
            tp <- intersect(de_compare_list,  edgeR_deg$Gene)
            fp <- setdiff( edgeR_deg$Gene, de_compare_list)
            fn <- setdiff(de_compare_list,  edgeR_deg$Gene)
            tpr <- length(tp)/length(de_compare_list)
            comparison_de_df <- rbind(comparison_de_df, c(files[iCT], groupComparison[iGC], "edgeR", iRep, iCO, pv, length(de_compare_list), length(tp), length(fp), length(fn), tpr, length(unique(edgeR_deg$Gene))))
            
            #Add the LFC
            edgeR_lfc <- edgeR_deg %>% filter(Gene %in% LFC_analysis_df$ID) %>% dplyr::select(Gene, logFC) %>% dplyr::rename(LFC=logFC)
            notFoundID <- data.frame(Gene=setdiff(LFC_analysis_df$ID, edgeR_lfc$Gene))
            notFoundID$LFC <- 0
            edgeR_lfc_toAdd <- rbind(edgeR_lfc, notFoundID) %>% arrange(Gene) %>% dplyr::select(LFC)
            LFC_analysis_df$edgeR_lfc_toAdd <- edgeR_lfc_toAdd
            colnames(LFC_analysis_df)[which(colnames(LFC_analysis_df)=="edgeR_lfc_toAdd")] <- sprintf("%s*%s*edgeR*%s*%s*%s", files[iCT], groupComparison[iGC], iRep, iCO, pv)
            
            #NOISeq
            noiseq_deg <- degenes(mynoiseqbio, q = 1-pv, M = NULL)
            noiseq_deg <- subset(noiseq_deg, abs(log2FC) >= 1)
            de_list[[sprintf("%s*%s*NOISeq*%s*%s*%s", files[iCT], groupComparison[iGC], iRep, iCO, pv)]] <- unique(rownames(noiseq_deg))
            
            #Add NOISeq results to the Comparison result dataframe
            tp <- intersect(de_compare_list,  rownames(noiseq_deg))
            fp <- setdiff(rownames(noiseq_deg), de_compare_list)
            fn <- setdiff(de_compare_list, rownames(noiseq_deg))
            tpr <- length(tp)/length(de_compare_list)
            comparison_de_df <- rbind(comparison_de_df, c(files[iCT], groupComparison[iGC], "NOISeq", iRep, iCO, pv, length(de_compare_list), length(tp), length(fp), length(fn), tpr, length(unique(rownames(noiseq_deg)))))
            
            #Add the LFC
            noiseq_lfc <- noiseq_deg %>% rownames_to_column(var="ID") %>% filter(ID %in% LFC_analysis_df$ID) %>% dplyr::select(ID, log2FC) %>% dplyr::rename(LFC=log2FC)
            notFoundID <- data.frame(ID=setdiff(LFC_analysis_df$ID, rownames(noiseq_deg)))
            notFoundID$LFC <- 0
            noiseq_lfc_toAdd <- rbind(noiseq_lfc, notFoundID) %>% arrange(ID) %>% dplyr::select(LFC)
            LFC_analysis_df$noiseq_lfc_toAdd <- noiseq_lfc_toAdd
            colnames(LFC_analysis_df)[which(colnames(LFC_analysis_df)=="noiseq_lfc_toAdd")] <- sprintf("%s*%s*NOISeq*%s*%s*%s", files[iCT], groupComparison[iGC], iRep, iCO, pv)
          }
        }
      }
    }
  }
}

colnames(comparison_de_df) <- c("Table", "Comparison", "Tool", "Replicate", "CutOff", "AdjPvalue", "TruthDE", "TruePositive", "FalsePositive", "FalseNegative", "TruePositiveRate", "TotalDEs")

#Convert to numeric: Replicate - FalsePositive - TruePositiveRate - TotalDEs
i <- c(4, 9, 11, 12)
comparison_de_df [ , i] <- apply(comparison_de_df [ , i], 2, function(x) as.numeric(as.character(x)))
comparison_de_df$CutOff <- factor(comparison_de_df$CutOff, levels = unique (as.numeric(as.character(comparison_de_df$CutOff)))) #Set the right order for the CutOff factors
comparison_de_df$AdjPvalue <- factor(comparison_de_df$AdjPvalue, levels = unique (as.numeric(as.character(comparison_de_df$AdjPvalue)))) #Set the right order for the AdjPvalue factors


##### ##### GGPLOT Series
#Colours
color_pvalue <- c("gray", "black")
color_table <- c("red", "blue", "green", "purple", "orange", "brown", "aquamarine", "yellow", "cyan")
color_total <- c(color_pvalue, color_table) #In case you want to descriminate the min and the max value of the AdjPvalue

#Pvalues
pv <- as.numeric(as.character(unique(comparison_de_df$AdjPvalue)))
min_pvalue <- min (pv)
max_pvalue <- max (pv)

###PLOT -> Number of Identified Differential Expressed Elements
for (iCom in unique(comparison_de_df$Comparison)) 
{
  ggplot_list[[sprintf("%s - %s - Identified DE elements", patternFile, iCom)]] <- comparison_de_df %>% 
    filter(Comparison==iCom) %>%   #Filter for Group Comparison and AdjPvalue
    ggplot(aes(x=Replicate, y=TotalDEs, shape=CutOff, color=AdjPvalue, group = interaction(Table))) + 
    geom_point() +
    facet_grid(cols = vars(Table), rows = vars(Tool)) +
    scale_colour_manual(values=color_table) +   #Colour for both tables and the Pvalues (min & max)
    scale_shape_manual (values=c(0:19)) +
    labs(title=sprintf("%s - Number of Identified DE elements", iCom)) +
    geom_hline(yintercept=as.numeric(unique(comparison_de_df$TruthDE)[1]), linetype="dashed", color = "black", size=0.5)
}

#Plot -> True Positive Rate vs. Number of False Positive 
for (iCom in unique(comparison_de_df$Comparison))    #Loop through the group comparisons
{
  ggplot_list[[sprintf("%s*%s*Number of Identified DE elements",patternFile, iCom)]] <- comparison_de_df %>% 
    filter(Comparison==iCom & AdjPvalue %in% c(min_pvalue, max_pvalue)) %>%   #Filter for Group Comparison and AdjPvalue
    ggplot(aes(x=FalsePositive, y=TruePositiveRate, shape=CutOff, color=Table, group = interaction(Table, AdjPvalue))) + 
    geom_point() +
    geom_line (aes(color=AdjPvalue)) +
    facet_grid(cols = vars(Replicate), rows = vars(Tool)) +
    scale_colour_manual(values=color_total) +   #Colour for both tables and the Pvalues (min & max)
    scale_shape_manual (values=c(0:19)) +
    theme(panel.spacing = unit(1.2, "lines")) +
    labs(x="N* False Positive", y="True Positive Rate", title=sprintf("%s - Number of Identified DE elements", iCom))
  
  ggplot_list[[sprintf("%s*%s*AxisX_log10",patternFile, iCom)]] <- comparison_de_df %>% 
    filter(Comparison==iCom & AdjPvalue %in% c(min_pvalue, max_pvalue)) %>%   #Filter for Group Comparison and AdjPvalue
    ggplot(aes(x=FalsePositive, y=TruePositiveRate, shape=CutOff, color=Table, group = interaction(Table, AdjPvalue))) + 
    geom_point() +
    geom_line (aes(color=AdjPvalue)) +
    facet_grid(cols = vars(Replicate), rows = vars(Tool)) +
    scale_colour_manual(values=color_total) +   #Colour for both tables and the Pvalues (min & max)
    scale_shape_manual (values=c(0:19)) +
    scale_x_continuous(trans = "log10") +
    theme(panel.spacing = unit(1.2, "lines")) +
    labs(x="N* False Positive", y="True Positive Rate")
}

###PLOT -> LFC Correlation Analysis
#1 -> Get Spearman correlation -> Absolute value of LFC
LFC_analysis_df <- LFC_analysis_df %>% column_to_rownames("ID")

LFC_analysis_df_2 <- as.data.frame(LFC_analysis_df)
LFC_analysis_df_2 [ , i] <- apply(LFC_analysis_df_2 [ , c(1:ncol(LFC_analysis_df_2))], 2, function(x) as.numeric(as.character(x)))
LFC_analysis_df_2[] <- lapply(LFC_analysis_df_2, abs)
cor_LFC_df <- cor(LFC_analysis_df_2, method = "spearman")
cor_truth_row <- subset(cor_LFC_df, rownames(cor_LFC_df) == "Truth")

#Found the Correlation between the truth and the other settings
cor_LFC_truth_df <- data.frame(ID=colnames(cor_LFC_df), Correlation=unname(cor_truth_row[1,]))
cor_LFC_truth_df <- cor_LFC_truth_df[!(cor_LFC_truth_df$ID=="Truth"), ]
cor_LFC_truth_df <- cor_LFC_truth_df %>% separate(ID, c("Table", "Comparison", "Tool", "Replicate", "CutOff", "AdjPvalue"), "\\*")

#Convert to numeric: Replicate - Correlation
i <- c(4,7)
cor_LFC_truth_df [ , i] <- apply(cor_LFC_truth_df [ , i], 2, function(x) as.numeric(as.character(x)))
cor_LFC_truth_df$CutOff <- factor(cor_LFC_truth_df$CutOff, levels = unique (as.numeric(as.character(cor_LFC_truth_df$CutOff)))) #Set the right order for the CutOff factors
cor_LFC_truth_df$AdjPvalue <- factor(cor_LFC_truth_df$AdjPvalue, levels = unique (as.numeric(as.character(cor_LFC_truth_df$AdjPvalue)))) #Set the right order for the AdjPvalue factors

#Plot
for (iCom in unique(cor_LFC_truth_df$Comparison)) 
{
  ggplot_list[[sprintf("%s*%s*Truth DE Correlation ", patternFile, iCom)]] <- cor_LFC_truth_df %>% 
    filter(Comparison==iCom) %>%   #Filter for Group Comparison and AdjPvalue
    ggplot(aes(x=Replicate, y=Correlation, shape=CutOff, color=AdjPvalue, group = interaction(Table))) + 
    geom_point() +
    facet_grid(cols = vars(Table), rows = vars(Tool)) +
    scale_colour_manual(values=color_table) +   #Colour for both tables and the Pvalues (min & max)
    scale_shape_manual (values=c(0:19)) +
    labs(title=sprintf("%s - %s - Correlation between truth DE", patternFile, iCom))
}


##### Save Data & Plots #####
dir.create(file.path(count_table_dir, output_dir), showWarnings = FALSE)

write.table(comparison_de_df, file= sprintf("%s/%s/comparison_tpr_fp.txt", count_table_dir, output_dir), row.names=FALSE, sep="\t", quote = FALSE)
write.table(cor_LFC_truth_df, file= sprintf("%s/%s/comparison_truthDE_lfc.txt", count_table_dir, output_dir), row.names=FALSE, sep="\t", quote = FALSE)
save(ggplot_list,file=sprintf("%s/%s/plots.rda", count_table_dir, output_dir))

pdf(sprintf("%s/%s/plots.pdf", count_table_dir, output_dir), onefile=TRUE)
for (i in seq(length(ggplot_list))) {grid.arrange (ggplot_list[[i]])}
dev.off()
