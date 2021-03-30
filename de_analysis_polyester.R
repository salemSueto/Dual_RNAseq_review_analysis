#!/de_analysis_polyester.R

#The script identified differential expressed elements from the merged count tables obtain from the Polyester simulation script
#with DESeq2, edgeR, and NOISeq tools at different number of replicates, row sum cut-off, and adjusted p-values.

#The inputs are:
#1) dirPath -> Path to the directory where the files are saved (DON'T add the last "/" to the path).
#2) replicate -> vector list of number of replicates to apply the analysis (>1).
#3) cutOff -> numeric vector list used to filter for genes/isoforms which has a row-sum higher than the cut-off value.
#4) padj_FDR_value ->  numeric vector list of adjusted p-values used to filter for differential expressed genes/isoforms.
#5) log2FC_value -> single value for the Log2 Fold Change used to select differential expressed genes/isoforms.

#The outputs are saved in the same "polyester_merged_dir" directory
#1) PCA plot file called "ggplots.pdf" which shows the results of the differential expressed analysis when compared to the correct list one.
#2) Text file called "TP_FP_FN.txt" which shows information concerning true positive, false positive, etc..
#3) Text file called "pca_res.txt" which shows the Principal Component Analysis (PCA) results.

#Load Library
library(tidyverse)  #Collection of R packages designed for data science: ggplot2, dplyr, tidyr, readr, purrr, tibble, stringr, forcats.
library(DESeq2)     #Differential gene expression analysis based on the negative binomial distribution
library(edgeR)      #Empirical Analysis of Digital Gene Expression Data in R      
library(NOISeq)     #Exploratory analysis and differential expression for RNA-seq data
library(factoextra) #Easy Multivariate Data Analyses and Elegant Visualization
library(gridExtra)  #Provides functions to work with "grid" graphics, notably to arrange multiple grid-based plots on a page, and draw tables.

rm(list = ls()) #Remove all saved objects

####### USER's Input ######
polyester_merged_dir <- "Hsapiens_Spneumoniae_SingleStrand_SE/polyester_merged_sim" #Path to the merged simulation directory obtained from Polyester
replicate <- c(2,10) #vector list of number of replicates to apply the analysis (>1).
cutOff <- c(0, 5000) #Filter for genes/isoforms which has a row-sum higher than the cut-off value.
padj_FDR_value <- c(0.05, 0.005, 5e-04, 5e-05, 5e-06) #Adjusted p-values used to filter for differential expressed genes/isoforms.
log2FC_value <- 1 #Single value for the Log2 Fold Change used to select differential expressed genes/isoforms.
###########################

###Global Object
de_list <- list()
de_typeSpecies_list <- list()
ggplot_list <- list()

truePositive_de <- list()
falsePositive_de <- list()
falseNegative_de <- list()
truth_de_polyester <- list()

###Get all files in the polyester_merged_dir directory
files <- list.files(path=polyester_merged_dir, pattern="")
list_of_files <- list()

for (i in files) 
{list_of_files[[i]] <- read.table(sprintf("%s/%s", polyester_merged_dir, i), header=T, sep = "\t", quote="", fill=FALSE)}


############################# DEGs Identification #############################
count_table_list <- list_of_files[c(2,3)]

for (iCT in 1:length(count_table_list)) 
{
  #Tidy the Count Table
  nowTable <- count_table_list[[iCT]]
  rowID <- nowTable[,1]
  nowTable <- nowTable %>% dplyr::select (-c(colnames(nowTable)[1])) %>% mutate_if (is.character,as.numeric)
  rownames(nowTable) <- rowID
  
  for (r in replicate)
  {
    for (c in cutOff) 
    {
      print(sprintf("Directory: %s * Table: %s * Replicate: %s * CutOff: %s", polyester_merged_dir, files[iCT], r, c))
      
      #Select the right samples for Group1 based on the number of replicates
      group1Samples <- list_of_files[["sim_rep_info.txt"]] %>% filter(group==1)
      group1Samples <- sample (group1Samples$rep_id, r)  #Randomly select the samples in the group
      
      #Select the right samples for Group2 based on the number of replicates
      group2Samples <- list_of_files[["sim_rep_info.txt"]] %>% filter(group==2)
      group2Samples <- sample (group2Samples$rep_id, r)  #Randomly select the samples in the group
      
      #Select the samples selecteed and tidy the table
      compaDF <- nowTable %>% 
        dplyr::select(all_of(c(group1Samples, group2Samples))) %>% 
        rownames_to_column() %>%
        mutate(Sum=rowSums(dplyr::select(., -starts_with("rowname")))) %>% 
        filter(Sum >= c) %>% 
        dplyr::select(-Sum) %>%
        column_to_rownames()
      
      # Group Factors
      sampleGroups <- append (rep(1,r), rep(2,r))
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
        #DESeq2
        deseq2_deg <- res %>% filter(padj < pv, abs(log2FoldChange) >= log2FC_value)
        de_list[[sprintf("%s*DESeq2*%s*%s*%s", names(count_table_list)[iCT], r, c, pv)]] <- rownames(deseq2_deg)
        
        #edgeR
        edgeR_deg <- tibble::rownames_to_column(edgeR_res$table, "Gene") %>% filter (FDR < pv, abs(logFC) >= log2FC_value)
        de_list[[sprintf("%s*edgeR*%s*%s*%s",  names(count_table_list)[iCT], r, c, pv)]] <- edgeR_deg$Gene
        
        #NOISeq
        noiseq_deg <- degenes(mynoiseqbio, q = 1-pv, M = NULL)
        noiseq_deg <- subset(noiseq_deg, abs(log2FC) >= 1)
        de_list[[sprintf("%s*NOISeq*%s*%s*%s",  names(count_table_list)[iCT], r, c, pv)]] <- rownames(noiseq_deg)
      }
    }
  }
}


############################# True Positive - False Positive - False Negative Analysis #############################
nRows_tp_tn_df <-  length(de_list)*length(unique(list_of_files[["sim_tx_info.txt"]]$Species))
tp_fp_fn_df <- as.data.frame(matrix(data=NA, ncol=12, nrow=nRows_tp_tn_df))
colnames(tp_fp_fn_df) <- c("Table", "Tool", "Replicate", "CutOff", "AdjPvalue", "Type", "Species", "TruthDE", "TruePositive", "FalsePositive", "FalseNegative", "TruePositiveRate")
i_TP_TN_FN <- 1

#Differential Expressed Type (Gene/Transcript)
for (iType in c("gene", "isoform"))  
{
  #Get the DE which has the Type (gene/tx)
  de_type_list <- de_list [grep(pattern = iType, names(de_list))]
  
  #Loop through the species
  for (iSpecies in unique(list_of_files[["sim_tx_info.txt"]]$Species)) 
  {
    #Get all the transcripts for the species * Get the truth DE
    all_tx_species_df <- list_of_files[["sim_tx_info.txt"]] %>% filter(Species==iSpecies)
    truth_de <- c()
    
    if(iType=="gene") 
    {
      truth_de <- all_tx_species_df %>% filter(DEstatus==TRUE) %>% dplyr::select(GeneID)
      truth_de <- as.vector(truth_de$GeneID)
      
      #Get the Species DE for each DE identified
      for (iDE in 1:length(de_type_list))
      {
        iDE_analysis <- de_type_list[[iDE]]
        iDE_analysis_type_species <- all_tx_species_df %>% filter(GeneID %in% iDE_analysis)
        iDE_analysis_type_species <- unique(iDE_analysis_type_species$GeneID)
        
        tp <- intersect(truth_de, iDE_analysis_type_species)
        fp <- setdiff(iDE_analysis_type_species, truth_de)
        fn <- setdiff(truth_de, iDE_analysis_type_species)
        tpr <- length(tp)/length(truth_de)
        
        #Table - Tool - Replicate - CutOff - AdjPvalue - Type - Species - TruthDE - TruePositive - FalsePositive - FalseNegative - TruePositiveRate
        iDE_analysis_info <- c(str_split(names(de_type_list)[iDE], pattern="\\*")[[1]], iType, iSpecies, length(truth_de), length(tp), length(fp), length(fn), tpr)
        tp_fp_fn_df[i_TP_TN_FN,] <- iDE_analysis_info
        i_TP_TN_FN <- i_TP_TN_FN+1
        
        de_typeSpecies_list[[paste(iDE_analysis_info[1:7], collapse = "*")]] <- iDE_analysis_type_species
      }
      truth_de_polyester[[sprintf("%s*%s", iType, iSpecies)]] <- truth_de
    }
    
    if(iType=="isoform") 
    {
      truth_de <- all_tx_species_df %>% filter(DEstatus==TRUE) %>% dplyr::select(TxID)
      truth_de <- as.vector(truth_de$TxID)
      
      #Get the Species DE for each DE identified
      for (iDE in 1:length(de_type_list))
      {
        iDE_analysis <- de_type_list[[iDE]]
        iDE_analysis_type_species <- all_tx_species_df %>% filter(TxID %in% iDE_analysis)
        iDE_analysis_type_species <- unique(iDE_analysis_type_species$TxID)
        
        tp <- intersect(truth_de, iDE_analysis_type_species)
        fp <- setdiff(iDE_analysis_type_species, truth_de)
        fn <- setdiff(truth_de, iDE_analysis_type_species)
        tpr <- length(tp)/length(truth_de)
        
        #Table - Tool - Replicate - CutOff - AdjPvalue - Type - Species - TruthDE - TruePositive - FalsePositive - FalseNegative - TruePositiveRate
        iDE_analysis_info <- c(str_split(names(de_type_list)[iDE], pattern="\\*")[[1]], iType, iSpecies, length(truth_de), length(tp), length(fp), length(fn), tpr)
        tp_fp_fn_df[i_TP_TN_FN,] <- iDE_analysis_info
        i_TP_TN_FN <- i_TP_TN_FN+1
        
        de_typeSpecies_list[[paste(iDE_analysis_info[1:7], collapse = "*")]] <- iDE_analysis_type_species
        
        #Save TruePositive - FalsePositive - FalseNegative DE
        truePositive_de[[paste(iDE_analysis_info[1:7], collapse = "*")]] <- tp
        falsePositive_de[[paste(iDE_analysis_info[1:7], collapse = "*")]] <- fp
        falseNegative_de[[paste(iDE_analysis_info[1:7], collapse = "*")]] <- fn
      }
      
      #Save TruePositive - FalsePositive - FalseNegative DE
      truth_de_polyester[[sprintf("%s*%s", iType, iSpecies)]] <- truth_de
    }
  }
}

#Convert to numeric certain columns
i <- c(3, 8, 10, 12) #Replicate - TruthDE - FalsePositive - TruePositiveRate
tp_fp_fn_df [ , i] <- apply(tp_fp_fn_df [ , i], 2, function(x) as.numeric(as.character(x)))
tp_fp_fn_df$CutOff <- factor(tp_fp_fn_df$CutOff, levels = unique (as.numeric(as.character(tp_fp_fn_df$CutOff)))) #Set the right order for the CutOff factors
tp_fp_fn_df$AdjPvalue <- factor(tp_fp_fn_df$AdjPvalue, levels = unique (as.numeric(as.character(tp_fp_fn_df$AdjPvalue)))) #Set the right order for the AdjPvalue factors


##### ##### GGPLOT
#Colours
color_pvalue <- c("gray", "black")
color_table <- c("red", "blue", "green", "purple", "orange", "brown", "aquamarine", "yellow", "cyan")
color_total <- c(color_pvalue, color_table) #In case you want to descriminate the min and the max value of the AdjPvalue

#Pvalues
pv <- as.numeric(as.character(unique(tp_fp_fn_df$AdjPvalue)))
min_pvalue <- min (pv)
max_pvalue <- max (pv)

for (sp in unique(tp_fp_fn_df$Species)) #Loop through the different species
{
  #Plot with the Min and MAx value of the AdjPvalue
  plotName <- paste(sp, " * ", "Adjusted P-value: ", min_pvalue, " - ", max_pvalue, sep = "")
  
  ggplot_list[[plotName]] <- tp_fp_fn_df %>% 
    filter(Species==sp & AdjPvalue %in% c(min_pvalue, max_pvalue)) %>%   #Filter for Species and AdjPvalue
    ggplot(aes(x=FalsePositive, y=TruePositiveRate, shape=CutOff, color=Table, group = interaction(Table, AdjPvalue))) + 
    geom_point() +
    geom_line (aes(color=AdjPvalue)) +
    facet_grid(cols = vars(Replicate), rows = vars(Tool)) +
    scale_colour_manual(values=color_total) +   #Colour for both tables and the Pvalues (min & max)
    scale_shape_manual (values=c(0:19)) +
    #scale_x_continuous(trans = "log10") +
    theme(panel.spacing = unit(1.2, "lines"), plot.title = element_text(hjust=0.5, size = 6)) +
    theme() +
    labs(x="N* False Positive", y="True Positive Rate", title = plotName)
  
  #Loop through the ADjPvalue
  for (pv in unique(tp_fp_fn_df$AdjPvalue)) 
  {
    plotName <- paste(sp, " * ", "Adjusted P-value: ", pv, sep = "")
    ggplot_list[[plotName]] <- tp_fp_fn_df %>% 
      filter(Species==sp & AdjPvalue==pv) %>%   #Filter for Species and AdjPvalue
      ggplot(aes(x=FalsePositive, y=TruePositiveRate, shape=CutOff, color=Table, group = interaction(Table, AdjPvalue))) + 
      geom_point() +
      geom_line () +
      facet_grid(cols = vars(Replicate), rows = vars(Tool)) +
      scale_colour_manual(values=color_table) +   #Colour for the different tables
      scale_shape_manual (values=c(0:19)) +
      #scale_x_continuous(trans = "log10") +
      theme(panel.spacing = unit(1.2, "lines"), plot.title = element_text(hjust=0.5, size = 6)) +
      #theme(legend.position = "none") +
      labs(x="N* False Positive", y="True Positive Rate", title = plotName)
    
    #Loop through CutOff
    for (co in unique(tp_fp_fn_df$CutOff)) 
    {
      plotName <- paste(sp, " * ", "Adjusted P-value: ", pv, " * ", "CutOff: ", co, sep = "")
      
      ggplot_list[[plotName]] <- tp_fp_fn_df %>% 
        filter (Species==sp & AdjPvalue==pv & CutOff==co) %>%  #Filter for species, pvalues, and cutoff values
        ggplot(aes(x=FalsePositive, y=TruePositiveRate, shape=Table, color=Table)) + 
        geom_point(alpha=1) +
        facet_grid(cols = vars(Replicate), rows = vars(Tool)) +
        scale_colour_manual(values=color_table) +   #Colour for the different tables
        scale_shape_manual (values=c(0:19)) +       #Shape types
        scale_x_continuous(trans = "log10") +
        theme(panel.spacing = unit(1.2, "lines"), plot.title = element_text(hjust=0.5, size = 6)) +
        #theme(legend.position = "none") +
        labs(x="N* False Positive", y="True Positive Rate", title = plotName)
    }
  }
}


############################# PCA Analysis  #############################
iType <- "gene"
iSpecies <- "Homo_sapiens.GRCh38"

for (iType in c("gene", "isoform"))  #Differential Expressed Type (Gene/Transcript)
{
  #Loop through the species
  for (iSpecies in unique(list_of_files[["sim_tx_info.txt"]]$Species)) 
  {
    #Get the correct list of DE based on the type and species
    de_select <- de_typeSpecies_list[intersect(grep(pattern=iType, names(de_typeSpecies_list)), grep(pattern=iSpecies, names(de_typeSpecies_list)))]
    de_select[["Truth*Truth*Truth*Truth*Truth*Truth*Truth"]] <- truth_de_polyester[[intersect(grep(pattern=iType, names(truth_de_polyester)), grep(pattern=iSpecies, names(truth_de_polyester)))]]
    
    #Get all the DE from the de_select
    de_type_element <- unname(unlist(de_select))
    de_type_element <- unique (de_type_element)
    
    #Create Presence/Absent dataframe
    de_matrix_01 <- as.data.frame(matrix(data=0, ncol=length(de_select), nrow=length(de_type_element)))
    colnames(de_matrix_01) <- names(de_select)
    rownames(de_matrix_01) <- de_type_element
    
    for (iCol in 1:length(de_select)) 
    {for (iRow in de_select[iCol]) {de_matrix_01[iRow,iCol] <- 1}}
    
    #Create Group dataframe to divide the information
    groupInfo <- data.frame(ID=names(de_select))
    groupInfo <- groupInfo %>% separate(col=ID, into = c("Table", "Tool", "Replicate", "CutOff", "AdjPvalue", "Type", "Species"), sep = "\\*")
    
    #PCA
    de_matrix_01 <- t(de_matrix_01)
    de_matrix_01 <- de_matrix_01[,which(apply(de_matrix_01, 2, var) != 0)] #Remove constant/zero columns to unit variance
    
    #res.pca <- prcomp(de_matrix_01, scale = TRUE)
    res.pca <- prcomp(de_matrix_01)
    res.pca_df <- as.data.frame(res.pca$x)
    
    percentage <- round(res.pca$sdev / sum(res.pca$sdev) * 100, 2)
    percentage <- paste( colnames(res.pca_df), " (", paste( as.character(percentage), "%", ")", sep=""), sep = "" )
    
    if (ncol(res.pca_df) < 2) {res.pca_df$PC2 <- 0} #In case there is only PC1 -> add PC2 equal to 0
    
    res.pca_df <- cbind(res.pca_df, groupInfo[,1:7])
    res.pca_df$AllInfo <- names(de_select)
    
    #Factor's Order
    res.pca_df$CutOff <- factor(res.pca_df$CutOff, levels = unique (as.character(res.pca_df$CutOff))) #Set the right order for the CutOff factors
    res.pca_df$AdjPvalue <- factor(res.pca_df$AdjPvalue, levels = unique (as.character(res.pca_df$AdjPvalue))) #Set the right order for the AdjPvalue factors
    res.pca_df$Replicate <- factor(res.pca_df$Replicate, levels = unique (as.character(res.pca_df$Replicate))) #Set the right order for the AdjPvalue factors
    
    set.seed(33)
    ggplot_list[[sprintf("PCA: %s * %s", iSpecies, iType)]] <- res.pca_df %>%  
      filter(Table != "Truth") %>% 
      ggplot(aes(x=PC1, y=PC2, color=CutOff, shape=AdjPvalue)) +
      geom_point(size=2) +
      annotate("text", x = res.pca_df$PC1[nrow(res.pca_df)], y = res.pca_df$PC2[nrow(res.pca_df)], colour = "black", label = "Truth") +
      facet_grid(cols = vars(Replicate), rows = vars(Tool)) +
      scale_shape_manual (values=c(0:19)) +       #Shape types
      theme(panel.spacing = unit(1.2, "lines"), plot.title = element_text(hjust=0.5, size = 6)) +
      xlab(percentage[1]) + ylab(percentage[2]) +
      labs(title = sprintf("PCA: %s * %s", iSpecies, iType))
  }
}


#Save all plot to one file
pdf(paste(polyester_merged_dir, "/ggplots.pdf", sep = ""), onefile=TRUE)
for (i in seq(length(ggplot_list))) {grid.arrange (ggplot_list[[i]])}
dev.off()

#Save dataframe
write.table(res.pca_df, file= sprintf("%s/pca_res.txt", polyester_merged_dir), row.names=FALSE, sep="\t", quote = FALSE)
write.table(tp_fp_fn_df, file= sprintf("%s/TP_FP_FN.txt", polyester_merged_dir), row.names=FALSE, sep="\t", quote = FALSE)
