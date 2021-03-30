#!/de_enrichment_countTables_comparison.R

#The script takes a list of count tables and it identifies the differential expressed terms with DESeq2 for each group comparison.
#Enrichment analysis is carried out for each group comparison and their similarity is showed in a heatmap for all count tables.

#The inputs are:
#1)count_table_dir -> list of path to the directory where the count tables are present (#ID #Sample1  #SampleN).
#2)patternFile -> Pattern used to select the right count tables.
#3)sample_group -> Path to the file where the sample's group are saved as two columns called "Sample" and "Group".
#4)groupComparison -> List of group comparison written as "Group1_Group2".
#5)cutOff -> Filter for genes/isoforms which has a row-sum higher than the cut-off value.
#6)padj_FDR_value -> Adjusted p-values used to filter for differential expressed genes/isoforms.
#7)log2FC_value -> Single value for the Log2 Fold Change used to select differential expressed genes/isoforms.
#8)gProfiler_supported_org -> List of Scientific Name for g:Profiler supported organisms.
#9)gProfiler_non_supported_org -> List of path to GMT custom annotation files for non-supported organisms.
#10)cutoff -> REVIGO's cut-off values (Allowed values: "0.90" "0.70" "0.50" "0.40")
#11)isPValue -> Is the REVIGO's input numbers p-values? (Allowed values: "yes"  "no")
#12)whatIsBetter -> In case of some other quantity than the pvalues where the "allowed" value is better. Allowed values: "higher" "lower" "absolute" "abs_log"
#13)goSizes -> Select a database with GO term sizes (Default is "0" which specifies "whole Uniprot").
#14)measure -> Select a semantic similarity measure to use. Allowed values: "RESNIK" "LIN" "SIMREL" "JIANG"

#The outputs are:
#1)DESeq2_de_res.txt -> DESeq2 result saved for all tables and group comparisons.
#2)A directory is created for each species used during the enrichment step with g:Profiler (e.i. "Homo sapiens").
#A sub-directory is created for each group comparison (e.i. "Homo sapiens/Case_Control"). And, the following files are saved:
#gProfiler_enrichment_res.csv -> g:Profiler enrichment analysis results filtered for Gene Ontology terms and KEGG.
#GOTERM_revigo.csv -> REVIGO's results for the specified GO terms (GO:BP, GO:CC, GO:BP).
#GOTERM_revigoCluster.csv -> clusters are found for the GO terms found after REVIGO analysis.
#GOTERM_revigoClusterHeatmap.csv -> table with ggplot heatmap information for the specified GO term.
#KEGG_heatmap.csv -> table with KEGG's heatmap input file for heatmap creation.
#heatmapPlots.pdf -> One single PDF file where the heatmap plots are saved

#Load Library
library(tidyverse)  #Collection of R packages designed for data science: ggplot2, dplyr, tidyr, readr, purrr, tibble, stringr, forcats.
library(DESeq2)     #Differential gene expression analysis based on the negative binomial distribution
library(gprofiler2) #Gene list functional enrichment analysis and namespace conversion
library(rvest)      #Wrappers around the 'xml2' and 'httr' packages to make it easy to download, then manipulate, HTML and XML.
library(tidytext)   #Text mining tasks, and plot generation are easier and more effective with tools already in wide use: dplyr, broom, tidyr, and ggplot2.
library(gridExtra)  #Provides the ability to arrange multiple grid-based plots on a page, and draw tables.
library(stringi)    #A multitude of language processing tools: pattern searching, concatenation, sorting, padding, wrapping, and many more.
library(scales)     #Graphical scales map data to aesthetics, and provide methods for automatically determining breaks and labels for axes and legends.
library(GO.db)      #Used to ....

rm(list = ls()) #Remove all saved objects

####### USER's Input ######
count_table_dir <- "/Users/salem_sueto/Desktop/Dual_RNAseq_review_methods/Aprianto/6_count_tables" #Path to the directory where the count tables are saved.
patternFile <- "host"   #Pattern used to select the right count tables
sample_group <- "/Users/salem_sueto/Desktop/Dual_RNAseq_review_methods/Aprianto/metadata/sample_groups.txt" #Path to the file where the sample's group are saved as two columns called "Sample" and "Group".

#groupComparison <- c("Pathogen4h_PathogenControl", "Pathogen24h_PathogenControl", "Pathogen48h_PathogenControl")
#groupComparison <- c("Host4h_HostControl", "Host24h_HostControl", "Host48h_HostControl")
#groupComparison <- c("InfBlood_InfHeart",  "InfBlood_Biofilm", "InfBlood_Planktonic", "InfHeart_Planktonic", "Biofilm_Planktonic")
#groupComparison <- c("cp30_cp60", "cp30_cp120", "cp60_cp120", "wt30_wt60", "wt30_wt120", "wt60_wt120", "wt30_cp30", "wt60_cp60", "wt120_cp120")  #List of group comparison: "Group1_Group2".
groupComparison <- c("cp30_cp60")
#groupComparison <- c("IP_PBS")  #List of group comparison: "Group1_Group2".
#groupComparison <- c("IP_E25", "IP_E37", "IP_S25", "IP_S37")  #List of group comparison: "Group1_Group2".
comparison_path <- "/Users/salem_sueto/Desktop/Dual_RNAseq_review_methods/Aprianto/metadata/host_paper_comparisons.txt" #Path to the file where the list of DEGs are listed for comparison for each group comparison

#DESeq2 input
cutOff <- 0                                               #Filter for genes/isoforms which has a row-sum higher than the cut-off value.
padj_FDR_value <- 0.05                                    #Adjusted p-values used to filter for differential expressed genes/isoforms.
log2FC_value <- 1                                        #Single value for the Log2 Fold Change used to select differential expressed genes/isoforms.

#g:Profiler's supported organims -> (https://biit.cs.ut.ee/gprofiler/page/organism-list)
gProfiler_supported_org <- c("Homo sapiens")                #List of Scientific Name for g:Profiler supported organisms.
#gProfiler_non_supported_org <- c("/Users/salem_sueto/Desktop/Dual_RNAseq_review_methods/00_GMT_Custom_Annotation/Mycobacterium_tuberculosis_strain_ATCC_25618_H37Rv.gmt")  #List of path to GMT custom annotation files for non-supported organisms.
gProfiler_non_supported_org <- c("")

#REVIGO's Input format
cutoff <- "0.40"                                            #Allowed values: "0.90" "0.70" "0.50" "0.40" 
isPValue <- "yes"                                           #Allowed values: "yes"  "no"
whatIsBetter <- "higher"                                    #Allowed values: "higher" "lower" "absolute" "abs_log"
goSizes <- "0"                                              #("whole UniProt")
measure <- "SIMREL"                                         #Allowed values: "RESNIK" "LIN" "SIMREL" "JIANG"
###########################

###Global Objects
baseurl <- "http://revigo.irb.hr/"                          #REVIGO website link (Do not change)
list_of_files <- list()
de_list <- list()
deg_plot <- list()
plot_all_list <- list()

#Create directory where to save the outputs
deseq2_dir_name <- sprintf("DESeq2_cutoff%s_FDR%s_LFC%s", cutOff, padj_FDR_value, log2FC_value)
dir.create(file.path(count_table_dir, deseq2_dir_name), showWarnings = FALSE)

###Load the files
sample_group_df <- read.table(sample_group, header=TRUE, sep = "\t")  #Get the sample and its group membership file
comparison_df <- read.table(comparison_path, header=TRUE, sep = "\t")

#Get the count tables
files <- list.files(path=count_table_dir, pattern=patternFile)
for (i in files) {list_of_files[[i]] <- read.table(sprintf("%s/%s", count_table_dir, i), header=TRUE, sep = "\t")}


### Differential expressed analysis ###
deseq2_res_df <- as.data.frame(matrix(data=NA, nrow = 0, ncol = 9))
deg_summary <- as.data.frame(matrix(data=NA, nrow = 0, ncol = 7))
deg_paper_summary <- as.data.frame(matrix(data=NA, nrow = 0, ncol = 3))

for (iCT in 1:length(list_of_files)) #Loop through the count tables
{
  #Tidy the Count Table
  nowTable <- list_of_files[[iCT]] %>% drop_na()
  rowID <- nowTable[,1]
  nowTable <- nowTable %>% dplyr::select (-c(colnames(nowTable)[1])) %>% mutate_if (is.character,as.numeric)
  rownames(nowTable) <- rowID
  
  for (iGC in 1:length(groupComparison))   #Loop through the group comparisons
  {
    print(sprintf("Table: %s * Group Comparison: %s", files[iCT], groupComparison[iGC]))
    
    #Get the info for the two groups
    groupInfo <- str_split(groupComparison[iGC], pattern="_")[[1]]
    
    if(length(groupInfo)==2)
    {
      #Select the right samples for Group1 based on the number of replicates
      group1Samples <- sample_group_df %>% filter(Group==groupInfo[1])
      group1Samples <- group1Samples$Sample
      
      #Select the right samples for Group2 based on the number of replicates
      group2Samples <- sample_group_df %>% filter(Group==groupInfo[2])
      group2Samples <- group2Samples$Sample
      
      #Select the samples selecteed and tidy the table
      compaDF <- nowTable %>% 
        dplyr::select(all_of(c(group1Samples, group2Samples))) %>% 
        rownames_to_column() %>%
        mutate(Sum=rowSums(dplyr::select(., -starts_with("rowname")))) %>% 
        filter(Sum >= cutOff) %>% 
        dplyr::select(-Sum) %>%
        column_to_rownames()
      
      # Group Factors
      sampleGroups <- append (rep(1, length(group1Samples)), rep(2, length(group2Samples)))
      deseq2_factors <- factor(sampleGroups)             #DESeq2's Factors
      
      #------- DESeq2 ------#
      matrixCountTable <- round (as.matrix(compaDF))
      coldata <- data.frame(row.names=colnames(matrixCountTable), deseq2_factors)
      dds <- DESeqDataSetFromMatrix(countData=matrixCountTable, colData=coldata, design=~deseq2_factors)
      dds <- DESeq(dds)
      res <- results(dds)
      res <- as.data.frame(res)
      
      deseq2_deg <- c()
      deseq2_deg <- res %>% filter(padj < padj_FDR_value, abs(log2FoldChange) >= log2FC_value)
      de_list[[sprintf("%s * %s", files[iCT], groupComparison[iGC])]] <- rownames(deseq2_deg)
      
      #DESeq2 info summary
      deseq2_deg$CountTable <- files[iCT]
      deseq2_deg$Comparison <- groupComparison[iGC]
      deseq2_deg <- deseq2_deg %>% rownames_to_column(var = "ID")
      deseq2_res_df <- rbind(deseq2_res_df, deseq2_deg)
      
      #DESeq2 DEGs Comparison summary
      paper_deg <- comparison_df [,grep(groupComparison[iGC], colnames(comparison_df))]
      paper_deg <- paper_deg[paper_deg != ""]
      
      tp <- intersect(paper_deg, deseq2_deg$ID)
      fp <- setdiff(deseq2_deg$ID, paper_deg)
      fn <- setdiff(paper_deg, deseq2_deg$ID)
      tpr <- length(tp)/length(paper_deg)
      deg_summary <- rbind(deg_summary, c(files[iCT], groupComparison[iGC], nrow(deseq2_deg), length(tp), tpr, length(fp), length(fn)))
      
      #Comparison DE elements summary
      deg_paper_summary <- rbind(deg_paper_summary, c("Paper", groupComparison[iGC], length(paper_deg)))
      
      #Add the comparison list of DEGs for this Group Comparison
      de_list[[sprintf("Paper * %s", groupComparison[iGC])]] <- paper_deg[paper_deg != ""]
    }
  }
}

#Save the DESeq2 results
write.table(deseq2_res_df, file=sprintf("%s/%s/DEGs_%s_info.txt", count_table_dir, deseq2_dir_name, sprintf("%s_results", patternFile), patternFile), row.names=FALSE, sep="\t", quote = FALSE)

colnames(deg_summary) <- c("Table", "Comparison", "DEGs", "Intersect", "IntersectRate", "DistinctPipeline", "DistinctPaper")
deg_summary <- deg_summary %>% distinct()
write.table(deg_summary, file=sprintf("%s/%s/DEGs_%s_summary.txt", count_table_dir, deseq2_dir_name, sprintf("%s_results", patternFile), patternFile), row.names=FALSE, sep="\t", quote = FALSE)

colnames(deg_paper_summary) <- c("Table", "Comparison", "DEGs")

#Save DE elements plots
deg_number_summary <- rbind(deg_summary[,1:3], deg_paper_summary)
deg_number_summary$Comparison <- factor(deg_number_summary$Comparison, levels=rev(groupComparison))

deg_plot[["Number of DEs 1"]] <- deg_number_summary%>%
  mutate_at(3, as.integer) %>%
  ggplot(aes(x=Table, y=Comparison, fill=DEGs)) +
  geom_raster(aes(fill = DEGs)) +
  scale_fill_continuous(high = "firebrick1", low = "ghostwhite", na.value="gray90") +
  geom_text(size=4, aes(label=DEGs), fontface = "bold") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 10, face = "bold")) +
  theme(axis.text.y = element_text(size=10, face="bold")) +
  labs(title="Number of DEGs 1", x="", y="") +
  coord_fixed(ratio=0.3)

deg_plot[["Number of DEs 2"]] <- deg_number_summary%>%
  mutate_at(3, as.integer) %>%
  ggplot(aes(x=Table, y=Comparison, fill=DEGs)) +
  geom_raster(aes(fill = DEGs)) +
  scale_fill_continuous(high = "firebrick1", low = "ghostwhite", na.value="gray90") +
  geom_text(size=4, aes(label=DEGs), fontface = "bold") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 10, face = "bold")) +
  theme(axis.text.y = element_text(size=10, face="bold")) +
  labs(title="Number of DEGs 2", x="", y="") +
  coord_fixed(ratio=0.3) +
  coord_flip()

deg_summary$Comparison <- factor(deg_summary$Comparison, levels=groupComparison)

deg_plot[["DEs_TPR_FP"]] <- deg_summary %>%
  mutate_at(3:7, as.numeric) %>%
  ggplot(aes(x=DistinctPipeline, y=IntersectRate, shape=Table, color=Comparison)) + 
  geom_point(size=4) +
  #scale_y_continuous(limits=c(0, 1), breaks = seq(0,1,0.2)) +
  labs(x="Number Distinct Pipeline DE elements", y="Intersect Rate", title="Number of Identified DE elements")

##### Save Heatmap plots #####
save(deg_plot, file = sprintf("%s/%s/DEGs_%s_summary_plots.RData", count_table_dir, deseq2_dir_name, sprintf("%s_results", patternFile), patternFile))

pdf(sprintf("%s/%s/DEGs_%s_summary_plots.pdf", count_table_dir, deseq2_dir_name, sprintf("%s_results", patternFile), patternFile), onefile=TRUE)
for (i in seq(length(deg_plot))) {grid.arrange (deg_plot[[i]])}
dev.off()


####### Enrichment Analysis g:Profiler & REVIGO ######
list_of_org <- c(gProfiler_supported_org, gProfiler_non_supported_org)
list_of_org <- list_of_org[list_of_org != ""]

#Loop through the group comparisons
for (iGC in groupComparison)
{
  #Get the elements from de_list which has the same group comparison
  de_gc <- de_list[grep(iGC, names(de_list))]
  
  #Loop through each species in the list_of_org vector list
  for (iSpecies in list_of_org)
  {
    ###Enrichment Analysis -> g:Profiler
    gProfiler_Res <- NULL
    gmt_CA_presence <- FALSE #Use to find whether the code is using a custom annotation for the species
    
    #Check if the iSpecies is a Scientific name and supported by g:Profiler
    if (length(grep(".gmt", iSpecies))==0)
    {
      gPr_iSpecies <- paste(substr (tolower(iSpecies), 1, 1), str_split (tolower(iSpecies), pattern = " ")[[1]][2], sep = "")
      
      tryCatch (
        {
          gProfiler_Res <- gost(de_gc, organism = gPr_iSpecies, ordered_query = FALSE,
                                multi_query = FALSE, significant = TRUE, exclude_iea = TRUE,
                                measure_underrepresentation = FALSE, evcodes = FALSE,
                                user_threshold = 0.05, correction_method = "gSCS",
                                domain_scope = "annotated", custom_bg = NULL,
                                numeric_ns = "", sources = NULL)
        },
        error = function(e) 
        {
          message(e)
          return(gProfiler_Res <- NULL)
        }
      )
    }
    
    #Check if the iSpecies is a custom annotation file
    if (length(grep(".gmt", iSpecies))>0)
    {
      gmt_CA_presence <- TRUE
      gmt_file <- upload_GMT_file(iSpecies)
    
      tryCatch (
        {
          gProfiler_Res <- gost(de_gc, organism = gmt_file, ordered_query = FALSE,
                                multi_query = FALSE, significant = TRUE, exclude_iea = TRUE,
                                measure_underrepresentation = FALSE, evcodes = FALSE,
                                user_threshold = 0.05, correction_method = "gSCS",
                                domain_scope = "annotated", custom_bg = NULL,
                                numeric_ns = "", sources = NULL)
        },
        error = function(e) 
        {
          message(e)
          return(gProfiler_Res <- NULL)
        }
      )
    }
    
    #Analyse the g:Profiler enrichment analysis
    if (is.null(gProfiler_Res) == FALSE)
    {
      #Plot List
      enrichPlot <- list()
      
      #Create Directories
      dir.create(file.path(count_table_dir, tail (str_split(iSpecies, pattern = "/")[[1]],1)), showWarnings = FALSE)  #Create directory for the species
      dir.create(file.path(sprintf("%s/%s", count_table_dir, tail (str_split(iSpecies, pattern = "/")[[1]],1)), iGC), showWarnings = FALSE)  #Create directory for the Group Comparison
      res_directory <- sprintf("%s/%s/%s", count_table_dir, tail (str_split(iSpecies, pattern = "/")[[1]],1), iGC)
      
      #Save the g:Profiler enrichment analysis results for GO and KEGG terms
      enrich_go_kegg <- gProfiler_Res$result[,c(1,3,9,10,11)]
      
      #In case of custom annotation find the source of the Gene Ontology terms
      if (gmt_CA_presence==TRUE)
      {
        enrich_go_kegg$species <- enrich_go_kegg$source
        enrich_go_kegg$source <- sprintf("GO:%s", unname (Ontology(enrich_go_kegg$term_id)))
      }
      
      enrich_go_kegg <- enrich_go_kegg %>% filter(source %in% c("GO:BP", "GO:CC", "GO:MF", "KEGG"))
      write.table(enrich_go_kegg, file = sprintf("%s/gProfiler_enrichment_res.txt", res_directory), row.names=FALSE, sep="\t", quote = FALSE)
      
      ##### REVIGO Analysis #####
      revigo_session <- html_session(baseurl)
      revigo_form <- html_form(revigo_session)[[1]] 
      
      # Loop through the GO terms
      goTerms <- unique (enrich_go_kegg$source)
      goTerms <- goTerms[goTerms != "KEGG"]
      
      for (go in goTerms) 
      {
        #Number of Cluster
        nCluster <- 40
        
        #Gene Ontology
        filterGO <- enrich_go_kegg %>% filter (source==go)
        goList <- paste(unique(filterGO$term_id), collapse="\n")
        
        filled_form <- set_values(revigo_form, 'goList'= goList, 
                                  'cutoff'= cutoff, 'isPValue'= isPValue, 
                                  'whatIsBetter'= whatIsBetter, 'goSizes'= goSizes, 'measure'= measure)
        
        #Get the Revigo Summary Table & Save it as a CSV file
        result_page <- submit_form(revigo_session, filled_form, submit='startRevigo')
        GORevigo <- result_page %>% follow_link("Export results to text table (CSV)")
        httr::content(GORevigo$response, as="text") %>% write_lines("revigoSummaryOnline.csv") 
        
        #Get the REVIGO summary result from the saved file
        revigoSummary <- read.csv2("revigoSummaryOnline.csv", header=TRUE, sep=",")
        revigoSummary$frequency <- gsub('%', '', revigoSummary$frequency)
        revigoSummary$frequency <- as.numeric( as.character(revigoSummary$frequency) );
        
        #Get only the GO terms represented in the "semantic space"
        uniqueGO <- revigoSummary %>% filter(plot_X != "null")
        iNum <- c(4,5,6,7,8) # Change the column into numeric: plot_X - plot_Y - plot_size - uniqueness - dispensability
        uniqueGO [ , iNum] <- apply(uniqueGO [ , iNum], 2, function(x) as.numeric(as.character(x)))
        
        # Add "HeadGO" Column which shows which GO REVIGO term has been put as the main GO for that group of terms
        revigoSummary$HeadGO <- NA
        headGO <- "NA"
        
        for (i in 1:nrow(revigoSummary)) 
        {
          if(revigoSummary$plot_X[i]!="null") {headGO <- revigoSummary$term_ID[i]}
          revigoSummary$HeadGO[i] <- headGO
        }
        write.table(enrich_go_kegg, file = sprintf("%s/%s_revigo.txt", res_directory, go), row.names=FALSE, sep="\t", quote = FALSE) #Save the REVIGO result
        
        ##### GENE ONTOLOGY REVIGO CLUSTER ANALYSIS
        
        #Find the clusters present in the REVIGO results based on their sematic space
        clusterDF <- uniqueGO %>% remove_rownames %>% column_to_rownames(var="term_ID")
        
        if (nCluster >= nrow(clusterDF)) #In case the nCluster is equal/higher than the number of GO term elements
        {
          if (nrow(clusterDF)==1) {nCluster <- nrow(clusterDF)}   #In case the GO terms are equal to 1
          if (nrow(clusterDF)>1) {nCluster <- nrow(clusterDF)-1}  #In case the GO terms are higher than 1
        }
        
        set.seed(100)
        clusters <- kmeans(clusterDF[,c(3,4)], nCluster)
        clusterDF$Cluster <- clusters$cluster
        
        # Group elements from the same cluster
        revigoCluster <- data.frame(matrix(ncol=10, nrow=nCluster))
        colnames(revigoCluster) <- c("Cluster", "PlotX", "PlotY", "RevigoGOs", "RevigoRep", "RevigoDescription", "AllGOs", "AllWords", "AllDescription", "FinalDescription")
        
        for (i in sort(unique(clusterDF$Cluster)))
        {
          filterGO_cluster <- clusterDF %>% filter(Cluster==i)
          
          # REVIGO GOs
          revigoCluster$Cluster[i] <- i
          revigoCluster$PlotX[i] <- mean(filterGO_cluster$plot_X)
          revigoCluster$PlotY[i] <- mean(filterGO_cluster$plot_Y)
          revigoCluster$RevigoGOs[i] <- paste (rownames(filterGO_cluster), collapse=" ")
          revigoCluster$RevigoRep[i] <- filterGO_cluster[which.max(filterGO_cluster$representative),][1,1] #Term with the max value of "representative of the cluster
          
          # REVIGO Cluster Description -> Concagenate the HEADs' term description of the cluster - Split - Delete generic terms (i.e. "of", "a", "an")
          wordOcc <- paste(filterGO_cluster$description, collapse=" ") %>% str_split(pattern = " ") 
          wordOcc <- wordOcc[[1]]
          indexWords <- !(wordOcc %in% stop_words$word)
          wordOcc <- as.data.frame(table (wordOcc[indexWords])) %>% arrange(desc (Freq))
          wordOcc <- as.vector(wordOcc$Var1[1:5])
          revigoCluster$RevigoDescription[i] <- paste (wordOcc[!is.na(wordOcc)], collapse = " ")
          
          #All GOs Columns
          filterAllGO <- revigoSummary %>% filter(HeadGO %in% rownames(filterGO_cluster))
          
          revigoCluster$AllGOs[i] <- paste (filterAllGO$term_ID, collapse=" ")
          revigoCluster$AllWords[i] <- paste(filterAllGO$description, collapse=" ")
          
          wordOcc <- paste(filterAllGO$description, collapse=" ") %>% str_split(pattern = " ") 
          wordOcc <- wordOcc[[1]]
          indexWords <- !(wordOcc %in% stop_words$word)
          wordOcc <- as.data.frame(table (wordOcc[indexWords])) %>% arrange(desc (Freq))
          wordOcc <- as.vector(wordOcc$Var1[1:5])
          revigoCluster$AllDescription[i] <- paste (wordOcc[!is.na(wordOcc)], collapse = " ")
          
          #Final Description of the Cluster
          revigoCluster$FinalDescription[i] <- paste(revigoCluster$RevigoRep[i], revigoCluster$AllDescription[i], nrow(filterAllGO), sep = " * ")
        }
        write.table(revigoCluster, file = sprintf("%s/%s_revigoCluster.txt", res_directory, go), row.names=FALSE, sep="\t", quote = FALSE) #Save the REVIGO_Cluster result
        
        #Get the lowest pvalue of the group comparison's GO terms and the cluster combination
        heatmapGO <- data.frame(matrix(ncol = length (unique(filterGO$query)), nrow = nrow(revigoCluster)))
        colnames(heatmapGO) <- unique(filterGO$query)
        rownames(heatmapGO) <- revigoCluster$FinalDescription
        
        for (n in 1:nrow(revigoCluster)) 
        {
          allGO <- str_split(revigoCluster$AllGOs[n], " ") [[1]] #Get all GOs 
          filterGOGroup <- filterGO %>% filter(term_id %in% allGO)  #Get all the rows which have those GOs
          
          uniqueGroup <- unique (filterGOGroup$query)
          for (k in uniqueGroup) 
          {
            pvaluesList <- filterGOGroup %>% filter(query == k)
            heatmapGO[n, which(colnames(heatmapGO)==k)] <- min(pvaluesList$p_value)
          }
        }
        
        #Create the Heatmap plot
        GroupName <- c()
        for (i in 1:nrow(heatmapGO)) 
        {
          groupPresence <- colnames(heatmapGO)[sapply(heatmapGO[i,], function(x) any(!is.na(x)))]
          groupPresence <- paste(groupPresence, collapse = " * ")
          GroupName <- append(GroupName, groupPresence)
        }
        
        heatmapGO$GroupName <- GroupName
        heatmapGO$NGroups <- length (unique(filterGO$query)) - rowSums(is.na(heatmapGO))
        heatmapGO <- heatmapGO %>% arrange(NGroups, GroupName) %>% dplyr::select(-c(NGroups, GroupName))
        
        heatmapGO <- tibble::rownames_to_column(heatmapGO, "Description")
        heatmapGO_melt <- reshape2::melt(heatmapGO,id.vars=c("Description"))
        heatmapGO_melt$Description <- factor(heatmapGO_melt$Description, levels = heatmapGO$Description)
        
        #Save data
        write.table(heatmapGO_melt, file = sprintf("%s/%s_revigoClusterHeatmap.txt", res_directory, go), row.names=FALSE, sep="\t", quote = FALSE)   #Save the REVIGO_Cluster_Heatmap result
        
        plotName <- sprintf("%s_heatmap", go)
        enrichPlot[[plotName]] <- heatmapGO_melt %>% 
          ggplot(aes(x=variable, y=factor(Description), fill=value)) +
          geom_raster(aes(fill = value)) +
          scale_fill_continuous(high = "ghostwhite", low = "firebrick1", na.value="gray90") +
          geom_text(size=3, aes(label=scientific(value, digits = 3))) +
          theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 10, face = "bold")) +
          theme(axis.text.y = element_text(size=10, face="bold")) +
          coord_fixed(ratio=ncol(heatmapGO)/nrow(heatmapGO)) +
          ggtitle(plotName)
        
        ##### Delete the "revigoSummaryOnline.csv" file
        if (file.exists("revigoSummaryOnline.csv")) {file.remove("revigoSummaryOnline.csv")}
      }
      
      ##### KEGG Heatmap creation #####
      filterKEGG <- filter(enrich_go_kegg, source=="KEGG")
      
      #Get KEGG heatmap in case there are any KEGG enriched term
      if (nrow(filterKEGG)>0)
      {
        heatmapKEGG <- data.frame(matrix(ncol = length (unique(filterKEGG$query)), nrow = length (unique (filterKEGG$term_id))))
        colnames(heatmapKEGG) <- unique(filterKEGG$query)
        rownames(heatmapKEGG) <- unique (filterKEGG$term_name)
        
        for (i in 1:nrow(filterKEGG)) 
        {heatmapKEGG[filterKEGG$term_name[i], filterKEGG$query[i]] <- filterKEGG$p_value[i]}
        write.table(heatmapKEGG, file = sprintf("%s/KEGG_heatmap.txt", res_directory), row.names=FALSE, sep="\t", quote = FALSE)   #Save the gProfiler_Heatmap result
        
        #Create the Heatmap plot
        GroupName <- c()
        for (i in 1:nrow(heatmapKEGG)) 
        {
          groupPresence <- colnames(heatmapKEGG)[sapply(heatmapKEGG[i,], function(x) any(!is.na(x)))]
          groupPresence <- paste(groupPresence, collapse = " * ")
          GroupName <- append(GroupName, groupPresence)
        }
        
        heatmapKEGG$GroupName <- GroupName
        heatmapKEGG$NGroups <- length (unique(filterGO$query)) - rowSums(is.na(heatmapKEGG))
        heatmapKEGG <- heatmapKEGG %>% arrange(NGroups, GroupName) %>% dplyr::select(-c(NGroups, GroupName))
        
        heatmapKEGG <- tibble::rownames_to_column(heatmapKEGG, "Description")
        heatmapKEGG_melt <- reshape2::melt(heatmapKEGG,id.vars=c("Description"))
        heatmapKEGG_melt$Description <- factor(heatmapKEGG_melt$Description, levels = heatmapKEGG$Description)
        heatValues <- heatmapKEGG_melt$value [!heatmapKEGG_melt$value %in% c(NA)]
        
        enrichPlot[["KEGG"]] <- heatmapKEGG_melt %>% 
          ggplot(aes(x=variable, y=Description, fill=value)) +
          geom_raster(aes(fill = value)) +
          scale_fill_continuous(high = "ghostwhite", low = "firebrick1", na.value="gray90") +
          geom_text(size=3, aes(label=scientific(value, digits = 3))) +
          theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 10, face = "bold")) +
          theme(axis.text.y = element_text(size=10, face="bold")) +
          coord_fixed(ratio=ncol(heatmapGO)/nrow(heatmapGO)) +
          ggtitle("KEGG_heatmap")
      }

      ##### Save Heatmap plots #####
      plot_all_list[[iGC]] <- enrichPlot
      save(enrichPlot, file = sprintf("%s/enrichment_summary_plots.RData", res_directory))
      
      pdf(sprintf("%s/heatmapPlots.pdf", res_directory), onefile=TRUE)
      for (i in seq(length(enrichPlot))) {grid.arrange (enrichPlot[[i]])}
      dev.off()
    }
  }
}
