#!/gene_transcript_info.R

#The script takes as input the path to directories where the transcriptome (.cdna.all.fa) and the annotation (*.gtf) files of the species are saved. 
#Then, it outputs two text files which displays gene-transcript relationship and the gene/transcript length; 
#whereas, one pdf file which shows the relation between the two input types. 
#The script creates a directory called "txome_annotation" inside the input-directory where the output files are saved.

#The input are:
#1) dir_list -> vector list of path to the directory where the transcriptome and the annotation files are saved.

#The output are:
#1) tx_gene_relation.txt -> find the gene and its transcript isoforms. The columnns are TranscriptID - GeneID - GTF - CDNA
#The GTF and CDNA columns shows whether the transcript is found in the file.
#2) tx_gene_length.txt -> find the length  for each elements. 
#Transcript length by checking the length from the transcriptome file. 
#Whereas, gene length is the sum of all the gene's transcript isoforms length.
#3) cdna_gtf_plots.pdf -> pdf file which shows the transcriptome and annotation venn groups.

#Load Library
library(rtracklayer)    #R interface to genome annotation files and the UCSC genome browser
library(tidyverse)      #Collection of R packages designed for data science: ggplot2, dplyr, tidyr, readr, purrr, tibble, stringr, forcats.
library(Biostrings)     #Fast manipulation of large biological sequences or sets of sequences.
library(ggVennDiagram)  #A 'ggplot2' Implement of Venn Diagram
library(rlang)          #Functions for Base Types and Core R and 'Tidyverse' Features
library(ggpubr)         #'ggplot2' Based Publication Ready Plots


rm(list = ls()) #Remove all saved objects

#User's INPUT
dir_list <- c("Hsapiens_Spneumoniae_SingleStrand_SE/reference/spneumoniae", "Hsapiens_Spneumoniae_SingleStrand_SE/reference/hsapiens")
######

#Functions (https://github.com/lambdamoses/blog) - (https://fromsystosys.netlify.app/2020/01/31/comparing-ensembl-gtf-and-cdna/)

plot_bar <- function(df, col_fill, name) {
  ggplot(df, aes(fct_reorder(gene_biotype, n, .fun = sum), n,
                 fill = .data[[col_fill]])) +
    geom_bar(position = position_dodge2(preserve = "single", width = 0.9), stat = "identity") +
    geom_text(aes(label = n, color = .data[[col_fill]]), hjust = -0.05,
              position = position_dodge2(preserve = "single", width = 0.9)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    scale_fill_discrete(name = name) +
    scale_color_discrete(name = name) +
    coord_flip() +
    theme(axis.title = element_blank(), legend.margin = margin(t = 14))
}

plot_bar_patch <- function(df, n_split, col_fill, name, title) {
  df <- df %>% 
    left_join(df %>% group_by(gene_biotype) %>% summarize(tot = sum(n)), by = "gene_biotype")
  df <- df %>% 
    mutate(bin = cut(tot, scales::breaks_log(n = n_split + 1)(range(df$tot))))
  dfs <- split(df, fct_rev(df$bin))
  plts <- map(dfs, plot_bar, col_fill = col_fill, name = name)
  hts <- rev(unname(table(df$bin)))
  hts <- hts + 8/min(hts)
  p <- ggarrange(plotlist = plts, ncol = 1, heights = hts,
                 common.legend = TRUE, legend = "top", align = "v")
  annotate_figure(p, fig.lab = title, fig.lab.pos = "top.left", fig.lab.size = 14)
}

##### Get CDNA and GTF information #####
for (iDir in dir_list) 
{
  #Check that the GTF and the CDNA files are from the same organism
  gtf_file_path <- sprintf("%s/%s", iDir, list.files(path=iDir, pattern=".gtf")[1])
  gtf_filename <- paste(head(str_split(basename(gtf_file_path), pattern = "\\.")[[1]], -2), collapse = ".")
  
  cdna_file_path <- sprintf("%s/%s", iDir, list.files(path=iDir, pattern=".cdna.all.fa")[1])
  cdna_filename <- str_replace(basename(cdna_file_path), pattern = ".cdna.all.fa", "")
  
  if (gtf_filename == cdna_filename)
  {
    #Create output directory "txome_annotation" where to save all the output files
    dir.create(file.path(iDir, "txome_annotation"))
    
    #Plot List
    gtf_cdna_plot_list <- list()
    
    ### GTF Analysis ###
    
    #Load the GTF file
    gtf <- rtracklayer::import(gtf_file_path)
    gtf_df <- as.data.frame(gtf)
    
    #Filter rows and columns
    gtf_filter_df <- gtf_df %>% 
      dplyr::filter(type=="transcript") %>%                                         #Filter for only transcript
      dplyr::select(transcript_id, gene_id, width) %>%                              #Select the columns
      dplyr::rename(TranscriptID=transcript_id, GeneID=gene_id, Length=width) %>%   #Rename the columns
      distinct(TranscriptID, GeneID, .keep_all = TRUE)                              #Remove eventual duplicates
    
    #Remove the substring "transcript:" from the transcript_id column, especially in bacteria annotation file.
    if (grepl(pattern = "transcript:", gtf_filter_df$TranscriptID[1], fixed=T)==TRUE)
    {gtf_filter_df$TranscriptID <- str_replace(gtf_filter_df$TranscriptID, "transcript:", "")}
    
    #Output: Gene and its isoforms
    gtf_rel <- gtf_filter_df[,c("TranscriptID", "GeneID")]
    
    #Output: Gene and Transcript length (ID - Length)
    gene_length <- gtf_filter_df %>% group_by(GeneID) %>% summarise(Length = sum(Length))  #Gene Length -> sum of its isoforms length
    tx_length <- gtf_filter_df[,c(1,3)] #Transcript Length
    colnames(gene_length) <- c("ID", "Length")
    colnames(tx_length) <- c("ID", "Length")
    gtf_length <- rbind(tx_length, gene_length)
    
    
    ### CDNA Analysis ###
    
    #Load the transcriptome file
    cdna_fa <- readDNAStringSet(cdna_file_path) #Read the cdna fasta file
    
    #Create dataframe for the cdna
    TranscriptID <- str_replace (sapply (str_split (names (cdna_fa), pattern=" "), `[`, 1), "\\..*", "") #Get Transcript ID for each sequence
    GeneID <- str_replace (sapply (str_split (sapply (str_split (names (cdna_fa), pattern=" gene:"), `[`, 2), pattern=" "), `[`, 1), "\\..*", "")
    Length <- width(cdna_fa)
    Sequence <- unname (as.character(cdna_fa))
    cdna_df <- as.data.frame(cbind(TranscriptID, GeneID, Length, Sequence))
    cdna_df$Length <- as.numeric(cdna_df$Length)
    
    #Output
    cdna_rel <- cdna_df %>% dplyr::select(TranscriptID, GeneID) %>% distinct(TranscriptID, GeneID, .keep_all = TRUE)
    
    cdna_gene_length <- cdna_df %>% group_by(GeneID) %>% summarise(Length = sum(Length))
    cdna_tx_length <- cdna_df[,c(1,3)]
    colnames(cdna_gene_length) <- c("ID", "Length")
    colnames(cdna_tx_length) <- c("ID", "Length")
    cdna_length <- rbind(cdna_gene_length, cdna_tx_length)
    
    
    ##### GTF - CDNA transcript comparison #####
    
    ###Venn Diagram
    cdna_tx <- cdna_df$TranscriptID #CDNA transcript IDs
    gtf_tx <- unique(gtf_df$transcript_id) %>% str_replace(pattern="transcript:", replacement="") #GTF transcript IDs
    gtf_tx <- gtf_tx[!is.na(gtf_tx)]
    
    #TX Venn groups
    inter_tx <- intersect(cdna_tx, gtf_tx)
    gtf_only_tx <- setdiff(gtf_tx, cdna_tx)
    cdna_only_tx <- setdiff(cdna_tx, gtf_tx)
    
    tx_venn <- list(cDNA=cdna_tx, GTF=gtf_tx)
    gtf_cdna_plot_list[["Venn_Transcript"]] <- ggVennDiagram(tx_venn) + 
      scale_fill_gradient(low="yellow",high = "red") + labs (title="Transcript Venn Diagram", size = 18)
    
    #Gene Venn groups
    inter_gene <- intersect(unique(cdna_df$GeneID), unique(gtf_df$gene_id))
    cdna_gene <- setdiff(unique(cdna_df$GeneID), (unique(gtf_df$gene_id)))
    gtf_gene <- setdiff(unique(gtf_df$gene_id), unique(cdna_df$GeneID))
    
    gtf_cdna_plot_list[["Venn_Gene"]] <- ggVennDiagram(list (cDNA= unique(cdna_df$GeneID), GTF= unique(gtf_df$gene_id))) + 
      scale_fill_gradient(low="yellow",high = "red") + 
      labs (title="Gene Venn Diagram", size = 18)
    
    
    ###GTF-only Trasncripts analysis
    gtf_meta <- as.data.frame(gtf[gtf$type == "transcript"])
    gtf_meta <- gtf_meta %>% mutate(gtf_only = !transcript_id %in% cdna_tx)
    
    #Add the "gene_biotype" column in case it is not present in gtf_meta object (e.i. bacterial) and set it equal to column "type"
    if (length (grep(pattern="gene_biotype", colnames(gtf_meta))) == 0) {gtf_meta$gene_biotype <- "protein_coding"} 
    
    n_txs_gtf <- gtf_meta %>% group_by(gtf_only, gene_biotype) %>% summarise(n = n())
    if (length (unique (n_txs_gtf$gene_biotype)) > 10) {nDiv <- 3}
    if (length (unique (n_txs_gtf$gene_biotype)) < 10) {nDiv <- length (unique (n_txs_gtf$gene_biotype))}
    
    gtf_cdna_plot_list[["GTF_only_number"]] <- plot_bar_patch(n_txs_gtf, nDiv, "gtf_only", "GTF only", "Number of transcripts in each gene biotype in GTF file")
    
    p <- ggplot(gtf_meta, aes(fct_reorder(gene_biotype, gtf_only, .fun = mean), fill = gtf_only)) +
      geom_bar(position = "fill", alpha = 0.5) +
      scale_y_continuous(expand = expansion(mult = c(0, 0))) +
      scale_fill_discrete(name = "GTF only") +
      labs(x = "gene biotype", y = "proportion") +
      coord_flip() + theme_bw() + labs(x = "") +
      theme(legend.position = "top", legend.justification = c(0,0.5), legend.margin = margin(t = 14)) +
      labs (title="Proportion of GTF only transcripts in each gene biotype in GTF", size = 18)
    gtf_cdna_plot_list[["GTF_only_proportion"]] <- p 
    
    
    ###CDNA only Trasncripts analysis
    cdna_meta <- tibble (transcript_id = cdna_tx,
                         cr = str_extract(names(cdna_fa), "(?<=((chromosome)|(scaffold)):GRCh38:).*?(?=\\s)"),
                         gene_biotype = str_extract(names(cdna_fa), "(?<=gene_biotype:).*?(?=\\s)"),
                         gene_id = str_extract(names(cdna_fa), "(?<=gene:).*?(?=\\.)"),
                         gene_symbol = str_extract(names(cdna_fa), "(?<=gene_symbol:).*?(?=\\s)"),
                         cdna_only = !transcript_id %in% gtf_tx) %>% 
      separate(cr, into = c("seqnames", "start", "end", "strand"), sep = ":") %>% 
      mutate(start = as.integer(start),
             end = as.integer(end),
             strand = case_when(
               strand == "1" ~ "+",
               strand == "-1" ~ "-",
               TRUE ~ "*"
             ))
    
    n_txs_cdna <- cdna_meta %>% group_by(cdna_only, gene_biotype) %>% summarise(n = n())
    if (length (unique (n_txs_cdna$gene_biotype)) > 10) {nDiv <- 3}
    if (length (unique (n_txs_cdna$gene_biotype)) < 10) {nDiv <- length (unique (n_txs_gtf$gene_biotype))}
    
    gtf_cdna_plot_list[["CDNA_only_number"]] <- plot_bar_patch(n_txs_cdna, nDiv, col_fill = "cdna_only", name = "cDNA only", title = "Number of transcripts in each gene biotype in cDNA fasta")
    
    p <- ggplot(cdna_meta, aes(fct_reorder(gene_biotype, cdna_only, .fun = mean), fill = cdna_only)) +
      geom_bar(position = "fill", alpha = 0.5) +
      scale_y_continuous(expand = expansion(mult = c(0, 0))) +
      scale_fill_discrete(name = "cDNA only") +
      labs(x = "gene biotype", y = "proportion") +
      coord_flip() + theme_bw() + labs(x = "") +
      theme(legend.position = "top", legend.justification = c(0,0.5), legend.margin = margin(t = 14)) +
      labs (title="Proportion of cDNA only transcripts in each gene biotype in cDNA fasta", size = 18)
    gtf_cdna_plot_list[["CDNA_only_proportion"]] <- p
    
    
    #Output: GENE - ISOFORM relation
    gene_isoform_rel <- rbind(gtf_rel, cdna_rel) %>% distinct(TranscriptID, GeneID, .keep_all = TRUE)
    gene_isoform_rel <- gene_isoform_rel %>% mutate(GTF = TranscriptID %in% gtf_tx)
    gene_isoform_rel <- gene_isoform_rel %>% mutate(CDNA = TranscriptID %in% cdna_tx)
    write.table(gene_isoform_rel, file= sprintf("%s/txome_annotation/tx_gene_relation.txt", iDir), row.names=F, quote=F, sep="\t")
    
    #OUTPUT: Transcript length
    inter_gtf_cdna <- cdna_length %>% filter(ID %in% inter_tx)
    gtf_only <- gtf_length %>% filter (ID %in% gtf_only_tx)
    cdna_only <- cdna_length %>% filter (ID %in% cdna_only_tx)
    tx_length <- rbind(inter_gtf_cdna, gtf_only, cdna_only)
    
    #OUTPUT: Gene length
    inter_gtf_cdna <- cdna_length %>% filter(ID %in% inter_gene)
    gtf_only <- gtf_length %>% filter (ID %in% gtf_gene)
    cdna_only <- cdna_length %>% filter (ID %in% cdna_gene)
    gene_length <- rbind(inter_gtf_cdna, gtf_only, cdna_only)
    
    #Id - Length
    id_length <- rbind(tx_length, gene_length)
    write.table(id_length, file= sprintf("%s/txome_annotation/tx_gene_length.txt", iDir), row.names=F, quote=F, sep="\t")

    #Save plot in the list 
    pdf (sprintf ("%s/txome_annotation/cdna_gtf_plots.pdf", iDir))
    for (i in 1:length(gtf_cdna_plot_list)) {print(gtf_cdna_plot_list[[i]])}
    dev.off()
    
    
  } else {print(sprintf("The GTF & the CDNA files in the %s directory are not from the same organism./nPlease ensure the filenames are the un-compressed files from Ensembl.", iDir))}
}
