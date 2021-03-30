#!/gmt_custom_annotation_creation.R

#The script creates a custom annotation for the selected organism.
#It is advisable to use a custom annotation when the species in question is non-supported one.
#The custom annotation is made from Uniprot.

#Load Library
library(UniProt.ws)
library(tidyverse)

rm(list = ls()) #Remove all saved objects

###USER's INPUT###
taxonID_list <- c(488223, 189423, 487214)
outDir <- "." #Directory where to save the GMT custom annotation files. DO NOT add the last "/" to the directory's path.
#################

###Global Function -> GMT files
write.gmt=function(obj,filename)
{
  conn=file(filename,'w')
  for(i in 1:length(obj))
  {
    geneList <- paste(obj[[i]]$entry, collapse = "\t")
    cat(obj[[i]]$head, obj[[i]]$desc, geneList, file=conn, sep='\t')
    cat('\n',file=conn)
  }
  close(conn)
  return(invisible())
}

###Loop through the taxonID_list 
for (iTx in taxonID_list) 
{
  orgSelected <- UniProt.ws(taxId=iTx)   #Upload the selected organism from the Uniprot database
  
  #Check that the selected organism exist in the Uniprot database
  if(is.null(orgSelected) == FALSE)
  {
    #Name of the selected species
    iTx_name <- str_replace_all (str_replace_all(species(orgSelected), "[()]", ""), " ", "_")
    
    #Select the information from the Uniprot database
    if(interactive()) {key <- keys(orgSelected, keytype="UNIPROTKB")}
    column <- c("GO", "GO-ID", "KEGG", "GENES", "FUNCTION", "SCORE", "EXISTENCE")
    
    res <- UniProt.ws::select(orgSelected, keys=key, columns=column)
    res$GeneID <- res$KEGG
    goID <- c()
    goName <- c()
    goGenes <- c()
    
    #Tidy the downloaded information
    for (nrow in 1:nrow(res)) 
    {
      #GO & Function Columns
      if (!is.na(res$`GO-ID`[nrow])) 
      {
        #GO Columns
        goSplit <- strsplit(res$GO[nrow], "; ") [[1]]
        nWords <- lengths (str_match_all(goSplit, "\\S+" ))
        
        goID <- append(goID, word(goSplit, -1))
        goName <- append(goName, word(goSplit, start=1, end=nWords-1))
        
        goSplit <- word(goSplit, start=1, end=nWords-1)
        res$GO[nrow] <- paste(goSplit, collapse = '; ')
        
        #Function Column
        res$FUNCTION[nrow] <- strsplit(res$FUNCTION[nrow], ": ")[[1]][2]
      }
      
      #GeneID Column
      res$GeneID[nrow] <- strsplit(res$GeneID[nrow], ":")[[1]][2]
      
      #Genes Column
      nWords <- lengths (str_match_all(res$GENES[nrow], "\\S+" ))
      if (nWords > 1)
      {
        goSplit <- strsplit(res$GENES[nrow], " ") [[1]]
        goSplit <- goSplit[-nWords]
        res$GENES[nrow] <- paste(goSplit, collapse = ' ')
      }
      
      #Score Column
      res$SCORE[nrow] <- gsub (" out of ", "/", res$SCORE[1])
    }
    
    #Get the list of Genes for each GOID
    goID <- str_sub(goID, 2, -2)
    maxGOGenes <- 0
    
    if (length(goID) > 0)
    {
      for (nrow in goID) 
      {
        filterGO <- res [str_detect(res$`GO-ID`, nrow, negate = F), "GeneID"]
        filterGO <- filterGO[!is.na(filterGO)]
        if (length(filterGO) > maxGOGenes) {maxGOGenes <- length(filterGO)}
        goGenes <- append(goGenes, paste(filterGO, collapse = "; "))
      }
      
      #Create the GMT dataframe & GMT list
      GMT_species <- as.data.frame(matrix(NA, nrow = length(goID), ncol = 2 + maxGOGenes))
      GMT_species[,1] <- goID
      GMT_species[,2] <- goName
      
      GMT_speciesList <- list()
      for (geneList in 1:length(goGenes)) 
      {
        #List
        head <- goID[geneList]
        desc <- goName[geneList]
        len <- length(strsplit(goGenes[geneList], "; ")[[1]])
        entry <- strsplit(goGenes[geneList], "; ")[[1]]
        
        goSingleList <- list(head = head, desc = desc, len = len, entry = entry)
        GMT_speciesList[[head]] <- goSingleList
        
        #GMT dataframe
        splitGene <- strsplit(goGenes[geneList], "; ")[[1]]
        for (geneElement in 1:length(splitGene))
        {GMT_species[geneList,2+geneElement] <- splitGene[geneElement]}
      }
      
      GMT_species <- GMT_species %>% distinct(V1, .keep_all = TRUE)
      GMT_species <- GMT_species %>% replace(., is.na(.), "")
      
      write.gmt(GMT_speciesList, sprintf("%s/%s.gmt", outDir, iTx_name))
    }else{sprintf("The TaxonID: %s corresponding to %s has no Gene Ontology! No GMT custom annotation will be created", iTx, species(orgSelected))}
  }
}
