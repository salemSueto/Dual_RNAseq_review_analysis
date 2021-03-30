#!/taxonID_search.R

#The script is used to identify the taxon ID of the species for whom the custom annotation will be used for. 
#Assign the scientific name of the species to the "org" character vector. 
#And save the taxon ID of the selected organisms.

#Load Library
library(UniProt.ws) #R Interface to UniProt Web Services

###USER's Input###
org <- "Streptococcus pneumoniae"
##################

availableUniprotSpecies(pattern=org, n=Inf) #Check for all instances of the query in the Uniprot database
