#!/dual_RNAseq_search_pubmed.R

#The script queries the Pubmed databbase with a list of search terms for a specific range of years.
#The script finds out how many queries are found when searching the Pubmed (2008-2019) for the following terms: 
#"Dual RNA-seq", "Parallel RNA-seq", and "Simultaneous RNA-seq".
#The user can add additional columns to the dataframe resulted from the search.

#The user's input are:
#1)terms -> vector list which contains the query terms
#2)years -> vector list which contains the years for which to select the results.
#3)Additional Insert section -> create all the vectors and add them to the "user_new_column" dataframe
#The vectors' length have to be the same as the length of the year's vector (e.i. 13 for 2008:2019 search).
#outputDir -> path to the directory where to save the geom line plot called "dual_RNAseq_query_res.pdf". 
#Use "." for the current directory. Whereas, do not add the last "/" for the output directory

#The plot can be found in the object called "query_plot".

#Load Library
library(rentrez)    #Interface to the NCBI's 'EUtils' API, allowing search in databases like "GenBank", and "PubMed".
library(reshape)    #Flexibly restructure and aggregate data using just two functions: melt and 'dcast' (or 'acast')
library(tidyverse)  #Collection of R packages designed for data science: ggplot2, dplyr, tidyr, readr, purrr, tibble, stringr, forcats.

rm(list=ls(all=TRUE)) #Remove all object saved in the system

###### USER's INPUT ######
terms <- c ("\"Dual\" [All Fields] AND \"RNA-seq\" [All Fields]", "\"Parallel\" [All Fields] AND \"RNA-seq\" [All Fields]", "\"Simultaneous\" [All Fields] AND \"RNA-seq\" [All Fields]")
years <- 2008:2019
outputDir <- "."   #Use "." for the current directory. Whereas, do not add the last "/" for the output directory

#Additional Insert
Manual <- c(0,0,1,2,3,10,14,15,21,15,14,50)
Reviews <- c(0,0,0,2,1,0,0,5,5,7,2,8)
user_new_column <- data.frame(Manual, Reviews)
#########################

######PUBMED search
resultDF <- data.frame(matrix(years))

# Function to find the number of papers per year
papers_by_year <- function(year, search_term)
{return(sapply(year, function(y) entrez_search(db="pubmed",term=search_term, mindate=y, maxdate=y, retmax=0)$count))}

for (i in 1:length(terms))
{
  trend_data <- sapply(terms[i], function(t) papers_by_year(years, t))
  resultDF <- cbind(resultDF, trend_data)
}

resultDF <- cbind(resultDF, user_new_column)

#PLOT
melted <- melt(resultDF, id.vars="matrix.years.")
Y_breaks <- seq(0, 185, by=5)

query_plot <- ggplot(data=melted, aes(x=matrix.years., y=value, group=variable, colour=variable)) + 
  geom_line(size=1) + 
  theme(axis.text.x = element_text(angle = 0, vjust = 1)) +
  scale_x_continuous("Years", labels = as.character(resultDF$matrix.years.), breaks = resultDF$matrix.years.) +
  scale_y_continuous("N* Papers", labels = as.character(Y_breaks), breaks = Y_breaks)

pdf(file = sprintf("%s/dual_RNAseq_query_res.pdf", outputDir))
query_plot
dev.off()
