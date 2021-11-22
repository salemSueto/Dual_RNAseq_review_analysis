#Comparison of dual RNA-seq studies of host-pathogen interaction

The following markdown displays the steps used to analysed a typical dual RNA-seq samples. 

##Table of Contents

1.	[Dual RNA-seq studies discovery](#section1.)
2.	[Dual RNA-seq study selection](#section2.)
	1.	[Reference](#section2.1.)
3.	[Preparation](#section3.)
	1. [Project directory creation](#section3.1.)
	2. [Tool](#section3.2.)
	3. [FastQ download](#section3.3.)
		1. [Aspera connect app download](#section3.3.1.)
		2. [Choose the dataset](#section3.3.2.)
		3. [Download the data](#section3.3.3.)
	4. [Trimmomatic Adapters](#section3.4.)
	5. [STAR: genome index generation](#section3.5.)
	6. [RSEM: genome index generation](#section3.6.)
	7. [Salmon: transcriptome index generation](#section3.7.)
	8. [Kallisto: transcriptome index generation](#section3.8.)
	9. [Scotty](#section3.9.)
4.	[R scripts](#section4.)
	1.	[R package](#section4.1.)
	2.	[Dual RNA-seq Pubmed search](#section4.2.)
	3.	[Gene-Transcript relation and length](#section4.3.)
	4.	[Polyester: dual RNA-seq simulation](#section4.4.)
	5.	[Merge HTSeq-count tables](#section4.5.)
	6.	[Merge RSEM tables](#section4.6.)
	7.	[Merge Salmon tables](#section4.7.)
	8.	[Merge Kallisto tables](#section4.8.)
	9.	[TPM Spearman correlation](#section4.9.)
	10.	[Scotty heatmap](#section4.10.)
	11.	[GMT custom annotation for un-supported organisms](#section4.11.)
		1.	[Taxon ID search](#section4.11.1.)
		2. 	[Custom annotation creation](#section4.11.2.)
	12.	[Polyester's simulated data differential expressed analysis](#section4.12.)
	13.	[Count tables differential expressed analysis](#section4.13.)
	14.	[Count tables differential expressed and enrichmment analysis comparison](#section4.14.) 
	15.	[R session Info](#section4.15.)
5.	[Polyester dual RNA-seq simulated data](#section5.)
	1.	[RNA-seq simulation data creation](#section5.1.)
	2.	[Merge fasta file](#section5.2.)
		1.	[Single-end sample](#section5.2.1.)
		2.	[Paired-end sample](#section5.2.2.)
	3.	[Merged data information analysis](#section5.1.)
6.	[Dual RNA-seq pipeline](#section6.)
	1.	[Single-end reads: Mapping & Quantification](#section6.1.)
	2.	[Paired-end reads: Mapping & Quantification](#section6.2.)
	3.	[Merge count tables](#section6.3.)
		1.	[HTSeq-count count table merge](#section6.3.1.)
		2.	[RSEM count table merge](#section6.3.2.)
		3.	[Salmon count table merge](#section6.3.3.)
		4.	[Kallisto count table merge](#section6.3.4.)
	4.	[Move count table](#section6.4.)
	5.	[TPM Spearman correlation](#section6.5.)
	6.	[Scotty analysis](#section6.6.)
	7.	[Custom annotation creation](#section6.7.)
	8.	[Polyester merged count table analysis](#section6.8.)
	9.	[Count tables DE and enrichment comparison analysis](#section6.9.)
		
		
---
[![](https://mermaid.ink/img/eyJjb2RlIjoiZ3JhcGggVERcbiAgQVtEdWFsIFJOQS1zZXEgc3R1ZGllc10gLS0-QltTdHVkeSBzZWxlY3Rpb25dXG4gICAgQiAtLT5DW1ByZXBhcmF0aW9uXVxuICAgIEMgLS0-RFtSIHNjcmlwdHNdXG4gICAgQyAtLT5FW1JhdyBEYXRhXVxuICAgIEUgLS0-RltNYXBwaW5nXVxuICAgIEYgLS0-R1tRdWFudGlmaWNhdGlvbl1cbiAgICBHIC0tPkhbQW5hbHlzaXNdXG4gICAgSCAtLT5JW1RQTSBTcGVhcm1hbiBjb3JyZWxhdGlvbl1cbiAgICBIIC0tPkpbU2NvdHR5IGFuYWx5c2lzXVxuICAgIEggLS0-S1tEaWZmZXJlbnRpYWwgZXhwcmVzc2VkIGFuYWx5c2lzXVxuICAgIEsgLS0-TFtFbnJpY2htZW50IGFuYWx5c2lzXVxuICAgIEMgLS0-REREW0V0Yy4uLl1cblxuICAgICAgICAgICAgIiwibWVybWFpZCI6eyJ0aGVtZSI6ImRlZmF1bHQiLCJ0aGVtZVZhcmlhYmxlcyI6eyJiYWNrZ3JvdW5kIjoid2hpdGUiLCJwcmltYXJ5Q29sb3IiOiIjRUNFQ0ZGIiwic2Vjb25kYXJ5Q29sb3IiOiIjZmZmZmRlIiwidGVydGlhcnlDb2xvciI6ImhzbCg4MCwgMTAwJSwgOTYuMjc0NTA5ODAzOSUpIiwicHJpbWFyeUJvcmRlckNvbG9yIjoiaHNsKDI0MCwgNjAlLCA4Ni4yNzQ1MDk4MDM5JSkiLCJzZWNvbmRhcnlCb3JkZXJDb2xvciI6ImhzbCg2MCwgNjAlLCA4My41Mjk0MTE3NjQ3JSkiLCJ0ZXJ0aWFyeUJvcmRlckNvbG9yIjoiaHNsKDgwLCA2MCUsIDg2LjI3NDUwOTgwMzklKSIsInByaW1hcnlUZXh0Q29sb3IiOiIjMTMxMzAwIiwic2Vjb25kYXJ5VGV4dENvbG9yIjoiIzAwMDAyMSIsInRlcnRpYXJ5VGV4dENvbG9yIjoicmdiKDkuNTAwMDAwMDAwMSwgOS41MDAwMDAwMDAxLCA5LjUwMDAwMDAwMDEpIiwibGluZUNvbG9yIjoiIzMzMzMzMyIsInRleHRDb2xvciI6IiMzMzMiLCJtYWluQmtnIjoiI0VDRUNGRiIsInNlY29uZEJrZyI6IiNmZmZmZGUiLCJib3JkZXIxIjoiIzkzNzBEQiIsImJvcmRlcjIiOiIjYWFhYTMzIiwiYXJyb3doZWFkQ29sb3IiOiIjMzMzMzMzIiwiZm9udEZhbWlseSI6IlwidHJlYnVjaGV0IG1zXCIsIHZlcmRhbmEsIGFyaWFsIiwiZm9udFNpemUiOiIxNnB4IiwibGFiZWxCYWNrZ3JvdW5kIjoiI2U4ZThlOCIsIm5vZGVCa2ciOiIjRUNFQ0ZGIiwibm9kZUJvcmRlciI6IiM5MzcwREIiLCJjbHVzdGVyQmtnIjoiI2ZmZmZkZSIsImNsdXN0ZXJCb3JkZXIiOiIjYWFhYTMzIiwiZGVmYXVsdExpbmtDb2xvciI6IiMzMzMzMzMiLCJ0aXRsZUNvbG9yIjoiIzMzMyIsImVkZ2VMYWJlbEJhY2tncm91bmQiOiIjZThlOGU4IiwiYWN0b3JCb3JkZXIiOiJoc2woMjU5LjYyNjE2ODIyNDMsIDU5Ljc3NjUzNjMxMjglLCA4Ny45MDE5NjA3ODQzJSkiLCJhY3RvckJrZyI6IiNFQ0VDRkYiLCJhY3RvclRleHRDb2xvciI6ImJsYWNrIiwiYWN0b3JMaW5lQ29sb3IiOiJncmV5Iiwic2lnbmFsQ29sb3IiOiIjMzMzIiwic2lnbmFsVGV4dENvbG9yIjoiIzMzMyIsImxhYmVsQm94QmtnQ29sb3IiOiIjRUNFQ0ZGIiwibGFiZWxCb3hCb3JkZXJDb2xvciI6ImhzbCgyNTkuNjI2MTY4MjI0MywgNTkuNzc2NTM2MzEyOCUsIDg3LjkwMTk2MDc4NDMlKSIsImxhYmVsVGV4dENvbG9yIjoiYmxhY2siLCJsb29wVGV4dENvbG9yIjoiYmxhY2siLCJub3RlQm9yZGVyQ29sb3IiOiIjYWFhYTMzIiwibm90ZUJrZ0NvbG9yIjoiI2ZmZjVhZCIsIm5vdGVUZXh0Q29sb3IiOiJibGFjayIsImFjdGl2YXRpb25Cb3JkZXJDb2xvciI6IiM2NjYiLCJhY3RpdmF0aW9uQmtnQ29sb3IiOiIjZjRmNGY0Iiwic2VxdWVuY2VOdW1iZXJDb2xvciI6IndoaXRlIiwic2VjdGlvbkJrZ0NvbG9yIjoicmdiYSgxMDIsIDEwMiwgMjU1LCAwLjQ5KSIsImFsdFNlY3Rpb25Ca2dDb2xvciI6IndoaXRlIiwic2VjdGlvbkJrZ0NvbG9yMiI6IiNmZmY0MDAiLCJ0YXNrQm9yZGVyQ29sb3IiOiIjNTM0ZmJjIiwidGFza0JrZ0NvbG9yIjoiIzhhOTBkZCIsInRhc2tUZXh0TGlnaHRDb2xvciI6IndoaXRlIiwidGFza1RleHRDb2xvciI6IndoaXRlIiwidGFza1RleHREYXJrQ29sb3IiOiJibGFjayIsInRhc2tUZXh0T3V0c2lkZUNvbG9yIjoiYmxhY2siLCJ0YXNrVGV4dENsaWNrYWJsZUNvbG9yIjoiIzAwMzE2MyIsImFjdGl2ZVRhc2tCb3JkZXJDb2xvciI6IiM1MzRmYmMiLCJhY3RpdmVUYXNrQmtnQ29sb3IiOiIjYmZjN2ZmIiwiZ3JpZENvbG9yIjoibGlnaHRncmV5IiwiZG9uZVRhc2tCa2dDb2xvciI6ImxpZ2h0Z3JleSIsImRvbmVUYXNrQm9yZGVyQ29sb3IiOiJncmV5IiwiY3JpdEJvcmRlckNvbG9yIjoiI2ZmODg4OCIsImNyaXRCa2dDb2xvciI6InJlZCIsInRvZGF5TGluZUNvbG9yIjoicmVkIiwibGFiZWxDb2xvciI6ImJsYWNrIiwiZXJyb3JCa2dDb2xvciI6IiM1NTIyMjIiLCJlcnJvclRleHRDb2xvciI6IiM1NTIyMjIiLCJjbGFzc1RleHQiOiIjMTMxMzAwIiwiZmlsbFR5cGUwIjoiI0VDRUNGRiIsImZpbGxUeXBlMSI6IiNmZmZmZGUiLCJmaWxsVHlwZTIiOiJoc2woMzA0LCAxMDAlLCA5Ni4yNzQ1MDk4MDM5JSkiLCJmaWxsVHlwZTMiOiJoc2woMTI0LCAxMDAlLCA5My41Mjk0MTE3NjQ3JSkiLCJmaWxsVHlwZTQiOiJoc2woMTc2LCAxMDAlLCA5Ni4yNzQ1MDk4MDM5JSkiLCJmaWxsVHlwZTUiOiJoc2woLTQsIDEwMCUsIDkzLjUyOTQxMTc2NDclKSIsImZpbGxUeXBlNiI6ImhzbCg4LCAxMDAlLCA5Ni4yNzQ1MDk4MDM5JSkiLCJmaWxsVHlwZTciOiJoc2woMTg4LCAxMDAlLCA5My41Mjk0MTE3NjQ3JSkifX0sInVwZGF0ZUVkaXRvciI6ZmFsc2V9)](https://mermaid-js.github.io/mermaid-live-editor/#/edit/eyJjb2RlIjoiZ3JhcGggVERcbiAgQVtEdWFsIFJOQS1zZXEgc3R1ZGllc10gLS0-QltTdHVkeSBzZWxlY3Rpb25dXG4gICAgQiAtLT5DW1ByZXBhcmF0aW9uXVxuICAgIEMgLS0-RFtSIHNjcmlwdHNdXG4gICAgQyAtLT5FW1JhdyBEYXRhXVxuICAgIEUgLS0-RltNYXBwaW5nXVxuICAgIEYgLS0-R1tRdWFudGlmaWNhdGlvbl1cbiAgICBHIC0tPkhbQW5hbHlzaXNdXG4gICAgSCAtLT5JW1RQTSBTcGVhcm1hbiBjb3JyZWxhdGlvbl1cbiAgICBIIC0tPkpbU2NvdHR5IGFuYWx5c2lzXVxuICAgIEggLS0-S1tEaWZmZXJlbnRpYWwgZXhwcmVzc2VkIGFuYWx5c2lzXVxuICAgIEsgLS0-TFtFbnJpY2htZW50IGFuYWx5c2lzXVxuICAgIEMgLS0-REREW0V0Yy4uLl1cblxuICAgICAgICAgICAgIiwibWVybWFpZCI6eyJ0aGVtZSI6ImRlZmF1bHQiLCJ0aGVtZVZhcmlhYmxlcyI6eyJiYWNrZ3JvdW5kIjoid2hpdGUiLCJwcmltYXJ5Q29sb3IiOiIjRUNFQ0ZGIiwic2Vjb25kYXJ5Q29sb3IiOiIjZmZmZmRlIiwidGVydGlhcnlDb2xvciI6ImhzbCg4MCwgMTAwJSwgOTYuMjc0NTA5ODAzOSUpIiwicHJpbWFyeUJvcmRlckNvbG9yIjoiaHNsKDI0MCwgNjAlLCA4Ni4yNzQ1MDk4MDM5JSkiLCJzZWNvbmRhcnlCb3JkZXJDb2xvciI6ImhzbCg2MCwgNjAlLCA4My41Mjk0MTE3NjQ3JSkiLCJ0ZXJ0aWFyeUJvcmRlckNvbG9yIjoiaHNsKDgwLCA2MCUsIDg2LjI3NDUwOTgwMzklKSIsInByaW1hcnlUZXh0Q29sb3IiOiIjMTMxMzAwIiwic2Vjb25kYXJ5VGV4dENvbG9yIjoiIzAwMDAyMSIsInRlcnRpYXJ5VGV4dENvbG9yIjoicmdiKDkuNTAwMDAwMDAwMSwgOS41MDAwMDAwMDAxLCA5LjUwMDAwMDAwMDEpIiwibGluZUNvbG9yIjoiIzMzMzMzMyIsInRleHRDb2xvciI6IiMzMzMiLCJtYWluQmtnIjoiI0VDRUNGRiIsInNlY29uZEJrZyI6IiNmZmZmZGUiLCJib3JkZXIxIjoiIzkzNzBEQiIsImJvcmRlcjIiOiIjYWFhYTMzIiwiYXJyb3doZWFkQ29sb3IiOiIjMzMzMzMzIiwiZm9udEZhbWlseSI6IlwidHJlYnVjaGV0IG1zXCIsIHZlcmRhbmEsIGFyaWFsIiwiZm9udFNpemUiOiIxNnB4IiwibGFiZWxCYWNrZ3JvdW5kIjoiI2U4ZThlOCIsIm5vZGVCa2ciOiIjRUNFQ0ZGIiwibm9kZUJvcmRlciI6IiM5MzcwREIiLCJjbHVzdGVyQmtnIjoiI2ZmZmZkZSIsImNsdXN0ZXJCb3JkZXIiOiIjYWFhYTMzIiwiZGVmYXVsdExpbmtDb2xvciI6IiMzMzMzMzMiLCJ0aXRsZUNvbG9yIjoiIzMzMyIsImVkZ2VMYWJlbEJhY2tncm91bmQiOiIjZThlOGU4IiwiYWN0b3JCb3JkZXIiOiJoc2woMjU5LjYyNjE2ODIyNDMsIDU5Ljc3NjUzNjMxMjglLCA4Ny45MDE5NjA3ODQzJSkiLCJhY3RvckJrZyI6IiNFQ0VDRkYiLCJhY3RvclRleHRDb2xvciI6ImJsYWNrIiwiYWN0b3JMaW5lQ29sb3IiOiJncmV5Iiwic2lnbmFsQ29sb3IiOiIjMzMzIiwic2lnbmFsVGV4dENvbG9yIjoiIzMzMyIsImxhYmVsQm94QmtnQ29sb3IiOiIjRUNFQ0ZGIiwibGFiZWxCb3hCb3JkZXJDb2xvciI6ImhzbCgyNTkuNjI2MTY4MjI0MywgNTkuNzc2NTM2MzEyOCUsIDg3LjkwMTk2MDc4NDMlKSIsImxhYmVsVGV4dENvbG9yIjoiYmxhY2siLCJsb29wVGV4dENvbG9yIjoiYmxhY2siLCJub3RlQm9yZGVyQ29sb3IiOiIjYWFhYTMzIiwibm90ZUJrZ0NvbG9yIjoiI2ZmZjVhZCIsIm5vdGVUZXh0Q29sb3IiOiJibGFjayIsImFjdGl2YXRpb25Cb3JkZXJDb2xvciI6IiM2NjYiLCJhY3RpdmF0aW9uQmtnQ29sb3IiOiIjZjRmNGY0Iiwic2VxdWVuY2VOdW1iZXJDb2xvciI6IndoaXRlIiwic2VjdGlvbkJrZ0NvbG9yIjoicmdiYSgxMDIsIDEwMiwgMjU1LCAwLjQ5KSIsImFsdFNlY3Rpb25Ca2dDb2xvciI6IndoaXRlIiwic2VjdGlvbkJrZ0NvbG9yMiI6IiNmZmY0MDAiLCJ0YXNrQm9yZGVyQ29sb3IiOiIjNTM0ZmJjIiwidGFza0JrZ0NvbG9yIjoiIzhhOTBkZCIsInRhc2tUZXh0TGlnaHRDb2xvciI6IndoaXRlIiwidGFza1RleHRDb2xvciI6IndoaXRlIiwidGFza1RleHREYXJrQ29sb3IiOiJibGFjayIsInRhc2tUZXh0T3V0c2lkZUNvbG9yIjoiYmxhY2siLCJ0YXNrVGV4dENsaWNrYWJsZUNvbG9yIjoiIzAwMzE2MyIsImFjdGl2ZVRhc2tCb3JkZXJDb2xvciI6IiM1MzRmYmMiLCJhY3RpdmVUYXNrQmtnQ29sb3IiOiIjYmZjN2ZmIiwiZ3JpZENvbG9yIjoibGlnaHRncmV5IiwiZG9uZVRhc2tCa2dDb2xvciI6ImxpZ2h0Z3JleSIsImRvbmVUYXNrQm9yZGVyQ29sb3IiOiJncmV5IiwiY3JpdEJvcmRlckNvbG9yIjoiI2ZmODg4OCIsImNyaXRCa2dDb2xvciI6InJlZCIsInRvZGF5TGluZUNvbG9yIjoicmVkIiwibGFiZWxDb2xvciI6ImJsYWNrIiwiZXJyb3JCa2dDb2xvciI6IiM1NTIyMjIiLCJlcnJvclRleHRDb2xvciI6IiM1NTIyMjIiLCJjbGFzc1RleHQiOiIjMTMxMzAwIiwiZmlsbFR5cGUwIjoiI0VDRUNGRiIsImZpbGxUeXBlMSI6IiNmZmZmZGUiLCJmaWxsVHlwZTIiOiJoc2woMzA0LCAxMDAlLCA5Ni4yNzQ1MDk4MDM5JSkiLCJmaWxsVHlwZTMiOiJoc2woMTI0LCAxMDAlLCA5My41Mjk0MTE3NjQ3JSkiLCJmaWxsVHlwZTQiOiJoc2woMTc2LCAxMDAlLCA5Ni4yNzQ1MDk4MDM5JSkiLCJmaWxsVHlwZTUiOiJoc2woLTQsIDEwMCUsIDkzLjUyOTQxMTc2NDclKSIsImZpbGxUeXBlNiI6ImhzbCg4LCAxMDAlLCA5Ni4yNzQ1MDk4MDM5JSkiLCJmaWxsVHlwZTciOiJoc2woMTg4LCAxMDAlLCA5My41Mjk0MTE3NjQ3JSkifX0sInVwZGF0ZUVkaXRvciI6ZmFsc2V9)
---		
		

<h2 id="section1.">1.	Dual RNA-seq studies discovery</h2>

Dual RNA-seq is a term used to identify studies where the RNA-seq platform is used to study the transcriptome of two interacting species without prior physical separation but a silico one at the transcript level. Nonetheless, several terms have been used for these studies such as "parallel" and "simultaneous". These three terms were used to query the Pubmed database to identify the potential dual RNA-seq studies. The results were scrutinized to filter for actual papers and reviews which use the dual RNA-seq strategy.  

The "dual\_search\_pubmed.R" script, saved in the "dual\_RNAseq/Rscript" directory, was used to find the number of queries for the three searched terms. And, the actual dual RNA-seq (Manual) and reviews were plotted together.

The user's inputs are displayed:

```
terms <- c ("\"Dual\" [All Fields] AND \"RNA-seq\" [All Fields]", "\"Parallel\" [All Fields] AND \"RNA-seq\" [All Fields]", "\"Simultaneous\" [All Fields] AND \"RNA-seq\" [All Fields]")
years <- 2008:2019
outputDir <- ""

#Additional Insert
Manual <- c(0,0,1,2,3,10,14,15,21,15,14,50)
Reviews <- c(0,0,0,2,1,0,0,5,5,7,2,8)
user_new_column <- data.frame(Manual, Reviews)
```

The length of inputs such as "Manual" and "Reviews" has to be the same as the length of the "years" object.


<h2 id="section2.">2.	Dual RNA-seq study selection</h2>

The manual selected dual RNA-seq papers were filtered for there parameters:

1. Human or mouse cells used as host.
2. Bacteria as pathogen. 
3. The gene count table is available for both the host and the pathogen.

| Paper                                                                                                         | Host         | Pathogen     | Host:Human/Mouse | Pathogen:Bacteria | Raw Count Table | 
|---------------------------------------------------------------------------------------------------------------|--------------|--------------|------------------|-------------------|-----------------| 
| [Kawahara et al. 2012](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3490861/)                                 | Rice         |  Fungus      | No               | No                | No              | 
| [Choi et al. 2014](http://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0002905)                 | Mosquito     |  Roundworm   | No               | No                | No              | 
| [Lowe et al. 2014](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0103098)                  | Rapeseed     |  Fungus      | No               | No                | No              | 
| [Asai et al. 2014](http://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1004443)            | Thale cress  |  Fungus      | No               | No                | No              | 
| [Hayden et al. 2014](https://link.springer.com/article/10.1007/s11295-014-0698-0)                             | Tanoak       |  Fungus      | No               | No                | No              | 
| [Rosani et al. 2015](https://onlinelibrary.wiley.com/doi/full/10.1111/1462-2920.12706)                        | Oyster       |  Virus       | No               | No                | No              | 
| [Field et al. 2015](http://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1005168)           | Bat          |  Fungus      | No               | No                | No              | 
| [Handa et al. 2015](https://academic.oup.com/pcp/article/56/8/1490/1845470#84725704)                          | Legume       |  Fungus      | No               | No                | No              | 
| [Meyer et al. 2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4773608/)                                    | Eucalypts    |  Fungus      | No               | No                | No              | 
| [Petitot et al. 2016](https://www.ncbi.nlm.nih.gov/m/pubmed/26610268/)                                        | Rice         |  Worm        | No               | No                | No              | 
| [Zhang et al. 2016](https://www.sciencedirect.com/science/article/pii/S1931524415004491?via%3Dihub)           | Zebrafish    |  Human       | No               | No                | No              | 
| [Zhou et al. 2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4865485/)                                     | Pig          |  Protozoa    | No               | No                | No              | 
| [Musungu et al. 2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5116468/)                                  | Corn         |  Fungus      | No               | No                | No              | 
| [Meyer et al. 2016](https://www.frontiersin.org/articles/10.3389/fpls.2016.00191/full)                        | Eucalypts    |  Fungus      | No               | No                | No              | 
| [Buddenborg et al. 2017](http://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0005984)           | Snail        |  Flatworm    | No               | No                | No              | 
| [Burns et al. 2017](https://elifesciences.org/articles/22054)                                                 | Salamander   |  Green Algae | No               | No                | No              | 
| [Daly et al. 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5476202/)                                     | Fungus       |  Fungus      | No               | No                | No              | 
| [Dong et al. 2017](https://www.ncbi.nlm.nih.gov/pubmed/28730600)                                              | Moth         |  Fungus      | No               | No                | No              | 
| [He et al. 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5462986/)                                       | Bacteria     |  Fungus      | No               | No                | No              | 
| [Hennessy et al. 2017](https://bmcresnotes.biomedcentral.com/articles/10.1186/s13104-017-2704-8)              | Fungus       |  Fungus      | No               | No                | No              | 
| [Huang et al. 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5329020/)                                    | Rice         |  Fungus      | No               | No                | No              | 
| [Ling et al. 2017](https://www.nature.com/articles/s41598-017-03563-6)                                        | Melon        |  Roundworm   | No               | No                | No              | 
| [Oeser et al. 2017](https://europepmc.org/abstract/pmc/pmc5379732)                                            | Rye          |  Fungus      | No               | No                | No              | 
| [Reeder et al. 2017](https://www.ncbi.nlm.nih.gov/m/pubmed/28614673/)                                         | Bat          |  Fungus      | No               | No                | No              | 
| [Tomada et al. 2017](http://onlinelibrary.wiley.com/doi/10.1111/1462-2920.13861/full)                         | Bacteria     |  Fungus      | No               | No                | No              | 
| [Vainio et al. 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5809736/)                                   | Fungus       |  Virus       | No               | No                | No              | 
| [Moolhuijzen et al. 2018](https://bmcresnotes.biomedcentral.com/articles/10.1186/s13104-018-3993-2)           | Wheat        |  Fungus      | No               | No                | No              | 
| [Field et al. 2018](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.14827)                               | Bat          |  Fungus      | No               | No                | No              | 
| [Varaldi et al. 2018](https://pubmed.ncbi.nlm.nih.gov/29806574/)                                              | Wasp         |  Virius      | No               | No                | No              | 
| [Wang et al. 2018](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6221353/)                                     | Wheat        |  Fungus      | No               | No                | No              | 
| [Li et al. 2019](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-019-5902-z)                    | Banana       |  Fungus      | No               | No                | No              | 
| [Becker et al. 2019](https://apsjournals.apsnet.org/doi/full/10.1094/MPMI-01-19-0028-R)                       | Rapeseed     |  Fungus      | No               | No                | No              | 
| [Bai et al. 2019](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6520846/)                                      | Snail        |  Virus       | No               | No                | No              | 
| [Kovalchuk et al. 2019](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6318961/)                                | Spruce       |  Fungus      | No               | No                | No              | 
| [Mateus et al. 2019](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6474227/)                                   | Cassava      |  Fungus      | No               | No                | No              | 
| [Rosenthal et al. 2011](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3146290/)                                | Bacteria     |  Bacteria    | No               | Yes               | No              | 
| [Yazawa et al. 2013](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3640049/)                                   | Sorghum      |  Bacteria    | No               | Yes               | No              | 
| [Vannucci et al. 2013](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-14-421)               | Pig          |  Bacteria    | No               | Yes               | No              | 
| [Camilios-Neto et al. 2014](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-15-378)          | Wheat        |  Bacteria    | No               | Yes               | No              | 
| [Teixeira et al. 2014](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4277218/)                                 | Cocoa Tree   |  Bacteria    | No               | Yes               | No              | 
| [Roux et al. 2014](https://onlinelibrary.wiley.com/doi/epdf/10.1111/tpj.12442)                                | Legume       |  Bacteria    | No               | Yes               | No              | 
| [Zuluaga et al. 2015](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-015-1460-1)               | Potato       |  Bacteria    | No               | Yes               | No              | 
| [Dutton et al. 2016](https://www.ncbi.nlm.nih.gov/m/pubmed/26042999/)                                         | Fungus       |  Bacteria    | No               | Yes               | No              | 
| [Paungfoo-Lonhienne et al. 2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5116747/)                       | Sugarcane    |  Bacteria    | No               | Yes               | No              | 
| [Bing et al. 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5489720/)                                     | Tse Tse fly  |  Bacteria    | No               | Yes               | No              | 
| [Grote et al. 2017](http://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0005357)                | Roundworm    |  Bacteria    | No               | Yes               | No              | 
| [Hendrickson et al. 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5329018/)                              | Bacteria     |  Bacteria    | No               | Yes               | No              | 
| [Munoz et al. 2017](https://academic.oup.com/gbe/article/9/9/2276/4101901)                                    | Tsetse fly   |  Bacteria    | No               | Yes               | No              | 
| [Mutha et al. 2018](https://onlinelibrary.wiley.com/doi/full/10.1111/omi.12248)                               | Bacteria     |  Bacteria    | No               | Yes               | No              | 
| [Kumar et al. 2018](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6262203/)                                    | Zebrafish    |  Bacteria    | No               | Yes               | No              | 
| [Valenzuela-Miranda et al. 2018](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6277808/)                       | Salmon       |  Bacteria    | No               | Yes               | No              | 
| [Miller et al. 2018](https://onlinelibrary.wiley.com/doi/full/10.1111/omi.12238)                              | Bacteria     |  Bacteria    | No               | Yes               | No              | 
| [Miller et al. 2018](https://onlinelibrary.wiley.com/doi/full/10.1111/omi.12238)                              | Bacteria     |  Bacteria    | No               | Yes               | No              | 
| [Liao et al. 2019](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6469767/)                                     | Rice         |  Bacteria    | No               | Yes               | No              | 
| [Huang et al. 2019](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6509204/)                                    | Fish         |  Bacteria    | No               | Yes               | No              | 
| [Mutha et al. 2019](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6529473/)                                    | Bacteria     |  Bacteria    | No               | Yes               | No              | 
| [Rubio et al. 2019](https://www.pnas.org/content/early/2019/06/19/1905747116.long)                            | Mollusca     |  Bacteria    | No               | Yes               | No              | 
| [Sun et al. 2019](https://pubmed.ncbi.nlm.nih.gov/30708058/)                                                  | Fish         |  Bacteria    | No               | Yes               | No              | 
| [Tang et al. 2019](https://pubmed.ncbi.nlm.nih.gov/31129188/)                                                 | Fish         |  Bacteria    | No               | Yes               | No              | 
| [Zhang et al. 2019](https://pubmed.ncbi.nlm.nih.gov/30776540/)                                                | Fish         |  Bacteria    | No               | Yes               | No              | 
| [Perazzoli et al. 2016](https://www.ncbi.nlm.nih.gov/pubmed/27784266)                                         | Bacteria(s)  | Bacteria(s)  | No               | Yes*              | No              | 
| [Tierney et al. 2012](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3299011/)                                  | Mouse        |  Fungus      | Yes              | No                | No              | 
| [Lisnic et al. 2013](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3784481/)                                   | Mouse        |  Virus       | Yes              | No                | No              | 
| [Pittman et al. 2014](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-15-806)                | Mouse        |  Protozoa    | Yes              | No                | No              | 
| [Yamagishi et al. 2014](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4158759/)                                | Human        |  Protozoa    | Yes              | No                | No              | 
| [Abate et al. 2015](http://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1005158)           | Human        |  Virus       | Yes              | No                | No              | 
| [Dillon et al. 2015](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-015-2237-2)                | Mouse        |  Leishmania  | Yes              | No                | No              | 
| [Park et al. 2015](https://www.ncbi.nlm.nih.gov/m/pubmed/26576844/)                                           | Mouse        |  Virus       | Yes              | No                | No              | 
| [Bruno et al. 2015](http://mbio.asm.org/content/6/2/e00182-15.full.pdf+html)                                  | Mouse        |  Fungus      | Yes              | No                | No              | 
| [Liu et al. 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4417116/)                                      | Human/Mouse  |  Fungus      | Yes              | No                | No              | 
| [Wilk et al. 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4557482/)                                     | Mouse        |  Virus       | Yes              | No                | No              | 
| [Wilk et al. 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4557482/pdf/12864_2015_Article_1867.pdf)      | Mouse        |  Virus       | Yes              | No                | No              | 
| [Fernandes et al. 2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4959658/)                                | Human        |  Leishmania  | Yes              | No                | No              | 
| [Ehret et al. 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5584376/)                                    | Mouse        |  Protozoa    | Yes              | No                | No              | 
| [Wesolowska-Andersen et al. 2017](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1140-8) | Human        |  Virus       | Yes              | No                | No              | 
| [Niemiec et al. 2018](https://pure.qub.ac.uk/portal/files/134673291/s12864_017_4097_4.pdf)                    | Human        |  Fungus      | Yes              | No                | No              | 
| [Fabozzi et al. 2018](https://jvi.asm.org/content/92/17/e00518-18.long)                                       | Human        |  Virus       | Yes              | No                | No              | 
| [Petrucelli et al.2018](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6070946/)                                | Human        |  Fungus      | Yes              | No                | No              | 
| [Lee et al. 2018](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6326353/)                                      | Human        |  Malaria     | Yes              | No                | No              | 
| [LaMonte et al. 2019](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6353872/)                                  | Human        |  Protozoa    | Yes              | No                | No              | 
| [Shadap et al. 2019](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6370631/)                                   | Mouse        |  Leishmania  | Yes              | No                | No              | 
| [Humphrys et al. 2013](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0080597)              | Human        |  Bacteria    | Yes              | Yes               | No              | 
| [Mavromatis et al. 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4950338/)                               | Mouse        |  Bacteria    | Yes              | Yes               | No              | 
| [Avraham et al. 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4578813/)                                  | Mouse        |  Bacteria    | Yes              | Yes               | No              | 
| [Rienksma et al. 2015](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-014-1197-2)              | Human        |  Bacteria    | Yes              | Yes               | No              | 
| [Baddal et al. 2016](https://www.ncbi.nlm.nih.gov/pubmed/26578681)                                            | Human        |  Bacteria    | Yes              | Yes               | No              | 
| [Damron et al. 2016](https://www.nature.com/articles/srep39172)                                               | Mouse        |  Bacteria    | Yes              | Yes               | No              | 
| [Westermann et al. 2016](https://www.nature.com/articles/nature16547)                                         | Human        |  Bacteria    | Yes              | Yes               | No              | 
| [Thanert et al. 2017](https://www.ncbi.nlm.nih.gov/pubmed/28155859)                                           | Mouse        |  Bacteria    | Yes              | Yes               | No              | 
| [Li et al. 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5596872/)                                       | Mouse        |  Bacteria    | Yes              | Yes               | No              | 
| [Aprianto et al. 2016](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1054-5)            | Human        |  Bacteria    | Yes              | Yes               | Yes             | 
| [Nuss et al. 2017](http://www.pnas.org/content/114/5/E791)                                                    | Mouse        |  Bacteria    | Yes              | Yes               | Yes             | 
| [Shenoy et al. 2017](http://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1006582)          | Mouse        |  Bacteria    | Yes              | Yes               | Yes             | 
| [Zimmermann et al. 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5566787/)                               | Human        |  Bacteria    | Yes              | Yes               | Yes             | 
| [Montoya et al. 2019](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6508871/)                                  | Human        |  Bacteria    | Yes              | Yes               | Yes*            | 
| [Griesenauer et al. 2019](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6581864/)                              | Human        |  Bacteria    | Yes              | Yes               | Yes*            | 
| [Perez-Losada et al. 2015](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0131819)          | Human        |  Bacteria(s) | Yes              | Yes*              | No              | 



<h3 id="section2.1.">2.1.	Reference</h3>

The following table shows the dual RNA-seq papers which passed the three filtering filters and the species used in such studies. 


| Study                                                                                                | RNAseq Data                                                         | BioProject - SRA | Organism                             | Role     | Genome     | Annotation    | Transcriptome | GenBank         | Download                                                                                                       | Path                                  | 
|------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------|------------------|--------------------------------------|----------|------------|---------------|---------------|-----------------|----------------------------------------------------------------------------------------------------------------|---------------------------------------| 
| [Aprianto et al. 2016](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1054-5)   | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE79595)  | SRP072326        | Homo sapiens                         | Host     | GRCh38     | GRCh38.97     | GRCh38        | GRCh38          | [Ensembl](https://www.ensembl.org/Homo_sapiens/Info/Index)                                                     | reference/hsapiens                    | 
| [Aprianto et al. 2016](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1054-5)   | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE79595)  | SRP072326        | Streptococcus pneumoniae D39         | Pathogen | ASM1436v1  | ASM1436v1.37  | ASM1436v1     | GCA_000014365.2 | [Ensembl Bacteria](https://bacteria.ensembl.org/Streptococcus_pneumoniae_d39/Info/Index)                       | reference/spneumoniae_D39             | 
| [Zimmermann et al. 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5566787/)                      | [ENA](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-5287/)  | ERP020415        | Homo sapiens                         | Host     | GRCh38     | GRCh38.97     | GRCh38        | GRCh38          | [Ensembl](https://www.ensembl.org/Homo_sapiens/Info/Index)                                                     | reference/hsapiens                    | 
| [Zimmermann et al. 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5566787/)                      | [ENA](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-5287/)  | ERP020415        | Mycobacterium tuberculosis H37Rv     | Pathogen | ASM19595v2 | ASM19595v2.37 | ASM19595v2    | GCA_000195955.2 | [Ensembl Bacteria](https://bacteria.ensembl.org/Mycobacterium_tuberculosis_h37rv/Info/Index)                   | reference/mtuberculosis_H37Rv         | 
| [Nuss et al. 2017](http://www.pnas.org/content/114/5/E791)                                           | [ENA](https://www.ebi.ac.uk/ena/data/search?query=PRJEB14242)       | PRJEB14242       | Mus musculus                         | Host     | GRCm38     | GRCm38.97     | GRCm38        | GRCm38          | [Ensembl](https://www.ensembl.org/Mus_musculus/Info/Index)                                                     | reference/mmusculus                   | 
| [Nuss et al. 2017](http://www.pnas.org/content/114/5/E791)                                           | [ENA](https://www.ebi.ac.uk/ena/data/search?query=PRJEB14242)       | PRJEB14242       | Yersinia pseudotuberculosis IP 32953 | Pathogen | ASM4736v1  | ASM4736v1.37  | ASM4736v1     | NC_006153       | [Ensembl Bacteria](https://bacteria.ensembl.org/Yersinia_pseudotuberculosis_ip_32953_gca_000834295/Info/Index) | reference/ypseudotuberculosis\_IP_32953 | 
| [Shenoy et al. 2017](http://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1006582) | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE86118)  | PRJNA340293      | Mus musculus                         | Host     | GRCm38     | GRCm38.97     | GRCm38        | GRCm38          | [Ensembl](https://www.ensembl.org/Mus_musculus/Info/Index)                                                     | reference/mmusculus                   | 
| [Shenoy et al. 2017](http://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1006582) | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE86118)  | PRJNA340293      | Streptococcus pneumoniae TIGR4       | Pathogen | ASM688v1   | ASM688v1.37   | ASM688v1      | NZ_AKBW01000001 | [Ensembl Bacteria](https://bacteria.ensembl.org/Streptococcus_pneumoniae_tigr4/Info/Index)                     | reference/spneumoniae_TIGR4           | 
| [Montoya et al. 2019](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6508871/)                         | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125943) | PRJNA518047      | Homo sapiens                         | Host     | GRCh38     | GRCh38.97     | GRCh38        | GRCh38          | [Ensembl](https://www.ensembl.org/Homo_sapiens/Info/Index)                                                     | reference/hsapiens                    | 
| [Montoya et al. 2019](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6508871/)                         | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125943) | PRJNA518047      | Mycobacterium leprae TN              | Pathogen | ASM19585v1 | ASM19585v1.48 | ASM19585v1    | NC_002677       | [Ensembl Bacteria](https://bacteria.ensembl.org/Mycobacterium_leprae_tn/Info/Index)                            | reference/mleprae_TN                  | 
| [Griesenauer et al. 2019](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6581864/)                     | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130901) | PRJNA541925      | Homo sapiens                         | Host     | GRCh38     | GRCh38.97     | GRCh38        | GRCh38          | [Ensembl](https://www.ensembl.org/Homo_sapiens/Info/Index)                                                     | reference/hsapiens                    | 
| [Griesenauer et al. 2019](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6581864/)                     | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130901) | PRJNA541925      | Haemophilus ducreyi 35000HP          | Pathogen | ASM794v1   | ASM794v1.48   | ASM794v1      | NC_002940       | [Ensembl Bacteria](https://bacteria.ensembl.org/Haemophilus_ducreyi_35000hp/Info/Index)                        | reference/hducreyi_35000HP            | 


<h2 id="section3.">3.	Preparation</h2>

The present section shows all the preparation that it is needed to conduct the analysis showcased in the methodology part of the document.


<h3 id="section3.1.">3.1.	Project directory creation </h3>

The following directories are created to help organise the dual RNA-seq study. The "dual\_RNAseq" is the main directory where all the information are saved for an efficient analysis such as references, tools, and R scripts.

```bash
$	mkdir dual_RNAseq
$	cd dual_RNAseq
$	mkdir tool
$	mkdir tool/adapters
$	mkdir reference
$	mkdir metadata
$	mkdir custom_annotation
$	mkdir Rscript
```

* The tool directory has all the downloaded version of the tools in case a package manager is not used. Futhermore, the sub-directory "adapters" has the sequence fasta files of the adapters used during the sequencing step.

* The referene directory has all the reference files for the selected species such as genome fasta file, genome annotation, and the species' transcriptomics file. Each species has its own sub-directory (e.i. hsapiens). Furthermore, it also has its own sub-directories to store star, rsem, salmon, kallisto indexes.

```
#Example for the Homo sapiens' reference directory
$	mkdir reference/hsapiens
$	mkdir reference/hsapiens/star
$	mkdir reference/hsapiens/rsem
$	mkdir reference/hsapiens/salmon
$	mkdir reference/hsapiens/kallisto
```

* The metadata directory stores information such gene and transcript relation, gene and transcript length, etc..

* The custom\_annotation directory stores the custom annotation of species which are not supported by the g:Profiler tool used for enrichment analysis.

* The Rscript directory stores all the R scripts used in this methodology.

When a new dual RNA-seq study is started, a new directory will be created (FirstAuthor_Year) which will be used to save all data related to that study.


<h3 id="section3.2.">3.2.	Tool </h3>

The Conda package manager was used to install all the command-line tools, except the BLAST tool. Conda installation can be found [here](https://bioconda.github.io/user/install.html). In case of conflict between the user's system and Conda, the command-line tools can be manually downloaded and installed. R packages were installed and uploaded directly from R, check [3.10.1. Install R packages](#section3.10.1.).

| Tool          | Version   | Type                | Date Accessed   | Manual                                                                                                          | Installation                            | Source                                                                  | 
|---------------|-----------|---------------------|-----------------|-----------------------------------------------------------------------------------------------------------------|-----------------------------------------|-------------------------------------------------------------------------| 
|  IBM Aspera Connect     | 3.10.1       |  App/Command-Line       |  NA             |  [Documentation](https://www.ibm.com/support/knowledgecenter/SSXMX3_3.10/connect_user_osx/guide.html)                                                              |  [Download](https://www.ibm.com/aspera/connect/)     |  NA                             | 
|  Samtools     | 1.9       |  Command-Line       |  NA             |  [Manual](http://www.htslib.org/doc/samtools.html)                                                              |  conda install -c bioconda samtools     |  [Website](http://www.htslib.org/download/)                             | 
|  FastQC       |  0.11.8   |  Command-Line       |  NA             |  [Manual](https://dnacore.missouri.edu/PDF/FastQC_Manual.pdf)                                                   |  conda install -c bioconda fastqc       |  [Website](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)  | 
|  MultiQC      | 1.2       |  Command-Line       |  NA             |  [Tutorial](https://multiqc.info/docs/)                                                                         |  conda install -c bioconda multiqc      |  [Github](https://github.com/ewels/MultiQC)                             | 
|  Trimmomatic  | 0.39      |  Command-Line       |  NA             |  [Manual](http://www.usadellab.org/cms/?page=trimmomatic)                                                       |  conda install -c bioconda trimmomatic  |  [Website](http://www.usadellab.org/cms/?page=trimmomatic)              | 
|  STAR         |  2.7.1a   |  Command-Line       |  NA             |  [Manual](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf)          |  conda install -c bioconda star         |   [Github](https://github.com/alexdobin/STAR)                           | 
|  HTSeq-count  | 0.6       |  Command-Line       |  NA             |  [Manual](https://htseq.readthedocs.io/en/release_0.11.1/count.html)                                            |  conda install -c bioconda htseq        |  [Website](https://pypi.org/project/HTSeq/)                             | 
|  RSEM         |  1.3.3    |  Command-Line       |  NA             |  [Manual](https://deweylab.github.io/RSEM/README.html)                                                          |  conda install -c bioconda rsem         |  [Github](https://github.com/deweylab/RSEM)                             | 
|  Salmon       | 0.14.1    |  Command-Line       |  NA             |  [Tutoral](https://combine-lab.github.io/salmon/getting_started/)                                               | conda install -c bioconda salmon        | [Github](https://github.com/COMBINE-lab/salmon)                         | 
| Kallisto      | 0.46.1    |  Command-Line       | NA              | [Manual](https://pachterlab.github.io/kallisto/manual)                                                          | conda install -c bioconda kallisto      | [Patcher Lab](https://pachterlab.github.io/kallisto/download.html)      | 
|  Scotty       |  NA       |  Matlab             |  Apr-20         |  [Help](http://scotty.genetics.utah.edu/help.html)                                                              |  NA                                     |  [Github](https://github.com/mbusby/Scotty)                             | 
|  DESeq2       |  1.24.0   |  R                  |  NA             |  [Vignette](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)               |  NA                                     |  NA                                                                     | 
|  edgeR        |  3.30.3   |  R                  |  NA             |  [Manual](https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf)      |  NA                                     |  NA                                                                     | 
|  NOISeq       |  2.31.0   |  R                  |  NA             |  [Bioconductor](https://www.bioconductor.org/packages/release/bioc/vignettes/NOISeq/inst/doc/NOISeq.pdf)        |  NA                                     |  NA                                                                     | 
|  g:Profiler2  |  0.2.0    |  R                  |  NA             |  [Vignette](https://cran.r-project.org/web/packages/gprofiler2/vignettes/gprofiler2.html)                       |  NA                                     |  NA                                                                     | 
|  REVIGO       |  NA       |  Website            |  May-20         |  [Tool Website](http://revigo.irb.hr/index.jsp?error=expired)                                                   |  NA                                     |  NA                                                                     | 


<h3 id="section3.3.">3.3.	FastQ download </h3>

The Aspera client was used to download the fastq files by following the [Biostars thread] (https://www.biostars.org/p/325010/). The raw data is saved in the directory called "dual\_RNAseq/FirstAuthor\_Year/0\_raw\_data".

####3.3.1.	Aspera connect app download

* [Download] (https://www.ibm.com/aspera/connect/) the IBM Aspera Connect app. And, install the program. The following executable files will be in their default locations:

```
#Linux
$HOME/.aspera/connect/bin/ascp #the executable
$HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh #Need later

Mac
$HOME/Applications/AsperaConnect.app/Contents/Resources/ascp #the executable
$HOME/Applications/AsperaConnect.app/Contents/Resources/asperaweb_id_dsa.openssh #Need later
```

* Add the folder with the ascp executable to your PATH.

####3.3.2.	Choose the dataset
* Access the [SRA-Explorer] (https://sra-explorer.info/) website. 
* Enter the accession number of the desired study and confirm. 
* Select the desired samples and add them to the cart.
* Click on the saved datasets button on the top right corner of the screen.
* Click on "FastQ Downloads" and then on "Aspera commands for downloading FastQ files" tabs.
* Select the right operation system (Linux/OSX) and click on copy/download.
* Download the samples metadata from the "Full Metadata" tab.

####3.3.3.	Download the data

Use the following script to download the fastq files.

```bash
## Either by a simple loop:
while read LIST; do
$LIST; done < download.txt

## or by using GNU parallel to have things parallelized:
cat download.txt | parallel "{}"
```

<h3 id="section3.4.">3.4.	Trimmomatic Adapters </h3>

[Download](https://github.com/timflutre/trimmomatic/tree/master/adapters) the adapter sequence files and save then in the dual\_RNAseq/tool/adapters directory.


<h3 id="section3.5.">3.5.	STAR: genome reference index generation </h3>

The following scripts can be used for the two type of organisms (Eukaryote and Bacteria) for the creation of STAR genome references indexes. The "sjdbOverhang" is an integer equal to ReadLength-1; whereas, the "genomeSAindexNbases" flag is calculated as min(14, log2(GenomeLength)/2 - 1).

```bash
#Host: Eukaryote Multicellular Organism -> (Homo sapiens/Mus musculus)
$	star --runMode genomeGenerate --genomeDir rootDirectory/../reference/species/star --genomeFastaFiles *dna.primary_assembly.fa --sjdbGTFfile *gtf --sjdbOverhang 74 -p 8

#Pathogen: Bacteria unicellular Organism (e.i. Streptococcus pneumoniae)
$	star --runMode genomeGenerate --genomeDir rootDirectory/../reference/species/star --genomeFastaFiles *dna.toplevel.fa --sjdbGTFfile *gtf --genomeSAindexNbases 9 --sjdbOverhang 74
```

<h3 id="section3.6.">3.6.	RSEM: genome reference index generation </h3>

The following scripts can be used for the two type of organisms (Eukaryote and Bacteria) for the creation of RSEM genome references indexes.

```bash
#Host & Pathogen (Eukaryote & Bacteria)
$	rsem-prepare-reference --gtf *gtf *dna.primary_assembly.fa rootDirectory/../reference/species/rsem/index
```

<h3 id="section3.7.">3.7.	Salmon: transcriptome reference index generation </h3>

The following scripts can be used for the two type of organisms (Eukaryote and Bacteria) for the creation of Salmon genome references indexes. The value of k at 31 work well for reads of 75bp or longer, but a shoter value of k may improve sensitivity even more when using selective alignment

```bash
#Host & Pathogen (Eukaryote & Bacteria)
$	salmon index -t *cdna.all.fa -i reference/species/salmon/ -k 27
```

<h3 id="section3.8.">3.8.	Kallisto: transcriptome reference index generation </h3>

The following scripts can be used for the two type of organisms (Eukaryote and Bacteria) for the creation of Kallisto genome references indexes.  The value of k at 31 work well for reads of 75bp or longer, but a shoter value of k may improve sensitivity even more when using selective alignment (default: 31, max value: 31).

```bash
#Host & Pathogen (Eukaryote & Bacteria)
$	kallisto index * cdna.all.fa -i reference/species/kallisto/index -k 27
```

<h3 id="section3.9.">3.9.	Scotty </h3>

The Scotty files were written using an older version of Matlab.  At some point Matlab depricated one of the functions used which will cause running errors on newer versions (Matlab 2013+). It is imperative to change the function "RandStream.setDefaultStream(s)" into "RandStream.setGlobalStream(s)" in the following files: getProbabilitiesFromFit (line 23) and getProbSequencedMonteCarloByReads (line 7).

Furthermore, in order to save the tables used to create the Power Optimization and Poisson Noise plots, the following lines were added to the file "getOptimizationCharts.m" after line 80 which showcase the line: "readsSX=readsS/10^6;".

```
csvwrite('powerOptimization.csv', pwrsCalc)
csvwrite('poissonNoise.csv', pwrBiases)
csvwrite('seqDepth.csv', readsSX)
csvwrite('nReps.csv', nReps)
```



<h2 id="section4.">4.	R scripts </h2>

The following R scripts are used during the pipeline and they are saved in the "dual\_RNAseq/Rscripts" directory.


<h3 id="section4.1.">4.1.	R package </h3>


| Program     | Version  | Purpose                                                                                 | Installation                       | Citation                                                                                                                                                                                                                            | 
|-------------|----------|-----------------------------------------------------------------------------------------|------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------| 
| Biostrings  | 2.56.0   | Manipulation of biological strings.                                                     | BiocManager::install("Biostrings") | H. Pagès, P. Aboyoun, R. Gentleman and S. DebRoy (2020). Biostrings: Efficient manipulation of biological strings. R package version 2.56.0.                                                                                        | 
| DESeq2      | 1.28.1   | Differential expression                                                                 | BiocManager::install("DESeq2")     | Love, M.I., Huber, W., Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2 Genome Biology 15(12):550 (2014)                                                                                  | 
| edgeR       | 3.30.3   | Differential expression                                                                 | BiocManager::install("edgeR")      | Robinson MD, McCarthy DJ and Smyth GK (2010). edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics 26, 139-140                                                        | 
| NOISeq      | 2.31.0   | Differential expression                                                                 | BiocManager::install("NOISeq")     | Tarazona, S., Furio-Tari, P., Turra, D., Di Pietro, A., Nueda, M.J., Ferrer, A., & Conesa, A. (2015). Data quality aware analysis of differential expression in RNA-seq with NOISeq R/Bioc package. Nucleic Acids Research.         | 
| polyester   | 1.24.0   | RNA-seq read simulation                                                                 | BiocManager::install("polyester")  | Alyssa C. Frazee, Andrew E. Jaffe, Rory Kirchner and Jeffrey T. Leek (2020). polyester: Simulate RNA-seq reads. R package version 1.24.0.                                                                                           | 
| tximport    | 1.16.1   | Import and summarize transcript-level estimates for transcript- and gene-level analysis | BiocManager::install("tximport")   | Charlotte Soneson, Michael I. Love, Mark D. Robinson (2015): Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences. F1000Research                                                             | 
| UniProt.ws  | 2.28.0   | R Interface to UniProt Web Services                                                     | BiocManager::install("UniProt.ws") | Marc Carlson and Bioconductor Package Maintainer (2020). UniProt.ws: R Interface to UniProt Web Services. R package version 2.28.0.                                                                                                 | 
| factoextra  | 1.0.7    | Multivariate Data Analyses and Elegant Visualization                                    | install.packages("factoextra")     | Alboukadel Kassambara and Fabian Mundt (2020). factoextra: Extract and Visualize the Results of Multivariate Data Analyses. R package version 1.0.7. https://CRAN.R-project.org/package=factoextra                                  | 
| forestmangr | 0.9.2    | Process forest inventory data                                                           | install.packages("forestmangr")    | Sollano Rabelo Braga, Marcio Leles Romarco de Oliveira and Eric Bastos Gorgens (2020). forestmangr: Forest Mensuration and Management. R package version 0.9.2. https://CRAN.R-project.org/package=forestmangr                      | 
| ggforce     | 0.3.2    | Aid in visual data investigations.                                                      | install.packages("ggforce")        | Thomas Lin Pedersen (2020). ggforce: Accelerating 'ggplot2'. R package version 0.3.2. https://CRAN.R-project.org/package=ggforce                                                                                                    | 
| gprofiler2  | 2_0.2.0  | Gene list functional enrichment analysis                                                | install.packages("gprofiler2")     | Liis Kolberg and Uku Raudvere (2020). gprofiler2: Interface to the 'g:Profiler' Toolset. R package version 0.2.0. https://CRAN.R-project.org/package=gprofiler2                                                                     | 
| gridExtra   | 2.3      | Arrange multiple grid-based plots on a page                                             | install.packages("gridExtra")      | Baptiste Auguie (2017). gridExtra: Miscellaneous Functions for "Grid" Graphics. R package version 2.3. https://CRAN.R-project.org/package=gridExtra                                                                                 | 
| rentrez     | 1.2.3    | Interface to the NCBI's 'EUtils' API                                                    | install.packages("rentrez")        | Winter, D. J. (2017) rentrez: an R package for the NCBI eUtils API The R Journal 9(2):520-526                                                                                                                                       | 
| reshape     | 0.8.8    | Flexibly restructure and aggregate data                                                 | install.packages("reshape")        | H. Wickham. Reshaping data with the reshape package. Journal of Statistical Software, 21(12), 2007.                                                                                                                                 | 
| rvest       | 0.3.6    | Download, then manipulate, HTML and XML                                                 | install.packages("rvest")          | Hadley Wickham (2020). rvest: Easily Harvest (Scrape) Web Pages. R package version 0.3.6. https://CRAN.R-project.org/package=rvest                                                                                                  | 
| scales      | 1.1.1    | Graphical scales map data to aesthetics                                                 | install.packages("scales")         | Hadley Wickham and Dana Seidel (2020). scales: Scale Functions for Visualization. R package version 1.1.1. https://CRAN.R-project.org/package=scales                                                                                | 
| stringi     | 1.5.3    | Language processing tools                                                               | install.packages("stringi")        | Gagolewski M. and others (2020). R package stringi: Character string processing facilities. http://www.gagolewski.com/software/stringi/.                                                                                            | 
| tidytext    | 0.2.6    | Text mining tasks and plot generation.                                                  | install.packages("tidytext")       | Silge J, Robinson D (2016). "tidytext: Text Mining and Analysis Using Tidy Data Principles in R." _JOSS_, *1*(3). doi: 10.21105/joss.00037 (URL:https://doi.org/10.21105/joss.00037), <URL: http://dx.doi.org/10.21105/joss.00037>. | 
| tidyverse   | 1.3.0    | Collection of R packages designed for data science                                      | install.packages("tidyverse")      | Wickham et al., (2019). Welcome to the tidyverse. Journal of Open Source Software, 4(43), 1686, https://doi.org/10.21105/joss.01686                                                                                                 | 
| ggVennDiagram   | 0.3    | A 'ggplot2' Implement of Venn Diagram                                      | install.packages("ggVennDiagram")      | Chun-Hui Gao (2019). ggVennDiagram: A 'ggplot2' Implement of Venn Diagram. R package version 0.3. https://CRAN.R-project.org/package=ggVennDiagram                                                                                                 | 
| rlang   | 0.4.9    | Functions for Base Types and Core R and 'Tidyverse' Features                                      | install.packages("rlang")      | Lionel Henry and Hadley Wickham (2020). rlang: Functions for Base Types and Core R and 'Tidyverse' Features. R package version 0.4.9. https://CRAN.R-project.org/package=rlang                                                                                                 | 
| ggpubr   | 0.4.0    | 'ggplot2' Based Publication Ready Plots                                   | install.packages("ggpubr")      | Alboukadel Kassambara (2020). ggpubr: 'ggplot2' Based Publication Ready Plots. R package version 0.4.0. https://CRAN.R-project.org/package=ggpubr                                                                                                 | 



<h3 id="section4.2.">4.2.	Dual RNA-seq Pubmed search </h3>

The "dual\_search\_pubmed.R" script queries the Pubmed databbase with a list of search terms for a specific range of years. The user can add additional columns to the dataframe resulted from the search. 

The user's input are:

* terms -> vector list which contains the query terms
* years -> vector list which contains the years for which to select the results.
* Additional Insert section -> create new vectors and add them to the "user\_new\_column" dataframe. The vectors' length have to be the same as the length of the year's vector (e.i. 12 for 2008:2019 search). 
* outputDir -> path to the directory where to save the geom line plot called "dual\_RNAseq\_query\_res.pdf".

Please check the RStudio plot in case the saved plot has some issues or bad aestetics. The plot can be found in the object called "query\_plot".


```R
#!/dual_search_pubmed.R

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
```

<h3 id="section4.3.">4.3.	Gene-Transcript relation and length </h3>

The gene\_transcript\_info.R script takes as input the path to directories where the transcriptome (.cdna.all.fa) and the annotation (*.gtf) files of the species are saved. Then, it outputs two text files which displays gene-transcript relationship and the gene/transcript length; whereas, one pdf file which shows the relation between the two input types. The script creates a directory called "txome_annotation" inside the directory the inputs are saved.

The input are:

* dir_list -> vector list of path to the directory where the transcriptome and the annotation files are saved.

The output are:

* tx_gene_relation.txt -> find the gene and its transcript isoforms (TranscriptID - GeneID - GTF - CDNA). The GTF and CDNA columns shows whether the transcript is found in the file.
* tx_gene_length.txt -> find the length  for each elements. Transcript length by checking the length from the transcriptome file. Whereas, gene length is the sum of all the gene's transcript isoforms length.
* cdna_gtf_plots.pdf -> pdf file which shows the transcriptome and annotation venn groups.

```R
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
```

<h3 id="section4.4.">4.4.	Polyester: dual RNA-seq simulation </h3>

The polyester\_sim.R script is used to create simulated RNA-seq fasta data. The simulation is performed on two or more organims. Then, the script merges all the simulated reads, sample by sample, to simulate a typical dual RNA-seq study.

The input are:

1. txome\_list -> list of path to the the transcriptome fasta files.
2. gene\_selection -> path to the file with the list of genes used to filter the transcriptome. The file contains one column for each transcriptome with header having the same name of the transcriptome. Each column has the list of genes which is used for the filtering step.
3. polyester\_settings -> path the file which has all the settings for Polyester to simulate data. The rownames has the same name of the transcriptome fasta files. The file has the following columns:
	* DE\_TXs (Number of differentially expressed transcripts)
	* Coverage (average coverage for all transcript)
	* Min\_FC (minimum value of fold change for differentially expressed transcripts)
	* Max\_FC (maximum value of fold change for differentially expressed transcripts) 
	* Read\_Length (length of the simulated reads)
	* Paired\_End (TRUE for paired-end; FALSE for single-end)
	* Strand\_Specific (TRUE for strand specific, FALSE for unstranded)
	* Replicate\_Group1 (number of replicates for group1)
	* Replicate\_Group2 (number of replicates for group2)
4. outputDir -> name of the directory where to save all the output files of the script. 

The script has different outputs and saved inside the output-directory (outputCir):

1. One output directory for each transcriptome file where:
	* List of simulated compressed fasta file equal to the sum of replicates for both groups.
	* id\_length.txt -> shows the gene/transcript ID and its length. The transcript's length is obtained from the transcriptome fasta file. And the gene's length is the sum of all the gene's isoforms length.
	* sim\_rep\_info.txt -> shows the sample's group membership.
	* sim\_tx\_info.txt -> shows the fold change and its DE status for each transcript simulated.
	* sim\_counts\_matrix.rda -> shows the correct count table for all the simulated transcripts in rda format.
	* sim\_isoform\_count.txt -> shows the correct transcript count for the simulated samples.
	* sim\_gene\_count.txt -> shows the correct gene count for the simulated samples.

2. The Polyester merged simulation (polyester\_merged\_sim) directory has the information from each simulated data is merged together.
	* id\_length.txt -> merges all the gene/transcript's length (1B) files together.
	* sim\_count\_gene.txt -> merges the gene count for each simulated transcriptome.
	* sim\_count\_isoform.txt -> merges the transcript count for each simulated transcriptome.
	* sim\_count\_percentage.txt -> shows the percentage of each transcriptome based on the total number of reads.
	* sim\_count\_summary.txt -> shows the total number of reads for each simulated transcriptome.
	* sim\_de\_gene.txt -> shows the differential expressed genes for each simulated transcriptome.
	* sim\_de\_tx.txt -> shows the differential expressed transcripts for each simulated transcriptome.
	* sim\_rep\_info.txt -> shows the sample's group membership for each simulated transcriptome.
	* sim\_tx\_info.txt -> shows the fold change and its DE status for each simulated transcriptome.

```R
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

#The script has different outputs and saved inside the output-directory (outputCir):

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
```

<h3 id="section4.5.">4.5.	Merge HTSeq-count tables </h3>

The htseq\_merge\_table.R script merges the gene count tables obtained from STAR in-built HTSeq-count, during the mapping step with the flag "--quantMode GeneCounts" in each directory. The outputs are count table for each of the three sequencing types.

The input are:

1. pathFiles -> list of directories where the files reside
2. patternFiles -> one pattern used to select the right files in the directory

The output files are saved in the same directory where the single count table files reside:

1. Count tables for each strand type (unstranded, strand, and reverse strand).
2. summary\_count.txt -> text file which has the summary info for all three strand types:
number of unmapped, number of multimapping, number of no features, number of ambiguous, and number of mapped reads.
3. One pdf file of the summary plots information.

```R
#!/htseq_merge_table.R

#!/htseq_merge_table.R
#The script takes the gene counts obtained from STAR in-built HTSeq-count, during the mapping step with the flag "--quantMode GeneCounts", 
#and merges them together to form one file for each of the three sequencing type.

#The inputs are:
#1)pathFiles -> list of directories where the files reside
#2)patternFiles -> one pattern used to select the right files in the directory

#The output are:
#1)Count tables for each strand type (unstranded, strand, and reverse strand).
#2)summary\_count.txt -> text file which has the summary info for all three strand types:
#number of unmapped, number of multimapping, number of no features, number of ambiguous, and number of mapped reads.
#3)One pdf file of the summary plots information.

#Load Library
library (tidyverse) #Collection of R packages designed for data science: ggplot2, dplyr, tidyr, readr, purrr, tibble, stringr, forcats.
library(gridExtra) #Provides the ability to arrange multiple grid-based plots on a page, and draw tables.

rm(list = ls()) #Remove all saved objects

####### USER's Input ######
pathFiles <- c("Hsapiens_Spneumoniae_SingleStrand_SE/3_star_spneumoniae/htseq") #Get list of all files in directory
patternFiles <- "sample" #Pattern used to select the right files
###########################

#Loop through each directory
for (iDir in pathFiles) 
{
  #Get all files which has the patternFiles character
  files <- list.files(path=iDir, pattern=patternFiles)
  list_of_files <- list()
  
  if (length(files) > 0)
  {
    # Save all files
    for (i in files) 
    {list_of_files[[i]] <- read.table(paste(iDir, "/", i, sep = ""), header=FALSE)}
    
    #Summary Info Dataframe
    summary_plot_df <- as.data.frame (matrix(nrow = 0, ncol = 4))
    colnames(summary_plot_df) <- c("ID", "variable", "value", "MapType")
    
    #Extract count for each sequencing type
    selectColID <- c("unstranded", "strand_yes", "strand_rev")
    
    for (selectCol in c(2,3,4)) #2: unstranded - #3: 1st read strand aligned with RNA (htseq-count -s yes) - #4: 1st read strand aligned with RNA (htseq-count -s reverse)
    {
      # Select the column for each files
      merge_htseq_df <- as.data.frame (matrix(data = NA, nrow = nrow(list_of_files[[1]])-4, ncol = length(list_of_files)))
      summary_htseq_df <- as.data.frame (matrix(data = NA, nrow = 4, ncol = length(list_of_files)))
      
      colnames(merge_htseq_df) <- names(list_of_files)
      colnames(summary_htseq_df) <- names(list_of_files)
      
      rownames(merge_htseq_df) <- list_of_files[[1]]$V1 [-(1:4)]
      rownames(summary_htseq_df) <- list_of_files[[1]]$V1 [1:4]
      
      for (i in 1:length(list_of_files)) 
      {
        sFile <- list_of_files[[i]] %>% dplyr::select (contains(as.character(selectCol)))
        summary_htseq_df[,i] <- sFile[1:4, 1]
        merge_htseq_df [,i] <- tail(sFile, -4)
      }
      
      merge_htseq_colsums <- colSums(merge_htseq_df)
      merge_htseq_df <- merge_htseq_df %>% rownames_to_column(var = "ID")
      colnames(merge_htseq_df) <- sapply (str_split(colnames(merge_htseq_df), "\\."), `[`, 1)
      
      summary_htseq_df <- rbind (summary_htseq_df, N_mapped_Quantified = merge_htseq_colsums) %>% rownames_to_column(var = "ID")
      summary_htseq_df$MapType <- selectColID[selectCol-1]
      colnames(summary_htseq_df) <- sapply (str_split(colnames(summary_htseq_df), "\\."), `[`, 1)
      
      temp_summary <- reshape2::melt(summary_htseq_df,id.vars=c("ID"))
      temp_summary$MapType <- selectColID[selectCol-1]
      summary_plot_df <- rbind(summary_plot_df, temp_summary)
      
      fileName <- paste("htseq-count", selectColID[selectCol-1], sep = "_")
      write.table(merge_htseq_df, paste(iDir, "/", fileName, "_count.txt", sep = ""), row.names=F, sep="\t", quote=F)
    }
    
    ####PLOT SUMMARY INFO #####
    summary_plot_df <- summary_plot_df %>% filter(value!=MapType)
    write.table(summary_plot_df, paste(iDir, "/summary_count.txt", sep = ""), row.names=F, sep="\t", quote=F)
    
    plots_list <- list()
    for (i in unique(summary_plot_df$ID))
    {
      plots_list[[paste("Variable:", i)]] <- summary_plot_df %>%
        filter(ID==i) %>%
        ggplot(aes(x=variable, y=value, fill=MapType)) + 
        geom_bar(position = "dodge", stat="identity") +
        labs(title=paste("Variable:", i)) + 
        theme(axis.text.x = element_text(angle=65, hjust=1.0))
    }
    
    pdf(paste(iDir, "/htseq_summary_plots.pdf", sep = ""), onefile=TRUE)
    for (i in seq(length(plots_list))) {grid.arrange (plots_list[[i]])}
    dev.off()
  }
}
```

<h3 id="section4.6.">4.6.	Merge RSEM tables </h3>

The rsem\_merge\_table.R script takes RSEM outputs from the transcriptome alignment of STAR and merges them together to form one file at the gene and transcript level. The user's input is pathFiles which is the list of directories where the RSEM output files reside.

The outputs are two txt files for each directory:

1. rsem\_gene\_count.txt -> merges all the gene count.
2. rsem\_isoform\_count.txt -> merges all the transcripts count.

```R
#!/rsem_merge_table.R

#The script takes RSEM outputs from the transcriptome alignment of STAR and merges them together to form one file at the gene and transcript level.

#The user's input is pathFiles which is the list of directories where the RSEM output files reside.

#The outputs are two txt files for each directory:
#1)"rsem_gene_count.txt" -> merges all the gene count
#1)"rsem_isoform_count.txt" -> merges all the transcripts count

#Load Library
library (tidyverse) #Collection of R packages designed for data science: ggplot2, dplyr, tidyr, readr, purrr, tibble, stringr, forcats.

#Remove all saved objects
rm(list = ls())

################
#USER's INPUT
pathFiles <- c("Hsapiens_Spneumoniae_SingleStrand_SE/2_star_hsapiens/rsem") #Get list of all files in directory
################

#Get the samples
for (iDir in pathFiles) #Loop through the directories
{
  for (iPat in c("gene", "isoform")) #Loop through the patterns
  {
    # Save all files
    files <- list.files(path=iDir, pattern=iPat)
    list_of_files <- list()
    
    if (length(files)>0)
    {
      # Save all files
      for (i in files) 
      {list_of_files[[i]] <- read.table(paste (iDir, "/", i, sep = ""), header=TRUE)}
      
      # Select the "expected_count" column for each files
      merge_rsem_df <- as.data.frame (matrix(data = NA, nrow = nrow(list_of_files[[1]]), ncol = length(list_of_files)+4))
      merge_rsem_df[,c(1,2,3,4)] <- list_of_files[[1]][,1:4]
      for (i in 1:length(list_of_files)) {merge_rsem_df [,i+4] <- list_of_files[[i]]$expected_count}
      merge_rsem_df <- merge_rsem_df %>% dplyr::select(-c(2,3,4))
      colnames(merge_rsem_df) <- c("ID", sapply(str_split(names(list_of_files), "\\."), `[`, 1))
      write.table (merge_rsem_df, file = sprintf("%s/rsem_%s_count.txt", iDir, iPat), sep="\t", row.names=F, quote=F)
    }
  }
}
```

<h3 id="section4.7.">4.7.	Merge Salmon tables </h3>

The salmon\_merge\_table.R script merges the Salmon reference-free alignments outputs into one file for both gene and transcript level. 

The input are:

1. fileDirectory -> list of directories where the Salmon's outputs are saved.
2. gene\_tx -> path to a file which has atleast two columns called "GeneID" and "TranscriptID".

The output files are saved in the same directory of the input files:

1. salmon\_gene\_count.txt" -> merges all the gene count
2. salmon\_isoform\_count.txt" -> merges all the transcripts count

```R
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
```

<h3 id="section4.8.">4.8.	Merge Kallisto tables </h3>

The kallisto\_merge\_table.R script merges the Kallisto reference-free alignments outputs into one file for both gene and transcript level. 

The input are:

1. fileDirectory -> list of directories where the Kallisto's outputs are saved.
2. gene\_tx -> path to a file which has atleast two columns called "GeneID" and "TranscriptID".

The two output files are saved in the same directory of the input files:

1. kallisto\_gene\_count.txt" -> merges all the gene count
2. kallisto\_isoform\_count.txt" -> merges all the transcripts count

```R

#!/kallisto_merge_table.R

#The script merges the Kallisto reference-free alignments at the gene and the transcript level.

#The script has 2 inputs:
#1)fileDirectory -> list of directories where the Kallisto's outputs are saved.
#2) gene_tx -> path to a file which has atleast two columns called "GeneID" and "TranscriptID".

#The script outputs two txt files for each directory:
#1)"kallisto\_gene\_count.txt" -> merges all the gene count
#2)"kallisto\_isoform\_count.txt" -> merges all the transcripts count

#Load Library
library(tidyverse)  #Collection of R packages designed for data science: ggplot2, dplyr, tidyr, readr, purrr, tibble, stringr, forcats.
library(tximport) #Import and summarize transcript-level estimates for transcript- and gene-level analysis

rm(list=ls(all=TRUE))

##### ##### ##### User's Input ##### ##### ##### 
fileDirectory <- c("Hsapiens_Spneumoniae_SingleStrand_SE/5_kallisto")  #Directory where the Kallisto results are saved
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
  #Get all the abundance.tsv sample names 
  files <- list.dirs(iDir, recursive=FALSE)
  files <- paste0(files, "/abundance.tsv")
  dirNames <- str_split(list.dirs(iDir, recursive=FALSE), pattern = "/")
  dirNames <- unlist (lapply(dirNames, function(x) x[length(x)]))
  names(files) <- dirNames
  
  #Remove the transcript version from the transcript ID column from the abundance.tsv files
  for (iTxt in files) 
  {
    iTxT_df <- read.csv2(iTxt, sep="\t", header = TRUE)
    iTxT_df$target_id <- sapply (str_split(iTxT_df$target_id, "\\."), `[`, 1)
    write.table(iTxT_df, file = unname(iTxt), row.names=F, sep="\t", quote=F)
  }
  
  #Get Gene Count Table
  txi_kallisto <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)
  count <- as.data.frame(txi_kallisto$counts) %>% rownames_to_column("ID")
  write.table(count, file = sprintf("%s/kallisto_gene.txt", iDir), row.names=F, sep="\t", quote=F)
  
  ##### Get Transcript Count Table
  txi_kallisto <- tximport(files, type = "kallisto", txOut = TRUE)
  count <- as.data.frame(txi_kallisto$counts) %>% rownames_to_column("ID")
  count$ID <- str_replace (count$ID, "\\..*", "")
  write.table(count, file = sprintf("%s/kallisto_isoform.txt", iDir), row.names=F, sep="\t", quote=F)
}
```

<h3 id="section4.9.">4.9.	TPM Spearman correlation</h3>

The tpm\_spearman\_cor.R script takes a series of count tables from a directory. First, it calculates for each gene/tx its TPM value. Then, it calculate the Spearman correlation among each count table duo combination. Finally, it plots a heatmap of the Spearman correlation.


The inputs are:

1. pathFiles -> list of directories where the count tables reside. The count tables, from the same directory, need to have the same headers and in the same order. And, the first column need to have the gene/transcript ID.
2. patternFiles -> select the right count tables, in the directory, based on the filename pattern.
3. gene\_tx\_Length -> file path to the file where the gene/transcript and its length are saved (#ID - #Length).
4. mean\_read\_seq -> mean of the fragment length sequenced.

Two output files are created for each directory-pattern combination and they are saved in the same directory as the input files:

1. pattern\_tpm\_spearman\_df.txt -> text file that shows the TPM Spearman correlation 
2. pattern\_tpm\_spearman\_plot.pdf -> heatmap plot of the TPM correlation.

```R
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
pathFiles <- c("Hsapiens_Spneumoniae_SingleStrand_SE/6_count_table/hsapiens")  #Directories where the count tables reside.
patternFiles <- c("gene", "isoform")           #Select the right count tables based on a specific pattern
ID_Length <- "Hsapiens_Spneumoniae_SingleStrand_SE/metadata/Hsapiens_Spneumoniae_length.txt" #Filepath with the gene/tx and its length (#1 ID - #2 Length)
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
  pdf(sprintf("%s/TPM_Spearman_Result/tpm_spearman_plot.txt", iDir), onefile=TRUE)
  tpm_spearman_plots
  dev.off()
}
```

<h3 id="section4.10.">4.10.	Scotty heatmap </h3>

The "scotty\_plot.R" script is used to take Scotty's analysis result and create plots which summarises the finding of all the analysis into few plots.

The inputs are the following:

* scotty\_dir -> directory's path where the Scotty results are saved.
Due to the high number of outputs from one Scotty analysis, the files are arranged in the following scheme: Main Directory -> Tool -> Species -> ReadType -> Comparison -> Scotty's result. The "Comparison" directory name has the following scheme: Group1\_Group2, with the underscore between the groups' name.

```
📦Scotty (Main Directory)
 ┣ 📦STAR-HTSeq (The tools used to generate the count table)
 	┣ 📦Hsapiens (The species used to generate the count table)
 	 	┣ 📦Gene (The reads' type present in the table as Gene/Isoform)
 	 		 	┣ 📦wildType_control (Group1_Group2)
 	 		 	┃ ┗ 📜powerPlot.csv
	 		  	┃ ┗ 📜powerBias.csv
	 		  	┃ ┗ 📜nReps.csv
	 		  	┃ ┗ 📜seqDepth.cs 	
	 		  	┃ ┗ 📜...		 	
```

* sampleGroupInfo -> path to the count table information file which has the following columns:
	1. Tool -> the column shows the list of tools used to obtain the count table.
	2. Species -> the column shows the list of species whose reads are in the count table.
	3. Type -> the column shows the type of reads used in the count table (Gene/Isoform).
	4. Sample -> the coumn shows the sample ID.
	5. Group -> the column shows the sample's group name.
	6. Count -> the the column shows the quantified number of reads for the sample.

> **The directories name, from scotty\_dir, and their counterpart in the "sampleGroupInfo" file need to have the same values.**

```
Example for sampleGroupInfo

Tool	Species	Type	Sample	Group	Count
STAR-HTSeq	Host	Gene	SRR3291465	wt30	13155211
STAR-HTSeq	Host	Gene	SRR3291466	wt30	7559415
STAR-HTSeq	Pathogen	Gene	SRR3291465	wt30	49792869
STAR-HTSeq	Pathogen	Gene	SRR3291466	wt30	42519919
STAR-HTSeq	HostPathogen	Gene	SRR3291465	wt30	62948080
STAR-HTSeq	HostPathogen	Gene	SRR3291466	wt30	50079334
...
```

* workCondition_nReps -> Filter Scotty's table results for the selected number of replicates.
* workCondition_seqDepth -> Filter Scotty's table results for the selected sequencing depth.
* outDir -> directory's path where to save the output files. 

The output is one pdf file called "scotty_plot.pdf" saved in the specified output-directory.

```R
#!/scotty_plot.R

#The script takes Scotty results and create heatmap plots for the Power Plot and Power Bias analysis.

#The inputs are:
#1)scotty_dir -> directory's path where the Scotty results are saved.
#Due to the high number of outputs from one Scotty analysis, the files are arranged in the following way:
#Create one main Directory called for example "Scotty". And create sub-directories based on the following scheme:
#Main Directory -> Tool -> Species -> ReadType -> Comparison -> Scotty's result
#The "Comparison" directory name has the following scheme: Group1_Group2.
#Example: Scotty -> STAR-HTSeq -> Hsapiens -> Gene -> wildType_control -> Scotty's result files

#2)sampleGroupInfo -> path to the file. The file has the following columns:
#Tool -> the column shows the list of tools used to obtain the count table used for Scotty analysis.
#Species -> the column shows the list of species whose reads are in the count table used for Scotty analysis.
#Type -> the column shows the type of reads used in the count table used for Scotty analysis such as Gene or Isoform.
#Sample -> the coumn shows the sample name.
#Group -> the column shows the group name the sample is part of.
#Count -> the the column shows the quantified number of reads.

#**** The directories name and their counterpart in the "sampleGroupInfo" file need to have the same values. ****

#3)workCondition_nReps -> Filter Scotty's table results for the selected number of replicates.
#4)workCondition_seqDepth -> Filter Scotty's table results for the selected sequencing depth.
#5)outDir -> directory path where to save the output files. 

#The output is one pdf file called "scotty_plot.pdf" saved in the specified ouput-directory. 

#Load Library
library(tidyverse) #Collection of R packages designed for data science: ggplot2, dplyr, tidyr, readr, purrr, tibble, stringr, forcats.
library(ggforce) #Aid in visual data investigations.

rm(list=ls(all=TRUE))

####### USER's Input ######
scotty_dir <- "Scotty"
sampleGroupInfo <- "Scotty/scotty_infos.txt"
workCondition_nReps <- 2        #Filter for the number of replicates for the two group comparions (from 2 to 10)
workCondition_seqDepth <- 44.6  #Filter for the sequencing depth of the reads expressed in Milion per reads.
outDir <- "Scotty"
###########################

###Global Objects
sampleGroupInfo_df <- read.table(sampleGroupInfo, header=TRUE, sep = "\t") #Get the quantification count for each sample group
plot_list <- list()  #Save all the plots

###Get all files in the Scotty directory and in all its sub-directories
powerPlot_list <- list.files (path=scotty_dir, pattern = "powerPlot.csv", recursive = TRUE)
powerBias_list <- list.files (path=scotty_dir, pattern = "powerBias.csv", recursive = TRUE)
nReps_list <- list.files (path=scotty_dir, pattern = "nReps.csv", recursive = TRUE)
seqDepth_list <- list.files (path=scotty_dir, pattern = "seqDepth.csv", recursive = TRUE)

###Get all information based on the working condition
nInfos <- str_count(powerPlot_list[1], "/")   #Number of Information based on the counnt of directories present by "/"
powerPlotBias_df <- as.data.frame(matrix(data = NA, nrow = length(powerPlot_list), ncol = nInfos + 3))
colnames(powerPlotBias_df) <- c("Tool", "Species", "Type", "Comparison", "PowerPlot", "PowerBias", "Count")

#Check that all the 4 files are present in each directory
if (var(c(length(powerPlot_list), length(powerBias_list), length(nReps_list), length(seqDepth_list))) == 0)
{
  for (i in 1:length(powerPlot_list)) 
  {
    #Get Power Plot & Power Bias Tables
    now_PowerPlot <- read.table(sprintf("%s/%s", scotty_dir, powerPlot_list[i]), header=FALSE, sep = ",") #Power Plot table
    now_PowerBias <- read.table(sprintf("%s/%s", scotty_dir, powerBias_list[i]), header=FALSE, sep = ",") #Power Bias table
    now_nReps <- read.table(sprintf("%s/%s", scotty_dir, nReps_list[i]), header=FALSE, sep = ",")         #Replicate table
    now_seqDepth <- read.table(sprintf("%s/%s", scotty_dir, seqDepth_list[i]), header=FALSE, sep = ",")   #Sequencing Depth table
    
    #Get the working conditions
    nowTableName <- str_replace_all(powerPlot_list[i], pattern = "/", "564752xuqyto4e488688098")
    nowTableName <- head (str_split (nowTableName, pattern = "564752xuqyto4e488688098")[[1]], -1) #Tool #Species #Gene/Isoform #Comparison
    getColumnIndex <- grep (as.character(round(workCondition_seqDepth, -1)), now_seqDepth[1,])
    
    #Add to the main table -> powerPlotBias_df
    powerPlotBias_df[i, c(1:length(nowTableName))] <- nowTableName  #Insert the general information
    powerPlotBias_df$PowerPlot[i] <- now_PowerPlot[workCondition_nReps, getColumnIndex] #PowerPlot
    powerPlotBias_df$PowerBias[i] <- now_PowerBias[workCondition_nReps, getColumnIndex]   #PowerBias
    
    groupComparison <- str_split(nowTableName[4], pattern = "_")[[1]]
    filter_comparison_df <- sampleGroupInfo_df %>% filter(Tool==nowTableName[1] & Species==nowTableName[2] & Type==nowTableName[3] & Group %in% groupComparison)
    powerPlotBias_df$Count[i] <- sum (filter_comparison_df$Count) / nrow(filter_comparison_df)
  }
} else(print("The four files (powerPlot.csv, powerBias.csv, nReps.csv, seqDepth.csv) are not present together in all directories."))


###Create Heatmap plot -> Tools used against the group comparisons
for (iSpecies in unique(powerPlotBias_df$Species)) 
{
  plot_list [[sprintf("%s - PowerPlot", iSpecies)]] <- powerPlotBias_df %>% 
    filter (Species==iSpecies) %>% 
    mutate_at (c(ncol(powerPlotBias_df), ncol(powerPlotBias_df)-1), as.numeric) %>%
    ggplot(aes(x=Tool, y=Comparison, fill=PowerPlot)) + 
    geom_tile() + 
    geom_text(aes(label = round(PowerPlot, 4))) +
    scale_fill_gradient(low="blue", high="red") +
    coord_equal() +
    labs(x="Tools", y="Group Comparisons", title = sprintf("%s - PowerPlot", iSpecies)) +
    theme(axis.text.x = element_text(size=9, face="bold", colour = "black"), axis.text.y = element_text(size=11, face="bold", colour = "black"))
  
  plot_list [[sprintf("%s - PowerBias", iSpecies)]] <- powerPlotBias_df %>% 
    filter (Species==iSpecies) %>% 
    mutate_at (c(ncol(powerPlotBias_df), ncol(powerPlotBias_df)-1), as.numeric) %>%
    ggplot(aes(x=Tool, y=Comparison, fill=PowerBias)) + 
    geom_tile() + 
    geom_text(aes(label = round(PowerBias, 4))) +
    scale_fill_gradient(low="blue", high="red") +
    coord_equal() +
    labs(x="Tools", y="Group Comparisons", title = sprintf("%s - PowerBias", iSpecies)) +
    theme(axis.text.x = element_text(size=9, face="bold", colour = "black"), axis.text.y = element_text(size=11, face="bold", colour = "black"))
}

###Plot -> PowerPlot Vs. PowerBias
plot_list [["PowerPlot - PowerBias"]] <- powerPlotBias_df %>%
  ggplot(aes(x=PowerPlot, y=PowerBias, shape=Comparison, color=Count)) +
  geom_point() +
  facet_grid(cols = vars(Tool), rows = vars(Species)) +
  scale_shape_manual (values=c(0:19)) +      #Shape types
  scale_color_gradient(low="blue", high="red") +
  labs(x="PowerPlot", y="PowerBias", title = "PowerPlot - PowerBias") +
  geom_mark_ellipse(aes(label=Type, group=Type))

###Save plots in one file
pdf(sprintf("%s/scotty_plot.pdf", outDir), onefile=TRUE)
plot_list
dev.off()
```

<h3 id="section4.11.">4.11.	GMT custom annotation for non-supported organisms </h3>

The enrichment analysis is carried out with the g:Profiler tool, but not all organisms are supported, as showcased [here](https://biit.cs.ut.ee/gprofiler/page/organism-list). Uniprot database is used to create GMT custom annotation for non-supported organims. Two R scripts are used. 

<h4 id="section4.11.1.">4.11.1.	Taxon ID search </h4>

The "taxonID\_search.R" script is used to identify the taxon ID of the species for whom the custom annotation will be used for. Assign the scientific name of the species to the "org" character vector. And save the taxon ID of the selected organisms.

```R
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
```

<h4 id="section4.11.2.">4.11.2.	Custom annotation creation </h4>

The "gmt\_custom\_annotation\_creation.R" script is used to create custom annotation for non-supported organisms. 

The input are:

1. taxonID_list -> list of taxon IDs (check "taxonID\_search.R" script).
2. outDir -> the path of the directory where the custom annotation will be saved. Do not add the last "/" to the path.

```R
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
```

<h3 id="section4.12.">4.12.	Polyester's simulated data differential expressed analysis </h3>

The "de\_analysis\_polyester.R" script analyses the merged count tables obtain from the [dual RNA-seq simulation](#section3.10.4.) with DESeq2, edgeR, and NOISeq tools at different number of replicates, row sum cut-off, and adjusted p-values.

The input are:

1. dirPath -> Path to the directory where the Polyester's merged files are saved (DON'T add the last "/" to the path).
2. replicate -> numeric vector list of number of replicates to apply the analysis (>1).
3. cutOff -> numeric vector list used to filter for genes/isoforms which has a row-sum higher than the cut-off value.
4. padj\_FDR\_value ->  numeric vector list of adjusted p-values used to filter for differential expressed genes/isoforms.
5. log2FC\_value -> single value for the Log2 Fold Change used to select differential expressed genes/isoforms.

The output are:

1. ggplots.pdf -> plot pdf file which shows the comparison results.
2. TP\_FP\_FN.txt -> shows information concerning true positive, false positive, etc..
3. pca\_res.txt -> shows the Principal Component Analysis (PCA) results.


```R
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
```

<h3 id="section4.13.">4.13.	Count tables differential expressed analysis </h3>

The de\_countTables\_comparison.R script takes a group of count tables and identifies its differential expressed (DE) elements with DESeq2, edgeR, and NOISeq tools at different values of cut-off, replicates, and adjusted p-values. Then, it compares the identified DE elements to a list furnished by the user. Finally, it creates a plot where it shows the True Positive Rate (TPR) and the number of False Positives of the comparison.

The input are:

1. count\_table\_dir -> Path to the directory where the count tables are saved (DO NOT add the last "/" on the path).
2. patternFile -> Pattern used to select the right count tables.
3. compare\_de -> Path to the file where the list of genes/isoforms is saved for the comparison. The header of the columns has the same name of the group comparison.
4. groupComparison -> List of group comparison written as "Group1_Group2".
5. sample\_group -> Path to the file where the sample's group are saved as two columns called "Sample" and "Group".
6. selectReplicate -> List of values for the number of replicates to use per group during the differential expressed analysis.
7. cutOff -> List of cut-off values for the filtering row (genes/isoforms) step which has a row-sum higher than the cut-off value.
8. padj_FDR_value -> List of adjusted p-values used to filter for differential expressed genes/isoforms.
9. log2FC_value -> Single value for the Log2 Fold Change used to select differential expressed genes/isoforms.

The outputs are:

1. comparison_res.txt -> table with the comparison results with the columns:
	* Table (table's name) 
	* Comparison (comparison name) 
	* Tool (tool used to find DE) 
	* Replicate (number of replicates) 
	* CutOff (cut-off value) 
	* AdjPvalue (adjusted p-value)
	* TruthDE (number of DE from the user's input) 
	* TruePositive (number of DE elements in common between the user's input and the identified DE)
	* FalsePositive (number of DE elements found in the identified DE but absent in the user's input)
	* FalseNegative (number of DE elements found in the user's input but absent in the identified DE)
	* TruePositiveRate (number of True Positive divided by the user's input DE)

2. One PDF plot file for each group comparison.


```R
#!/de_countTables_comparison.R

#The script takes a group of count tables and identifies its differential expressed (DE) elements with DESeq2, edgeR, and NOISeq at different values of cut-off, replicates, and adjusted p-values.
#Then, it compares the identified DE elements to a list furnished by the user.
#Finally, it creates a plot where it shows the True positive rate and the False positive of the comparison.

#The inputs are:
#1)count_table_dir -> Path to the directory where the count tables are saved. DO NOT add the last "/" on the path.
#2)patternFile -> Pattern used to select the right count tables.
#3)compare_de -> Path to the file where the list of genes/isoforms is saved for the comparison. The header of the columns has the same name of the group comparison.
#4)groupComparison -> List of group comparison written as "Group1_Group2".
#5)sample_group -> Path to the file where the sample's group are saved as two columns called "Sample" and "Group".
#6)selectReplicate -> List of values for the number of replicates to use per group during the differential expressed analysis.
#7)cutOff -> List of cut-off values for the filtering row (genes/isoforms) step which has a row-sum higher than the cut-off value.
#8)padj_FDR_value -> List of adjusted p-values used to filter for differential expressed genes/isoforms.
#9)log2FC_value -> Single value for the Log2 Fold Change used to select differential expressed genes/isoforms.

#The outputs are:
#1)comparison_res.txt -> table with the comparison results with the columns:
#Table (table's name) #Comparison (comparison name) #Tool (tool used to find DE) #Replicate (number of replicates) #CutOff (cut-off value) #AdjPvalue (adjusted p-value)
#TruthDE (number of DE from the user's input) 
#TruePositive (number of DE elements in common between the user's input and the identified DE)
#FalsePositive (number of DE elements found in the identified DE but absent in the user's input)
#FalseNegative (number of DE elements found in the user's input but absent in the identified DE)
#TruePositiveRate (number of True Positive divided by the user's input DE)
#2)One PDF plot file for each group comparison.

#Load Library
library(tidyverse)  #Collection of R packages designed for data science: ggplot2, dplyr, tidyr, readr, purrr, tibble, stringr, forcats.
library(DESeq2)     #Differential gene expression analysis based on the negative binomial distribution
library(edgeR)      #Empirical Analysis of Digital Gene Expression Data in R      
library(NOISeq)     #Exploratory analysis and differential expression for RNA-seq data
library(gridExtra)  #Provides functions to work with "grid" graphics, notably to arrange multiple grid-based plots on a page, and draw tables.

rm(list = ls()) #Remove all saved objects

####### USER's Input ######
count_table_dir <- "Polyester_Count_Tables"                   #Path to the directory where the count tables are saved.
patternFile <- "isoforms_ct"                                  #Pattern used to select the right count tables from the directory
compare_de <- "Polyester_Count_Tables/truth_de_isoforms.txt"  #Path to the file where the list of genes/isoforms is saved for the comparison. Column name as Group1_Group2.
sample_group <- "Polyester_Count_Tables/sample_group.txt"     #Path to the file where the sample's group are saved as two columns called "Sample" and "Group".

groupComparison <- c("Case_Control")                          #List of group comparison written as "Group1_Group2".
selectReplicate <- c(2, 10)                                #The highest number of replicates to use per group during the differential expressed analysis
cutOff <- c(0, 3000)                                     #Filter for genes/isoforms which has a row-sum higher than the cut-off value.
padj_FDR_value <- c(0.05, 5e-05)          #Adjusted p-values used to filter for differential expressed genes/isoforms.
log2FC_value <- 1                                             #Single value for the Log2 Fold Change used to select differential expressed genes/isoforms.
###########################

###Global Objects
list_of_files <- list()
de_list <- list()
ggplot_list <- list()

compare_de_df <- read.table(compare_de, header=TRUE, sep = "\t")      #Get the compare DE list file
sample_group_df <- read.table(sample_group, header=TRUE, sep = "\t")  #Get the sample and its group membership file
comparison_de_df <- as.data.frame(matrix(data=NA, ncol=11, nrow=0))   #Comparison Table Result

#Get the count tables
files <- list.files(path=count_table_dir, pattern=patternFile)
for (i in files) {list_of_files[[i]] <- read.table(sprintf("%s/%s", count_table_dir, i), header=TRUE, sep = "\t")}

### Differential expressed analysis ###
for (iCT in 1:length(list_of_files)) #Loop through the count tables
{
  #Tidy the Count Table
  nowTable <- list_of_files[[iCT]]
  rowID <- nowTable[,1]
  nowTable <- nowTable %>% dplyr::select (-c(colnames(nowTable)[1])) %>% mutate_if (is.character,as.numeric)
  rownames(nowTable) <- rowID
  
  for (iGC in 1:length(groupComparison))   #Loop through the group comparisons
  {
    #Get the list of DE to compare to
    de_compare_list <- compare_de_df[, grep(pattern=groupComparison[iGC], names(compare_de_df))]
    
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
            de_list[[sprintf("%s*%s*DESeq2*%s*%s*%s", files[iCT], groupComparison[iGC], iRep, iCO, pv)]] <- rownames(deseq2_deg)
            
            #Add DESeq2 results to the Comparison result dataframe
            tp <- intersect(de_compare_list, rownames(deseq2_deg))
            fp <- setdiff(rownames(deseq2_deg), de_compare_list)
            fn <- setdiff(de_compare_list, rownames(deseq2_deg))
            tpr <- length(tp)/length(de_compare_list)
            comparison_de_df <- rbind(comparison_de_df, c(files[iCT], groupComparison[iGC], "DESeq2", iRep, iCO, pv, length(de_compare_list), length(tp), length(fp), length(fn), tpr))
            
            ###edgeR
            edgeR_deg <- tibble::rownames_to_column(edgeR_res$table, "Gene") %>% filter (FDR < pv, abs(logFC) >= log2FC_value)
            de_list[[sprintf("%s*%s*edgeR*%s*%s*%s", files[iCT], groupComparison[iGC], iRep, iCO, pv)]] <- edgeR_deg$Gene
            
            #Add edgeR results to the Comparison result dataframe
            tp <- intersect(de_compare_list,  edgeR_deg$Gene)
            fp <- setdiff( edgeR_deg$Gene, de_compare_list)
            fn <- setdiff(de_compare_list,  edgeR_deg$Gene)
            tpr <- length(tp)/length(de_compare_list)
            comparison_de_df <- rbind(comparison_de_df, c(files[iCT], groupComparison[iGC], "edgeR", iRep, iCO, pv, length(de_compare_list), length(tp), length(fp), length(fn), tpr))
            
            #NOISeq
            noiseq_deg <- degenes(mynoiseqbio, q = 1-pv, M = NULL)
            noiseq_deg <- subset(noiseq_deg, abs(log2FC) >= 1)
            de_list[[sprintf("%s*%s*NOISeq*%s*%s*%s", files[iCT], groupComparison[iGC], iRep, iCO, pv)]] <- rownames(noiseq_deg)
            
            #Add NOISeq results to the Comparison result dataframe
            tp <- intersect(de_compare_list,  rownames(noiseq_deg))
            fp <- setdiff(rownames(noiseq_deg), de_compare_list)
            fn <- setdiff(de_compare_list, rownames(noiseq_deg))
            tpr <- length(tp)/length(de_compare_list)
            comparison_de_df <- rbind(comparison_de_df, c(files[iCT], groupComparison[iGC], "NOISeq", iRep, iCO, pv, length(de_compare_list), length(tp), length(fp), length(fn), tpr))
          }
        }
      }
    }
  }
}

colnames(comparison_de_df) <- c("Table", "Comparison", "Tool", "Replicate", "CutOff", "AdjPvalue", "TruthDE", "TruePositive", "FalsePositive", "FalseNegative", "TruePositiveRate")

#Convert to numeric: Replicate - FalsePositive - TruePositiveRate 
i <- c(4, 9, 11)
comparison_de_df [ , i] <- apply(comparison_de_df [ , i], 2, function(x) as.numeric(as.character(x)))
comparison_de_df$CutOff <- factor(comparison_de_df$CutOff, levels = unique (as.numeric(as.character(comparison_de_df$CutOff)))) #Set the right order for the CutOff factors
comparison_de_df$AdjPvalue <- factor(comparison_de_df$AdjPvalue, levels = unique (as.numeric(as.character(comparison_de_df$AdjPvalue)))) #Set the right order for the AdjPvalue factors


##### ##### GGPLOT
#Colours
color_pvalue <- c("gray", "black")
color_table <- c("red", "blue", "green", "purple", "orange", "brown", "aquamarine", "yellow", "cyan")
color_total <- c(color_pvalue, color_table) #In case you want to descriminate the min and the max value of the AdjPvalue

#Pvalues
pv <- as.numeric(as.character(unique(comparison_de_df$AdjPvalue)))
min_pvalue <- min (pv)
max_pvalue <- max (pv)

for (iGC in 1:length(groupComparison))   #Loop through the group comparisons
{
  ggplot_list[[sprintf("%s", groupComparison[iGC])]] <- comparison_de_df %>% 
    filter(Comparison==groupComparison[iGC] & AdjPvalue %in% c(min_pvalue, max_pvalue)) %>%   #Filter for Group Comparison and AdjPvalue
    ggplot(aes(x=FalsePositive, y=TruePositiveRate, shape=CutOff, color=Table, group = interaction(Table, AdjPvalue))) + 
    geom_point() +
    geom_line (aes(color=AdjPvalue)) +
    facet_grid(cols = vars(Replicate), rows = vars(Tool)) +
    scale_colour_manual(values=color_total) +   #Colour for both tables and the Pvalues (min & max)
    scale_shape_manual (values=c(0:19)) +
    theme(panel.spacing = unit(1.2, "lines")) +
    labs(x="N* False Positive", y="True Positive Rate")
  
  ggplot_list[[sprintf("%s*AxisX_log10", groupComparison[iGC])]] <- comparison_de_df %>% 
    filter(Comparison==groupComparison[iGC] & AdjPvalue %in% c(min_pvalue, max_pvalue)) %>%   #Filter for Group Comparison and AdjPvalue
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

#Save plots
pdf(sprintf("%s/%s_plot.pdf",count_table_dir, groupComparison[iGC]), onefile=TRUE)
for (i in seq(length(ggplot_list))) {grid.arrange (ggplot_list[[i]])}
dev.off()

#Save Comparison result table
write.table(comparison_de_df, file= sprintf("%s/comparison_res.txt", count_table_dir), row.names=FALSE, sep="\t", quote = FALSE)
```

<h3 id="section4.14.">4.14.	Count tables differential expressed and enrichmment analysis comparison </h3>

The de\_enrichment\_countTables\_comparison.R script takes a list of count tables and it identifies the differential expressed terms with DESeq2 for each group comparison.
Enrichment analysis is carried out for each group comparison and their similarity is showed in a heatmap for all count tables.

The input are:

1. count_table_dir -> list of path to the directory where the count tables are present (#ID #Sample1  #SampleN).
2. patternFile -> Pattern used to select the right count tables.
3. sample_group -> Path to the file where the sample's group are saved as two columns called "Sample" and "Group".
4. groupComparison -> List of group comparison written as "Group1_Group2".
5. cutOff -> Filter for genes/isoforms which has a row-sum higher than the cut-off value.
6. padj_FDR_value -> Adjusted p-values used to filter for differential expressed genes/isoforms.
7. log2FC_value -> Single value for the Log2 Fold Change used to select differential expressed genes/isoforms.
8. gProfiler_supported_org -> List of Scientific Name for g:Profiler supported organisms.
9. gProfiler_non_supported_org -> List of path to GMT custom annotation files for non-supported organisms.
10. baseurl -> REVIGO website link (Do not change)
11. cutoff -> REVIGO's cut-off values (Allowed values: "0.90" - "0.70" - "0.50" - "0.40")
12. isPValue -> Is the REVIGO's input numbers p-values? (Allowed values: "yes" - "no")
13. whatIsBetter -> In case of some other quantity than the pvalues where the "allowed" value is better. Allowed values: "higher" - "lower" - "absolute" - "abs_log"
14. goSizes -> Select a database with GO term sizes (Default is "0" which specifies "whole Uniprot").
15. measure -> Select a semantic similarity measure to use. Allowed values: "RESNIK" - "LIN" - "SIMREL" - "JIANG"

The outputs are:

1. DESeq2_de_res.txt -> DESeq2 result saved for all tables and group comparisons.
2. A directory is created for each species used during the enrichment step with g:Profiler (e.i. "Homo sapiens"). A sub-directory is created for each group comparison (e.i. "Homo sapiens/Case_Control"). And, the following files are saved:
	* gProfiler_enrichment_res.txt -> g:Profiler enrichment analysis results.
	* GOTERM_revigo.txt -> REVIGO's results for the specified GO terms (GO:BP, GO:CC, GO:BP).
	* GOTERM_revigoCluster.txt -> clusters are found for the GO terms found after REVIGO analysis.
	* GOTERM_revigoClusterHeatmap.txt -> ggplot heatmap information for the specified GO term.
	* KEGG_heatmap.txt -> KEGG's heatmap input file for heatmap creation.
	* heatmapPlots.pdf -> One single PDF file where the heatmap plots are saved

```
📦Homo sapiens (Species)
 ┣ 📦Case_Control (Group Comparison)
    ┃ ┗ 📜gProfiler_enrichment_res.txt
    ┃ ┗ 📜GOTERM_revigo.txt
    ┃ ┗ 📜GOTERM_revigoCluster.txt
    ┃ ┗ 📜GOTERM_revigoClusterHeatmap.txt    
    ┃ ┗ 📜KEGG_heatmap.txt
    ┃ ┗ 📜heatmapPlots.pdf
```

```R

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
#10)baseurl -> REVIGO website link (Do not change)
#11)cutoff -> REVIGO's cut-off values (Allowed values: "0.90" "0.70" "0.50" "0.40")
#12)isPValue -> Is the REVIGO's input numbers p-values? (Allowed values: "yes"  "no")
#13)whatIsBetter -> In case of some other quantity than the pvalues where the "allowed" value is better. Allowed values: "higher" "lower" "absolute" "abs_log"
#14)goSizes -> Select a database with GO term sizes (Default is "0" which specifies "whole Uniprot").
#15)measure -> Select a semantic similarity measure to use. Allowed values: "RESNIK" "LIN" "SIMREL" "JIANG"

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

rm(list = ls()) #Remove all saved objects

####### USER's Input ######
count_table_dir <- "Hsapiens_Spneumoniae_SingleStrand_SE/6_count_table/spneumoniae"                 #Path to the directory where the count tables are saved.
patternFile <- "gene"                                #Pattern used to select the right count tables
sample_group <- "Hsapiens_Spneumoniae_SingleStrand_SE/metadata/sim_rep_info.txt"   #Path to the file where the sample's group are saved as two columns called "Sample" and "Group".
groupComparison <- c("1_2")                        #List of group comparison written as "Group1_Group2".

#DESeq2 input
cutOff <- 500                                               #Filter for genes/isoforms which has a row-sum higher than the cut-off value.
padj_FDR_value <- 0.0005                                    #Adjusted p-values used to filter for differential expressed genes/isoforms.
log2FC_value <- 1                                           #Single value for the Log2 Fold Change used to select differential expressed genes/isoforms.

#g:Profiler's supported organims -> (https://biit.cs.ut.ee/gprofiler/page/organism-list)
gProfiler_supported_org <- c("Homo sapiens")                #List of Scientific Name for g:Profiler supported organisms.
gProfiler_non_supported_org <- c("")  #List of path to GMT custom annotation files for non-supported organisms.

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

###Load the files
sample_group_df <- read.table(sample_group, header=TRUE, sep = "\t")  #Get the sample and its group membership file

#Get the count tables
files <- list.files(path=count_table_dir, pattern=patternFile)
for (i in files) {list_of_files[[i]] <- read.table(sprintf("%s/%s", count_table_dir, i), header=TRUE, sep = "\t")}


### Differential expressed analysis ###
deseq2_res_df <- as.data.frame(matrix(data=NA, nrow = 0, ncol = 9))

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
    groupInfo <- str_split(groupComparison, pattern="_")[[1]]
    
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
      
      deseq2_deg <- res %>% filter(padj < padj_FDR_value, abs(log2FoldChange) >= log2FC_value)
      de_list[[sprintf("%s * %s", files[iCT], groupComparison[iGC])]] <- rownames(deseq2_deg)
      
      deseq2_deg$CountTable <- files[iCT]
      deseq2_deg$Comparison <- groupComparison[iGC]
      deseq2_deg <- deseq2_deg %>% rownames_to_column(var = "ID")
      deseq2_res_df <- rbind(deseq2_res_df, deseq2_deg)
    }
  }
}

#Save the DESeq2 results
write.table(deseq2_res_df, file=sprintf("%s/DESeq2_%s_result.txt", count_table_dir, patternFile), row.names=FALSE, sep="\t", quote = FALSE)


####### Enrichment Analysis g:Profiler & REVIGO ######
list_of_org <- c(gProfiler_supported_org, gProfiler_non_supported_org)
list_of_org <- list_of_org[list_of_org != ""]

#Loop through the group comparisons
iGC <- groupComparison
iSpecies <- list_of_org[2]

for (iGC in groupComparison) 
{
  #Get the elements from de_list which has the same group comparison
  de_gc <- de_list[grep(iGC, names(de_list))]
  
  #Loop through each species in the list_of_org vector list
  for (iSpecies in list_of_org)
  {
    ###Enrichment Analysis -> g:Profiler
    gProfiler_Res <- NULL
    
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
      gmt_file <- upload_GMT_file(iSpecies)
      #go_bacteria <- data.frame(matrix(ncol = 8, nrow = 0))
      #colnames(go_bacteria) <- c("p_value","term_size","query_size","intersection_size","term_id","source", "term_name","Group")
      
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
      dir.create(file.path(count_table_dir, iSpecies), showWarnings = FALSE)  #Create directory for the species
      dir.create(file.path(sprintf("%s/%s", count_table_dir, iSpecies), iGC), showWarnings = FALSE)  #Create directory for the Group Comparison
      res_directory <- sprintf("%s/%s/%s", count_table_dir, iSpecies, iGC)
      
      #Save the g:Profiler enrichment analysis results for GO and KEGG terms
      enrich_go_kegg <- gProfiler_Res$result[,c(1,3,9,10,11)] %>% filter(source %in% c("GO:BP", "GO:CC", "GO:MF", "KEGG"))
      write.table(enrich_go_kegg, file = sprintf("%s/gProfiler_enrichment_res.txt", res_directory), row.names=FALSE, sep="\t", quote = FALSE)
      
      ##### REVIGO Analysis #####
      revigo_session <- html_session(baseurl)
      revigo_form <- html_form(revigo_session)[[1]] 
      
      # Loop through the GO terms
      for (go in c("GO:BP", "GO:CC", "GO:MF")) 
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
        
        # GENE ONTOLOGY REVIGO CLUSTER ANALYSIS
        #Find the clusters present in the REVIGO results based on their sematic space
        clusterDF <- uniqueGO %>% remove_rownames %>% column_to_rownames(var="term_ID")
        if (nCluster >= nrow(clusterDF)) {nCluster <- nrow(clusterDF)-1}
        
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
        write.table(heatmapGO, file = sprintf("%s/%s_revigoClusterHeatmap.txt", res_directory, go), row.names=FALSE, sep="\t", quote = FALSE)   #Save the REVIGO_Cluster_Heatmap result
        
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
      heatmapKEGG <- data.frame(matrix(ncol = length (unique(filterKEGG$query)), nrow = length (unique (filterKEGG$term_id))))
      colnames(heatmapKEGG) <- unique(filterKEGG$query)
      rownames(heatmapKEGG) <- unique (filterKEGG$term_name)
      
      for (i in 1:nrow(filterKEGG)) 
      {heatmapKEGG[filterKEGG$term_name[i], filterKEGG$query[i]] <- filterKEGG$p_value[i]}
      write.csv(heatmapKEGG, file = sprintf("%s/KEGG_heatmap.txt", res_directory), row.names=FALSE, sep="\t", quote = FALSE)   #Save the gProfiler_Heatmap result
      
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
      
      ##### Save Heatmap plots #####
      pdf(sprintf("%s/heatmapPlots.pdf", res_directory), onefile=TRUE)
      for (i in seq(length(enrichPlot))) {grid.arrange (enrichPlot[[i]])}
      dev.off()
    }
  }
}
```

<h3 id="section4.15.">4.15.	R session Info </h3>

The R session information is showcased down-below.

```
> sessionInfo()
R version 4.0.2 (2020-06-22)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Catalina 10.15.7

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib

locale:
[1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8

attached base packages:
 [1] grid      splines   stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] scales_1.1.1                stringi_1.5.3               tidytext_0.2.6              rvest_0.3.6                 xml2_1.3.2                  gprofiler2_0.2.0           
 [7] ggforce_0.3.2               ggVennDiagram_0.3           rlang_0.4.9                 ggpubr_0.4.0                VennDiagram_1.6.20          futile.logger_1.4.3        
[13] GenomicFeatures_1.40.1      AnnotationDbi_1.50.3        rtracklayer_1.48.0          forestmangr_0.9.2           tximport_1.16.1             gridExtra_2.3              
[19] factoextra_1.0.7            NOISeq_2.31.0               Matrix_1.3-0                edgeR_3.30.3                limma_3.44.3                DESeq2_1.28.1              
[25] SummarizedExperiment_1.18.2 DelayedArray_0.14.1         matrixStats_0.57.0          Biobase_2.48.0              GenomicRanges_1.40.0        GenomeInfoDb_1.24.2        
[31] forcats_0.5.0               stringr_1.4.0               dplyr_1.0.2                 purrr_0.3.4                 readr_1.4.0                 tidyr_1.1.2                
[37] tibble_3.0.4                ggplot2_3.3.3               tidyverse_1.3.0             polyester_1.24.0            Biostrings_2.56.0           XVector_0.28.0             
[43] IRanges_2.22.2              S4Vectors_0.26.1            BiocGenerics_0.34.0        

loaded via a namespace (and not attached):
  [1] readxl_1.3.1             backports_1.2.1          BiocFileCache_1.12.1     selectr_0.4-2            plyr_1.8.6               lazyeval_0.2.2           BiocParallel_1.22.0     
  [8] SnowballC_0.7.0          digest_0.6.27            htmltools_0.5.0          fansi_0.4.1              magrittr_2.0.1           memoise_1.1.0            openxlsx_4.2.3          
 [15] annotate_1.66.0          modelr_0.1.8             askpass_1.1              prettyunits_1.1.1        colorspace_2.0-0         blob_1.2.1               rappdirs_0.3.1          
 [22] ggrepel_0.9.0            haven_2.3.1              crayon_1.3.4             RCurl_1.98-1.2           jsonlite_1.7.2           genefilter_1.70.0        survival_3.2-7          
 [29] glue_1.4.2               polyclip_1.10-0          gtable_0.3.0             zlibbioc_1.34.0          car_3.0-10               Rhdf5lib_1.10.1          abind_1.4-5             
 [36] futile.options_1.0.1     DBI_1.1.0                rstatix_0.6.0            Rcpp_1.0.5               viridisLite_0.3.0        xtable_1.8-4             progress_1.2.2          
 [43] units_0.6-7              foreign_0.8-81           bit_4.0.4                htmlwidgets_1.5.3        httr_1.4.2               RColorBrewer_1.1-2       logspline_2.1.16        
 [50] ellipsis_0.3.1           pkgconfig_2.0.3          XML_3.99-0.5             farver_2.0.3             dbplyr_2.0.0             locfit_1.5-9.4           utf8_1.1.4              
 [57] tidyselect_1.1.0         labeling_0.4.2           reshape2_1.4.4           munsell_0.5.0            cellranger_1.1.0         tools_4.0.2              cli_2.2.0               
 [64] generics_0.1.0           RSQLite_2.2.1            broom_0.7.3              yaml_2.2.1               bit64_4.0.5              fs_1.5.0                 zip_2.1.1               
 [71] formatR_1.7              biomaRt_2.44.4           tokenizers_0.2.1         compiler_4.0.2           rstudioapi_0.13          plotly_4.9.2.2           curl_4.3                
 [78] e1071_1.7-4              ggsignif_0.6.0           reprex_0.3.0             tweenr_1.0.1             geneplotter_1.66.0       lattice_0.20-41          classInt_0.4-3          
 [85] vctrs_0.3.6              pillar_1.4.7             lifecycle_0.2.0          BiocManager_1.30.10      data.table_1.13.2        cowplot_1.1.1            bitops_1.0-6            
 [92] R6_2.5.0                 KernSmooth_2.23-18       rio_0.5.16               janeaustenr_0.1.5        lambda.r_1.2.4           MASS_7.3-53              assertthat_0.2.1        
 [99] rhdf5_2.32.4             openssl_1.4.3            withr_2.3.0              GenomicAlignments_1.24.0 Rsamtools_2.4.0          GenomeInfoDbData_1.2.3   hms_0.5.3               
[106] class_7.3-17             carData_3.0-4            sf_0.9-6                 lubridate_1.7.9.2  
```

<h2 id="section5.">5.	Polyester: dual RNA-seq simulated data </h2>

RNA-seq data was simulated with Polyester for one host species (human) and one bacterial pathogen (Streptococcus pneumoniae D39).

<h3 id="section5.1.">5.1.	RNA-seq simulation data creation </h3>

The polyester\_sim.R script was used to create simulated RNA-seq for both species. The inputs are: 

* The "txome\_list" vector character with the filenames of the transcriptome files.

```
txome_list <- c("reference/hsapiens/Homo_sapiens.GRCh38.cdna.all.fa", "reference/spneumoniae/Streptococcus_pneumoniae_d39.ASM1436v1.cdna.all.fa")
```

* The "gene\_selection" object has the path to the file with the list of genes for the transcripts filtering step. The headers of the file have the same name of the transcriptome fasta filename and the list of genes are saved column-wise.

```
head head gene_selection_file.txt

Homo_sapiens.GRCh38.cdna.all.fa	Streptococcus_pneumoniae_d39.ASM1436v1.cdna.all.fa
ENSG00000137154	SPD_1318
ENSG00000156508	SPD_0636
ENSG00000197956	SPD_0420
ENSG00000202198	SPD_1823
```

* The "polyester\_settings" object has path to the file which has the setting values for the Polyester simulation.

```
cat polyester_setting_input.txt

	DE_TXs	Coverage	Min_FC	Max_FC	Read_Length	Paired_End	Strand_Specific	Replicate_Group1	Replicate_Group2
Homo_sapiens.GRCh38.cdna.all.fa	10	10	2	20	75	FALSE	TRUE	10	10
Streptococcus_pneumoniae_d39_gca_000014365.ASM1436v2.cdna.all.fa	10	3	2	20	75	FALSE	TRUE	10	10
```

* The "outputDir" object has the path to the directory where the outputs will be saved. The directory will be saved in case it does not exist.

```
outputDir <- "Polyester_Hs_Sp"
```

<h3 id="section5.2.">5.2.	Merge fasta file </h3>

The simulated fasta file samples are merged together depending on the sample name.


<h4 id="section5.2.1.">5.2.1.	Single-end sample </h4>	

Polyester creates single-end samples with the following pattern: sample_01.fasta.gz. In order to simulate dual RNA-seq data, samples with the same name in different species directories will be merged together. The following bash script is used and it needs the following information which is highlighted under the "USER INPUT" description. The script creates a new directory, given by the user in the script, where the merged fasta files are saved.

```bash
#!/bin/bash

#USER INPUT -> Save each simulated directory in a variable
d1="Homo_sapiens.GRCh38.cdna.all.fa.selectedTXs_cov_22_rep_4_4"
d2="Streptococcus_pneumoniae_d39.ASM1436v1.cdna.all.fa.selectedTXs_cov_200_rep_4_4"

#USER INPUT -> Save all directory's variable in a list
dirList=($d1 $d2)

#USER INPUT -> Save all sample's pattern
s1="sample_01"
s2="sample_02"
s3="sample_03"
s4="sample_04"

#USER INPUT -> Save all sample's pattern variables in a list
sampleList=($s1 $s2 $s3 $s4)

#USER INPUT -> create new directory where to save the merged fasta file
dirOut="Polyester_dual_RNAseq"

#Create new directories
mkdir $dirOut

#Merge the fasta file with the same sample name in the different directories
for s in "${sampleList[@]}"      ### Sample Loop ###
do
    declare -a sampleList=()
    for d in "${dirList[@]}" ### Directory loop ###
    do
    sampleList+=($(ls $d/$s*))
    done
    mergeCommand="cat "${sampleList[*]}" >> $dirOut/$s.fasta.gz"
    eval "$mergeCommand"
done
```

<h4 id="section5.2.2.">5.2.2.	Paired-end sample </h4>

Polyester creates paired-end samples with the following pattern: sample_01_1.fasta.gz and sample_01_2.fasta.gz. In order to simulate dual RNA-seq data, samples with the same name in different directories will be merged together. The following bash script is used. The script needs the following information which is highlighted under the "USER INPUT" description. The script creates a new directory where merged samples for the forward reads and the reverse reads are saved.


```bash
#!/bin/bash

#USER INPUT -> Save each simulated directory in a variable
d1="Homo_sapiens.GRCh38.cdna.all.fa.selectedTXs_cov_22_rep_4_4"
d2="Streptococcus_pneumoniae_d39.ASM1436v1.cdna.all.fa.selectedTXs_cov_200_rep_4_4"

#USER INPUT -> Save all directory's variable in the dirList object list
dirList=($d1 $d2)

#USER INPUT -> Save sample's pattern for the forward and reverse file
s1="sample_01"
s2="sample_02"
s3="sample_03"
s4="sample_04"

#USER INPUT -> Save all sample's pattern variables in a list
sampleList=($s1 $s2 $s3 $s4)

#USER INPUT -> create new directory where to save the merged fasta file
dirOut="polyester_merged_RNAseq"

#Create new directory
mkdir $dirOut

#Merge the fasta files with the same sample name for both the forward and reward
for s in "${sampleList[@]}"      ### Sample Loop ###
do
    declare -a forwardList=()
    declare -a reverseList=()
    for d in "${dirList[@]}" ### Directory loop ###
    do
        forwardList+=($(ls $d/$s*1*))
        reverseList+=($(ls $d/$s*2*))
    done
    forwardCommand="cat "${forwardList[*]}" >> $dirOut/$s.1.fasta.gz"
    reverseCommand="cat "${reverseList[*]}" >> $dirOut/$s.2.fasta.gz"
    eval "$forwardCommand"
    eval "$reverseCommand"
done
```

<h3 id="section5.3.">5.3.	Merged data information analysis </h3>

The [de\_analysis\_polyester.R](#section4.12.) script is used only for the Polyester simulated merged information directory. The script identifies differential expressed elements from the merged count tables obtain from the Polyester simulation script with DESeq2, edgeR, and NOISeq tools at different number of replicates, row sum cut-off, and adjusted p-values. The input are:

* polyester_merged_dir -> directory where the Polyester's merged files are saved 
* replicate -> numeric vector list of number of replicates to apply the analysis (>1).
* cutOff -> numeric vector list used to filter rows which has a row-sum higher value.
* padj\_FDR\_value -> numeric vector list of adjusted p-values.
* log2FC_value -> Log2 Fold Change value used to select differential expressed genes/isoforms.

```
polyester_merged_dir <- "Hsapiens_Spneumoniae_SingleStrand_SE/polyester_merged_sim"
replicate <- c(2,10)
cutOff <- c(0, 5000)
padj_FDR_value <- c(0.05, 0.005, 5e-04, 5e-05, 5e-06)
log2FC_value <- 1
```


<h2 id="section6.">6.	Dual RNA-seq pipeline </h2>

The first step of the analysis is the creation of directories, as showcased below. 

```
📦Author_Year
 ┣ 📦0_raw_data
 	┣ 📦qc
 ┣ 📦1_trimmed_data
 	┣ 📦qc
 ┣ 📦2_star_host
 	┣ 📦htseq
 	┣ 📦rsem
 ┣ 📦3_star_pathogen
 	┣ 📦htseq
 	┣ 📦rsem
 ┣ 📦4_salmon_host
 ┣ 📦5_salmon_pathogen
 ┣ 📦6_kallisto_host
 ┣ 📦7_kallisto_pathogen
 ┣ 📦8_count_table
 	┣ 📦host
 	┣ 📦pathogen
```

The project's main directory is called Author_Year. The raw fastq files are saved in the 0\_raw\_data directory. The raw data can be downloaded by following the [FastQ download pipeline](#section3.3.).

```bash
#!/bin/bash

###User input
mainDir="Aprianto_2016"

###Create directories 
$	mkdir $mainDir
$	cd $mainDirs
$	mkdir 0_raw_data
$	mkdir 0_raw_data/qc
$	mkdir 1_trimmed_data
$	mkdir 1_trimmed_data/qc
$	mkdir 2_star_host
$	mkdir 2_star_host/htseq
$	mkdir 2_star_host/rsem
$	mkdir 3_star_pathogen
$	mkdir 3_star_pathogen/htseq
$	mkdir 3_star_pathogen/rsem
$	mkdir 4_salmon_host
$	mkdir 5_salmon_pathogen
$	mkdir 6_kallisto_host
$	mkdir 7_kallisto_pathogen
$	mkdir 8_count_table
$	mkdir 8_count_table/host
$	mkdir 8_count_table/pathogen
```

Both single-end and paired-end pipelines has the same steps:

1. Pre-processing step: The trimming step is performed on the raw data. Whereas, quality control is performed after and before the trimming.
2. STAR mapping step: Align the trimmed reads to the reference genome by the STAR mapper. The HTSeq-count gene quantification is done during the mapping.
3. RSEM quantification step: gene and transcript level.
4. Salmon reference-free quantification step: gene and transcript level.
5. Kallisto reference-free quantification step: gene and transcript level.
6. Analysis on the count tables obtained. 


<h3 id="section6.1.">6.1. Single-end reads </h3>

First, compressed the raw data (e.i. fasta.gz, fastq.gz), saved in the 0_raw_data directory. And run the following bash script.

```bash
#!/bin/bash

#####User input#####
#Tool's path
fastqc="fastqc"
multiqc="multiqc"
trimmomatic="trimmomatic"
star="star"
rsem-calculate-expression="rsem-calculate-expression"
salmon="salmon"
kallisto="kallisto"

#General
ngs_format=".fastq.gz" #RNA seq reads format
hostRef="rootDirectory/../reference/hsapiens" #Host reference directory
pathogenRef="rootDirectory/../reference/spneumoniae_D39" #Pathogen reference directory

#Pre-processing: trimmomatic
adapter="rootDirectory/../TruSeq3-SE.fa" #Adapter fasta file absolute path
illumina_clip="2:30:10" #Cut adapter and other illumina-specific sequences from the read
sw="5:20" #Sliding window & Cut within the window once the average quality falls below the threshold.
minlen="50" #Drop the read if it is below a specified length

#RSEM options
fp="1.0" #RSEM's forward_prob setting
avReadLength="75"

#Kallisto
avReadLength="75" #Average Read Length
sdReadLen="20" #Estimated standard deviation of fragment length
##################


#####Pre-processing#####
#Raw Data Quality Control
$fastqc 0_raw_data/*$ngs_format -o 0_raw_data/qc
$multiqc 0_raw_data/qc -o 0_raw_data/qc

#Raw Data Trimming
cd 0_raw_data
for f in *$ngs_format
do
	$trimmomatic SE $f ../1_trimmed_data/${f%%$ngs_format}"_trim.$ngs_format" ILLUMINACLIP:$adapter:$illumina_clip SLIDINGWINDOW:$sw MINLEN:$minlen;
done
cd ../

#Trimmed Data Quality Control
$fastqc 1_trimmed_data/*$ngs_extension -o 1_trimmed_data/qc
$multiqc 1_trimmed_data/qc -o 1_trimmed_data/qc


#####STAR mapping#####
#Host
cd 1_trimmed_data
for f in *$ngs_extension
do 
	$star --genomeDir $hostRef/star --readFilesIn $f --quantMode TranscriptomeSAM GeneCounts --outFilterMismatchNoverLmax 0.05 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate    --outFileNamePrefix ../2_star_host/${f%%$ngs_extension}"."; 
done
cd ../
mv 2_star_host/*ReadsPerGene* 2_star_host/htseq

#Pathogen
cd 1_trimmed_data

for f in *$ngs_extension
do 
	$star --genomeDir $pathogenRef/star --readFilesIn $f --readFilesCommand zcat --alignIntronMax 1 --outSAMtype BAM SortedByCoordinate --outFilterMismatchNoverLmax 0.05 --quantMode TranscriptomeSAM GeneCounts --outFileNamePrefix ../3_star_pathogen/${f%%$ngs_extension}"."; 
done
cd ../
mv 3_star_pathogen/*ReadsPerGene* 3_star_pathogen/htseq


#####RSEM quantification#####
#Host
cd  2_star_host
for f in *toTranscriptome*
do 
	$rsem-calculate-expression --bam --no-bam-output --estimate-rspd --forward-prob $fp $f $hostRef/rsem/index rsem/${f%%.toTranscriptome.out.bam};
done
cd ..

#Pathogen
cd 3_star_pathogen
for f in *toTranscriptome*
do 
	$rsem-calculate-expression --bam --no-bam-output --estimate-rspd --forward-prob $fp $f $pathogenRef/rsem/index rsem/${f%%.toTranscriptome.out.bam};
done
cd ..


#####Salmon#####
#Host
cd 1_trimmed_data
for f *$ngs_extension
do
	$salmon quant -i $hostRef/salmon -l A -r $f -p 8 --validateMappings -o ../4_salmon_host/${f%%$ngs_extension} --numBootstraps 10;
done
cd ..

#Pathogen
cd 1_trimmed_data
for f in *$ngs_extension
do
	$salmon quant -i $pathogenRef/salmon -l A -r $f -p 8 --validateMappings -o ../5_salmon_pathogen/${f%%$ngs_extension} --numBootstraps 10;
done
cd ../


#####Kallisto#####
#Host
cd 1_trimmed_data
for f in *$ngs_extension 
do
	$kallisto quant -i $hostRef/kallisto/index -l $avReadLength -o ../6_kallisto_host/${f%%$ngs_extension} -b 10 --single -s $sdReadLen $f;
done
cd ../

#Pathogen
cd 1_trimmed_data
for f in *$ngs_extension
do
	$kallisto quant -i $pathogenRef/kallisto/index -l $avReadLength -o ../7_kallisto_pathogen/${f%%$ngs_extension} -b 10 --single  -s $sdReadLen $f;
done
cd ..
```

<h3 id="section6.2.">6.2. Paired-end reads </h3>

Paired-end samples have the same filename pattern. For example, the files "sample\_01\_2.fastq.gz" and  "sample\_01\_2.fastq.gz" are paired-end samples. The sub-string "sample\_01" identifies the paired ones; and, the use of "1" and "2" differentiates between the forward and the reverse, respectively. Due to the paired-end samples, the user needs to give the list of IDs used to identify the paired-end samples such as "sample\_01".

```bash
#!/bin/bash

#####User input#####

#Tool's path
fastqc="fastqc"
multiqc="multiqc"
trimmomatic="trimmomatic"
star="star"
rsem-calculate-expression="rsem-calculate-expression"
salmon="salmon"
kallisto="kallisto"

#General information
ngs_format=".fastq.gz" #RNA seq reads format
hostRef="rootDirectory/../reference/hsapiens"
pathogenRef="rootDirectory/../reference/spneumoniae_D39"

#Pre-processing: trimmomatic
adapter="rootDirectory/../TruSeq3-SE.fa" #Adapter fasta file absolute path
illumina_clip="2:30:10" #Cut adapter and other illumina-specific sequences from the read
sw="5:20" #Sliding window & Cut within the window once the average quality falls below the threshold.
minlen="50" #Drop the read if it is below a specified length

#RSEM options
fp="1.0" #RSEM's forward_prob setting
avReadLength="75"

#Kallisto
avReadLength="75" #Average Read Length
sdReadLen="20" #Estimated standard deviation of fragment length

#Sample
sampleList=("sample_01" "sample_02" "sample_03" "sample_04" "sample_05" "sample_06" "sample_07" "sample_08" "sample_09" "sample_10")
##################


#####Pre-processing#####
#Raw Data Quality Control
$fastqc 0_raw_data/*$ngs_format -o 0_raw_data/qc
$multiqc 0_raw_data/qc -o 0_raw_data/qc

#Raw Data Trimming
cd 0_raw_data
for s in "${sampleList[@]}"
do
	$trimmomatic PE -p 8 $s*1* $s*2* $s*1"_paired$ngs_format" $s*1"_unpaired$ngs_format" $s*2"_paired"$ngs_format $s*2"_unpaired"$ngs_format ILLUMINACLIP:$adapter:$illumina_clip SLIDINGWINDOW:$sliding_window MINLEN:$minlen;
done
cd ../
mv 0_raw_data/*paired* 1_trimmed_data

#Trimmed Data Quality Control
$fastqc 1_trimmed_data/*_paired* -o 1_trimmed_data/qc
$multiqc 1_trimmed_data/qc -o 1_trimmed_data/qc


#####STAR mapping#####
#Host
cd 1_trimmed_data
for s in "${sampleList[@]}"
do
	$star --genomeDir $hostRef --readFilesIn $s*1_paired* $s*2_paired --quantMode TranscriptomeSAM GeneCounts --outFilterMismatchNoverLmax 0.05 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ../2_star_host/$s".";
done
cd ../
mv 2_star_host/*ReadsPerGene* 2_star_host/htseq

#Pathogen
cd 1_trimmed_data
for s in "${sampleList[@]}"
do 
	$star --genomeDir $pathogenRef --readFilesIn $s*1_paired* $s*2_paired* --quantMode TranscriptomeSAM GeneCounts --outFilterMismatchNoverLmax 0.05 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --alignIntronMax 1 --outFileNamePrefix ../3_star_pathogen/$s".";
done
cd ../
mv 3_star_pathogen/*ReadsPerGene* 3_star_pathogen/htseq


#####RSEM quantification#####
#Host
cd 2_star_host
for s in *toTranscriptome*
do
	s0=${s%%Aligned.toTranscriptome.out.bam}
	$rsem-calculate-expression -paired-end --forward-prob $fp -p 8 --bam --no-bam-output --estimate-rspd $s $hostRef/rsem/index rsem/$s0; 
done
cd ../

#Pathogen
cd 2_star_pathogen
for s in *toTranscriptome*
do
	s0=${s%%.Aligned.toTranscriptome.out.bam}
	$rsem-calculate-expression --paired-end --forward-prob $fp -p 8 --bam --no-bam-output --estimate-rspd $f $pathogenRef/rsem/index rsem/$s0; 
done
cd ../


#####Salmon#####
#Host
cd 1_trimmed_data
for s in "${sampleList[@]}"
do
	$salmon quant -i $hostRef/salmon -l A -1 $s*1_paired* -2 $s*2_paired* -p 8 --validateMappings -o ../4_salmon_host/$s --numBootstraps 10;
done
cd ../

#Pathogen
cd 1_trimmed_data
for s in "${sampleList[@]}"
do
	$salmon quant -i $pathogenRef/salmon -l A -1 $s*1_paired* -2 $s*2_paired* -p 8 --validateMappings -o ../5_salmon_pathogen/$s --numBootstraps 10;
done
cd ../


#####Kallisto#####
#Host
cd 1_trimmed_data
for s in "${sampleList[@]}"
do
    kallisto quant -i $hostRef/index -o ../6_kallisto_host/$s -b 10 $s*1_paired* $s*2_paired*;
done
cd ../

#Pathogen
cd 1_trimmed_data
for s in "${sampleList[@]}"
do
    $kallisto quant -i $pathogenRef/index -o ../7_kallisto_pathogen/$s -b 10 $s*1_paired* $s*2_paired*;
done
cd ../
```

<h3 id="section6.3.">6.3. Merge count tables </h3>

The following R scripts were used to merge the count tables obtained from the same quantification tool.

<h4 id="section6.3.1.">6.3.1. HTSeq-count count table merge </h4>

The [htseq\_merge\_table.R](#section4.5.) script merges HTSeq-counts count table. The scripts outputs three text files corresponding to the three type of RNA sequencing: un-stranded, stranded, and reverse stranded. The input are:

* Host & Pathogen

```
pathFiles <- c("rootDirectory/../Author_Year/2_star_host/htseq", "rootDirectory/../Author_Year/3_star_pathogen/htseq")
patternFiles <- "sample"
```

<h4 id="section6.3.2.">6.3.2. RSEM count table merge </h4>

The [rsem\_merge\_table.R](#section4.6.) script merges RSEM count tables. The input are:

* Host & Pathogen

```
pathFiles <- c("rootDirectory/../Author_Year/3_star_host/rsem", "rootDirectory/../Author_Year/3_star_pathogen/rsem")
```

<h4 id="section6.3.3.">6.3.3. Salmon count table merge </h4>

The [salmon\_merge\_table.R](#section4.7.) script merges Salmon count tables. Check [4.3. Gene-Transcript relation and length](#section4.3.) to create the file with the gene and transcript's length used for the gene_tx object.

* Host input

```
fileDirectory <- "rootDirectory/../Author_Year/4_salmon_host"
gene_tx <- "rootDirectory/../reference/host/txome_annotation/tx_gene_length.txt"
```

* Pathogen input

```
fileDirectory <- "rootDirectory/../Author_Year/5_salmon_pathogen"
gene_tx <- "rootDirectory/../reference/pathogen/txome_annotation/tx_gene_length.txt"
```

<h4 id="section6.3.4.">6.3.4. Kallisto count table merge </h4>

The [kallisto\_merge\_table.R](#section4.8.) script merges Kallisto count tables. Check [4.3. Gene-Transcript relation and length](#section4.3.) to create the file with the gene and transcript's length used for the gene_tx object.

* Host input

```
fileDirectory <- "rootDirectory/../Author_Year/6_kallisto_host"
gene_tx <- "rootDirectory/../reference/host/txome_annotation/tx_gene_length.txt"
```

* Pathogen input

```
fileDirectory <- "rootDirectory/../Author_Year/7_kallisto_pathogen"
gene_tx <- "rootDirectory/../reference/pathogen/txome_annotation/tx_gene_length.txt"
```

<h3 id="section6.4.">6.4. Move count table </h3>

Copy and rename the merged count tables to the right 8\_count\_table directory.

```bash
#!/bin/bash

#Host merged tables -> Gene & Isoform
cp 2_star_host/htseq/*strand_yes* 8_count_table/host/STAR-HTSeq_host_gene.txt
cp 2_star_host/rsem/*gene*txt 8_count_table/host/RSEM_host_gene.txt
cp 2_star_host/rsem/*isoform*txt 8_count_table/host/RSEM_host_isoform.txt
cp 4_salmon_host/*gene*txt 8_count_table/host/Salmon_host_gene.txt
cp 4_salmon_host/* isoform*txt 8_count_table/host/Salmon_host_ isoform.txt
cp 6_kallisto_host/*gene*txt 8_count_table/host/Kallisto_host_gene.txt
cp 6_kallisto_host/*isoform*txt 8_count_table/host/Kallisto_host_isoform.txt

#Pathogen merged tables -> Gene
cp 3_star_pathogen/htseq/*strand_yes* 8_count_table/pathogen/STAR-HTSeq_pathogen_gene.txt
cp 3_star_pathogen/rsem/*gene*txt 8_count_table/pathogen/RSEM_pathogen_gene.txt
cp 3_star_pathogen/rsem/*isoform*txt 8_count_table/pathogen/RSEM_pathogen_isoform.txt
cp 5_salmon_pathogen/*gene*txt 8_count_table/pathogen/Salmon_pathogen_gene.txt
cp 5_salmon_pathogen/*isoform*txt 8_count_table/pathogen/Salmon_pathogen_isoform.txt
cp 7_kallisto_pathogen/*gene*txt 8_count_table/pathogen/Kallisto_pathogen_gene.txt
cp 7_kallisto_pathogen/*isoform*txt 8_count_table/pathogen/Kallisto_pathogen_isoform.txt
```

<h3 id="section6.5.">6.5. TPM Spearman correlation </h3>

The [tpm\_spearman\_cor.R](#section4.9.) script calculates the TPM values for each count table and calculate the Spearman correlation between the tables. The TPM Spearman correlation is calculated between count table that have the same read type (Gene/Isoform) and come from the same species (#Host & Pathogen). The input are:

1. pathFiles -> directory's path where the count tables reside. The count tables need to have the same headers and in the same order. And, the first column need to have the gene/transcript ID.
2. patternFiles -> select the right count tables based on the pattern.
3. ID\_Length -> path to the file where the gene/transcript and its length are saved (#1 ID - #2 Length). Such file can be created by using [4.3. Gene-Transcript relation and length](#section4.3.) script.
4. mean\_read\_seq -> mean of the fragment length sequenced.

* Host

```
pathFiles <- "rootDirectory/../Author_Year/8_count_table/host"
patternFiles <- c("gene", "isoform")
ID_Length <- "rootDirectory/../reference/hsapiens/txome_annotation/tx_gene_length.txt"
mean_read_seq <- 75
```

* Pathogen

```
pathFiles <- "rootDirectory/../Author_Year/8_count_table/pathogen"
patternFiles <- c("gene", "isoform")
ID_Length <- "rootDirectory/../reference/hsapiens/txome_annotation/tx_gene_length.txt"
mean_read_seq <- 75
```

<h3 id="section6.6.">6.6. Scotty analysis </h3>

The Scotty can be executed only in the Matlab software:

* Open Matlab.
* Change Matlab's working directory to the Scotty directory.
* Create a text file (e.i. scotty.txt) as the count table where the gene/transcript ID is the first column; followed by samples of Group 1, then followed by samples of Group 2.
* Save the scotty.txt file in the Scotty directory.

* Run the following command into Matlab (check the ReadMe.pdf file for additional information) for the desired group comparisons.

```
scottyEstimate('./scotty.txt', '2', '2', '1990634', '2.0', '0.05', '10', '0', '0', '0', 'Inf', '10', '10000000', '100000000', '50', '50', '87', 'Scotty');
```

* Save the files as discussed in [4.10. Scotty heatmap](#section4.10.).

```
📦Scotty (Main Directory)
 ┣ 📦STAR-HTSeq (The tools used to generate the count table)
    ┣ 📦Hsapiens (The species used to generate the count table)
        ┣ 📦Gene (The reads' type present in the table as Gene/Isoforms)
                ┣ 📦wildType_control (Group1_Group2)
                ┃ ┗ 📜powerPlot.csv
                ┃ ┗ 📜powerBias.csv
                ┃ ┗ 📜nReps.csv
                ┃ ┗ 📜seqDepth.cs    
                ┃ ┗ 📜...   
```

The [scotty\_plot.R](#section4.10.) script takes Scotty's analysis result and create plots which summarises the finding of all the analysis into few plots.


<h3 id="section6.7.">6.7. Custom annotation creation </h3>

The g:Profiler is used to find enriched terms from the list of differential expressed elements and the species under study. Unfortunately, not all species are supported but it is possible to use custom annotation file to overcome such limitation. There are two steps : 

* Find the species' taxon ID by using its scientific name from the [taxonID\_search.R](#section4.11.1.) script.

```
org <- "Streptococcus pneumoniae"
```

* Create the custom annotation by using [gmt\_custom\_annotation\_creation.R](#section4.11.2.) script. The input are a list of taxon IDs (taxonID_list) and a directory to save the output files (outDir).

```
taxonID_list <- c(488223, 189423, 487214)
outDir <- "rootDirectory/../dual_RNAseq/custom_annotation"
```

<h3 id="section6.8.">6.8. Polyester merged count table analysis </h3>

The [de\_countTables\_comparison.R](#section4.13.) script takes a group of count tables and identifies its differential expressed (DE) elements with DESeq2, edgeR, and NOISeq tools at different values of cut-off, replicates, and adjusted p-values. Then, it compares the identified DE elements to a list furnished by the user. The input are:

* count_table_dir -> path to the directory where the count tables are saved.
* patternFile -> pattern used to select the right count tables.
* compare_de -> path to the file where the list of genes/isoforms is saved for the comparison. The header of the columns has the same name of the group comparison as Group1_Group2.
* groupComparison -> list of group comparison written as "Group1_Group2".
* sample_group -> path to the file where the sample's group are saved (column: "Sample" and "Group").
* selectReplicate -> list of values for the number of replicates to use per group during the differential expressed analysis.
* cutOff -> list of cut-off values for filtering rows which has a row-sum higher than the value.
* padjFDRvalue -> list of adjusted p-values.
* log2FC_value -> single value for the Log2 Fold Change.

The example below compare host count tables with isoform reads and compare them with the user's own list of truth differential expressed isoforms. In the case of Polyester simulation is obtained from the "sim_de_tx.txt" file; whereas, for comparison between published paper and the present pipeline we used the author's list of differential expressed elements.

```
count_table_dir <- "rootDirectory/../Author_Year/8_count_table/host/"
patternFile <- "isoform"
compare_de <- "truth_de_isoforms.txt"
groupComparison <- c("Case_Control")
sample_group <- "Polyester_Count_Tables/sample_group.txt"
selectReplicate <- c(2, 10)
cutOff <- c(0, 3000)
padj_FDR_value <- c(0.05, 5e-05)
log2FC_value <- 1
```

<h3 id="section6.9.">6.9. Count tables DE and enrichment comparison analysis </h3>

The [de\_enrichment\_countTables\_comparison.R](#section4.14.) script takes a list of count tables and it identifies the differential expressed terms with DESeq2 for each group comparison. The enrichment analysis is carried out for each group comparison and their similarity is showed in a heatmap for all count tables. The inputs are:

* count_table_dir -> path to the directory where the count tables are present
* patternFile -> Pattern used to select the right count tables.
* sample_group -> Path to the file where the sample's group are saved as two columns called "Sample" and "Group".
* groupComparison -> List of group comparison written as "Group1_Group2".
* cutOff -> Filter for genes/isoforms which has a row-sum higher than the cut-off value.
* padj_FDR_value -> Adjusted p-values used to filter for differential expressed genes/isoforms.
* log2FC_value -> Single value for the Log2 Fold Change used to select differential expressed genes/isoforms.
* gProfiler_supported_org -> List of Scientific Name for g:Profiler supported organisms.
* gProfiler_non_supported_org -> List of path to GMT custom annotation files for non-supported organisms.
* cutoff -> REVIGO's cut-off values (Allowed values: "0.90" "0.70" "0.50" "0.40")
* isPValue -> Is the REVIGO's input numbers p-values? (Allowed values: "yes"  "no")
* whatIsBetter -> In case of some other quantity than the pvalues where the "allowed" value is better. Allowed values: "higher" "lower" "absolute" "abs_log"
* goSizes -> Select a database with GO term sizes (Default is "0" which specifies "whole Uniprot").
* measure -> Select a semantic similarity measure to use. Allowed values: "RESNIK" "LIN" "SIMREL" "JIANG"

```
count_table_dir <- "rootDirectory/../Author_Year/Hsapiens_Spneumoniae_SingleStrand_SE/6_count_table/spneumoniae"
patternFile <- "gene"
sample_group <- "Hsapiens_Spneumoniae_SingleStrand_SE/metadata/sim_rep_info.txt"
groupComparison <- c("1_2")

#DESeq2 input
cutOff <- 500
padj_FDR_value <- 0.0005
log2FC_value <- 1

#g:Profiler input
gProfiler_supported_org <- c("Homo sapiens")
gProfiler_non_supported_org <- c("")

#REVIGO's Input format
cutoff <- "0.40"
isPValue <- "yes"
whatIsBetter <- "higher"
goSizes <- "0"
measure <- "SIMREL"
```



