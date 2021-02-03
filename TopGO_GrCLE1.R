#GO analyses using TopGO in GrCLE1-1 unique genes with cluster 1 DEGs from 1289 DEGs as example
# Also functional enrichment analysis using TopGO.
# TopGO does not provide clusters, and therefore the functional network is built using only the gene-term sets

setwd("/Users/ma2292/Documents/cornell/WangLab/CLE_project/GO_analyses")
if (!requireNamespace("BiocManager", quietly=TRUE))
  + install.packages("BiocManager")
BiocManager::install()
if (!requireNamespace("BiocManager", quietly=TRUE))
  + install.packages("BiocManager")
BiocManager::install("topGO")

#installing biomaRt for ID conversion and GO ID
library("topGO")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")
library(biomaRt)
ensembl <- useMart("plants_mart", dataset = "stuberosum_eg_gene", host = 'plants.ensembl.org')
genesToGO <- biomaRt::getBM(attributes = c("ensembl_gene_id", "go_id"), mart = ensembl)
#convert from table format to list format
genesToGO <- by(genesToGO$go_id, genesToGO$ensembl_gene_id, function(x) as.character(x))
head(genesToGO)
#remove blank entries
genesToGO <- genesToGO[genesToGO$go_id !="",]
attributes <- listAttributes(ensembl) #check attributes for the next step and filter
features <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", "ensembl_peptide_id", "ensembl_exon_id", "description", "chromosome_name"), filters = "chromosome_name", values = 1, mart = ensembl)
all.genes <- sort(unique(as.character(genesToGO$ensembl_gene_id))) #select only unique genes from all the total genes or the gene universe
int.genes <- read.table("cluster1_uniqueGrCLE1.txt", sep="\t", header=T)
int.genes <- as.character(int.genes$GeneID)
int.genes <- factor(as.integer(all.genes %in% int.genes))
names(int.genes) <- all.genes
go.obj <- new("topGOdata", description = "Cluster 1 GrCLE1-1 uniqueDEGs", ontology="BP", allGenes = int.genes, annot = annFUN.gene2GO, gene2GO = genesToGO)
go.obj # check your GO results
sg <- sigGenes(go.obj)
str(sg) #see internal structure
numSigGenes(go.obj) #significant no. of genes 

# perform enrichment testing
resultFisher <- runTest(go.obj, algorithm = "weight01", statistic = "fisher")
resultsKS.elim <- runTest(go.obj, algorithm = "elim", statistic = "ks")
# for visualisation
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Rgraphviz")

showSigOfNodes(go.obj, score(resultFisher), firstSigNodes = 5, useInfo = "all")
printGraph(go.obj, resultFisher, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)




