#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

#BiocManager::install("GenomicFeatures")
#BiocManager::install("GenomicRanges")
#BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
#BiocManager::install("annotatr")
## para sacar secuencia : BSgenome.Hsapiens.UCSC.hg38
## biomart


library('rentrez')
library("GenomicFeatures")
library("GenomicRanges")
library('stringr')
library('annotatr')
library("TxDb.Hsapiens.UCSC.hg38.knownGene")

getCrossRef <- function(ID, dbFrom, dbTo){
	crossi <- entrez_link(dbfrom=dbFrom, id=ID, db=dbTo)
	return(crossi)
}

#Pasarle los MIM numbers de morbidmap.txt (OMIM_data)
getDataFromOMIM <- function(omimID){
	termS <- paste(omimID, '[MIM]', sep= '')
	searchClinVar <- entrez_search(db="clinvar", term=termS, retmax= 150)
	

}

getClinVarData <- function(){


}


cleaningClinVarData <- function(data){
	preData <- data[,-15]
	names(data) <- c('Name', 'Gene(s)', 'Protein_change', 'Condition(s)', 'Clinical_significance_(Last-reviewed)', 'Review_status', 'Accession', 'GRCh37Chromosome', 'GRCh37Location', 'GRCh38Chromosome', 'GRCh38Location', 'VariationID', 'dbSNP_ID', 'Canonical_SPDI')
	clinCat <- data[,'Clinical_significance_(Last-reviewed)']
	cleanClinCat <- as.vector(unlist(sapply(clinCat, gsub, pattern= "\\(.*?\\)", replacement= '')))

}

crossiOC <- entrez_link(dbfrom= 'omim', id= '615787', db= 'clinvar')
crossiCO <- entrez_link(dbfrom= 'clinvar', id= '225696', db= 'omim')

###############################################
#Despues de quitarle primera linea, resultado a "mendelian"
data <- read.table('clinvarData.txt', sep = '\t')

names(data) <- c('Name', 'Gene(s)', 'Protein_change', 'Condition(s)', 'Clinical_significance_(Last-reviewed)', 'Review_status', 'Accession', 'GRCh37Chromosome', 'GRCh37Location', 'GRCh38Chromosome', 'GRCh38Location', 'VariationID', 'dbSNP_ID', 'Canonical_SPDI')

clinCat <- data[,'Clinical_significance_(Last-reviewed)']

#Pathogenic ...
#likely pathogenic
cleanClinCat <- as.vector(unlist(sapply(clinCat, gsub, pattern= "\\(.*?\\)", replacement= '')))

table(cleanClinCat) #??

# no assertion criteria provided 
patosClinCatComple <-  cleanClinCat[grep('Likely pathogenic', cleanClinCat)]
patosClinCat <- cleanClinCat[grep('Pathogenic', cleanClinCat)]

patosClinWhole <- c(patosClinCat, patosClinCatComple)

table(patosClinWhole)

patosData1 <-  data[grep('Pathogenic', cleanClinCat),]
patosData2 <-  data[grep('Likely pathogenic', cleanClinCat),]

patosData <- rbind(patosData1, patosData2)
reviewData <- patosData[,'Review_status']
cleanPatosData <- patosData[-grep('no assertion criteria provided', reviewData),]

varID <- patosData$VariationID
cosito <- patosData[1,]
locationCosito <- as.vector(str_split(cosito$GRCh38Location, pattern= ' - ', simplify=TRUE))

#-----

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

select(txdb, keys = "9636",columns= columns(txdb),  keytype="GENEID")

genesH38 <- genes(txdb)
exonsH38 <- exons(txdb)

genPruebaIndex <- match(x= '9636', table= genesH38$gene_id)

#-----









