#BiocManager::install("org.Hs.eg.db")

library('rentrez')
library("xml2")
library('XML')
library('methods')
library("GenomicFeatures")
library("GenomicRanges")
library('stringr')
library('annotatr')
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library('org.Hs.eg.db')
library(biomaRt)


######################################################### FUNCIONES #########################################################
#58067

getLocationWithFormat <- function(startS, stopS){
	diffe <- as.numeric(startS) - as.numeric(stopS)
	if (diffe == 0){
		return(startS)
	}
	else{
		locationF <- paste(startS, '-', stopS)
		return(locationF)
	}
}

##
getSnpDataFromClinvarID <- function(clinvarID){
	cross <- entrez_link(dbfrom='clinvar', id=clinvarID, db='snp')
	if (length(cross$link) == 0){
		return(NA)}
	else{
		crossID <- cross$links$clinvar_snp
		SNPsummary <- entrez_summary(db= 'snp', id= crossID)
		snp_id <- paste('rs', SNPsummary$snp_id, sep='')
		return(snp_id)}
}

getConcaCats <- function(veve){
	if (length(veve) > 1){
		vevo <- str_replace_all(string= paste(veve, collapse= '\t'), pattern= '\t', replacement= '|')
		return(vevo)}
	else{
		return(veve)}
}

#entrez_link(dbfrom='clinvar', id=clinvarID, db='all')
#fetcheo <- entrez_fetch(db= 'clinvar', id= clinvarID, rettype= 'docsum', retmode= 'xml')
#library('XML')
#library('methods')
#xmlParse(fetcheo)

getClinVarDataByID <- function(clinvarID){
	message(clinvarID)
	esummary <- entrez_summary(db="clinvar", id=clinvarID)
	genes <- esummary$genes
	symbolGenes <- getConcaCats(genes$symbol)
	locationData <- esummary$variation_set$variation_loc[[1]]
	grCh37Data <- locationData[locationData$assembly_name == 'GRCh37',]
	locaCh37 <- getLocationWithFormat(grCh37Data$start, grCh37Data$stop)
	if ('GRCh38' %in% locationData$assembly_name){
		grCh38Data <- locationData[locationData$assembly_name == 'GRCh38',]
		locaCh38 <- getLocationWithFormat(grCh38Data$start, grCh38Data$stop)}
	else {
		grCh38Data <- NA
		locaCh38 <- NA
	}
	dbSNP_ID <- getSnpDataFromClinvarID(clinvarID)
	conditions <- getConcaCats(esummary$trait_set$trait_name)
	alleleIDs <- getConcaCats(esummary$variation_set$measure_id)
	rowData <- data.frame(Name = esummary$title, Genes = symbolGenes, Protein_change = esummary$protein_change, Condition = conditions, Clinical_significance = esummary$clinical_significance$description, Review_status = esummary$clinical_significance$review_status, Accession = esummary$accession, GRCh37Chromosome= grCh37Data$chr, GRCh37Location= locaCh37, GRCh38Chromosome = grCh38Data$chr, GRCh38Location = locaCh38, VariationID = clinvarID, AlleleID=esummary$variation_set$measure_id, dbSNP_ID = dbSNP_ID, Canonical_SPDI= esummary$variation_set$canonical_spdi)
	return(rowData)
}

# Data de los genes ligados a el ID [symbol	geneid	strand	source]
getGenDataByClinvarID <- function(clinvarID){
	esummary <- entrez_summary(db="clinvar", id=clinvarID)
	return(esummary$genes)
}

#Pasarle los MIM numbers de morbidmap.txt (OMIM_data)
getDataFromOMIM <- function(omimID, retMax){
	termS <- paste(omimID, '[MIM]', sep= '')
	searchClinVar <- entrez_search(db="clinvar", term=termS, retmax= retMax)
	searchIds <- searchClinVar$ids
	dataFromOMIM <- data.frame() #??
	for (id in searchIds){
		tempData <- getClinVarDataByID(id)
		dataFromOMIM <- rbind(dataFromOMIM, tempData)}
	dirtyDataFromOMIM<- dataFromOMIM[rowSums(is.na(dataFromOMIM)) > 1, ]
	indexes <- as.numeric(row.names(dirtyDataFromOMIM))
	cleanDataFromOMIM <- dataFromOMIM[-indexes,]
	return(cleanDataFromOMIM)
}


# Limpia y filtra
cleaningClinVarData <- function(data, type){
	names(data) <- c('Name', 'Gene(s)', 'Protein_change', 'Condition(s)', 'Clinical_significance', 'Review_status', 'Accession', 'GRCh37Chromosome', 'GRCh37Location', 'GRCh38Chromosome', 'GRCh38Location', 'VariationID','AlleleID(s)', 'dbSNP_ID', 'Canonical_SPDI')
	clinCat <- data[,'Clinical_significance']
	if (type == 'file'){
		cleanClinCat <- as.vector(unlist(sapply(clinCat, gsub, pattern= "\\(.*?\\)", replacement= '')))
		data[,'Clinical_significance'] <- cleanClinCat
	}
	else{
		cleanClinCat <- clinCat}
	patosClinCatComple <-  cleanClinCat[grep('Likely pathogenic', cleanClinCat)]
	patosClinCat <- cleanClinCat[grep('Pathogenic', cleanClinCat)]
	patosClinWhole <- c(patosClinCat, patosClinCatComple)
	patosData1 <-  data[grep('Pathogenic', cleanClinCat),]
	patosData2 <-  data[grep('Likely pathogenic', cleanClinCat),]
	patosData <- rbind(patosData1, patosData2)
	cleanPatosData <- patosData[patosData$Review_status != 'no assertion criteria provided',]
	return(cleanPatosData)
}

crossiOC <- entrez_link(dbfrom= 'omim', id= '615787', db= 'clinvar')
crossiCO <- entrez_link(dbfrom= 'clinvar', id= '225696', db= 'omim')

######################################################## C. Principal ########################################################
#Despues de quitarle primera linea, resultado al query "mendelian"
data <- read.table('clinvarData.txt', sep = '\t')

#Limpieza y filtrado de datos (Clinical_significance y Review_status)
cleanPatosData <- cleaningClinVarData(data, type= 'file') 




# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

varID <- patosData$VariationID
cosito <- patosData[1,]
locationCosito <- as.vector(str_split(cosito$GRCh38Location, pattern= ' - ', simplify=TRUE))

#-----

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

select(txdb, keys = "9636",columns= columns(txdb),  keytype="GENEID")

genesH38 <- genes(txdb)
exonsH38 <- exons(txdb)

genPruebaIndex <- match(x= '9636', table= genesH38$gene_id)


#------
example <- cleanPatosData[10,]

locationExample <- as.numeric(as.vector(str_split(example$GRCh38Location, pattern= ' - ', simplify=TRUE)))

gr <- GRanges(
    seqnames = paste('chr',  example$GRCh38Chromosome, sep= ''),
    ranges = IRanges(start= locationExample[1], end= locationExample[2]),
    strand = '*',
    geneID = example[, 'Gene(s)'], dbSNPID= example$dbSNP_ID)

gogo <- GRanges(seqnames= 'chr1', ranges= IRanges(start= 10869, end= 11868), strand= '+')

annots = c( 'hg38_basicgenes', 'hg38_genes_intergenic', 'hg38_genes_intronexonboundaries')

annotations = build_annotations(genome = 'hg38', annotations = annots)

bolillo = annotate_regions(
    regions = gr,
    annotations = annotations,
    ignore.strand = TRUE,
    quiet = FALSE)

telera = annotate_regions(
    regions = gogo,
    annotations = annotations,
    ignore.strand = TRUE,
    quiet = FALSE)
teleraDf = data.frame(telera)

doggo <- findOverlaps(exonsH38, gr)
doggo <- findOverlaps(genesH38, gr)