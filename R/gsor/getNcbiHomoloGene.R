# 01/18/2018 developed by Yan Li at Center for Research Informatics, University of Chicago.
# 3 functions in this R script
# 1. getNcbiHomoloGene(): obtain the whole NCBI Homologene database
# 2. findHomo(genelist, genelistFormat, from, to, logSave =  FALSE): get homologous genes
# 3. homoMap(geneLists, speciesLists, genelistFormat): map homologous gene summary?
# ----------------
library(clusterProfiler)
library(curl)
library(tidyr)
library(dplyr)
# ----------------
# 1. 'getNcbiHomoloGene': This function script is used to access and get NCBI Homologene database for homologous genes analysis.
# currently this function does not work on gardner
getNcbiHomoloGene <- function(version){
  # ---
  # 1.0. specify downloading version number
  if (missing(version)) version = 'current' else version = paste('build', version, sep = '')
  # 1. specify downloading url path, homologen data, and TaxID/TaxName
  urlPath           = paste('ftp://ftp.ncbi.nlm.nih.gov/pub/HomoloGene', version, sep = '/')
  urlHomologeneData = paste(urlPath, 'homologene.data', sep = '/')
  urlTaxidTaxname   = paste(urlPath, 'build_inputs/taxid_taxname', sep = '/')
  if (version > 65 |  version == 'current') {
    urlReleaseNo    = paste(urlPath, 'RELEASE_NUMBER', sep = '/')
    print(sprintf('The used HomoloGene database is version %s', readLines(urlReleaseNo)))
  } else {
    print(sprintf('The used HomoloGene database is version %s', version))
  }
  # ---
  # 1.2.1 create temperary dir/files names to download 'homologene.data' and 'taxid_taxname'
  tmpDir            <- paste(getwd(), 'download', sep = '/')
  if(!dir.exists(tmpDir)) dir.create(tmpDir)
  tmpHomologeneData <- tempfile(pattern = 'homologene.', tmpdir = tmpDir, fileext = '.data')
  tmpTaxidTaxname   <- tempfile(pattern = 'taxid_taxname.', tmpdir = tmpDir)
  # 1.2.2 downloaed 'homologene.data' and 'taxid_taxname'
  curl_download(url = urlHomologeneData, destfile = tmpHomologeneData)
  curl_download(url = urlTaxidTaxname, destfile = tmpTaxidTaxname)
  # ---
  # 1.3. load 'homologene.data' and 'taxid_taxname' into R
  HomologeneData <- read.delim2(file = tmpHomologeneData, header = F, sep = '\t')
  TaxidTaxname   <- read.delim2(file = tmpTaxidTaxname, header = F, sep = '\t')
  # ---
  # 1.4. process HomologeneData with TaxidTaxname
  colNamesHomologeneData   = c('hgID', 'taxID', 'entrezGeneID', 'geneSymbol', 'protGI', 'protAccession')
  colNamesTaxidTaxname     = c('taxID', 'species')
  colnames(HomologeneData) <- colNamesHomologeneData
  # print( head(HomologeneData) )
  colnames(TaxidTaxname)   <- colNamesTaxidTaxname
  # print( head(TaxidTaxname) )
  # ---
  # 1.5. complete process, remove temperary download directory
  file.remove(c(tmpHomologeneData, tmpTaxidTaxname))
  file.remove(tmpDir)
  # ---
  # 1.6. save 'HomologeneData' & 'TaxidTaxname' into global environment
  HomologeneData <<- HomologeneData
  TaxidTaxname   <<- TaxidTaxname
  # 1.6.1 add short name for TaxidTaxname species
  TaxidTaxname$genome <- sapply(X = as.character(TaxidTaxname$species), FUN = function(x) paste(tolower(substring(strsplit(gsub("\\(|\\)", "", x), split = " ")[[1]], 1, 1)), collapse = '') )
  TaxidTaxname <- TaxidTaxname %>% dplyr::mutate(genome = replace(genome, species == 'Macaca mulatta', 'mam'))
  TaxidTaxname <<- TaxidTaxname
}
# ----------------
# 2. 'findHomo': This function script is used to identify homologous genes from one species to another species.
findHomo <- function(genelist, genelistFormat, from, to, logSave =  FALSE) {
  # # flags test for function
  # # genelist <- worm_2venn_degs$`C.elegans-CCM3_C.elegans-Kri1`
  # genelist <- mouse_2venn_degs$`Mouse-CCM3_Mouse-Kri1`
  # genelistFormat = 'geneSymbol'
  # # from           = 'ce'
  # from           = 'mm'
  # to             = 'hs'
  # ---
  # 2.0 set default genelistFormat into GeneSymbol, and convert them to lowercase
  if (missing(genelistFormat)) genelistFormat = 'genesymbol'
  genelistFormat <- tolower(as.character(genelistFormat))
  # ---
  # 2.1 convert input genelist into data frame with entrez id and print out how many genes lost due to no corresponding entrezID
  if (genelistFormat=='genesymbol') {
    entrezGenelists <- bitr(geneID = genelist, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = as.character(get.orgdb(ref.genome = as.character(from))))
  }else if (genelistFormat=='entrezid') {
    entrezGenelists <- genelist
  }
  noGenesMissingEntrezid <- length(genelist) - dim(entrezGenelists)[1]
  print(sprintf('out of %s input genelists, %s (%s%%) input genes has no corresponding entrezID, and the rest %s (%s%%) has corresponding entrezID', 
                length(genelist), noGenesMissingEntrezid, round(noGenesMissingEntrezid*100/length(genelist), digits = 2), 
                dim(entrezGenelists)[1], round(dim(entrezGenelists)[1]*100/length(genelist), digits = 2) ))
  # ---
  # 2.2 find out the corresponding taxid of 2 input 'from' and 'to' species for further homologous genes serch 
  fromTaxid     <- as.character(TaxidTaxname$taxID[match(as.character(from), TaxidTaxname$genome)])
  toTaxid       <- as.character(TaxidTaxname$taxID[match(as.character(to), TaxidTaxname$genome)])
  # ---
  # 2.3 subset/filter the whole 'HomologeneData' based on input 'from' and 'to' species
  fromHomologenData <- HomologeneData %>% filter(taxID == fromTaxid) %>% mutate_at('entrezGeneID', funs(as.character))  
  toHomologenData   <- HomologeneData %>% filter(taxID == toTaxid) %>% mutate_at('entrezGeneID', funs(as.character))    
  print(sprintf("%s genes out of %s input genelists with corresponding entrezID can be found in the HomologGene database based on input 'from' species: %s", 
                dim(fromHomologenData)[1], dim(entrezGenelists)[1], 
                as.character(dplyr::filter(TaxidTaxname, genome == as.character(from))$species)))
  # ---
  # 2.4 identify 1) subsetted corresponding 'HomologeneData' based on input 'genelist' and 'from' species, 
  #              2) remove no corresponding data in 'HomologeneData' from input 'genelist' and 'from' species,
  #              3) sort the subseted 'HomologeneData' with corresponding dta in 'HomologeneData' by 'hgID' 
  #              4) print out gene no. with corresponding data in 'HomologeneData' based on input 'genelist'
  #     Note, right side y(ENTREZID) is from the input genelists unique entrezID
  fromGenelistHomologenData     <- right_join(x = fromHomologenData, y = entrezGenelists, by = c("entrezGeneID" = "ENTREZID") ) 
  fromGenelistHomologenDataRmNa <- fromGenelistHomologenData %>% 
    filter(!is.na(taxID)) %>%
    arrange(hgID)  
  
  # print('---')
  print(sprintf("%s (%s%%) input genes out of %s input genelists from %s species with entrezID in NCBI HomologeneData database", 
                dim(fromGenelistHomologenDataRmNa)[1], 
                (round(dim(fromGenelistHomologenDataRmNa)[1]*100/dim(fromHomologenData)[1], digits = 2)), 
                dim(fromHomologenData)[1], 
                as.character(dplyr::filter(TaxidTaxname, genome == as.character(from))$species) ))
  # print('---')
  # ---
  # 2.5 identify 1) matched homologens from x to y based on hgID
  #              2) sorted them by gene symbol based on 'from' input 'genelist' alphabatically
  #              3) print out no. of genes from input 'genelist' has homologous genes from 'from' to 'to' species
  #     Note: 1) if one gene from input genelist has more than 1 homologeous genes in y, only the first match will be kept with distinct()
  #           2) not every gene from input list has a corresponding homologous genes in y, as shown in below with 'HomologenMatchPrep3'&'HomologenMatch3'
  # HomologenMatchPrep3   <- dplyr::left_join(x = fromGenelistHomologenDataRmNa, y = toHomologenData, by = c("hgID", "hgID") )
  # HomologenMatch3       <- distinct(HomologenMatchPrep3, SYMBOL, .keep_all = T) %>% na.omit() %>% arrange(desc(-hgID))
  HomologenMatchPrep   <- dplyr::inner_join(x = fromGenelistHomologenDataRmNa, y = toHomologenData, by = c("hgID", "hgID") ) 
  HomologenMatch       <- distinct(HomologenMatchPrep, SYMBOL, .keep_all = T) %>% arrange(desc(-hgID))
  
  HomologenMatchReturn <- dplyr::select(HomologenMatch, c('SYMBOL', 'geneSymbol.y'))
  colnames(HomologenMatchReturn) <- c('geneSymbol_from', 'geneSymbol_to')
  if (logSave) save(HomologenMatch, file = paste(getwd(), 'NcbiHomoloGeneSearchResults.Rdata', sep = '/'))
  print(sprintf("out of %s input genelist from %s species from NCBI HoologeneData database, %s (%s%%) has corresonding homologous genes in %s species, %s input genes has %s (%s%%) homologous genes in %s species", 
                dim(fromGenelistHomologenDataRmNa)[1], as.character(dplyr::filter(TaxidTaxname, genome == as.character(from))$species), 
                dim(HomologenMatch)[1], 
                round(dim(HomologenMatch)[1]*100/dim(fromGenelistHomologenDataRmNa)[1]), 
                dplyr::filter(TaxidTaxname, genome == as.character(to))$species,
                length(unique(HomologenMatch$entrezGeneID.x)),
                length(unique(HomologenMatch$entrezGeneID.y)), 
                round(length(unique(HomologenMatch$entrezGeneID.y))*100/dim(HomologenMatch)[1]), 
                dplyr::filter(TaxidTaxname, genome == as.character(to))$species))
  return(HomologenMatchReturn)
}
# ----------------
# 3. 'homoVenn': 
homoVenn <- function(speciesLists){
  HomologenDataSpecies <- lapply(speciesLists, function(x) as.data.frame(HomologeneData %>% filter(taxID == filter(TaxidTaxname, genome == x)$taxID) %>% mutate_at('entrezGeneID', funs(as.character))))
  names(HomologenDataSpecies) <- speciesLists
  HomoGeneSpecies  <- lapply(HomologenDataSpecies, function(x) unique(x$hgID))
  library(gplots)
  venn(HomoGeneSpecies)
}
# ----------------
# 3. 'homoMap': 
homoMap <- function(geneLists, speciesLists, genelistFormat){
  if (missing(genelistFormat)) genelistFormat = 'GeneSymbol'
  genelistFormat <- tolower(genelistFormat)
  # speciesLists must be species abbrevation of provided geneLists 
  names(geneLists) <- speciesLists
  # homoloGene databse species specified data
  HomologenDataSpecies <- lapply(speciesLists, function(x) as.data.frame(HomologeneData %>% filter(taxID == filter(TaxidTaxname, genome == x)$taxID) %>% mutate_at('entrezGeneID', funs(as.character))))
  names(HomologenDataSpecies) <- speciesLists
  
  # correpsonding entrezID of input geneLists
  if (genelistFormat=='genesymbol') {
    entrezGenelists <- lapply(1:length(geneLists), function(x) bitr(geneID = as.character(geneLists[[x]]), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = as.character(get.orgdb(ref.genome = as.character(names(geneLists)[x])))))
  }else if (genelistFormat=='entrezid') {
    entrezGenelists <- geneLists
  }
  # cooresponding homoloGene data of input geneLists with corresponding entrezID
  GenelistsHomologenData     <- lapply(1:length(geneLists), function(x) right_join(x = HomologenDataSpecies[[x]], y = entrezGenelists[[x]], by = c("entrezGeneID" = "ENTREZID") ))
  GenelistsHomologenDataRmNa <- lapply(1:length(geneLists), function(x) GenelistsHomologenData[[x]] %>% filter(!is.na(taxID)) %>% arrange(hgID))
  names(GenelistsHomologenDataRmNa) <- speciesLists
  # Summary
  resSummary <- data.frame('inputGenesNo'          = sapply(geneLists, length), 
                           'inputGenesWentrezIDNo' = sapply(entrezGenelists, function(x) length(x$ENTREZID)),
                           'speciesHomoloGeneDBNo' = sapply(HomologenDataSpecies, function(x) dim(x)[1]), 
                           'uniqEntrezGeneIDNo'    = sapply(GenelistsHomologenDataRmNa, function(x) length(unique(x$entrezGeneID))), 
                           'UniqHomoloGeneNo'      = sapply(GenelistsHomologenDataRmNa, function(x) length(unique(x$hgID))))
  return(list(summary = resSummary, HomologenDataSpecies = HomologenDataSpecies, mappedHomologenes =  GenelistsHomologenDataRmNa, speciesLists = speciesLists, homoloSpeciesDb = HomologenDataSpecies))
}
