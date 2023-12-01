#This analysis is based on clusterprofiler(CP).
## CP includes 3 major functions 1). enrichGo() for GO over-representation test (hypergeometric test);
## 2). groupGO() for groupGo (for gene clssicification);
## 3). gseGO() for GSEA, which use permutation test, user can set nPerm for number of permutations, 
## only gene set size in [minGSSize, maxGSSize] will be tested
## GO-terms Semnatic Similarity Measures: https://bioconductor.org/packages/release/bioc/html/GOSemSim.html
##------------------
## In this script, we are only focous on the enrichGO()
## 3 steps involved into this analysis
## step1: getGeneListEntrezID() in getEntrenzID.R - convert gene list into entrezID based on their corresponding reference genome
## step2: getOrgDB() in getOrgDB.R - load the corresponding GO db
## step3: 1) enrichGO() from clusterProfiler to conduct GO over-representation test (hypergeometrix test) 
##        it is enbedded in cp.enrich.go() in getEnrichGoRes.R with below 6 flags:
##        gsorInputGenes for step1, gsorInputGenesType for step1, refGenome for step2, GOpool for step3, resFnamePrefix for step3
##        if input is a list, then it will be clustered results with column 'Cluster', if input is a single gene list, the result is only enrichment results.
##        if GOpool is 'T, GO enrichment results has column 'ONTOLOGY', otherwise not.
## step4: GO analysis results visulization (dot plot): flag5&6 resFnamePrefix
##------------------
##start GO over-representation test for 3 different GO categories (BP, MF, CC)
# source('R/getEntrenzID.R') below 2 functions
#----------
getGeneListEntrezID <- function(inputGeneList, geneListsType, refGenome) {
  library(clusterProfiler)
  if (geneListsType == 'GeneSymbol') {
    geneIDres <- bitr(inputGeneList, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = as.character(getOrgDB(refGenome = as.character(refGenome))))
    # geneIDres <- bitr(inputGeneList, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
    geneListEntrezID <- geneIDres$ENTREZID
  } else if (geneListsType == 'ensembleID') {
    geneIDres <- bitr(inputGeneList, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = as.character(getOrgDB(refGenome = as.character(refGenome))))
    # geneIDres <- bitr(inputGeneList, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
    geneListEntrezID <- geneIDres$ENTREZID
  } else if (geneListsType == 'entrezID') {
    geneListEntrezID <- inputGeneList
  }
  return(geneListEntrezID)
}
#----------
##opposite from entrezId to geneSymbol
#----------
getGeneListGeneSymbol <- function(inputGeneList, geneListsType, refGenome) {
  library(clusterProfiler)
  if (geneListsType == 'entrezID') {
    geneSymbolRes <- bitr(inputGeneList, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = as.character(getOrgDB(refGenome = as.character(refGenome))))
    geneListGeneSymbol <- geneSymbolRes$SYMBOL
  } else if (geneListsType == 'ensembleID') {
    geneSymbolRes <- bitr(inputGeneList, fromType = "ENTREZID", toType = "ENSEMBL", OrgDb = as.character(getOrgDB(refGenome = as.character(refGenome))))
    geneListGeneSymbol <- geneSymbolRes$SYMBOL
  } else if (geneListsType == 'GeneSymbol') {
    geneListGeneSymbol <- inputGeneList
  }
  return(geneListGeneSymbol)
}
#----------
getOrgDB <- function(refGenome) {
  if (refGenome == 'human' | refGenome == 'hs') {
    library(org.Hs.eg.db)
    orgdb <- "org.Hs.eg.db"
  } else if (refGenome == 'mouse' | refGenome == 'mm') {
    library(org.Mm.eg.db)
    orgdb <- "org.Mm.eg.db"
  } else if (refGenome == 'rat' | refGenome == 'rn') {
    library(org.Rn.eg.db)
    orgdb <- "org.Rn.eg.db"
  } else if (refGenome == 'worm' | refGenome == 'ce') {
    library(org.Ce.eg.db)
    orgdb <- "org.Ce.eg.db"
  }
  return(orgdb)
}
#---------
getReactomeOrg <- function(refGenome) {
  if (refGenome == 'human' | refGenome == 'hs' ) {
    reactomeOrg <- 'human'
  } else if (refGenome == 'mouse' | refGenome == 'mm') {
    reactomeOrg <- 'mouse'
  } else if (refGenome == 'rat' | refGenome == 'rn') {
    reactomeOrg <- 'rat'
  } else if (refGenome == 'worm' | refGenome == 'ce') {
    reactomeOrg <- 'celegans'
  }
  return(reactomeOrg)
} 
## ---
getKeggOrg <- function(refGenome) {
  if (refGenome == 'human' | refGenome == 'hs' ) {
    keggOrg <- 'hsa'
  } else if (refGenome == 'mouse' | refGenome == 'mm') {
    keggOrg <- 'mmu'
  } else if (refGenome == 'rat' | refGenome == 'rn') {
    keggOrg <- 'rno'
  } else if (refGenome == 'worm' | refGenome == 'ce') {
    keggOrg <- 'cel'
  }
  return(keggOrg)
} 
## ---
get.t2g <- function(refGenome, category, subcategory=NULL) {
  print(sprintf("GSOR conducted on BROAD MSigDB category: %s.", category))
  library(msigdbr)
  if (refGenome == 'human' | refGenome == 'hs' ) {
    if (!is.null(subcategory)){
      m_t2g       <- msigdbr(species = "Homo sapiens", category = category, subcategory = subcategory) %>% dplyr::select(gs_name, entrez_gene) %>% as.data.frame()
    } else {
      m_t2g       <- msigdbr(species = "Homo sapiens", category = category) %>% dplyr::select(gs_name, entrez_gene) %>% as.data.frame()
    }
    
  } else if (refGenome == 'mouse' | refGenome == 'mm' ) {
    if (!is.null(subcategory)){
      m_t2g       <- msigdbr(species = "Mus musculus", category = category, subcategory = subcategory) %>% dplyr::select(gs_name, entrez_gene) %>% as.data.frame()
    } else {
      m_t2g       <- msigdbr(species = "Mus musculus", category = category) %>% dplyr::select(gs_name, entrez_gene) %>% as.data.frame()
    }
  } else {
    stop("HALLMARK analysis only can be conducted with human or mouse gene lists input.")
  }
  return(m_t2g)
}


## ---
cpORTest <- function(gsorInputGenes, gsorInputGenesType, refGenome, functionDB) {
  print('*********')
  # #step 1: get gene list entrezID 
  print(sprintf('Start to transform input %s genelist (%s) into entrez gene id', refGenome, gsorInputGenesType))
  ##pathwayGeneEntrezidInput is the result from this step, input for step3
  if (is.list(gsorInputGenes)) {
    ## get entrez ID in list based on input gene lists
    pathwayGeneEntrezidInput <- list()
    for (l in 1:length(gsorInputGenes)) {
      pathwayGeneEntrezidInput[[l]] <- getGeneListEntrezID(inputGeneList = gsorInputGenes[[l]], geneListsType = geneListsType, refGenome = refGenome)
    }
    if (is.null(names(gsorInputGenes))) {
      names(pathwayGeneEntrezidInput) <- paste('X', 1:length(gsorInputGenes), sep = '')
    }else {
      names(pathwayGeneEntrezidInput) <- names(gsorInputGenes)
    }
    ## end loop over through the gene lists to get entrez id
  } else {
    pathwayGeneEntrezidInput <- getGeneListEntrezID(gsorInputGenes, geneListsType, refGenome)
  }
  print('End to transform input genelist into entrez gene id')
  print('=-=-=-')
  # ---
  ##step 2.1: get enrichGO() database
  gsorDb                       <- getOrgDB(refGenome = as.character(refGenome))
  ##step 2.2: get enrichMKEGG(), enrichKEGG() corresponding KEGG supported organism listed in 'http://www.genome.jp/kegg/catalog/org_list.html'
  gsorKeggOrg                  <- getKeggOrg(refGenome = as.character(refGenome))
  ##step 2.3: get enrichPathway() corresponding orgranis
  reactomeOrg                  <- getReactomeOrg(refGenome = as.character(refGenome))
  # ---
  # #step 3: conduct gsor analysis based on different functional DB
  if ( functionDB == 'go' ) {
    # ---
    print('===')
    print(sprintf('START over representation for pooled GeneOntologies, including all 3 categories.'))
    if (is.list(gsorInputGenes)) {
      GOres                      <- compareCluster(geneCluster   = pathwayGeneEntrezidInput, 
                                                   fun           = 'enrichGO',
                                                   OrgDb         = gsorDb,
                                                   pvalueCutoff  = 1,
                                                   qvalueCutoff  = 1,
                                                   pAdjustMethod = 'BH',
                                                   ont           = 'ALL',
                                                   pool          = T,
                                                   readable      = T)
      ##GOres is a object with 4 slots, reformat the slot 'compareClusterResult' to have only qvalue>0.2
      orRes4plot                      <- GOres
      orRes4plot@compareClusterResult <- subset(orRes4plot@compareClusterResult, p.adjust <= 0.25)
      gsorCPres                       <- as.data.frame(GOres@compareClusterResult) 
      gsorCPres$geneSymbol            = gsorCPres$geneID
    } else {
      GOres                      <- enrichGO(gene          = pathwayGeneEntrezidInput, 
                                             OrgDb         = gsorDb, 
                                             # keytype       = "ENTREZID",
                                             ont           = 'ALL',
                                             pvalueCutoff  = 1, 
                                             qvalueCutoff  = 1, 
                                             pAdjustMethod = 'BH',
                                             pool          = T,
                                             readable      = T)
      orRes4plot                 <- GOres
      orRes4plot                 <- subset(orRes4plot, p.adjust <= 0.25)
      gsorCPres                  <- as.data.frame(GOres@result) # without @results slot also can directly get results, but other slot can check the analysis options if needed
      gsorCPres$geneSymbol       = gsorCPres$geneID
    }
    print(sprintf('END over representation for pooled GeneOntologies, including all 3 categories.'))
    print('***')
    # ---
  } else if ( functionDB == 'gobp' ) {
    # ---
    print('===')
    print(sprintf('START over representation for pooled GeneOntologies BP.'))
    if (is.list(gsorInputGenes)) {
      gobpRes                         <- compareCluster(geneCluster   = pathwayGeneEntrezidInput, 
                                                        fun           = 'enrichGO',
                                                        OrgDb         = gsorDb,
                                                        pvalueCutoff  = 1,
                                                        qvalueCutoff  = 1,
                                                        pAdjustMethod = 'BH',
                                                        ont           = 'BP',
                                                        pool          = T,
                                                        readable      = T)
      ##gobpRes is a object with 4 slots, reformat the slot 'compareClusterResult' to have only qvalue>0.2
      orRes4plot                      <- gobpRes
      orRes4plot@compareClusterResult <- subset(orRes4plot@compareClusterResult, p.adjust <= 0.25)
      gsorCPres                       <- as.data.frame(gobpRes@compareClusterResult) 
      gsorCPres$geneSymbol            = gsorCPres$geneID
    } else {
      gobpRes                    <- enrichGO(gene          = pathwayGeneEntrezidInput, 
                                             OrgDb         = gsorDb, 
                                             # keytype       = "ENTREZID",
                                             ont           = 'BP',
                                             pvalueCutoff  = 1, 
                                             qvalueCutoff  = 1, 
                                             pAdjustMethod = 'BH',
                                             pool          = T,
                                             readable      = T)
      orRes4plot                 <- gobpRes
      orRes4plot                 <- subset(gobpRes, p.adjust <= 0.25)
      gsorCPres                  <- as.data.frame(gobpRes@result) # without @results slot also can directly get results, but other slot can check the analysis options if needed
      gsorCPres$geneSymbol       = gsorCPres$geneID
    }
    print(sprintf('END over representation for pooled GeneOntologies BP.'))
    print('***')
    # ---
  } else if ( functionDB == 'gocc' ) {
    # ---
    print('===')
    print(sprintf('START over representation for pooled GeneOntologies CC.'))
    if (is.list(gsorInputGenes)) {
      goccRes                         <- compareCluster(geneCluster   = pathwayGeneEntrezidInput, 
                                                        fun           = 'enrichGO',
                                                        OrgDb         = gsorDb,
                                                        pvalueCutoff  = 1,
                                                        qvalueCutoff  = 1,
                                                        pAdjustMethod = 'BH',
                                                        ont           = 'CC',
                                                        pool          = T,
                                                        readable      = T)
      ##goccRes is a object with 4 slots, reformat the slot 'compareClusterResult' to have only qvalue>0.2
      orRes4plot                      <- goccRes
      orRes4plot@compareClusterResult <- subset(orRes4plot@compareClusterResult, p.adjust <= 0.25)
      gsorCPres                       <- as.data.frame(goccRes@compareClusterResult) 
      gsorCPres$geneSymbol            = gsorCPres$geneID
    } else {
      goccRes                    <- enrichGO(gene          = pathwayGeneEntrezidInput, 
                                             OrgDb         = gsorDb, 
                                             # keytype       = "ENTREZID",
                                             ont           = 'CC',
                                             pvalueCutoff  = 1, 
                                             qvalueCutoff  = 1, 
                                             pAdjustMethod = 'BH',
                                             pool          = T,
                                             readable      = T)
      orRes4plot                 <- goccRes
      orRes4plot                 <- subset(goccRes, p.adjust <= 0.25)
      gsorCPres                  <- as.data.frame(goccRes@result) # without @results slot also can directly get results, but other slot can check the analysis options if needed
      gsorCPres$geneSymbol       = gsorCPres$geneID
    }
    print(sprintf('END over representation for pooled GeneOntologies CC.'))
    print('***')
    # ---
  } else if ( functionDB == 'gomf' ) {
    # ---
    print('===')
    print(sprintf('START over representation for pooled GeneOntologies MF.'))
    if (is.list(gsorInputGenes)) {
      gomfRes                         <- compareCluster(geneCluster   = pathwayGeneEntrezidInput, 
                                                        fun           = 'enrichGO',
                                                        OrgDb         = gsorDb,
                                                        pvalueCutoff  = 1,
                                                        qvalueCutoff  = 1,
                                                        pAdjustMethod = 'BH',
                                                        ont           = 'MF',
                                                        pool          = T,
                                                        readable      = T)
      ##gomfRes is a object with 4 slots, reformat the slot 'compareClusterResult' to have only qvalue>0.2
      orRes4plot                      <- gomfRes
      orRes4plot@compareClusterResult <- subset(orRes4plot@compareClusterResult, p.adjust <= 0.25)
      gsorCPres                       <- as.data.frame(gomfRes@compareClusterResult) 
      gsorCPres$geneSymbol            = gsorCPres$geneID
    } else {
      gomfRes                    <- enrichGO(gene          = pathwayGeneEntrezidInput, 
                                             OrgDb         = gsorDb, 
                                             # keytype       = "ENTREZID",
                                             ont           = 'MF',
                                             pvalueCutoff  = 1, 
                                             qvalueCutoff  = 1, 
                                             pAdjustMethod = 'BH',
                                             pool          = T,
                                             readable      = T)
      orRes4plot                 <- gomfRes
      orRes4plot                 <- subset(gomfRes, p.adjust <= 0.25)
      gsorCPres                  <- as.data.frame(gomfRes@result) # without @results slot also can directly get results, but other slot can check the analysis options if needed
      gsorCPres$geneSymbol       = gsorCPres$geneID
    }
    print(sprintf('END over representation for pooled GeneOntologies MF.'))
    print('***')
    # ---
  } else if ( functionDB == 'kegg' ) {
    # ---
    print('===')
    print(sprintf('START over representation for kegg DB.'))
    print(sprintf("KEGG DB organism is %s", gsorKeggOrg))
    if (is.list(gsorInputGenes)) {
      keggORres                       <- compareCluster(geneCluster   = pathwayGeneEntrezidInput, 
                                                        fun           = 'enrichKEGG',
                                                        organism = gsorKeggOrg, keyType = 'kegg', 
                                                        minGSSize = 10, maxGSSize = 500, 
                                                        pvalueCutoff  = 1, 
                                                        qvalueCutoff  = 1, 
                                                        pAdjustMethod = 'BH')
      orRes4plot                      <- keggORres
      orRes4plot@compareClusterResult <- subset(orRes4plot@compareClusterResult, p.adjust <= 0.25)
      gsorCPres                       <- as.data.frame(keggORres@compareClusterResult) 
    } else {
      keggORres                 <- enrichKEGG(gene = pathwayGeneEntrezidInput, 
                                              organism = gsorKeggOrg, keyType = 'kegg', 
                                              minGSSize = 10, maxGSSize = 500, 
                                              pvalueCutoff  = 1, 
                                              qvalueCutoff  = 1, 
                                              pAdjustMethod = 'BH')
      orRes4plot                 <- keggORres
      orRes4plot                 <- subset(keggORres, p.adjust <= 0.25)
      gsorCPres                  <- keggORres@result
    }
    for (r in 1:dim(gsorCPres)[1]) {
      gsorCPres$geneSymbol[r] <- paste(getGeneListGeneSymbol(inputGeneList = strsplit(gsorCPres$geneID[r], split = '/')[[1]], geneListsType = 'entrezID', refGenome = refGenome), collapse = '/')
    }
    print(sprintf('END over representation for kegg DB.'))
    print('***')
    # ---
  } else if ( functionDB == 'mkegg' ) {
    # ---
    print('===')
    print(sprintf('START over representation for kegg module categories DB.'))
    if (is.list(gsorInputGenes)) {
      mkeggORres                      <- compareCluster(geneCluster   = pathwayGeneEntrezidInput, 
                                                        fun           = 'enrichMKEGG',
                                                        organism = gsorKeggOrg, keyType = 'kegg', 
                                                        minGSSize = 1, maxGSSize = 500, 
                                                        pvalueCutoff  = 1, 
                                                        qvalueCutoff  = 1, 
                                                        pAdjustMethod = 'BH')
      orRes4plot                      <- mkeggORres
      orRes4plot@compareClusterResult <- subset(orRes4plot@compareClusterResult, p.adjust <= 0.25)
      if (!is.null(mkeggORres)) {
        print('No MKEGG results can be obatined')
        gsorCPres                     <- NULL
      } else {
        gsorCPres                     <- as.data.frame(mkeggORres@compareClusterResult)
      }
    } else {
      mkeggORres                <- enrichMKEGG(gene = pathwayGeneEntrezidInput, 
                                               organism = gsorKeggOrg, keyType = 'kegg', 
                                               minGSSize = 1, maxGSSize = 500, 
                                               pvalueCutoff  = 1, 
                                               qvalueCutoff  = 1, 
                                               pAdjustMethod = 'BH')
      
      orRes4plot                 <- mkeggORres
      orRes4plot                 <- subset(mkeggORres, p.adjust <= 0.25)
      if (!is.null(mkeggORres)) {
        print('No MKEGG results can be obatined')
        gsorCPres                <- NULL
      } else {
        gsorCPres                <- mkeggORres@result
      }
    }
    if (!is.null(gsorCPres)) {
      for (r in 1:dim(mkeggORres@result)[1]) {
        gsorCPres$geneSymbol[r] <- paste(getGeneListGeneSymbol(inputGeneList = strsplit(gsorCPres$geneID[r], split = '/')[[1]], gsorInputGenesType = 'entrezID', refGenome = refGenome), collapse = '/')
      }
    }
    print(sprintf('END over representation for MKEGG module categories DB.'))
    print('***')
    # ---
  } else if ( functionDB == 'broadMSigDB' ) {
    # ---
    print('===')
    print(sprintf('START over representation for Broad MSigDB gene sets.'))
    if (is.list(gsorInputGenes)) {
      hmORres                      <- compareCluster(geneCluster   = pathwayGeneEntrezidInput, 
                                                     fun           = 'enricher',
                                                     TERM2GENE = get.t2g(refGenome = refGenome, category = broadMSigDB.category, subcategory = broadMSigDB.subcategory), 
                                                     minGSSize = 1, maxGSSize = 500, 
                                                     pvalueCutoff  = 1, 
                                                     qvalueCutoff  = 1, 
                                                     pAdjustMethod = 'BH')
      orRes4plot                      <- hmORres
      orRes4plot@compareClusterResult <- subset(orRes4plot@compareClusterResult, p.adjust <= 0.25)
      gsorCPres                       <- as.data.frame(hmORres@compareClusterResult)
    } else {
      hmORres                <- enricher(gene = pathwayGeneEntrezidInput, 
                                         TERM2GENE = get.t2g(refGenome = refGenome, category = broadMSigDB.category, subcategory = broadMSigDB.subcategory), 
                                         minGSSize = 1, maxGSSize = 500, 
                                         pvalueCutoff  = 1, 
                                         qvalueCutoff  = 1, 
                                         pAdjustMethod = 'BH')
      
      orRes4plot                 <- hmORres
      orRes4plot                 <- subset(hmORres, p.adjust <= 0.25)
      gsorCPres                  <- hmORres@result
    }
    for (r in 1:dim(gsorCPres)[1]) {
      gsorCPres$geneSymbol[r] <- suppressMessages(paste(getGeneListGeneSymbol(inputGeneList = strsplit(gsorCPres$geneID[r], split = '/')[[1]], geneListsType = 'entrezID', refGenome = refGenome), collapse = '/'))
      
    }
    print(sprintf('END over representation for Broad MSigDB gene module categories DB.'))
    print('***')
    # ---
  } else if ( functionDB == 'ncg' ) {
    library(DOSE)
    # ---
    print('===')
    print(sprintf('START over representation for NCG (Network of Cancer Genes database) DB.'))
    if (is.list(gsorInputGenes)) {
      ncgORres                        <- compareCluster(geneCluster   = pathwayGeneEntrezidInput, 
                                                        fun           = 'enrichNCG',
                                                        minGSSize = 10, maxGSSize = 500, 
                                                        pvalueCutoff  = 1, 
                                                        qvalueCutoff  = 1, 
                                                        pAdjustMethod = 'BH', 
                                                        readable = TRUE)
    } else {
      ncgORres                   <- enrichNCG(gene = pathwayGeneEntrezidInput, 
                                              minGSSize = 10, maxGSSize = 500, 
                                              pvalueCutoff  = 1, 
                                              qvalueCutoff  = 1, 
                                              pAdjustMethod = 'BH', 
                                              readable = TRUE)
    }
    if (is.null(ncgORres)) {
      print('No Network of Cancer genes identified from input gene list')
      gsorCPres                <- NULL
    } else {
      if (is.list(gsorInputGenes)) {
        orRes4plot                      <- ncgORres
        orRes4plot@compareClusterResult <- subset(orRes4plot@compareClusterResult, p.adjust <= 0.25)
        gsorCPres                       <- as.data.frame(ncgORres@compareClusterResult) 
      } else {
        orRes4plot                      <- ncgORres
        orRes4plot                      <- subset(ncgORres, p.adjust <= 0.25)
        gsorCPres                       <- ncgORres@result
      }
      gsorCPres$geneSymbol = gsorCPres$geneID
    }
    print(sprintf('END over representation for NCG (Network of Cancer Genes database) DB.'))
    print('***')
    # ---
  } else if ( functionDB == 'DisGeNET' ) {
    library(DOSE)
    # ---
    print('===')
    print(sprintf('START over representation for DisGeNET (http://www.disgenet.org/) DB.'))
    if (is.list(gsorInputGenes)) {
      dgnORres                   <- compareCluster(geneCluster   = pathwayGeneEntrezidInput, 
                                                   fun           = 'enrichDGN',
                                                   minGSSize = 10, maxGSSize = 500, 
                                                   pvalueCutoff  = 1, 
                                                   qvalueCutoff  = 1, 
                                                   pAdjustMethod = 'BH', 
                                                   readable = TRUE)
    } else {
      dgnORres                   <- enrichDGN(gene = pathwayGeneEntrezidInput, 
                                              minGSSize = 10, maxGSSize = 500, 
                                              pvalueCutoff  = 1, 
                                              qvalueCutoff  = 1, 
                                              pAdjustMethod = 'BH', 
                                              readable = TRUE)
    }
    if (is.null(dgnORres)) {
      print('no gene-disease associations obtained')
      gsorCPres                <- NULL
    } else {
      if (is.list(gsorInputGenes)) {
        orRes4plot                      <- dgnORres
        orRes4plot@compareClusterResult <- subset(orRes4plot@compareClusterResult, p.adjust <= 0.25)
        gsorCPres                       <- as.data.frame(dgnORres@compareClusterResult)
      } else {
        orRes4plot                      <- ncgORres
        orRes4plot                      <- subset(ncgORres, p.adjust <= 0.25)
        gsorCPres                       <- dgnORres@result
      }
      gsorCPres$geneSymbol              = gsorCPres$geneID
    }
    print(sprintf('END over representation for DisGeNET (http://www.disgenet.org/) DB.'))
    print('***')
    # ---
  } else if ( functionDB == 'reactome' ) {
    # ---
    library(ReactomePA)
    print('===')
    print(sprintf('START over representation for Reactome pathway.'))
    if (is.list(gsorInputGenes)) {
      reactomeORres              <- compareCluster(geneCluster   = pathwayGeneEntrezidInput, 
                                                   fun           = 'enrichPathway',
                                                   organism = reactomeOrg, 
                                                   pvalueCutoff = 1, 
                                                   pAdjustMethod = 'BH', 
                                                   qvalueCutoff  = 1, 
                                                   minGSSize = 10, maxGSSize = 500, 
                                                   readable = TRUE)
    } else {
      reactomeORres              <- enrichPathway(gene = pathwayGeneEntrezidInput, 
                                                  organism = reactomeOrg, 
                                                  pvalueCutoff = 1, 
                                                  pAdjustMethod = 'BH', 
                                                  qvalueCutoff  = 1, 
                                                  minGSSize = 10, maxGSSize = 500, 
                                                  readable = TRUE)
    }
    if (is.null(reactomeORres)) {
      print('no reactome pathway enrichment analysis results obtained')
      gsorCPres                <- NULL
    } else {
      if (is.list(gsorInputGenes)) {
        orRes4plot                      <- reactomeORres
        orRes4plot@compareClusterResult <- subset(reactomeORres@compareClusterResult, p.adjust <= 0.25)
        gsorCPres                       <- as.data.frame(reactomeORres@compareClusterResult)
      } else {
        orRes4plot                      <- reactomeORres
        orRes4plot                      <- subset(reactomeORres, p.adjust <= 0.25)
        gsorCPres                       <- reactomeORres@result
      }
      gsorCPres$geneSymbol     = gsorCPres$geneID
    }
    print(sprintf('END over representation for Reactome pathway.'))
    print('***')
    # ---
  }
  ####complete step 3.2  pool GOs analysis.
  print('**********')
  return(list( 'resFull' = gsorCPres, 'res4plot' = orRes4plot))
}
## ------
outputCpORTestRes <- function(cpORTestRes, resFnamePrefix, plotTitle) {
  if (missing(resFnamePrefix)) resFnamePrefix <- paste(getwd(), 'TEST', sep = '/')
  write.table(x = cpORTestRes$resFull, file = paste(resFnamePrefix, '.txt', sep = ''), quote = F, sep = '\t', row.names = F, col.names = T)
  ##make GO enrichment plot
  no.showcategory <- 10
  if ('Cluster' %in% colnames(cpORTestRes$resFull)) {
    print(sprintf('There are a total of %s %s in below clusters', dim(cpORTestRes$res4plot)[1], plotTitle))
    print('A total of enriched functions in each cluster is listed below')
    print(table(as.data.frame(cpORTestRes$resFull)$Cluster))
  }
  if (dim(cpORTestRes$res4plot)[1]!=0) {
    if (no.showcategory == 10) {
      pdf(file = paste(resFnamePrefix, 'top10_category.pdf', sep = '_'), width = 10, height = 8  )
    } else if (no.showcategory == 5) {
      pdf(file = paste(resFnamePrefix, 'top5_category.pdf', sep = '_'), width = round(length(cpORTestRes$res4plot@geneClusters)/1, digits = 0), height = round(0.75*length(cpORTestRes$res4plot@geneClusters), digits = 0)  )
    }
    print(dotplot(cpORTestRes$res4plot, showCategory = no.showcategory, title = plotTitle)  + theme(axis.text.x = element_text(angle = 90)))
    dev.off()
  }
}
###################

# pdf(file = paste(resFnamePrefix, 'full_category.pdf', sep = '_'), width=25, height=30 )
# print(dotplot(cpORreactomeRes$res4plot, showCategory = length(unique((cpORreactomeRes$res4plot@compareClusterResult)$ID)), title = 'FULL'))
# dev.off()
