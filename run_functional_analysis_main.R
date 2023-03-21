##--------------------------------------------------------------------------------------##
rm(list = ls())
## ---
suppressPackageStartupMessages(library(argparse))
currentDir          <- getwd()
## ---------------------------------------------------------------------------------------
parser              <- ArgumentParser()
## add parser argument options
parser$add_argument("--workingDirectory", nargs = 1,
                    type="character", help="Required: code working/execution directory",
                    metavar = '')

parser$add_argument("--repoDirectory", nargs = 1,
                    type="character", help="Required: GSOR bitbucket repo downloaded full path",
                    metavar = '')

parser$add_argument("--inputFnames", nargs = 1,
                    type="character", help="Required: GSOR input file name",
                    metavar = '')

parser$add_argument("--inputNames", nargs = 1,
                    type="character", help="Required: GSOR input file name corresponding short names.",
                    metavar = '')

parser$add_argument("--inputType", nargs = 1,
                    default = 'geneListsInput',
                    type="character", help="Required: GSOR input type inside input file name, options are DEG_OL, DEGs, CytoScapeOutputEssential, CytoScapeOutputModules, geneListsInput, scClusterPosGenes, cripsr_mle.",
                    metavar = '')

parser$add_argument("--topN", nargs = 1,
                    default = '0',
                    type="integer", help="If specified, not 0, will only use this top number of genes for functional enrichment analysis.",
                    metavar = '')

parser$add_argument("--outputFname", nargs = 1,
                    default = 'GSOR_analysis_results',
                    type="character", help="Required: GSOR analysis results directory file name./",
                    metavar = '')

parser$add_argument("--genome", nargs = 1,
                    default = 'human',
                    type="character", help="Required: reference genome, options are: human/hs, mouse/mm, rat/rn, worm/ce.",
                    metavar = '')

parser$add_argument("--genelistType", nargs = 1,
                    default = 'GeneSymbol',
                    type="character", help="Required: the input type of gene list in the 'inputFnames' option. options are GeneSymbol, ensembleID, or entrezID.",
                    metavar = '')

parser$add_argument("--resFnameSuffix", nargs = 1,
                    default = 'gsor_results',
                    type="character", help="Required: GSOR analysis results file name suffix.",
                    metavar = '')

parser$add_argument("--runGO", nargs = 1,
                    default = 'T',
                    type="character", help="Required: whether to turn on GO enrichment analysis.",
                    metavar = '')

parser$add_argument("--goPool", nargs = 1,
                    default = 'T',
                    type="character", help="Required: when run GO enrichment analysis, whether to run MF/CC/BP separately or together.",
                    metavar = '')

parser$add_argument("--runKegg", nargs = 1,
                    default = 'T',
                    type="character", help="Required: whether to turn on KEGG enrichment analysis.",
                    metavar = '')

parser$add_argument("--runReactome", nargs = 1,
                    default = 'T',
                    type="character", help="Required: whether to turn on Reactome enrichment analysis.",
                    metavar = '')

parser$add_argument("--runHumanDisease", nargs = 1,
                    default = 'F',
                    type="character", help="Required: whether to turn on Human Disease enrichment analysis.",
                    metavar = '')

parser$add_argument("--runNetworkCancerGenes", nargs = 1,
                    default = 'F',
                    type="character", help="Required: whether to turn on Network Cancer Genes enrichment analysis.",
                    metavar = '')

parser$add_argument("--runBroadMSigDB", nargs = 1,
                    default = 'T',
                    type="character", help="Required: whether to turn on Hallmark genes enrichment analysis.",
                    metavar = '')

parser$add_argument("--broadMSigDBCategory", nargs = 1,
                    default = 'H',
                    type="character", help="Default running on BroadMSigDB Hallmark category. If 'runBroadMSigDB' is on, you can specifiy the Broad MSigDB category here, options are H (hallmark gene sets), C1 (positional gene sets), C2 (curated gene sets), CGP (chemical and genetic pertubations), CP (cannonical pathways), C7 (immunologica signature gene sets).",
                    metavar = '')

## ------
args <- parser$parse_args()
print(args)
print('---===---')
## ----------------------------------------------------------------- ##
workDir                     <- as.character(args$workingDirectory)
currentProjFullPath         <- workDir
## ---
repoDownload                <- as.character(args$repoDirectory)
gsorSourceCodes             <- list.files(path = paste(repoDownload, 'R/gsor', sep = '/'), pattern = '*.R', full.names = T )
for (f in gsorSourceCodes) source(f)
## ---
cpEnrichGoInputFnames       <- as.list(as.character(paste(unlist(strsplit(args$inputFnames, split = ', ')), sep = ', ' )))
cpEnrichGoInputNames        <- as.character((paste(unlist(strsplit(args$inputNames, split = ', ')), sep = ', ' )))
if (length(cpEnrichGoInputFnames)!=length(cpEnrichGoInputNames)) stop("Provide corresponing --inputNames for --inputFnames")
names(cpEnrichGoInputFnames) <- cpEnrichGoInputNames
## debug
# print(str(cpEnrichGoInputFnames))
# print(cpEnrichGoInputFnames[[1]])
# print(sprintf("input type is %s", as.character(args$inputType)))
# print('0000000000')
## ---
cpEnrichGoInputType         <- as.character(args$inputType)  # options are 'DEGs', 'DEG_ORG', 'CytoScapeOutputEssential', 'CytoScapeOutputModules', 'geneListsInput', 'scClusterPosGenes'
cpEnrichGoDirName           <- as.character(args$outputFname)
refGenome4cpenrichGo        <- as.character(args$genome)
geneListsType               <- as.character(args$genelistType)
resFnameSuffix              <- as.character(args$resFnameSuffix)
runGO                       <- as.logical(args$runGO)
goPool                      <- as.logical(args$goPool)
runKegg                     <- as.logical(args$runKegg)
runReactome                 <- as.logical(args$runReactome)
runHumanDisease             <- as.logical(args$runHumanDisease)
runNetworkCancerGenes       <- as.logical(args$runNetworkCancerGenes)
runBroadMSigDB              <- as.logical(args$runBroadMSigDB)
broadMSigDB.category        <- as.character(args$broadMSigDBCategory)
topN                        <- as.numeric(args$topN)
## ----------------------------------------------------------------- ##
# print(sprintf("cpEnrichGoInputType=%s", cpEnrichGoInputType)) ##for debug
## ---------------------------------------------------------------------------------------------------------------- ##
## ----------MAIN----------MAIN----------MAIN----------MAIN----------MAIN----------MAIN----------MAIN----------MAIN ##
## Main Programs to run functional enrichment analysis with below options specified in running program.
# 0.1. setup enrichGoResSaveDir for results to be save
enrichGoResSaveDir  <- paste(workDir, cpEnrichGoDirName, sep = '/')
if (!dir.exists(enrichGoResSaveDir)) dir.create(enrichGoResSaveDir)
## ----------------------------------------------------------------- ##
print('=========')
print(sprintf('functional enrichment (GSOR) analysis results will be saved at %s', enrichGoResSaveDir))
print('------')
## ----------------------------------------------------------------- ##
print('Start step 1: input gene lists processing')
# 1.1 get input gene lists for cp.enrich.go() function
if (cpEnrichGoInputType == 'DEG_OL') {
  library(openxlsx)
  sheetNames <- getSheetNames(cpEnrichGoInputFnames)
  # sheetNames <- sheetNames(degFullResFname[[l]]) #gdata()
  detach("package:openxlsx", unload=TRUE)
  # -
  cpEnrInputGeneLists <- lapply(sheetNames, function(x) {
    resFull  <- read.xlsx(file = cpEnrichGoInputFnames, header = T, sheetName = x, stringsAsFactors = 'F')
    resGenes <- resFull[,1]
  })
  cpEnrichGoInputRename <- as.logical(F)
  if (cpEnrichGoInputRename) {
    names(cpEnrInputGeneLists) <- cpEnrichGoInputNames
  } else {
    names(cpEnrInputGeneLists) <- sheetNames
  }
  # -
} else if ( cpEnrichGoInputType == 'DEGs' ) {
  print(str(cpEnrichGoInputFnames))
  cpEnrInputGeneLists         <- lapply(cpEnrichGoInputFnames, function(x) {
    resFull                   <- read.delim(file = x, header = T, check.names = F, sep = '\t', row.names = 1)
    resGenes                  <- rownames(resFull)
    return(resGenes)
  })
  names(cpEnrInputGeneLists) <- names(cpEnrichGoInputFnames)
  # -
} else if ( cpEnrichGoInputType == 'CytoScapeOutputEssential' ) {
  # -
  cpEnrInputGeneLists <- list()
  for (i in 1:length(cpEnrichGoInputFnames)) {
    cpEnrInputGeneLists[[i]] <- getEssentialGlist(cytoResCsvFname = cpEnrichGoInputFnames[[i]], degreeCutoff = degreeCutoffVals[i], refGenome4cpenrichGo = refGenome4cpenrichGo)
  }
  names(cpEnrInputGeneLists) <- names(cpEnrichGoInputFnames)
  # -
} else if ( cpEnrichGoInputType == 'CytoScapeOutputModules' ) {
  # -
  cpEnrInputGeneLists <- list()
  for (i in 1:length(cpEnrichGoInputFnames)) {
    moduleGlists <- getModuleGlist(cytoResCsvFname = cpEnrichGoInputFnames[[i]], moduleCutoff = moduleCutoffVals[i])
    names(moduleGlists) <- paste(names(cpEnrichGoInputFnames)[i], names(moduleGlists), sep = ':')
    cpEnrInputGeneLists <- append(cpEnrInputGeneLists, moduleGlists)
  }
  # -
} else if ( cpEnrichGoInputType == 'geneListsInput' ) {
  # cpEnrInputGeneLists         <- cpEnrichGoInputFnames
  print(str(cpEnrichGoInputFnames))
  cpEnrInputGeneLists         <- lapply(cpEnrichGoInputFnames, function(x) {
    resFull                   <- read.delim(file = x, header = T, check.names = F, sep = '\t')
    if ('gene' %in% colnames(resFull)) {
      if (topN!=0) {
        resGenes                  <- resFull$gene[1:topN]
      } else {
        resGenes                  <- resFull$gene
      }
      
    } else {
      if (topN!=0) {
        resGenes                  <- resFull[1:topN,1]
      } else {
        resGenes                  <- resFull[,1]
      }
      
    }
    return(resGenes)
  })
  names(cpEnrInputGeneLists) <- names(cpEnrichGoInputFnames)
} else if ( cpEnrichGoInputType == 'scClusterPosGenes' ) {
  # scRNAseqClusterPosGenes     <- read.delim(file = cpEnrichGoInputFnames, header = T, sep = '\t', check.names = F, row.names = 1)
  # scRNAseqClusterRes          <- scRNAseqClusterPosGenes %>% dplyr::group_by(cluster)
  # scRNAseqClusterGroups       <- attr(scRNAseqClusterRes, "groups")
  # cpEnrInputGeneLists         <- list()
  # for (i in 1:length(scRNAseqClusterGroups$cluster)) {
  #   cpEnrInputGeneLists[[i]]  <- scRNAseqClusterPosGenes[scRNAseqClusterGroups$.rows[[i]],] %>% dplyr::pull(gene)
  # }
  # names(cpEnrInputGeneLists)  <- paste('cluster_', scRNAseqClusterGroups$cluster, sep = '')
  # ## ---
  for (i in 1:length(cpEnrichGoInputFnames)) {
    scRNAseqClusterPosGenes   <- read.delim(file = cpEnrichGoInputFnames[[i]], header = T, sep = '\t', check.names = F)
    if (colnames(scRNAseqClusterPosGenes)[1]=="") scRNAseqClusterPosGenes <- scRNAseqClusterPosGenes[,-1]
    scRNAseqClusterRes        <- scRNAseqClusterPosGenes %>% dplyr::group_by(cluster)
    scRNAseqClusterGroups     <- attr(scRNAseqClusterRes, "groups")
    
    cpEnrInputGeneListsPrep   <- list()
    for (j in 1:length(scRNAseqClusterGroups$cluster)) {
      if (topN!=0) {
        cpEnrInputGeneListsPrep[[j]]  <- scRNAseqClusterPosGenes[scRNAseqClusterGroups$.rows[[j]],] %>% dplyr::arrange(desc(avg_log2FC)) %>% dplyr::top_n(topN, avg_log2FC) %>% dplyr::pull(gene) 
      } else {
        cpEnrInputGeneListsPrep[[j]]  <- scRNAseqClusterPosGenes[scRNAseqClusterGroups$.rows[[j]],] %>% dplyr::pull(gene)
      }
      
    }
    names(cpEnrInputGeneListsPrep)  <- paste(names(cpEnrichGoInputFnames)[i], gsub('[ ]', '-', scRNAseqClusterGroups$cluster), sep = '_')
    ## -
    if (i ==1) {
      cpEnrInputGeneLists      = cpEnrInputGeneListsPrep
    } else {
      cpEnrInputGeneLists      = c(cpEnrInputGeneLists, cpEnrInputGeneListsPrep)
    }
  }
  # print(str(cpEnrInputGeneLists))
  # print('9999999999999999')
  ## -
} else if ( cpEnrichGoInputType == 'cripsr_mle' ) {
  # print('999999999')
  # print(sprintf("input file name is %s", cpEnrichGoInputFnames))
  scRNAseqClusterPosGenes     <- read.delim(file = as.character(cpEnrichGoInputFnames), header = T, sep = '\t', check.names = F, row.names = 1)
  scRNAseqClusterPosGenes$gene <- rownames(scRNAseqClusterPosGenes)
  scRNAseqClusterRes          <- scRNAseqClusterPosGenes %>% dplyr::group_by(cluster)
  scRNAseqClusterGroups       <- attr(scRNAseqClusterRes, "groups")
  cpEnrInputGeneLists         <- list()
  for (i in 1:length(scRNAseqClusterGroups$cluster)) {
    cpEnrInputGeneLists[[i]]  <- scRNAseqClusterPosGenes[scRNAseqClusterGroups$.rows[[i]],] %>% dplyr::pull(gene)
  }
  names(cpEnrInputGeneLists)  <- as.character(scRNAseqClusterGroups$cluster)
}
print('END step 1: input gene lists processing')
print('Input gene list looks like as below:')
str(cpEnrInputGeneLists)
print('=========')
## ----------------------------------------------------------------- ##
# 2.1. if 'runGO' on, run enrichment GO analysis with cp.enrich.go() function
print('------')
print('Start step 2: Start GSOR analysis')
# print(sprintf("running %s", runBroadMSigDB)) #for debug
if (runGO) {
  print('START Step 2.1 clustered input gene lists GO over-representation enrichment analysis')
  if (goPool) {
    cpORgoPoolRes <- cpORTest(gsorInputGenes     = cpEnrInputGeneLists,
                              # gsorInputGenes     = lapply(cpEnrInputGeneLists, '[', 1:10),
                              gsorInputGenesType = geneListsType,
                              refGenome          = refGenome4cpenrichGo,
                              functionDB         = 'go')
    resFnamePrefix = paste(enrichGoResSaveDir, '/', resFnameSuffix, '_cpEnrich_GoPool_clustered_res', sep = '')
    outputCpORTestRes(cpORTestRes = cpORgoPoolRes, resFnamePrefix = resFnamePrefix, plotTitle = 'GO enrichment')
  } else {
    cpORgoSepResBP <- cpORTest(gsorInputGenes     = cpEnrInputGeneLists,
                              gsorInputGenesType = geneListsType,
                              refGenome          = refGenome4cpenrichGo,
                              functionDB         = 'gobp')
    cpORgoSepResCC <- cpORTest(gsorInputGenes     = cpEnrInputGeneLists,
                              gsorInputGenesType = geneListsType,
                              refGenome          = refGenome4cpenrichGo,
                              functionDB         = 'gocc')
    cpORgoSepResMF <- cpORTest(gsorInputGenes     = cpEnrInputGeneLists,
                              gsorInputGenesType = geneListsType,
                              refGenome          = refGenome4cpenrichGo,
                              functionDB         = 'gomf')
    resFnamePrefix = paste(enrichGoResSaveDir, '/', resFnameSuffix, '_cpEnrich_GoSep_clustered_res', sep = '')
    outputCpORTestRes(cpORTestRes = cpORgoSepResBP, resFnamePrefix = paste(resFnamePrefix, '_BP', sep = ''), plotTitle = 'GO BP enrichment')
    outputCpORTestRes(cpORTestRes = cpORgoSepResCC, resFnamePrefix = paste(resFnamePrefix, '_CC', sep = ''), plotTitle = 'GO CC enrichment')
    outputCpORTestRes(cpORTestRes = cpORgoSepResMF, resFnamePrefix = paste(resFnamePrefix, '_MF', sep = ''), plotTitle = 'GO MF enrichment')
  }
  print('COPMPLETE Step 2.1 clustered input gene lists GO enrichment analysis')
  print('======')
}
## ---
## 2.2 if 'runReactome' on, use cpORTest() function to conduct reactome pathway over representation analysis
if (runReactome) {
  print('Start step 2.2 clusterProfiler over-representation on Reactome Pathway')
  cpORreactomeRes <-  cpORTest(gsorInputGenes     = cpEnrInputGeneLists,
                               gsorInputGenesType = geneListsType,
                               refGenome          = refGenome4cpenrichGo,
                               functionDB         = 'reactome')
  resFnamePrefix = paste(enrichGoResSaveDir, '/', resFnameSuffix, '_cpEnrich_reactome_clustered_res', sep = '')
  outputCpORTestRes(cpORTestRes = cpORreactomeRes, resFnamePrefix = resFnamePrefix, plotTitle = 'Reactome Enrichment')
  print('END step 2.2 clusterProfiler over-representation on Reactome Pathway')
}
## ---
## 2.3 if 'runKegg' on, use cpORTest() function to conduct KEGG pathway over representation analysis
if (runKegg) {
  print('Start step 2.3 clusterProfiler over-representation on KEGG Pathway')
  cpORkeggRes <-  cpORTest(gsorInputGenes     = cpEnrInputGeneLists,
                           gsorInputGenesType = geneListsType,
                           refGenome          = refGenome4cpenrichGo,
                           functionDB         = 'kegg')
  resFnamePrefix = paste(enrichGoResSaveDir, '/', resFnameSuffix, '_cpEnrich_kegg_clustered_res', sep = '')
  outputCpORTestRes(cpORTestRes = cpORkeggRes, resFnamePrefix = resFnamePrefix, plotTitle = 'KEGG Enrichment')
  print('END step 2.3 clusterProfiler over-representation on KEGG Pathway')
}
## 2.4 if 'runHumanDisease' on, use cpORTest() function to conduct DisGeNET pathway over representation analysis
if (runHumanDisease) {
  print('Start step 2.4 clusterProfiler over-representation on DisGeNET Pathway')
  cpORHumanDiseaseRes <-  cpORTest(gsorInputGenes     = cpEnrInputGeneLists,
                                   gsorInputGenesType = geneListsType,
                                   refGenome          = refGenome4cpenrichGo,
                                   functionDB         = 'DisGeNET')
  resFnamePrefix = paste(enrichGoResSaveDir, '/', resFnameSuffix, '_cpEnrich_DisGeNET_clustered_res', sep = '')
  outputCpORTestRes(cpORTestRes = cpORHumanDiseaseRes, resFnamePrefix = resFnamePrefix, plotTitle = 'DisGeNET Enrichment')
  print('END step 2.4 clusterProfiler over-representation on DisGeNET Pathway')
}
## ---
## 2.5 if 'ncg' on, use cpORTest() function to conduct NCG (national cancer genes) pathway over representation analysis
if (runNetworkCancerGenes) {
  print('Start step 2.5 clusterProfiler over-representation on NCG Pathway')
  cpORNCGRes <-  cpORTest(gsorInputGenes     = cpEnrInputGeneLists,
                          gsorInputGenesType = geneListsType,
                          refGenome          = refGenome4cpenrichGo,
                          functionDB         = 'ncg')
  resFnamePrefix = paste(enrichGoResSaveDir, '/', resFnameSuffix, '_cpEnrich_NCG_clustered_res', sep = '')
  outputCpORTestRes(cpORTestRes = cpORNCGRes, resFnamePrefix = resFnamePrefix, plotTitle = 'NCG Enrichment')
  print('END step 2.5 clusterProfiler over-representation on NCG Pathway')
}
## ---
## ---
## 2.6 if 'broadMSigDB' on, use cpORTest() function to conduct Broad MSigDB pathway over representation analysis
if (runBroadMSigDB) {
  # print("00000000000") ## for debug
  print(sprintf('Start step 2.6 hallmarker genes (%s) over-representation analysis', broadMSigDB.category))
  # print("00000000000") ## for debug
  cpHallmarkRes <-  cpORTest(gsorInputGenes     = cpEnrInputGeneLists,
                             gsorInputGenesType = geneListsType,
                             refGenome          = refGenome4cpenrichGo,
                             functionDB         = 'broadMSigDB')
  resFnamePrefix = paste(enrichGoResSaveDir, '/', resFnameSuffix, sprintf('_cpEnrich_broadMsigDB_category_%s_clustered_res', broadMSigDB.category), sep = '')
  outputCpORTestRes(cpORTestRes = cpHallmarkRes, resFnamePrefix = resFnamePrefix, plotTitle = 'hallmark Enrichment')
  print(sprintf('END step 2.6 hallmarker genes (%s) over-representation analysis', broadMSigDB.category))
}
## ---
print('END step 2: END GSOR analysis')
print('------')
## -----------END-----------END-----------END-----------END-----------END-----------END-----------END-----------END ##
## ---------------------------------------------------------------------------------------------------------------- ##











