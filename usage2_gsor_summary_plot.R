rm(list = ls())
## ---------
# runGarder                   = as.logical(T)
runGarder                   = as.logical(F)
if (runGarder) {
  # bitbucketGitPath <- '/gpfs/data/biocore-analysis/yli/functional_pathway_analysis/'
  bitbucketGitPath <- '/gpfs/data/bioinformatics/yli2/gsor_functional_enrichment_analysis'
} else {
  # bitbucketGitPath <- '/Volumes/yli/functional_pathway_analysis/'
  bitbucketGitPath <- '/Volumes/bioinformatics/yli2/gsor_functional_enrichment_analysis'
}
## -
## source all related codes in downloaded repo
gsorSourceCodes            <- list.files(path = paste(bitbucketGitPath, 'R/gsor', sep = '/'), pattern = '*.R', full.names = T )
for (f in gsorSourceCodes) source(f)
## ---------------------------------------------------------------------------- ##
workDir          = "/Volumes/yli/CRI-BIO-757-MED-ELengyel-yli-lszhu/pre-post_menopausal_analysis/FT_pre_post_menopausal_integration_results/functional_enrichment_results_preManuscript"
setwd(workDir)
## ---
# projPrep(runHPC = runGarder, workDirPath = workDir) ## after this workDir, sysDir will be updated
## ---------------------------------------------------------------------------- ##
enrichResDirName           <- paste(workDir, 'menstrual_proliferativeSecComp_DEGs_STs_gsor_GOsep/', sep = '/')
resFnamePrefix             <- 'menstrual_proSecComp_STs_DEGs'  
print(sprintf("If corresponding, %s items should be in the list, they are '%s'", 
              length(list.files(path = enrichResDirName, pattern = '*cpEnrich.*txt$', full.names = T)), 
              paste(basename(list.files(path = enrichResDirName, pattern = '*cpEnrich.*txt$', full.names = T)), collapse = '; ')) )
## if all 8 DB conducted, sort by order'DisGeNET, GoPool, GoSep_BP, GoSep_CC, GoSep_MF, kegg, NetworkCancerGenes, reactome'
## this one has 3 GO sep pathways
# ## ------
## 1. full plots of all cell types
toi.selection              <- as.logical(F)
cluster.selection          <- as.logical(F)
# cluster.sel                <- sprintf('Beta%s_p1', 1:6)
# cluster.sel.fnameSuffix    <- 'BetaOnly_panel1'
fdrOn                      <- as.logical(T)
padjValListInput           <- 0.05 ##order is bp, cc(0 sel1), mf
fulldotplotValInput        <- list(c(9.5, 40),
                                   c(8.5, 16),
                                   c(8, 18)) ##order BP, CC, MF
plotTop20 = as.logical(F)
xLevelsReorder             <- as.logical('T')
xLevelsReorderVals <- c('ST1_proliferative', 'ST1_secretory', 'ST2_proliferative', 'ST2_secretory', 'ST3_proliferative', 'ST3_secretory')
precise.match = as.logical('F')
axis.y.size = 14
## ------
# ## select refined GO BP, to be added------------------
# toi.selection              <- as.logical(T)
# toi.fnameSuffix            <- 'MJsel4'
# ## March 1 slack selection, only keep highlighted ones
# toi.sel <- c('epithelium migration', 
#              'ERK1 and ERK2 cascade', 
#              'regulation of immune effector process',
#              'response to steroid hormone', 
#              'MHC class II protein complex assembly', 
#              'peptide antigen assembly with MHC class II protein complex', 
#              'signal transduction by p53 class mediator', 
#              'DNA damage response, signal transduction by p53 class mediator', 
#              'antigen processing and presentation of exogenous peptide antigen via MHC class II', 
#              'developmental growth involved in morphogenesis', 
#              'microtubule-based movement', 'microtubule-based transport', 
#              'response to fibroblast growth factor', 
#              'cellular response to fibroblast growth factor stimulus',
#              'DNA recombination', 
#              'DNA replication', 
#              'DNA-dependent DNA replication maintenance of fidelity', 
#              'somatic cell DNA recombination', 
#              'chromatin remodeling at centromere', 
#              'cell cycle DNA replication')
# cluster.selection          <- as.logical(F)
# # cluster.sel                <- sprintf('Beta%s_p1', 1:6)
# # cluster.sel.fnameSuffix    <- 'BetaOnly_panel1'
# fdrOn                      <- as.logical(T)
# padjValListInput           <- 0.1 ##order is bp, cc(0 sel1), mf
# fulldotplotValInput        <- list(c(14, 10),
#                                    c(11, 10),
#                                    c(10.5, 7)) ##order BP, CC, MF
# plotTop20 = as.logical(F)
# # topEnrichDotplotValInput   <- list(c(16, 14), ##No enrichment, no use
# #                                    c(12, 16),
# #                                    c(18, 16))
# xLevelsReorder             <- as.logical('T')
# xLevelsReorderVals <- sprintf('OE_%s', c('SE1-pre', 'SE2-pre', 'SE3-pre', 'SE-post', 'SE-post-specific'))
# precise.match = as.logical('T')
# legend.position = 'right'
## ------------
source(file.path(sprintf('%s/run_summarizeEnrichGoRes_gsor.R', bitbucketGitPath)))

## ---------------------------------------------------------------------------- ##
