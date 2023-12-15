## This script could summarize enrichment analysis results under folder 'enrichResDirName' 
## with all full analysis results saved in '*cpEnrich.*txt$' files
## with respect to different enrichment analysis (including GO, kegg, reactome etc.)
## -------------------------------------------------------------------------------------##
if(!exists('ylimit')) ylimit = 0.05
if(!exists('xlimit')) xlimit = NULL
if (toi.selection) {
  if (cluster.selection) {
    if (exists("toi.fnameSuffix")) {
      ressummaryDir               <- paste(enrichResDirName, sprintf('GSOR_results_summary_toi%s_cluster%s', toi.fnameSuffix, cluster.sel.fnameSuffix), sep = '/')
    } else {
      ressummaryDir               <- paste(enrichResDirName, sprintf('GSOR_results_summary_toi_cluster%s', cluster.sel.fnameSuffix), sep = '/')
    }
  } else {
    if (exists("toi.fnameSuffix")) {
      ressummaryDir               <- paste(enrichResDirName, sprintf('GSOR_results_summary_toi%s', toi.fnameSuffix), sep = '/')
    } else {
      ressummaryDir               <- paste(enrichResDirName, 'GSOR_results_summary_toiSelection', sep = '/')
    }
    
  }
} else {
  if (cluster.selection) {
    ressummaryDir               <- paste(enrichResDirName, sprintf('GSOR_results_summary_cluster%s', cluster.sel.fnameSuffix), sep = '/')
  } else {
    ressummaryDir               <- paste(enrichResDirName, 'GSOR_results_summary', sep = '/')
  }
}
if (!dir.exists(ressummaryDir)) dir.create(ressummaryDir)
## --------------------
enrichResFnames             <- list.files(path = enrichResDirName, pattern = sprintf('*%s.*txt$', resFnamePrefix), full.names = T)
# # options(java.parameters = "-Xmx32000m")
# options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx8192m"))
# ---
## 1 read in 1) all 'cpEnrich*txt' files inside 'enrichResDirName'; 2) assoicated pVals for each enrichment result
cpORenrichmentResAll        <- lapply(enrichResFnames, function(x) read.delim(file = x, header = T, sep = '\t', stringsAsFactors = F) )
cpORenrichmentResAllNames   <- sapply(enrichResFnames, function(x) gsub(pattern = '.txt', replacement = '', x = basename(x)))
names(cpORenrichmentResAll) <- cpORenrichmentResAllNames
if (exists("enrichres.sel")) {
  cpORenrichmentResAll      <- cpORenrichmentResAll[enrichres.sel]
  cpORenrichmentResAllNames <- cpORenrichmentResAllNames[enrichres.sel]
}
## ---
## 2.whether to select terms of interest
if (toi.selection) {
  if (exists("toi.selection.fname")) {
    if (tools::file_ext(toi.selection.fname)=='xlsx') {
      toi.sel.org <- xlsx::read.xlsx(file = toi.selection.fname, sheetIndex = 1, header = F)
      toi.sel <- toi.sel.org$X1
    } else {
      toi.sel.org <- read.delim(file = toi.selection.fname, header = T, sep = '\t')
      toi.sel <- toi.sel.org$X1
    }
  } else {
    toi.sel <- toi.sel
  }
  if(precise.match) {
    cpORenrichmentResAll2 <- lapply(cpORenrichmentResAll, function(x) x[grep(paste(sprintf('^%s$', toi.sel), collapse = '|'), x$Description),])
  } else {
    cpORenrichmentResAll2 <- lapply(cpORenrichmentResAll, function(x) x[grep(paste(toi.sel, collapse = '|'), x$Description),])
  }
  
  no.tois <- data.frame('full' = sapply(cpORenrichmentResAll, function(x) dim(x)[1]), 'toi.sel' = sapply(cpORenrichmentResAll2, function(x) dim(x)[1]))
  print('-=-=-=-=-=-=-')
  print(no.tois)
  print('-=-=-=-=-=-=-')
  cpORenrichmentResAll <- cpORenrichmentResAll2
  } 
## ---
if (cluster.selection) {
  cpORenrichmentResAll3 <- lapply(cpORenrichmentResAll, function(x) x[grep(paste(cluster.sel, collapse = '|'), x$Cluster),])
  cpORenrichmentResAll  <- cpORenrichmentResAll3
}
## ---
if (length(padjValListInput) == 1) {
  padjValList               <- lapply(cpORenrichmentResAll, function(x) rep( x = padjValListInput, times = length(levels(factor(levels(factor(x$Cluster))))) ) )
} else {
  if (length(padjValListInput)!=length(cpORenrichmentResAllNames)) stop("please provide corresponding length list or only 1 item list")
  ## -
  padjValList               <- padjValListInput
  names(padjValList)        <- cpORenrichmentResAllNames
}
## -
if (length(fulldotplotValInput) == 1) {
  fulldotplotVal            <- lapply(cpORenrichmentResAll, function(x) fulldotplotValInput[[1]])
} else {
  if (length(fulldotplotValInput)!=length(cpORenrichmentResAllNames)) stop("please provide corresponding length list or only 1 item list")
  ## -
  fulldotplotVal            <- fulldotplotValInput
  names(fulldotplotVal)     <- cpORenrichmentResAllNames
}
## -
if (plotTop20|plotTop30) {
  if (length(topEnrichDotplotValInput) == 1) {
    topEnrichDotplotVal       <- lapply(cpORenrichmentResAll, function(x) topEnrichDotplotValInput[[1]])
  } else {
    if (length(topEnrichDotplotValInput)!=length(cpORenrichmentResAllNames)) stop("please provide corresponding length list or only 1 item list")
    topEnrichDotplotVal       <- topEnrichDotplotValInput
    names(topEnrichDotplotVal)<- cpORenrichmentResAllNames
  }
}

## ---
## summarize CP OR enrichment analysis results, loop over all *txt files inside 'enrichResDirName', 
## if all 7 DB conducted, sort by order'DisGeNET, GoPool, GoSep_BP, GoSep_CC, GoSep_MF, kegg, reactome'
for (i in 1:length(cpORenrichmentResAll)) {
  if (dim(cpORenrichmentResAll[[i]])[1]!=0) {
    print('=========')
    print(sprintf("%s. processing '%s.txt' results", i, names(cpORenrichmentResAll)[i]))
    # 1.1 summarize the overall enricheGO results
    dataCombDfFull            <- cpORenrichmentResAll[[i]]
    clusterNames              <- levels(factor(cpORenrichmentResAll[[i]]$Cluster))
    enrichNoFullPrep          <- data.frame(table(dataCombDfFull$Cluster))
    enrichNoFull              <- data.frame(enrichNoFullPrep$Freq, row.names = enrichNoFullPrep$Var1 )
    # -
    # 1.2 summarize the enriched analysis results on specified adjusted p-value
    padjVal                   <- padjValList[[i]]
    enrichNoFull$padjVal      <- padjVal[1:dim(enrichNoFull)[1]]
    
    if (fdrOn) {
      dataCombDfadjP            <- resLists2DfCluster(listRes = resDfCluster2list(cpORenrichmentResClusterDF = dataCombDfFull), adjpVal = enrichNoFull$padjVal)
    }else {
      dataCombDfadjP            <- resLists2DfCluster(listRes = resDfCluster2list(cpORenrichmentResClusterDF = dataCombDfFull), adjpVal = enrichNoFull$padjVal, fdrOn = as.logical(F))
    }
    
    if (dim(dataCombDfadjP)[1] == 0) {
      print('NO enriched functions identified')
    } else {
      enrichNoAdjpPrep        <- as.data.frame(table(dataCombDfadjP$Cluster))
      enrichNoAdjp            <- data.frame(enrichNoAdjpPrep$Freq, row.names = enrichNoAdjpPrep$Var1)
      if (grepl(pattern = 'GoPool', names(cpORenrichmentResAll)[i])) {
        # 1.3 summarize No. BP/CC/MF of enriched GO
        ontCatOnFull          <- as.data.frame(summarizeGoCat(dataCombDfFull = dataCombDfFull))
        ontCatOnAdjp          <- as.data.frame(summarizeGoCat(dataCombDfFull = dataCombDfadjP))
        ernchiNoSummary       <- dplyr.merge.all(enrichNoFull, ontCatOnFull, enrichNoAdjp, ontCatOnAdjp)
        colnames(ernchiNoSummary) <- gsub(pattern = '.Freq|.x|.y', x = names(ernchiNoSummary), replacement = '')
      } else {
        ernchiNoSummary       <- dplyr.merge.all(enrichNoFull, enrichNoAdjp)
        colnames(ernchiNoSummary) <-  gsub(pattern = '.Freq', x = names(ernchiNoSummary), replacement = '')
      }
      # -
      summyFname              <- paste(ressummaryDir, '/', names(cpORenrichmentResAll)[i], '_resSummary.txt', sep = '')
      write.table(x = ernchiNoSummary, file = summyFname, col.names = NA, row.names = T, sep = '\t', quote = F)
      ## ---
      # 1.2 output enriched GO at specified adjp values.
      if (length(unique(padjVal)) == 1) {
        if (fdrOn) {
          cpORenrichmentResFname  <- paste(ressummaryDir, '/', names(cpORenrichmentResAll)[i], '_adjustPval', gsub('[.]', '', unique(padjVal)), '.xlsx', sep = '')
          adjPdotPlotFname        <- paste(ressummaryDir, '/', names(cpORenrichmentResAll)[i], '_adjustPval', gsub('[.]', '', unique(padjVal)), '_dotplot.pdf', sep = '')
          adjPtopDotPlotFname     <- paste(ressummaryDir, '/', names(cpORenrichmentResAll)[i], '_adjustPval', gsub('[.]', '', unique(padjVal)), '_Top20_dotplot.pdf', sep = '')
          # topGOpieChartResFname <- paste(ressummaryDir, '/topGO_pieChart_', gsub('[.]', '', unique(padjVal)), '.xlsx', sep = '')
        } else {
          cpORenrichmentResFname  <- paste(ressummaryDir, '/', names(cpORenrichmentResAll)[i], '_nominalPval', gsub('[.]', '', unique(padjVal)), '.xlsx', sep = '')
          adjPdotPlotFname        <- paste(ressummaryDir, '/', names(cpORenrichmentResAll)[i], '_nominalPval', gsub('[.]', '', unique(padjVal)), '_dotplot.pdf', sep = '')
          adjPtopDotPlotFname     <- paste(ressummaryDir, '/', names(cpORenrichmentResAll)[i], '_nominalPval', gsub('[.]', '', unique(padjVal)), '_Top20_dotplot.pdf', sep = '')
          # topGOpieChartResFname <- paste(ressummaryDir, '/topGO_pieChart_', gsub('[.]', '', unique(padjVal)), '.xlsx', sep = '')
        }
        
      } else {
        if (fdrOn){
          cpORenrichmentResFname  <- paste(ressummaryDir, '/', names(cpORenrichmentResAll)[i], '_adjustPval', paste(unlist(gsub('[.]', '', padjVal)), collapse = '_'), '.xlsx', sep = '')
          adjPdotPlotFname        <- paste(ressummaryDir, '/', names(cpORenrichmentResAll)[i], '_adjustPval', paste(unlist(gsub('[.]', '', padjVal)), collapse = '_'), '_dotplot.pdf', sep = '')
          adjPtopDotPlotFname     <- paste(ressummaryDir, '/', names(cpORenrichmentResAll)[i], '_adjustPval', gsub('[.]', '', unique(padjVal)), '_Top20_dotplot.pdf', sep = '')
          # topGOpieChartResFname <- paste(ressummaryDir, '/topGO_pieChart_', paste(unlist(gsub('[.]', '', padjVal)), collapse = '_'), '.xlsx', sep = '')
        } else {
          cpORenrichmentResFname  <- paste(ressummaryDir, '/', names(cpORenrichmentResAll)[i], '_nominalPval', paste(unlist(gsub('[.]', '', padjVal)), collapse = '_'), '.xlsx', sep = '')
          adjPdotPlotFname        <- paste(ressummaryDir, '/', names(cpORenrichmentResAll)[i], '_nominalPval', paste(unlist(gsub('[.]', '', padjVal)), collapse = '_'), '_dotplot.pdf', sep = '')
          adjPtopDotPlotFname     <- paste(ressummaryDir, '/', names(cpORenrichmentResAll)[i], '_nominalPval', gsub('[.]', '', unique(padjVal)), '_Top20_dotplot.pdf', sep = '')
          # topGOpieChartResFname <- paste(ressummaryDir, '/topGO_pieChart_', paste(unlist(gsub('[.]', '', padjVal)), collapse = '_'), '.xlsx', sep = '')
        }
      }
      # -
      clusterLevels             <- levels(factor(dataCombDfadjP$Cluster)) 
      for (l in 1:length(clusterLevels)) {
        if (l ==1 ) {
          write.xlsx2(x = dataCombDfadjP %>% dplyr::filter(Cluster == clusterLevels[l]),
                      file = cpORenrichmentResFname, row.names = F, sheetName = gsub(pattern = '-|:', replacement = '', clusterLevels[l]))
          gc()
        } else {
          if (dim(dataCombDfadjP %>% dplyr::filter(Cluster == clusterLevels[l]))[1]!=0) {
            write.xlsx2(x = dataCombDfadjP %>% dplyr::filter(Cluster == clusterLevels[l]),
                        file = cpORenrichmentResFname, row.names = F, sheetName = gsub(pattern = '-|:', replacement = '', clusterLevels[l]), append = l > 1)
            gc()
          }
        }
        
      }
      ## -
      # ---
      # 3. make dotplot for specified adjPval
      if ( xLevelsReorder ) {
        plotlinkAll      <- makeDotPlot(res2plot = dataCombDfadjP, padjon = as.logical(fdrOn), xLevels = xLevelsReorderVals, ylimit = ylimit, xlimit = xlimit)
      } else {
        plotlinkAll      <- makeDotPlot(res2plot = dataCombDfadjP, padjon = as.logical(fdrOn), xLevels = clusterLevels, ylimit = ylimit, xlimit = xlimit)
      }
      
      ggsave(filename = adjPdotPlotFname, plot = plotlinkAll, device = 'pdf', width = fulldotplotVal[[i]][1], height = fulldotplotVal[[i]][2], limitsize =F)
      
      if (plotTop20|plotTop30) {
        # 3.2 top enrichment dotplot
        for (c in 1:length(clusterLevels)) {
          if (plotTop20) {
            topOrRes <- dplyr::filter(.data = dataCombDfadjP, Cluster == as.character(clusterLevels[c])) %>% dplyr::filter(between(row_number(p.adjust), 1, 20))
          }
          if (plotTop30) {
            topOrRes <- dplyr::filter(.data = dataCombDfadjP, Cluster == as.character(clusterLevels[c])) %>% dplyr::filter(between(row_number(p.adjust), 1, 30))
          }
          if ( c == 1 ) {
            dotplotInput <- topOrRes
          } else {
            dotplotInput <- rbind(dotplotInput, topOrRes)
          }
        }
        if ( xLevelsReorder ) {
          plotlinkAll <- makeDotPlot(res2plot = dotplotInput, padjon = as.logical(fdrOn), xLevels = xLevelsReorderVals, ylimit = ylimit, xlimit = xlimit)
        } else {
          plotlinkAll <- makeDotPlot(res2plot = dotplotInput, padjon = as.logical(fdrOn), xLevels = clusterLevels, ylimit = ylimit, xlimit = xlimit)
        }
        
        if (plotTop30) {
          ggsave(filename = sub("_Top20_dotplot", "_Top30_dotplot", adjPtopDotPlotFname), plot = plotlinkAll, device = 'pdf', width = topEnrichDotplotVal[[i]][1], height = topEnrichDotplotVal[[i]][2], limitsize =F)
        } else {
          ggsave(filename = adjPtopDotPlotFname, plot = plotlinkAll, device = 'pdf', width = topEnrichDotplotVal[[i]][1], height = topEnrichDotplotVal[[i]][2], limitsize =F)
        }
      }
      ## ---
    }
  } else {
    print('=========')
    print(sprintf("%s. processing '%s.txt' results", i, names(cpORenrichmentResAll)[i]))
    print("No erniched functions for this data")
  }
}
## -------------------------------------------------------------------------------------##
