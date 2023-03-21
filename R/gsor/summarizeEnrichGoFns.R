library(xlsx)
library(dplyr)
library(ggplot2)
library(gplots)
if (!exists("axis.y.size")) axis.y.size = 20
if (!exists("axis.x.size")) axis.x.size = 20
if (!exists("legend.position")) legend.position = 'bottom'
if (!exists("legend.size")) legend.size=14
if (!exists("legend.title.size")) legend.title.size= 15
  

trim <- function (x) gsub("^\\s+|\\s+$", "", x)
options(java.parameters = "-Xmx4000m")
jgc <- function()
{
  gc()
  .jcall("java/lang/System", method = "gc")
}    
## -
trim <- function (x) gsub("^\\s+|\\s+$", "", x)
## -
merge.all <- function(x, ..., by = "row.names") {
  L <- list(...)
  for (i in seq_along(L)) {
    x <- merge(x, L[[i]], by = by)
    rownames(x) <- x$Row.names
    x$Row.names <- NULL
  }
  return(x)
}
dplyr.merge.all <- function(x, ..., by = "row.names") {
  L <- list(...)
  for (i in seq_along(L)) {
    x2merge <- x %>% mutate('rowName' = rownames(x))
    L2merge <- L[[i]] %>% mutate('rowName' = rownames(L[[i]])) 
    x       <- dplyr::full_join(x = x2merge, y = L2merge, by = 'rowName')
    rownames(x) <- x$rowName
    x       <- x %>% dplyr::select(-rowName)
  }
  return(x)
}
## ---
## 0.2 define functions used for change the enrichment results in list to 
##     the combined DF with adding cluter info as the names() from enrichment list
resLists2DfCluster <- function(listRes, topNo, adjpVal, adjpValRoundup) {
  if (missing(adjpValRoundup)) adjpValRoundup <- as.logical(F)
  if (missing(topNo) & !missing(adjpVal)) {
    listResCluster <- lapply(1:length(listRes), function(x) {
      if (adjpValRoundup) {
        lisResSep    <- listRes[[x]] %>% dplyr::mutate('p.adjust.roundup' = round(p.adjust, digits = 1)) 
        listResSel   <- lisResSep%>% dplyr::filter(p.adjust.roundup<=adjpVal[x])
      } else {
        listResSel   <- listRes[[x]] %>% dplyr::filter(p.adjust<=adjpVal[x])
      }
      listResSel   <- dplyr::mutate(listResSel, 'Cluster' = names(listRes)[x])
      print(sprintf('No. of enrichGO results for \'%s\' res at FDR adjusted p-value of %s is %s', unique(listResSel$Cluster), adjpVal[x], dim(listResSel)[1] ))
      listResSel
    })
  } else if (!missing(topNo) & missing(adjpVal)) {
    listResCluster <- lapply(1:length(listRes), function(x) {
      listResDF    <- dplyr::mutate(listRes[[x]], 'Cluster' = names(listRes)[x])
      listResTop   <- listResDF %>% dplyr::filter(between(row_number(p.adjust), 1, topNo[x]))
      listResTop
    })
  } else if (missing(topNo) & missing(adjpVal)) {
    listResCluster <- lapply(1:length(listRes), function(x) {
      listResDF    <- dplyr::mutate(listRes[[x]], 'Cluster' = names(listRes)[x])
      print(sprintf('total No. enrichGO results for \'%s\' res is %s', unique(listResDF$Cluster), dim(listResDF)[1] ))
      listResDF
    })
  }
  resListComb <- dplyr::bind_rows(listResCluster)
  return(resListComb)
}
# ---
# 0.3 define function to make dotplot with combined Df obtained from resLists2DfCluster()
makeDotPlot <- function(res2plot, padjon, xLevels, yLevels, ylimit=0.05, xlimit=NULL) {
  heightval = length((res2plot$res)$Description) * 10
  res2plot$GeneRatio2 <- as.numeric(as.character(res2plot$Count))/as.numeric(as.character(sapply(res2plot$GeneRatio, function(x) strsplit(as.character(x), split = '/')[[1]][2]))) 
  
  if (missing(xLevels)) {
    xLevels <- levels(factor(res2plot$Cluster))
  } else {
    cluster.factor <- as.factor(res2plot$Cluster)
    levels(cluster.factor) <- xLevels
  }
  if (missing(yLevels)) {
    yLevels <- unique(factor(res2plot$Description))
  } else {
    Description.factor <- as.factor(res2plot$Description)
    levels(Description.factor) <- yLevels
  }
  # print(head(res2plot))
  # print(table(res2plot$Cluster))
  # print('000000')
  if(is.null(xlimit)) xlimit = max(res2plot$Count)
  print(sprintf("xlimit('max(res2plot$Count)') is %s; ylimit is %s", xlimit, ylimit))
  dotPlot <- ggplot(res2plot, aes(x=factor(Cluster, levels = as.character(xLevels)), y=factor(Description, levels = rev(yLevels) )))
  if (padjon) {
    dotPlot <- dotPlot + geom_point(aes(size = as.numeric(Count), colour = as.numeric(as.character(p.adjust))), show.legend = T ) + theme_bw()
  }else {
    dotPlot <- dotPlot + geom_point(aes(size = as.numeric(Count), colour = as.numeric(as.character(pvalue))), show.legend = T ) + theme_bw()
  }
  # print(xLevels)
  dotPlot <- dotPlot + scale_x_discrete(labels = as.character(xLevels), limits = as.character(xLevels))
  dotPlot <- dotPlot + xlab(xLevels)
  dotPlot <- dotPlot + labs(x = "", y = "")
  dotPlot <- dotPlot + theme(axis.text.y = element_text(colour = "black", size = axis.y.size))
  dotPlot <- dotPlot + theme(axis.text.x = element_text(colour = "black", size = axis.x.size, angle = 90, vjust = 0.5, hjust = 1))
  # dotPlot <- dotPlot + labs(title = "Selected over represented GOs") 
  dotPlot <- dotPlot + theme(plot.title = element_text(size = 20, hjust = 0.5))
  if (padjon) {
    # dotPlot <- dotPlot + scale_colour_gradient(name = "adjusted \n p-value\n")
    dotPlot <- dotPlot + scale_colour_gradient(name = "FDR adjusted\n p-value\n", limits= c(0, ylimit))
    # dotPlot <- dotPlot + scale_colour_gradient(name = "adjusted \n p-value\n",
    #                                      limits=c(0, round2(max(cpEnrichGoRes$p.adjust), pround.digits)), low="red", high="blue")
  }else {
    dotPlot <- dotPlot + scale_colour_gradient(name = "original \n p-value\n", limits= c(0, ylimit))
    # dotPlot <- dotPlot + scale_colour_gradient(name = "original \n p-value\n",
    #                                      limits=c(0, round2(max(cpEnrichGoRes$pvalue), pround.digits)), low="red", high="blue")
  }
  dotPlot <- dotPlot + scale_size(name = "Gene Count")
  # dotPlot <- dotPlot + scale_size(name = "Gene Ratio", range = c(1,10), breaks = c(0.25, 0.50, 0.75, 1))
  # dotPlot <- dotPlot + scale_size(name = "count", range = c(min(res.pcut$Count), max(res.pcut$Count)),
  #                     breaks=c(min(res.pcut$Count), 5, 10, max(res.pcut$Count)))
  dotPlot <- dotPlot + theme(plot.margin = unit(c(0.5,0.5,0.1,0.1), "cm") ) ##top, right, bottom, left margins
  dotPlot <- dotPlot + theme(legend.text=element_text(size=legend.size), legend.title = element_text(size=legend.title.size) )
  if (legend.position == 'bottom') {
    dotPlot <- dotPlot + theme(legend.direction = "horizontal", legend.position = legend.position, legend.box = "vertical", legend.spacing.y = unit(-0.2, "cm"))
    dotPlot <- dotPlot + guides(colour = guide_colorbar(order = 1, label.hjust = 1, label.theme = element_text(angle=90, size = legend.size) ), 
                                size   = guide_legend(order = 2, label.theme = element_text(size = legend.size, angle = 0) ))
  } else {
    dotPlot <- dotPlot + theme(legend.direction = "vertical", legend.position = legend.position, legend.box = "vertical", legend.spacing.y = unit(0, "cm"))
    dotPlot <- dotPlot + guides(colour = guide_colorbar(order = 1, label.hjust = 0.5, label.theme = element_text(angle=0, size = legend.size) ), 
                                size   = guide_legend(order = 2, label.theme = element_text(size = legend.size, angle = 0) ))
  }
  
  dotPlot
}

# ---
# 00. used for subset topGO and make corresponding pie chart, topNo is for BP, CC, and MF by order, e.g. (c 10, 3, 5)
subsetTopGo <- function(gseaRes, topNo) {
  # ---
  OntNo <- gseaRes %>% group_by(ONTOLOGY) %>% count(ONTOLOGY) 
  OntNo <- as.data.frame(OntNo) %>% mutate('per' = round(x = n*100/sum(n), digits = 0))
  # print(OntNo)
  gseaRes$`bgGeneNo` <- sapply(gseaRes$BgRatio, function(x) as.numeric(unlist(strsplit(x, split = '/'))[1]))
  gseaRes$`grNo`     <- sapply(gseaRes$GeneRatio, function(x) as.numeric(unlist(strsplit(x, split = '/'))[1]))
  gseaRes$`note`     <- sapply(1:length(gseaRes$GeneRatio), function(x) paste(unlist(strsplit(gseaRes$GeneRatio[x], split = '/'))[1], unlist(strsplit(gseaRes$BgRatio[x], split = '/'))[1], sep = '/') )
  # print(head(gseaRes[,c(1:8,12:14)]))
  # ---
  if (length(topNo)==1) {
    OntTop   <- gseaRes %>% group_by(ONTOLOGY) %>% dplyr::filter(between(row_number(p.adjust), 1, topNo)) %>% arrange(row_number(ONTOLOGY))
  } else {
    OntTopBP <- gseaRes %>% dplyr::filter(ONTOLOGY=='BP') %>% dplyr::filter(between(row_number(p.adjust), 1, topNo[1]))
    OntTopCC <- gseaRes %>% dplyr::filter(ONTOLOGY=='CC') %>% dplyr::filter(between(row_number(p.adjust), 1, topNo[2]))
    OntTopMF <- gseaRes %>% dplyr::filter(ONTOLOGY=='MF') %>% dplyr::filter(between(row_number(p.adjust), 1, topNo[3]))
    OntTop   <- dplyr::bind_rows(OntTopBP, OntTopCC, OntTopMF) 
  }
  
  # note: grouping doesn't change how the data looks (apart from listing, how it's grouped)
  # print(OntTop[,c(1:2,4:8,12:14)])
  # ---
  return(list(OntCategoryNo = OntNo, OntTop = OntTop))
}
# 1. make multilayer piechart
# 1.0 function to make multi-layer piechart
makeMultilayerPie <- function(gseaTop, pieChartName, pieWidth, pieHeight) {
  innerLayer     <- gseaTop[[1]]
  innerLayerPer  <- innerLayer$per
  outerLayerData <- as.data.frame(gseaTop[[2]]) %>% 
    group_by(ONTOLOGY) %>% 
    mutate('grNoPer'=round(x = grNo*100/sum(grNo), digits = 2)) %>% 
    arrange(row_number(ONTOLOGY))
  outerLayerColNo <- outerLayerData %>% group_by(ONTOLOGY) %>% count(ONTOLOGY)
  for ( i in 1:length(innerLayer$ONTOLOGY) ) {
    # ---
    catPer = round(x = innerLayer$per[i] * outerLayerData$grNoPer[outerLayerData$ONTOLOGY == as.character(innerLayer$ONTOLOGY[i])] /100, digits = 1)
    if (i ==1) {
      outerLayerPer <- catPer
    } else {
      outerLayerPer <- c(outerLayerPer, catPer)
    }
    # ---
  }
  # start to make the pie-chart
  if(length(gseaTop$OntCategoryNo$ONTOLOGY)==3) {
    pdf(file = pieChartName, width = pieWidth, height = pieHeight)
    library(plotrix)
    iniR=0.2 # initial radius
    colorsInner = list(BP='#800000',CC='green',MF='blue')
    # from outer circle to inner circle
    # 0 circle: blank
    pie(1, radius=iniR, init.angle=90, col=c('white'), border = NA, labels='')
    # outter circle: show top 30 GO cat
    library(dichromat)
    colfuncBP <-colorRampPalette(c("#B22222", "#E46A6A"))
    colfuncCC <-colorRampPalette(c("#006600", "#CCFFCC"))
    colfuncMF <-colorRampPalette(c("#0066CC", "#CCE5FF"))
    outterPie <- floating.pie(0, 0, outerLayerPer, radius=4*iniR, startpos=pi/2, 
                              col=c(colfuncBP(outerLayerColNo$n[outerLayerColNo$ONTOLOGY=='BP']), 
                                    colfuncCC(outerLayerColNo$n[outerLayerColNo$ONTOLOGY=='CC']), 
                                    colfuncMF(outerLayerColNo$n[outerLayerColNo$ONTOLOGY=='MF'])), border=NA)
    pie.labels(0, -0.02, outterPie, radius = 4.3*iniR, paste(outerLayerData$Description, outerLayerData$note, sep = ': '), cex = 1.5)
    # inner circle: show ontology generic
    floating.pie(0, 0, innerLayerPer, radius=2*iniR, startpos=pi/2, col=as.character(colorsInner[as.character(innerLayer$ONTOLOGY)]),border=NA)
    legend('topleft', names(colorsInner), col=as.character(colorsInner[names(colorsInner)]), pch=15, bty='n', ncol = 3, adj = c(0, 0.5), cex = 2)
    
    dev.off()
  } else {
    print('Can NOT generate multiple piechart, some GO category has 0 enriched GOs at specified adjusted p-val')
  }
  # Test on color chosen
  # # plot 10 gradient color starting from #990000 to #FFCCCC
  # colfunc<-colorRampPalette(c("#990000", "#E46A6A"))
  # plot(rep(1,3),col=colfunc(3),pch=19,cex=3)
}

# ---
# summarize the category number of BP/CC/MF of enricheGO results
summarizeGoCat <- function(dataCombDfFull) {
  clusterLevels <- names(table(dataCombDfFull$Cluster))
  for ( i in 1:length(clusterLevels)) {
    ontCatFull_i     <- as.data.frame(table((dataCombDfFull %>% dplyr::filter(Cluster == clusterLevels[i]) )$ONTOLOGY))
    ontCatFull_i     <- ontCatFull_i %>% dplyr::mutate_if(is.factor, as.character)
    if (i==1) {
      ontCatFullPrep <- ontCatFull_i
    } else {
      ontCatFullPrep <- dplyr::full_join(x = ontCatFullPrep, y = ontCatFull_i, by = 'Var1') 
    }
  }
  ontCatFull           <- t(ontCatFullPrep[,-1])
  colnames(ontCatFull) <- ontCatFullPrep$Var1
  rownames(ontCatFull) <- clusterLevels
  return(ontCatFull)
}

# ---
# subset int go
queryGoInt <- function(dataCombDfFull, goInts) {
  termIntPrep      <- as.character(goInts)
  termInt          <- goInts
  dataCombDfInt    <- dataCombDfFull %>% filter(grepl(paste(termInt, collapse = '|'), Description))
  dataCombDfIntRes <- distinct(dataCombDfInt, Cluster, ID, .keep_all = T)
  print(sprintf('GO of interest is: %s', paste(termInt, collapse = ', ')))
  print(table(dataCombDfIntRes$Cluster))
  return(dataCombDfIntRes)
}
# ---
## enrichment DF results with 'Cluster' column coverted into list with respect to each 'Cluster'
resDfCluster2list   <- function(cpORenrichmentResClusterDF) {
  clusterNames      <- levels(factor(cpORenrichmentResClusterDF$Cluster))
  resPrep           <- list()
  for (c in 1:length(clusterNames)) {
    resPrep[[c]]    <- cpORenrichmentResClusterDF %>% dplyr::filter(Cluster == as.character(clusterNames[c]))
  }
  names(resPrep)    <- as.character(clusterNames)
  return(resPrep)
}
## ---
## ---------