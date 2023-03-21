## fns used for ''
## --------------------------------------------------------------------------- ##
# resLists2DfCluster <- function(listRes, topNo, adjpVal, padjon) {
#   if (missing(padjon)) padjon <- as.logical(T) 
#   if (missing(topNo) & !missing(adjpVal)) {
#     listResCluster <- lapply(1:length(listRes), function(x) {
#       listRes    <- dplyr::mutate(listRes[[x]], 'Cluster' = names(listRes)[x])
#       if (padjon) {
#         listResSel <- listRes %>% dplyr::filter(p.adjust<=adjpVal[x])
#       } else {
#         listResSel <- listRes %>% dplyr::filter(pvalue<=adjpVal[x])
#       }
#     })
#   } else if (!missing(topNo) & missing(adjpVal)) {
#     listResCluster <- lapply(1:length(listRes), function(x) {
#       listRes    <- dplyr::mutate(listRes[[x]], 'Cluster' = names(listRes)[x])
#       if (padjon) {
#         listResTop <- listRes %>% dplyr::filter(between(row_number(p.adjust), 1, topNo))
#       } else {
#         listResTop <- listRes %>% dplyr::filter(between(row_number(pvalue), 1, topNo))
#       }
#       
#     })
#   } else if (missing(topNo) & missing(adjpVal)) {
#     listResCluster <- lapply(1:length(listRes), function(x) {
#       listRes    <- dplyr::mutate(listRes[[x]], 'Cluster' = names(listRes)[x])
#     })
#   }
#   resListComb <- dplyr::bind_rows(listResCluster)
#   resListComb
# }
# ---
# 0.2
makeDotPlot <- function(res2plot, padjon, high.limit = NULL, xLevelsOrder = NULL) {
  # heightval = length((res2plot$res)$Description) * 10
  res2plot$GeneRatio2 <- as.numeric(as.character(res2plot$Count))/as.numeric(as.character(sapply(res2plot$GeneRatio, function(x) strsplit(as.character(x), split = '/')[[1]][2]))) 
  if (is.null(xLevelsOrder)) {
    xLevels <- levels(factor(res2plot$Cluster))
  } else {
    xLevels <- xLevelsOrder
  }

  # yLevels <- unique(factor(res2plot$ID))
  yLevels <- unique(factor(res2plot$Description))
  if (padjon) {
    res2plot$log10P <- -log10(res2plot$p.adjust)
  } else {
    res2plot$log10P <- -log10(res2plot$pvalue)
  }
  print(head(res2plot))
  print(table(res2plot$Cluster))
  print('0000000000')
  if (length(unique(res2plot$Cluster))>1) {
    dotPlot <- ggplot(res2plot, aes(x=factor(Cluster, levels = as.character(xLevels)), y=factor(Description, levels = rev(yLevels) )))
  } else {
    dotPlot <- ggplot(res2plot, aes(x=log10P, y=factor(Description, levels = rev(yLevels) )))
  }
  
  if (padjon) {
    dotPlot <- dotPlot + geom_point(aes(size = as.numeric(Count), colour = as.numeric(as.character(p.adjust))), show.legend = T ) + theme_bw()
    if(is.null(high.limit)) high.limit = max(res2plot$p.adjust) 
  }else {
    dotPlot <- dotPlot + geom_point(aes(size = as.numeric(Count), colour = as.numeric(as.character(pvalue))), show.legend = T ) + theme_bw()
    if(is.null(high.limit)) high.limit = max(res2plot$pvalue) 
  }
  dotPlot <- dotPlot + scale_x_discrete('Cluster', labels = xLevelsOrder, breaks = Cluster)
  dotPlot <- dotPlot + theme(axis.text.y = element_text(colour = "black", size = 14))
  if (length(unique(res2plot$Cluster))>1) {
    dotPlot <- dotPlot + labs(x = "", y = "")
    dotPlot <- dotPlot + theme(axis.text.x = element_text(colour = "black", size = 14, angle = 45, hjust = 1))
  } else {
    if (padjon) {
      dotPlot <- dotPlot + labs(x = "-Log10 (FDR)", y = "")
    } else {
      dotPlot <- dotPlot + labs(x = "-Log10 (nominal p-value)", y = "")
    }
    dotPlot <- dotPlot + theme(axis.text.x = element_text(colour = "black", size = 14, angle = 0, hjust = 0))
  }
  
  # dotPlot <- dotPlot + labs(title = "Selected over represented GOs") 
  dotPlot <- dotPlot + theme(plot.title = element_text(size = 10, hjust = 0.5),
                             axis.title.x = element_text(size = 20))
  if (padjon) {
    # dotPlot <- dotPlot + scale_colour_gradient(name = "adjusted \n p-value\n")
    dotPlot <- dotPlot + scale_colour_gradient(name = "FDR",
                                         limits=c(0, high.limit), low="red", high="blue")
  }else {
    # dotPlot <- dotPlot + scale_colour_gradient(name = "original \n p-value\n")
    dotPlot <- dotPlot + scale_colour_gradient(name = "original \n p-value\n",
                                         limits=c(0, high.limit), low="red", high="blue")
  }
  dotPlot <- dotPlot + scale_size(name = "Enriched Gene Number")
  # dotPlot <- dotPlot + scale_size(name = "Gene Ratio", range = c(1,10), breaks = c(0.25, 0.50, 0.75, 1))
  # dotPlot <- dotPlot + scale_size(name = "count", range = c(min(res.pcut$Count), max(res.pcut$Count)),
  #                     breaks=c(min(res.pcut$Count), 5, 10, max(res.pcut$Count)))
  dotPlot <- dotPlot + theme(plot.margin = unit(c(0.5,0.01,0.01,0.01), "cm") ) ##top, right, bottom, left margins
  dotPlot <- dotPlot + theme(legend.text=element_text(size=14), legend.title = element_text(size=14) )
  # dotPlot <- dotPlot + theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "vertical", legend.spacing.y = unit(-0.2, "cm"))
  dotPlot <- dotPlot + guides(colour = guide_colorbar(order = 1, label.hjust = 0.5, label.theme = element_text(angle=0, size = 20) ), 
                              size   = guide_legend(order = 2, label.theme = element_text(size = 14, angle = 0) ))
  dotPlot <- dotPlot + theme(legend.title=element_text(size=16))
  dotPlot
}
## ---------------------------------------------------------------------------