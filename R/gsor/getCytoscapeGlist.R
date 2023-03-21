# This script is used to deal with cytoscap results to 1) extract essential gene lists; and 2) extract module gene lists.
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

getEssentialGlist <- function(cytoResCsvFname, degreeCutoff, refGenome4cpenrichGo) {
  cytoRes      <- read.csv(file = cytoResCsvFname, header = T, check.names = F, stringsAsFactors = F)
  print('---')
  print( table(cytoRes$edgeDegree) )
  print('***')
  cytoRes$edgeDegreeUpdate = as.numeric(sapply(strsplit(x = as.character(cytoRes$edgeDegree), split = '_'), '[[', 1))
  essGlistPrep <- cytoRes %>% dplyr::filter(edgeDegreeUpdate>=degreeCutoff)
  essGlistOrg  <- essGlistPrep %>% distinct(name) %>% pull()
  print(sprintf('%s essential genes from \'%s\' above %s edge degree', length(essGlistOrg), basename(cytoResCsvFname), as.numeric(degreeCutoff)))
  if (refGenome4cpenrichGo == 'mouse') {
    essGlist   <- firstup(tolower(essGlistOrg))
  } else if (refGenome4cpenrichGo == 'human') {
    essGlist   <- essGlistOrg
  } else if (refGenome4cpenrichGo == 'rat') {
    essGlist   <- tolower(essGlistOrg)
  } else {
    essGlist   <- essGlistOrg
  }
  return(essGlist)
}

getModuleGlist <- function(cytoResCsvFname, moduleCutoff) {
  cytoRes      <- read.csv(file = cytoResCsvFname, header = T, check.names = F, stringsAsFactors = F)
  print('---')
  # print( table(cytoRes$module) )
  # print('***')
  moduleGlistPrep <- cytoRes %>% dplyr::filter(module<=moduleCutoff)
  moduleGlistSep  <- moduleGlistPrep %>% group_by(module) %>% do(data = (.))
  names(moduleGlistSep$data) <- paste('Module', moduleGlistSep$module, sep = '_')
  print(sapply((moduleGlistSep$data), function(x) dim(x)[1]))
  essGlistOrg        <- lapply(moduleGlistSep$data, function(x) x %>% distinct(name) %>% pull() )
  # print(sprintf('%s essential genes from %s.', sapply(essGlistOrg, length), names(essGlistOrg)))
  if (refGenome4cpenrichGo == 'mouse') {
    essGlist   <- lapply(essGlistOrg, function(x) firstup(tolower(x)) )
  } else if (refGenome4cpenrichGo == 'human') {
    essGlist   <- essGlistOrg
  } else if (refGenome4cpenrichGo == 'rat') {
    essGlist   <- lapply(essGlistOrg, function(x) tolower(x) ) 
  } else {
    essGlist   <- essGlistOrg
  }
  print('***')
  return(essGlist)
}
