## step 2: cp.enrich.go() in getEnrichGoRes.R
## This step is used to obtain reference genome name and load reference genome to R
#----------
get.orgdb <- function(ref.genome) {
  if (ref.genome == 'human' | ref.genome == 'hs') {
    library(org.Hs.eg.db)
    orgdb <- "org.Hs.eg.db"
  } else if (ref.genome == 'mouse' | ref.genome == 'mm') {
    library(org.Mm.eg.db)
    orgdb <- "org.Mm.eg.db"
  } else if (ref.genome == 'rat' | ref.genome == 'rn') {
    library(org.Rn.eg.db)
    orgdb <- "org.Rn.eg.db"
  } else if (ref.genome == 'worm' | ref.genome == 'ce') {
    library(org.Ce.eg.db)
    orgdb <- "org.Ce.eg.db"
  }
  return(orgdb)
}
#----------