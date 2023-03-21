# This function transfer excel input spreadsheet into lists object in R
excelSheet2list <- function(excelFname, header) {
  library(openxlsx)
  sheetNames <- getSheetNames(excelFname)
  detach("package:openxlsx", unload=TRUE)
  sheetLists <- lapply(sheetNames, function(x) read.xlsx(file = excelFname, header = as.logical(header), sheetName = x, stringsAsFactors = 'F') )
  names(sheetLists) <- sheetNames
  return(sheetLists)
}
# listWrite2excel(list2write, fname2write): write list object into excel spreadsheet with defined fileName
listWrite2excel <- function(list2write, fname2write) {
  for (l in 1:length(list2write)) {
    sheetName    <- names(list2write)[l]
    write.xlsx2(list2write[[l]], fname2write, col.names = T, row.names = T, sheetName = as.character(sheetName), append = l > 1)
  }
}
