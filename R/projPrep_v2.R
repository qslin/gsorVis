# This script includes below functions:
# 0. projPrep(runHPC, workDirPath) either run on Gardner or local, where sysDir and workDir will be saved as global variable
# 'workDirPath' defines the full path (after '/group/bioinformatics/) without prefix on HPC or local 
# ------
# 0. funciton: projPrep(), where patternUpdate is consistent with sysDir (the program running enviroment)
#                        , where definded patternUpdate could be not exsit, if not exsit, the 'workDir' will be created
library(stringr)
projPrep <- function(runHPC, workDirPath, runHPCgroup) {
  print('Start to prepare the project input')
  if (missing(runHPCgroup)) runHPCgroup <- as.logical('T')
  if (runHPCgroup) {
    patternUpdate     <<- c('/group/bioinformatics', '/Volumes/bioinformatics/')
    # ------
    # 0. Define the global environment path on gardner for local/gardner running
    if (runHPC){
      sysDir          <<- '/group/bioinformatics'
    } else {
      sysDir          <<- "/Volumes/bioinformatics"
    }
  } else {
    patternUpdate     <<- c('/gpfs/data/biocore-analysis/yli/', '/Volumes/yli/')
    # ------
    # 0. Define the global environment path on gardner for local/gardner running
    if (runHPC){
      sysDir          <<- '/gpfs/data/biocore-analysis/yli/'
    } else {
      sysDir          <<- "/Volumes/yli/"
    }
  }
  workDirPathUpdata <- str_replace(string = as.character(workDirPath), pattern = paste(as.character(patternUpdate), collapse = '|'), replacement = '')
  workDir           <<- paste(sysDir, workDirPathUpdata, sep = '/')
  if (dir.exists(workDir)) {
    setwd(workDir)
  } else {
    stop("ERROR: the specified 'workDirPath' for 'workDir' does not exist")
  }
  print(sprintf('current work directory is: %s', workDir))
  print('---')
}
# 1. functions: projPathUpdate(): this function must be used after funciton projPrep(), where global environmental variable sysDir and workDir are defined
projPathUpdate <- function(projPath, runHPCgroup, runHPC) {
  if (missing(runHPCgroup)) runHPCgroup <- as.logical('T')
  if (missing(runHPC)) runHPC <- as.logical('F')
  if (runHPCgroup){
    patternUpdate     <<- c('/group/bioinformatics', '/Volumes/bioinformatics/')
    # ------
    # 0. Define the global environment path on gardner for local/gardner running
    if (runHPC){
      sysDir          <<- '/group/bioinformatics'
    } else {
      sysDir          <<- "/Volumes/bioinformatics"
    }
  } else {
    patternUpdate     <<- c('/gpfs/data/biocore-analysis/yli/', '/Volumes/yli/')
    # ------
    # 0. Define the global environment path on gardner for local/gardner running
    if (runHPC){
      sysDir          <<- '/gpfs/data/biocore-analysis/yli/'
    } else {
      sysDir          <<- "/Volumes/yli/"
    }
  }
  projPathUpdata <- str_replace(string = as.character(projPath), pattern = paste(as.character(patternUpdate), collapse = '|'), replacement = paste(as.character(sysDir), '/', sep = ''))
  return(projPathUpdata)
}
# ------
