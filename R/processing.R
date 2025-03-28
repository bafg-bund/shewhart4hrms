


processNewFiles <- function(filePaths) {
  
  koka <- ntsworkflow::loadReport(F, filePaths$report)
  rawPaths <- getNewMeasFiles(filePaths)
  
  addLineToLogFile(filePaths, "Processing: ")
  addLineToLogFile(filePaths, paste(rawPaths, collapse = ", "))
  
  koka$addRawFiles(F, rawPaths)
  koka <- tryToProcess(filePaths, koka)
  koka$clearAndSave(F, filePaths$report)
  
  writeResultsTable(filePaths)
  addLineToLogFile(filePaths, "Completed processing new files")
  invisible(TRUE)
}


getNewMeasFiles <- function(filePaths) {
  if (file.exists(filePaths$results)) {
    alreadyProcessed <- read.csv(filePaths$results)$samp
  } else {
    return(logical(0))
  }
  allFiles <- getDataFiles(filePaths)
  allFiles[!is.element(basename(allFiles), alreadyProcessed)]
}

tryToProcess <- function(filePaths, reportObj) {
  tryCatch(
    reportObj$process_all(),
    error = function(cnd) {
      addLineToLogFile(filePaths, "Processing failed", ERROR)
      addLineToLogFile(filePaths, conditionMessage(cnd), ERROR)
      warning(conditionMessage(cnd))
    }
  )
  reportObj
}

writeResultsTable <- function(filePaths) {
  reportObj <- loadReport(F, filePaths$report)
  resTable <- reportObj$ISresults
  resTable <- addMzRtData(filePaths, resTable)
  resTable <- addPeakWidth(resTable)
  resTable <- addMeasTime(filePaths, resTable)
  write.csv(
    resTable, 
    file = filePaths$results, 
    row.names = F
  )
}

addMeasTime <- function(filePaths, resTable) {
  fullPaths <- file.path(filePaths$mzxmlDataFiles, resTable$samp)
  modTimes <- as.integer(file.mtime(fullPaths))
  resTable$time <- as.POSIXct(modTimes, origin = "1970-01-01 00:00.00 UTC", tz = "Europe/Berlin")
  resTable
}

addPeakWidth <- function(resTable) {
  resTable$peak_width_min <- resTable$peak_end - resTable$peak_start
  resTable
}

addMzRtData <- function(filePaths, resTable) {
  resTable <- merge(resTable, makeIsTableForPol(filePaths))
  resTable$delta_mz_mDa <- abs(resTable$mz - resTable$mz_calc) * 1000 
  resTable$delta_rt_min <- abs(resTable$rt - resTable$rt_known)
  resTable
}

makeIsTableForPol <- function(filePaths) {
  isPol <- read.csv(filePaths$isTable)
  isPol$mass <- mapply(ntsworkflow::get_mass, isPol$formula, isPol$adduct, MoreArgs = list(charge = switch(filePaths$pol, pos = 1, neg = -1)))
  
  isPol <- isPol[, c("name", "mass", "rt")]
  colnames(isPol) <- c("IS", "mz_calc", "rt_known")
  isPol
}


