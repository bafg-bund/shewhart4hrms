
#' Start shewhart4hrms visualizer
#' 
#' @export
viewShewhart <- function(path) {
  appDir <- system.file("shiny", "shewhart4hrms-visualizer", package = "shewhart4hrms")
  if (appDir == "") {
    stop("Could not find app directory. Try re-installing 'shewhart4hrms'.", call. = FALSE)
  }
  globalwd <<- path
  shiny::runApp(appDir, display.mode = "normal", launch.browser = TRUE)
}

#' @export
newFilePaths <- function(rootPath, polarity) {
  x <- list()
  x[["mzxmlDataFiles"]] <- file.path(rootPath, "mzXML-data-files")
  x[["logFile"]] <- file.path(rootPath, "log.txt")
  x[["settings"]] <- file.path(rootPath, "settings", "processingSettings.yml")
  x[["warningLevels"]] <- file.path(rootPath, "settings", "warningLevels.csv")
  x[["rootDir"]] <- rootPath
  x[["report"]] <- file.path(rootPath, "results", glue::glue("shewhart-{polarity}.Report"))
  x[["results"]] <- file.path(rootPath, "results", glue::glue("{polarity}_results.csv"))
  x[["isTable"]] <- file.path(rootPath, "settings", glue::glue("is-table-{polarity}.csv"))
  x[["pol"]] <- polarity
  x[["fileNamePolPattern"]] <- getFileNamePolPattern(x, polarity)
  structure(x, class = c("FilePaths", "list"))
}

getFileNamePolPattern <- function(filePaths, polarity) {
  se <- readSettings(filePaths)
  switch(polarity, pos = se$file_pattern_pos, neg = se$file_pattern_neg)
}

#' @export
checkShewartDir <- function(shewhartPath, polarity = "pos") {
  testFilePaths <- newFilePaths(shewhartPath, polarity)
  if (!all(c("mzXML-data-files", "results", "settings") %in% list.dirs(shewhartPath, full.names = F)))
    stop("Incorrect directory in ", globalwd)
  
  if (!file.exists(testFilePaths$settings))
    stop("No settings file found in ", testFilePaths$settings)
  
  if (!file.exists(testFilePaths$isTable))
    stop("No is-table file found in ", testFilePaths$isTable)
}



checkFileNamesHavePol <- function(filePaths) {
  filesOk <- polInFileName(filePaths)
  if (!all(filesOk)) {
    stop(makePolPatternErrorText(filePaths, filesOk))
  }
}

#' @export
checkFileNamesHavePolShiny <- function(filePaths, prog) {
  filesOk <- polInFileName(filePaths)
  if (!any(filesOk)) {
    showNotification(
      makePolPatternErrorText(filePaths, filesOk), 
      type = "error"
    )
    prog$close()
  }
  filesHavePol <- polInFileName(filePaths)
  req(any(filesHavePol))
}

polInFileName <- function(filePaths) {
  grepl(filePaths$fileNamePolPattern, list.files(filePaths$mzxmlDataFiles))
}

makePolPatternErrorText <- function(filePaths, filesOk) {
  paste(
    "Some files must contain '",filePaths$fileNamePolPattern, "' in the name",
    "The following files are not valid:",
    paste(list.files(filePaths$mzxmlDataFiles)[!filesOk], collapse = ", ")
  )
}


readSettings <- function(filePaths) {
  yaml::read_yaml(filePaths$settings)
}

getDataFiles <- function(filePaths) {
  namePatternSetting <- filePaths$fileNamePolPattern 
  dataFiles <- list.files(
    filePaths$mzxmlDataFiles, 
    full.names = TRUE,
    pattern = paste0(namePatternSetting, "\\.mzX?ML$")
  )
}

addLineToLogFile <- function(filePaths, messageText, logLevel = INFO) {
  log_appender(appender_file(filePaths$logFile))
  log_level(logLevel, messageText)
}

#' @export
getResults <- function(filePaths) {
  results <- read.csv(filePaths$results)
  results$time <- as.POSIXct(results$time, origin = "1970-01-01 00:00.00 UTC", tz = "Europe/Berlin")
  results
}

#' @export
isInitiated <- function(filePaths) {
  if (file.exists(filePaths$results))
    nrow(read.csv(filePaths$results)) > 0 else FALSE
}

# Copyright 2025 Bundesanstalt für Gewässerkunde
# This file is part of shewhart4hrms