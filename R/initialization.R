



  
#' Create directory for shewhart4hrms
#' 
#' Creates a new shewhart4hrms directory tree and example files
#'
#' @param newDir Path to new directory
#'
#' @return TRUE if successful
#' @export
#'
#' @examples shewhart4hrms::newDirTree("~/linsensuppe")
newDirTree <- function(path) {
  if (!dir.exists(path))
    dir.create(path) else stop("You must create a new directory")
  
  dir.create(file.path(path, "mzXML-data-files"))
  dir.create(file.path(path, "results"))
  dir.create(file.path(path, "settings"))
  
  file.copy(
    system.file("example-settings", "is-table-pos.csv", package = "shewhart4hrms"),
    file.path(path, "settings")
  )
  
  file.copy(
    system.file("example-settings", "is-table-neg.csv", package = "shewhart4hrms"),
    file.path(path, "settings")
  )
  
  file.copy(
    system.file("example-settings", "processingSettings.yml", package = "shewhart4hrms"),
    file.path(path, "settings")
  )
  
  file.copy(
    system.file("example-settings", "warningLevels.csv", package = "shewhart4hrms"),
    file.path(path, "settings")
  )
  
  # TEMPORARY Place dummy DB in dir and add path to settings
  file.copy(
    system.file("example-settings", "CSL_olmesartan-d6.db", package = "shewhart4hrms"),
    file.path(path)
  )
  pathDb <- file.path(normalizePath(path, "/"), "CSL_olmesartan-d6.db")
  textToAdd <- c("# Dummy DB location (temporary)", paste0("db_path: ", pathDb))
  settingsPath <- file.path(path, "settings", "processingSettings.yml")
  cat(textToAdd, file = settingsPath, append = TRUE, sep = "\n")
  
  # place bat file in directory for starting the app
  
  pathToR <- normalizePath(file.path(R.home(), "bin", "R.exe"))
  shewhartDirForR <- normalizePath(path, winslash = "/")
  currentWorkDir <- normalizePath(getwd(), winslash = "\\")
  batPath <- file.path(path, sprintf("shewhart4hrms-%s.bat", basename(path)))
  cat(sprintf("cd /d \"%s\"\n\"%s\" -e \"shewhart4hrms::viewShewhart('%s')\"\npause", currentWorkDir, pathToR, shewhartDirForR), 
      file = batPath)
  
  # place readme file in the directory
  readme <- sprintf("
    ---- Using shewhart4hrms ----
    Step 1: Convert your data into mzXML using Proteowizard
    Step 2: Save data into the %s/mzXML-data-files directory
    Step 3: Open %s
    Step 4: Click 'Refresh' to process the new files
    ", shewhartDirForR, basename(batPath)
  )
  
  cat(readme, file = file.path(path, "README.txt"))
  
  message("Edit the IS tables in ", file.path(path, "settings"))
  message("Copy mzXML files to ", file.path(path, "mzXML-data-files"))
  
  invisible(TRUE)
}

#' @export
initializeResultsTable <- function(filePaths) {
  
  numFiles <- length(getDataFiles(filePaths))
  
  if (numFiles == 0) {
    stop(
      "No mzXML files with '", 
      filePaths$fileNamePolPattern,
      "' in name found in ",
      filePaths$mzxmlDataFiles
    )
  }

  processFirstFile(filePaths)
  
  invisible(TRUE)
}

processFirstFile <- function(filePaths) {  # pol <- "pos"
  settings <- readSettings(filePaths)
  
  message("Initiating ", filePaths$pol, " in ", filePaths$rootDir)  
  firstFile <- getFirstFile(filePaths)
  
  shewhart <- Report$new()
  shewhart$addDB(F, settings$db_path)   # Temp add dummy DB
  shewhart$addRawFiles(F, firstFile)
  shewhart$addIS(F, filePaths$isTable)
  shewhart$changeSettings("area_threshold", settings[[paste0("area_threshold_", filePaths$pol)]])
  shewhart$changeSettings("use_int_threshold","area")
  shewhart$changeSettings("EIC_extraction", settings[[paste0("eic_extraction_width_", filePaths$pol)]])
  shewhart$changeSettings("ISmztol", settings[[paste0("mz_tol_", filePaths$pol)]])
  shewhart$changeSettings("ISrttolm", settings[[paste0("rt_tol_", filePaths$pol)]])
  shewhart$changeSettings("pol", filePaths$pol)
  
  shewhart <- tryToProcess(filePaths, shewhart)
  
  # check that IS were found 
  if (nrow(shewhart$ISresults) == 0)
    stop("No IS found in ", filePaths$pol, " mode in file ", firstFile)
  
  shewhart$clearAndSave(F, filePaths$report)
  
  writeResultsTable(filePaths)
  
  message("Completed initialization ", filePaths$pol)
}

getFirstFile <- function(filePaths) {
  dataFiles <- getDataFiles(filePaths)
  fileDetails <- file.info(dataFiles)
  rownames(fileDetails[with(fileDetails, order(as.POSIXct(mtime))), ])[1]
}

# Copyright 2025 Bundesanstalt für Gewässerkunde
# This file is part of shewhart4hrms