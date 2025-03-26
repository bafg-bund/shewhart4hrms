



  
#' Create directory for shewhart4hrms
#' 
#' Creates a new shewhart4hrms directory tree and example files
#'
#' @param newDir Path to new directory
#'
#' @return TRUE if successful
#' @export
#'
#' @examples shewhart4hrms::create_dir("~/linsensuppe")
newDirTree <- function(path) {
  if (!dir.exists(path))
    dir.create(path) else stop("You must create a new directory")
  
  dir.create(file.path(path, "mzXML-data-files"))
  dir.create(file.path(path, "results"))
  dir.create(file.path(path, "settings"))
  
  # in the settings folder, place IS-Table examples
  file.copy(
    system.file("example-settings", "is-table-pos.csv", package = "shewhart4hrms"),
    file.path(path, "settings")
  )
  
  file.copy(
    system.file("example-settings", "is-table-neg.csv", package = "shewhart4hrms"),
    file.path(path, "settings")
  )
  
  # place an example settings.yml file
  
  file.copy(
    system.file("example-settings", "settings.yml", package = "shewhart4hrms"),
    file.path(path, "settings")
  )
  
  # TEMPORARY Place dummy DB in dir and add path to settings
  file.copy(
    system.file("example-settings", "CSL_olmesartan-d6.db", package = "shewhart4hrms"),
    file.path(path)
  )
  pathDb <- file.path(normalizePath(path, "/"), "CSL_olmesartan-d6.db")
  textToAdd <- c("# Dummy DB location (temporary)", paste0("db_path: ", pathDb))
  settingsPath <- file.path(path, "settings", "settings.yml")
  cat(textToAdd, file = settingsPath, append = TRUE, sep = "\n")
  
  # place bat file in directory for starting the app
  
  ptr <- normalizePath(file.path(R.home(), "bin", "R.exe"))
  pth <- normalizePath(path, winslash = "/")
  batpth <- file.path(path, sprintf("shewhart4hrms-%s.bat", basename(path)))
  cat(sprintf("\"%s\" -e \"shewhart4hrms::viewshewhart('%s')\"\npause", ptr, pth), 
      file = batpth)
  
  # place readme file in the directory
  readme <- sprintf("
    ---- Using shewhart4hrms ----
    Step 1: Convert your data into mzXML using Proteowizard
    Step 2: Save data into the %s/mzXML-data-files directory
    Step 3: Open %s
    Step 4: Click 'Refresh' to process the new files
    ", pth, basename(batpth)
  )
  
  cat(readme, file = file.path(path, "README.txt"))
  
  message("Edit the IS tables in ", file.path(path, "settings"))
  message("Copy mzXML files to ", file.path(path, "mzXML-data-files"))
  
  invisible(TRUE)
}



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
  
  # Temp add dummy DB
  shewhart$addDB(F, settings$db_path)
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



# Copyright 2020-2025 Bundesanstalt für Gewässerkunde
# This file is part of shewhart4hrms