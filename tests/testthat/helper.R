


setupTestDir <- function(nameTestDir, .local_envir) {
  localDir <- file.path(withr::local_tempdir(.local_envir = .local_envir), nameTestDir)
  suppressMessages(newDirTree(localDir))
  dataFiles <- list.files(test_path("fixtures", "hrmsDataFiles"), full.names = T)
  
  newDirForData <- file.path(localDir, "mzXML-data-files")
  for (f in dataFiles) {
    file.copy(f, newDirForData)
  }
  # For different modification time
  dataFilesNew <- list.files(newDirForData, full.names = T)
  dataFilesNew <- dataFilesNew[order(dataFilesNew)]
  measTimes <- stringr::str_extract(basename(dataFilesNew), "\\d{8}")
  mapply(
    function(file, newTime) Sys.setFileTime(file, newTime), 
    dataFilesNew,
    lubridate::ymd(measTimes)
  )
  localDir
}
