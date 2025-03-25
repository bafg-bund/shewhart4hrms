
test_that("Directory tree is created", {
  localDir <- file.path(withr::local_tempdir(), "test")
  suppressMessages(newDirTree(localDir))
  
  newDirs <- list.dirs(localDir)
  
  expect_length(newDirs, 4)
  expect_contains(newDirs, file.path(localDir, "mzXML-data-files"))

})

test_that("initiate produces the correct csv output files", {
  localDir <- setupTestDir("shewartTestDir", parent.frame())
  filePaths <- newFilePaths(localDir, "pos")
  
  allFiles <- list.files(localDir, r = T, f = T)
  resultsTable <- read.csv(filePaths$results)
  
  expect_length(allFiles, 10)
  expect_equal(resultsTable[1, "int_h"], 448)
  expect_equal(nrow(resultsTable), 1)
  
})
