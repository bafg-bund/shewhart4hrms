
test_that("Directory tree is created", {
  localDir <- file.path(withr::local_tempdir(), "test")
  suppressMessages(newDirTree(localDir))
  
  newDirs <- list.dirs(localDir)
  allFiles <- list.files(localDir, recursive = T)
  
  expect_length(newDirs, 4)
  expect_contains(newDirs, file.path(localDir, "mzXML-data-files"))
  expect_length(allFiles, 7)
  expect_contains(allFiles, "settings/warningLevels.csv")
  

})

test_that("initiate pos produces the correct csv output files", {
  localDir <- setupTestDir("shewartTestDir", parent.frame())
  filePaths <- newFilePaths(localDir, "pos")
  initializeResultsTable(filePaths)
  
  allFiles <- list.files(localDir, r = T, f = T)
  resultsTable <- read.csv(filePaths$results)
  
  expect_length(allFiles, 13)
  expect_equal(resultsTable[1, "intensity"], 448)
  expect_equal(nrow(resultsTable), 1)
  
})

test_that("initiate neg produces the correct csv output files", {
  localDir <- setupTestDir("shewartTestDir", parent.frame())
  filePaths <- newFilePaths(localDir, "neg")
  initializeResultsTable(filePaths)
  
  allFiles <- list.files(localDir, r = T, f = T)
  resultsTable <- read.csv(filePaths$results)
  
  expect_length(allFiles, 13)
  expect_equal(resultsTable[1, "intensity"], 15210)
  expect_equal(nrow(resultsTable), 1)
  
})

# Copyright 2025 Bundesanstalt für Gewässerkunde
# This file is part of shewhart4hrms
