
test_that("We can process a second file after initiation in pos mode", {
  localDir <- setupTestDir("shewartTestDir", parent.frame())
  filePaths <- newFilePaths(localDir, "pos")
  initializeResultsTable(filePaths)
  processNewFiles(filePaths)
  
  allFiles <- list.files(localDir, r = T, f = T)
  resultsTable <- read.csv(filePaths$results)
  
  expect_length(allFiles, 14)
  expect_equal(resultsTable[1, "intensity"], 448)
  expect_equal(resultsTable[2, "intensity"], 529)
  expect_equal(nrow(resultsTable), 2)
})

test_that("We can process a second file after initiation in neg mode", {
  localDir <- setupTestDir("shewartTestDir", parent.frame())
  filePaths <- newFilePaths(localDir, "neg")
  initializeResultsTable(filePaths)
  processNewFiles(filePaths)
  
  allFiles <- list.files(localDir, r = T, f = T)
  resultsTable <- read.csv(filePaths$results)
  
  expect_length(allFiles, 14)
  expect_equal(resultsTable[1, "intensity"], 15210)
  expect_equal(resultsTable[2, "intensity"], 15587)
  expect_equal(nrow(resultsTable), 2)
})

# Copyright 2025 Bundesanstalt für Gewässerkunde
# This file is part of shewhart4hrms
