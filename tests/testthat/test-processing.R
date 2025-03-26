
test_that("We can process a second file after initiation in pos mode", {
  localDir <- setupTestDir("shewartTestDir", parent.frame(), "pos")
  filePaths <- newFilePaths(localDir, "pos")
  processNewFiles(filePaths)
  
  allFiles <- list.files(localDir, r = T, f = T)
  resultsTable <- read.csv(filePaths$results)
  
  expect_length(allFiles, 11)
  expect_equal(resultsTable[1, "int_h"], 448)
  expect_equal(resultsTable[2, "int_h"], 529)
  expect_equal(nrow(resultsTable), 2)
})

test_that("We can process a second file after initiation in neg mode", {
  localDir <- setupTestDir("shewartTestDir", parent.frame(), "neg")
  filePaths <- newFilePaths(localDir, "neg")
  processNewFiles(filePaths)
  
  allFiles <- list.files(localDir, r = T, f = T)
  resultsTable <- read.csv(filePaths$results)
  
  expect_length(allFiles, 11)
  expect_equal(resultsTable[1, "int_h"], 15210)
  expect_equal(resultsTable[2, "int_h"], 15587)
  expect_equal(nrow(resultsTable), 2)
})
