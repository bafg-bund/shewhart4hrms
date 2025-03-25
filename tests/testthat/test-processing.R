
test_that("We can process a second file after initiation", {
  localDir <- setupTestDir("shewartTestDir", parent.frame())
  filePaths <- newFilePaths(localDir, "pos")
  processNewFiles(filePaths)
  
  allFiles <- list.files(localDir, r = T, f = T)
  resultsTable <- read.csv(filePaths$results)
  
  expect_length(allFiles, 11)
  expect_equal(resultsTable[1, "int_h"], 448)
  expect_equal(resultsTable[2, "int_h"], 529)
  expect_equal(nrow(resultsTable), 2)
})

