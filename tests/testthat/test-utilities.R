

test_that("You can tell a tree is initiated", {
  localDir <- setupTestDir("shewartTestDir", .local_envir = parent.frame())
  filePaths <- newFilePaths(localDir, polarity = "neg")
  initializeResultsTable(filePaths)
  expect_true(isInitiated(filePaths))
})
