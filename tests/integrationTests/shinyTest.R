
devtools::load_all(".")
source(test_path("helper.R"))
localDir <- setupTestDir("shinyTest", parent.frame())
filePaths <- newFilePaths(localDir, "pos")
initializeResultsTable(filePaths)
viewShewhart(localDir)

source(test_path("helper.R"))
localDir <- setupTestDir("shinyTest", parent.frame())
filePaths <- newFilePaths(localDir, "neg")
initializeResultsTable(filePaths)
viewShewhart(localDir)



source(test_path("helper.R"))
localDir <- setupTestDir("shinyTest", parent.frame())
viewShewhart(localDir)



rstudioapi::filesPaneNavigate(localDir)
