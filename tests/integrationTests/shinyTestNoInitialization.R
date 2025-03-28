
devtools::load_all(".")

source(test_path("helper.R"))
localDir <- setupTestDir("shinyTest", parent.frame())
viewShewhart(localDir)


