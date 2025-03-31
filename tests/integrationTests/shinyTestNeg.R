
devtools::load_all(".")


source(test_path("helper.R"))
localDir <- setupTestDir("shinyTest", parent.frame())
filePaths <- newFilePaths(localDir, "neg")
initializeResultsTable(filePaths)
viewShewhart(localDir)

# Copyright 2025 Bundesanstalt für Gewässerkunde
# This file is part of shewhart4hrms

