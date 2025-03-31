
devtools::load_all(".")

source(test_path("helper.R"))
localDir <- setupTestDir("shinyTest", parent.frame())
filePaths <- newFilePaths(localDir, "pos")
initializeResultsTable(filePaths)
file.copy(test_path("fixtures", "pos_resultsDummyDataset.csv"), filePaths$results, overwrite = T)
file.copy(test_path("fixtures", "warningLevels.csv"), filePaths$warningLevels, overwrite = T)

viewShewhart(localDir)

# Copyright 2025 Bundesanstalt für Gewässerkunde
# This file is part of shewhart4hrms
