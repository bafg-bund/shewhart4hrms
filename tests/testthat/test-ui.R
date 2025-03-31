

test_that("A trend plot can be made", {
  results <- read.csv(test_path("fixtures", "pos_resultsDummyDataset.csv"))
  results$time <- as.POSIXct(results$time, origin = "1970-01-01 00:00.00 UTC", tz = "Europe/Berlin")
  warningLevels <- read.csv(test_path("fixtures", "warningLevels.csv"))
  dates = c(lubridate::ymd(20230812, tz = "Europe/Berlin"), lubridate::ymd(20230814, tz = "Europe/Berlin"))
  testPlot <- makeTimeSeries(results = results, warningLevels = warningLevels, dates = dates, type = "int_h")
  expect_s3_class(testPlot, "ggplot")
  expect_length(testPlot$layers, 4)
  expect_equal(as.character(testPlot$facet$params$rows), "~name_is")
})

# Copyright 2025 Bundesanstalt für Gewässerkunde
# This file is part of shewhart4hrms