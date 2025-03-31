

makeTabSetPanel <- function(plotHeightSpec) {
  tabsetPanel(
    tabPanel(
      "intensity",
      column(8, plotOutput("intensityTimeSeries", height = plotHeightSpec)),
      column(4, plotOutput("intensityBellCurve", height = plotHeightSpec))
    ),
    tabPanel(
      "area",
      column(8, plotOutput("areaTimeSeries", height = plotHeightSpec)),
      column(4, plotOutput("areaBellCurve", height = plotHeightSpec))
    ),
    tabPanel(
      "delta_mz_mDa",
      column(8, plotOutput("delta_mz_mDaTimeSeries", height = plotHeightSpec)),
      column(4, plotOutput("delta_mz_mDaBellCurve", height = plotHeightSpec))
    ),
    tabPanel(
      "delta_rt_min", 
      column(8, plotOutput("delta_rt_minTimeSeries", height = plotHeightSpec)),
      column(4, plotOutput("delta_rt_minBellCurve", height = plotHeightSpec))
    ),
    tabPanel(
      "peak_width_min", 
      column(8, plotOutput("peak_width_minTimeSeries", height = plotHeightSpec)),
      column(4, plotOutput("peak_width_minBellCurve", height = plotHeightSpec))
    )
  )
}


makeTimeSeries <- function(results, warningLevels, dates, type) {
  
  warningLevels <- subset(warningLevels, parameter == type)  # warningLevels is already subsetted for polarity
  ggplot(results, aes(time, !!ensym(type), color = name_is)) + 
    geom_point(size = 3) + 
    geom_line(linewidth = 1) + 
    geom_hline(aes(yintercept = upper_warning), data = warningLevels, color = "red", alpha = 0.7, linewidth = 1, linetype = "dashed") +
    geom_hline(aes(yintercept = lower_warning), data = warningLevels, color = "red", alpha = 0.7, linewidth = 1, linetype = "dashed") +
    facet_grid(rows = vars(name_is), scales = "free") +
    scale_x_datetime() + 
    scale_colour_discrete() +
    coord_cartesian(xlim = dates) +
    theme_grey(18) +
    theme(legend.position = "none")
}


makeHistogram <- function(results, dates, type) {
  latestLines <- subset(results, time == dates[2])
  results <- results[results$time >= dates[1] & results$time <= dates[2], ]
  
  ggplot() +
    facet_grid(rows = vars(name_is), scales = "free") +
    geom_histogram(aes(x = !!ensym(type), fill = name_is), results, binwidth = computeBinwidth) +
    geom_vline(aes(xintercept = !!ensym(type)), latestLines, color = "black", alpha = 0.7, linewidth = 1, linetype = "dashed") +
    scale_colour_discrete() +
    coord_flip() +
    theme_grey(18) +
    ylab("Count (-)") +
    theme(legend.position = "none")
}

computeBinwidth <- function(x) {
  n <- length(x)
  range <- max(x) - min(x)
  k <- ceiling(log2(n) + 1) # Sturges' formula for number of bins
  range / k
}

requireFileExists <- function(path) {
  fileExits <- file.exists(path)
  validate(need(fileExits, paste("No results table found, check", path)))
  req(fileExits)
}
