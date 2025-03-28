




makeTrendPlot <- function(results, warningLevels, dates, type) {
  
  # warningLevels is for one polarity only
  
  warningLevels <- subset(warningLevels, parameter == type)
  ploti <- ggplot(results, aes_(quote(time), as.name(type), color = quote(IS))) + 
    geom_point(size = 5) + 
    geom_line(linewidth = 2) + 
    geom_hline(aes(yintercept = upper_warning), data = warningLevels, color = "red", alpha = 0.7, linewidth = 2, linetype = "dashed") +
    geom_hline(aes(yintercept = lower_warning), data = warningLevels, color = "red", alpha = 0.7, linewidth = 2, linetype = "dashed") +
    facet_grid(rows = vars(IS), scales = "free") +
    scale_x_datetime() + 
    scale_colour_discrete() +
    coord_cartesian(xlim = dates) +
    theme_grey(18) +
    theme(legend.position = "none")
  
  ploti
}


makeBellCurve <- function(df, dates, type, filePaths) {
  
  newDat <- df[df$time >= dates[1] & df$time <= dates[2], ]
  
  newPlot <- ggplot(newDat, aes_(as.name(type), color = quote(IS))) +
    facet_grid(rows = vars(IS), scales = "free") +
    scale_colour_discrete() +
    ylab("Density") +
    geom_density(linewidth = 2) +
    coord_flip() +
    theme_grey(18) +
    theme(legend.position = "none")
  
  newPlot
}


