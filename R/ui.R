




fixDates <- function(dates) {
  attr(dates, "tzone") <- "Europe/Berlin"
  c(dates[1]-10, dates[2]+10)
}

makeTrendPlot <- function(df, dates, type) {
  ploti <- ggplot(df, aes_(quote(time), as.name(type))) + 
    geom_point() + 
    geom_line() + 
    facet_grid(rows = vars(IS), scales = "free") +
    scale_x_datetime() + 
    coord_cartesian(xlim = dates) +
    theme_grey(18)
  
  ploti
}


makeBellCurve <- function(df, dates, type, filePaths) {
  
  newDat <- df[df$time >= dates[1] & df$time <= dates[2], ]
  
  newPlot <- ggplot(newDat, aes_(as.name(type))) +
    facet_grid(rows = vars(IS), scales = "free") +
    ylab("Density") +
    geom_density() +
    coord_flip() +
    theme_grey(18)
  
  newPlot
}

makeBellCurveOld <- function(df, dates, type, pol_i, filePaths) {
  dates <- fixDates(dates)
  se <- readSettings(filePaths)
  # minimale und maximale Intensitätsgrenzen
  POS_INT_RANGE_HEIGHT <- as.numeric(se$intensity_range_pos)
  NEG_INT_RANGE_HEIGHT <- as.numeric(se$intensity_range_neg)
  POS_INT_RANGE_AREA <- as.numeric(se$area_range_pos)
  NEG_INT_RANGE_AREA <- as.numeric(se$area_range_neg)
  
  # Grenzen Massenshift, Retentionszeitdifferenz und Peakbreite
  MAX_MASS_SHIFT_MDA <- se$max_mass_shift * 1000
  MAX_RT_SHIFT_MIN <- se$max_rt_shift
  MAX_PEAKWIDTH_MIN <- se$max_peakwidth
  
  newDat <- df[df$time >= dates[1] & df$time <= dates[2], ]
  
  types <- c("int_h", "int_a", "delta_mz_mDa", "delta_rt_min", "peak_width_min")
  means <- vapply(types, function(x) mean(df[, x], na.rm = T, trim = .15), numeric(1))
  stdevs <- vapply(types, function(x) sd(newDat[, x], na.rm = T), numeric(1))
  
  ploti <- ggplot(newDat, aes_(as.name(type))) + geom_density() + ylab("Density") +
    geom_vline(xintercept = newDat[which.max(newDat$time), type]) + 
    geom_vline(xintercept = means[type], color = "blue", alpha = .2) +
    geom_vline(xintercept = means[type] + stdevs[type], color = "red", alpha = .2) +
    geom_vline(xintercept = means[type] - stdevs[type], color = "red", alpha = .2) +
    annotate("text", x = means[type], y = Inf, label = "mean", color = "blue", alpha = .2, hjust = -.1,
             vjust = 1) +
    annotate("text", x = means[type] + stdevs[type], y = Inf, label = "1 sd", color = "red", alpha = .2, 
             hjust = -.1, vjust = 1) +
    annotate("text", x = means[type] - stdevs[type], y = Inf, label = "1 sd", color = "red", alpha = .2,
             hjust = -.1, vjust = 1) +
    annotate("text", x = newDat[which.max(newDat$time), type], y = Inf, label = "latest", 
             color = "black", hjust = -.1, vjust = 1) +
    theme_grey(18)
  
  ploti2 <- ggplot(newDat, aes_(as.name(type), color = quote(IS))) + ylab("Density") +
    geom_density() + theme_grey(18)
  # Warning levels for intensity height
  if (type == "int_h" && (!is.null(POS_INT_RANGE_HEIGHT) || !is.null(NEG_INT_RANGE_HEIGHT))) {
    switch(pol_i,
           pos = {
             ploti <- ploti + 
               geom_vline(xintercept = POS_INT_RANGE_HEIGHT[1], color = "darkred", 
                          alpha = .3) +
               annotate("text", POS_INT_RANGE_HEIGHT[1], Inf, label = "low warning", 
                        color = "darkred", alpha = .3, hjust = -.1, vjust = 1) + 
               geom_vline(xintercept = POS_INT_RANGE_HEIGHT[2], color = "darkred", 
                          alpha = .3) +
               annotate("text", POS_INT_RANGE_HEIGHT[2], Inf, label = "high warning", 
                        color = "darkred", alpha = .3, hjust = -.1, vjust = 1)
           },
           neg = {
             ploti <- ploti + geom_vline(xintercept = NEG_INT_RANGE_HEIGHT[1], color = "darkred", alpha = .3) +
               annotate("text", NEG_INT_RANGE_HEIGHT[1], Inf, label = "low warning", color = "darkred", alpha = .3, 
                        hjust = -.1, vjust = 1) + geom_vline(xintercept = NEG_INT_RANGE_HEIGHT[2], color = "darkred", alpha = .3) +
               annotate("text", NEG_INT_RANGE_HEIGHT[2], Inf, label = "high warning", color = "darkred", alpha = .3, 
                        hjust = -.1, vjust = 1)
           })
    
  } else if (type == "int_a" && (!is.null(POS_INT_RANGE_AREA) || !is.null(NEG_INT_RANGE_AREA))) {
    # Warning levels for intensity area
    
    switch(pol_i,
           pos = {
             ploti <- ploti + 
               geom_vline(xintercept = POS_INT_RANGE_AREA[1], color = "darkred", 
                          alpha = .3) +
               annotate("text", POS_INT_RANGE_AREA[1], Inf, label = "low warning", 
                        color = "darkred", alpha = .3, hjust = -.1, vjust = 1) + 
               geom_vline(xintercept = POS_INT_RANGE_AREA[2], color = "darkred", 
                          alpha = .3) +
               annotate("text", POS_INT_RANGE_AREA[2], Inf, label = "high warning", 
                        color = "darkred", alpha = .3, hjust = -.1, vjust = 1)
           },
           neg = {
             ploti <- ploti + geom_vline(xintercept = NEG_INT_RANGE_AREA[1], color = "darkred", alpha = .3) +
               annotate("text", NEG_INT_RANGE_AREA[1], Inf, label = "low warning", color = "darkred", alpha = .3, 
                        hjust = -.1, vjust = 1) + geom_vline(xintercept = NEG_INT_RANGE_AREA[2], color = "darkred", alpha = .3) +
               annotate("text", NEG_INT_RANGE_AREA[2], Inf, label = "high warning", color = "darkred", alpha = .3, 
                        hjust = -.1, vjust = 1)
           })
    
  } else if (type == "delta_mz_mDa" && !is.null(MAX_MASS_SHIFT_MDA)) {
    # warning level for mass shift in mDa
    
    ploti <- ploti + 
      geom_vline(xintercept = MAX_MASS_SHIFT_MDA, color = "darkred", 
                 alpha = .3) +
      annotate("text", MAX_MASS_SHIFT_MDA, Inf, label = "high warning", 
               color = "darkred", alpha = .3, hjust = -.1, vjust = 1)
  } else if (type == "delta_rt_min" && !is.null(MAX_RT_SHIFT_MIN)) {
    # warning level for RT shift in min
    
    ploti <- ploti + 
      geom_vline(xintercept = MAX_RT_SHIFT_MIN, color = "darkred", 
                 alpha = .3) +
      annotate("text", MAX_RT_SHIFT_MIN, Inf, label = "high warning", 
               color = "darkred", alpha = .3, hjust = -.1, vjust = 1)
  } else if (type == "Width_mDa" && !is.null(MAX_PEAKWIDTH_MIN)) {
    # warning level for peak width
    
    ploti <- ploti + 
      geom_vline(xintercept = MAX_PEAKWIDTH_MIN, color = "darkred", 
                 alpha = .3) +
      annotate("text", MAX_PEAKWIDTH_MIN, Inf, label = "high warning", 
               color = "darkred", alpha = .3, hjust = -.1, vjust = 1)
  } else {
    # no change
  }
  
  
  cowplot::plot_grid(ploti, ploti2, align = "h", nrow = 1)
}