
# Copyright 2020-2022 Bundesanstalt für Gewässerkunde
# This file is part of shewhart4hrms
# shewhart4hrms is free software: you can redistribute it and/or modify it under the 
# terms of the GNU General Public License as published by the Free Software 
# Foundation, either version 3 of the License, or (at your option) any 
# later version.
# 
# shewhart4hrms is distributed in the hope that it will be useful, but WITHOUT ANY 
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along 
# with shewhart4hrms. If not, see <https://www.gnu.org/licenses/>.


library(shiny)
library(ggplot2)

# Setup ####
#browser()
# Set wd to shewhart folder...
if (!exists("globalwd"))
  stop("Path to shewhart4hrms directory not found, use shewhart4hrms::create_dir
       to create a new folder")

if (!all(c("mzXML-data-files", "results", "settings") %in% list.dirs(globalwd, full.names = F)))
  stop("Incorrect directory in ", globalwd)

if (!file.exists(file.path(globalwd, "settings", "settings.yml")))
  stop("No settings file found in ", file.path(globalwd, "settings"))

# Define constants ####
se <- yaml::read_yaml(file.path(globalwd, "settings", "settings.yml"))
#se <- yaml::read_yaml("~/rprojects/shewhart4hrms/inst/example-settings/settings.yml")

# minimale und maximale Intensitätsgrenzen
POS_INT_RANGE_HEIGHT <- as.numeric(se$intensity_range_pos)
NEG_INT_RANGE_HEIGHT <- as.numeric(se$intensity_range_neg)
POS_INT_RANGE_AREA <- as.numeric(se$area_range_pos)
NEG_INT_RANGE_AREA <- as.numeric(se$area_range_neg)

# Grenzen Massenshift, Retentionszeitdifferenz und Peakbreite
MAX_MASS_SHIFT_MDA <- se$max_mass_shift * 1000
MAX_RT_SHIFT_MIN <- se$max_rt_shift
MAX_PEAKWIDTH_MIN <- se$max_peakwidth

# location of files
KOKAP <- file.path(globalwd, "results", "shewhart-pos.Report")
KOKAN <- file.path(globalwd, "results", "shewhart-neg.Report")
POS_RES <- file.path(globalwd, "results", "pos_res.csv")
NEG_RES <- file.path(globalwd, "results", "neg_res.csv")
POS_RES_GRAPH <- file.path(globalwd, "results", "pos_res_graph.csv")
NEG_RES_GRAPH <- file.path(globalwd, "results", "neg_res_graph.csv")
IS_TABLE_POS <- file.path(globalwd, "settings", "is-table-pos.csv")
IS_TABLE_NEG <- file.path(globalwd, "settings", "is-table-neg.csv")
MZXML_DATA_FILES <- file.path(globalwd, "mzXML-data-files")
LOGFILE <- file.path(globalwd, "log.txt")

if (!all(file.exists(POS_RES, POS_RES_GRAPH, IS_TABLE_POS)))
  stop("No data found, first run shewhart4hrms::intiate to process the first file.")

# Define functions ####
# get files and info
getDat <- function(graphCSV) {
  #browser()
  
  dat <- read.csv(graphCSV)
  dat$time <- as.POSIXct(dat$time, origin = "1970-01-01 00:00.00 UTC", tz = "Europe/Berlin")
  
  dat
}

fixDates <- function(dates) {
  attr(dates, "tzone") <- "Europe/Berlin"
  c(dates[1]-10, dates[2]+10)
}

makeTrendPlot <- function(df, dates, type) {
  ploti <- ggplot(df, aes_(quote(time), as.name(type))) + 
    geom_point() + geom_line() + facet_wrap(~ IS, scales = "free") +
    scale_x_datetime(limits = dates) + theme_grey(18)
  if (type == "int_h" || type == "int_a") 
    ploti <- ploti + scale_y_continuous(labels = function(x) format(x, digits = 1, scientific = TRUE))
  ploti
}

makeBellCurve <- function(df, dates, type, pol_i) {
  
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

# Load data ####

dat1 <- getDat(POS_RES_GRAPH)

# App ####
## ui ####
ui <- fluidPage(
  fluidRow(
    column(3, h2(paste("shewhart4hrms", basename(globalwd)))),
    column(6, sliderInput("dates", "", timeFormat = "%F", value = c(min(dat1$time), max(dat1$time)),
                     min = min(dat1$time), max = max(dat1$time), width = "100%")),
    column(
      1, 
      radioButtons("pol", "", choices = c("pos", "neg"), selected = "pos", inline = F),
      offset = 1
    ),
    column(1, style = "margin-top: 25px;", actionButton("newFiles", "Refresh"))
  ),
  fluidRow(column(12,
      tabsetPanel(
        tabPanel("Intensity_cps", 
                 fluidRow(plotOutput("inten")), fluidRow(plotOutput("intenHist"))),
        tabPanel("Intensity_area",
                 fluidRow(plotOutput("intenA")), fluidRow(plotOutput("intenAHist"))),
        tabPanel("Mass shift mDa", 
                 fluidRow(plotOutput("mDa")), fluidRow(plotOutput("mDaHist"))
                 ),
        tabPanel("RT shift min", 
                 fluidRow(plotOutput("RT")), fluidRow(plotOutput("RTHist"))
                 ),
        tabPanel("Width min", 
                 fluidRow(plotOutput("width")), fluidRow(plotOutput("widthHist"))
                 ))
    )
  )
)

## server ####
server <- function(input, output, session) {
  
  dat <- reactive({
    df <- read.csv(switch(input$pol, pos = POS_RES_GRAPH, neg = NEG_RES_GRAPH))
    df$time <- as.POSIXct(df$time, origin = "1970-01-01 00:00.00 UTC", tz = "Europe/Berlin")
    df
  })
  
  observeEvent(input$pol, {
    updateSliderInput(session, "dates", min = min(dat()$time), max = max(dat()$time),
                     value = c(min(dat()$time), max(dat()$time)))
  })
  
  observeEvent(input$newFiles, {
    
    prog <- Progress$new(session, min = 0, max = 3)
    prog$set(message = "Processing new IS files", detail = "looking for new files")
    #withProgress(message = 'Processing new IS files', value = 0, {
     
    # get files to process from caller

    
    # check that all files have either pos or neg in name
    filesOk <- grepl("pos|neg", list.files(MZXML_DATA_FILES))
    if (!all(filesOk)) {
      showNotification(
        paste(
          "All files must contain 'pos' or 'neg' in the name",
          "The following files are not valid:",
          paste(list.files(MZXML_DATA_FILES)[!filesOk], collapse = ", ")
        ), 
        type = "error"
      )
      prog$close()
    }
    req(all(filesOk))
    
    new_files <- function(pola) {
      polRes <- read.csv(switch(pola, pos = POS_RES, neg = NEG_RES))
      allFilesPol <- list.files(MZXML_DATA_FILES, pattern = pola, full.names = T)
      allFilesPol[!is.element(basename(allFilesPol), polRes$samp)]
    }
    
    process_new_files <- function(pola) {
      
      polRes <- read.csv(switch(pola, pos = POS_RES, neg = NEG_RES))
      koka <- switch(
        pola,
        pos = shewhart4hrms::loadReport(F, KOKAP),
        neg = shewhart4hrms::loadReport(F, KOKAN)
      )
      
      # Log info
      cat("\n", file = LOGFILE, append = T)
      cat(date(), file = LOGFILE, append = T)
      cat(" Processing: ", file = LOGFILE, append = T)
      cat(paste(new_files(pola), collapse = ", "), file = LOGFILE, append = T)
      cat("\n", file = LOGFILE, append = T)
      
      rawPaths <- new_files(pola)
      # check that these files are new
      if (all(is.element(basename(rawPaths), basename(koka$rawFiles)))) {
        warning("Error in finding new files, this should not happen, possibly due to corrupt files")
        return(FALSE)
      }
      # process files
      koka$addRawFiles(F, rawPaths)
      koka$process_all()
      koka$clearAndSave(F, switch(pola, pos = KOKAP, neg = KOKAN))
      
      # get new part of IS results
      IStab <- koka$ISresults
      IStab <- IStab[!is.element(IStab$samp, polRes$samp), ]
      new_mtime <- data.frame(
        samp = basename(rawPaths), 
        mtime = as.integer(file.mtime(rawPaths)), 
        stringsAsFactors = F
      )
      IStab <- merge(IStab, new_mtime, by = "samp")
      polRes <- rbind(polRes, IStab)
      # add this to results table
      write.csv(polRes, file = switch(pola, pos = POS_RES, neg = NEG_RES), row.names = F)
      
      # compute delta mz and delta rt
      isPol <- switch(pola, pos = read.csv(IS_TABLE_POS), neg = read.csv(IS_TABLE_NEG))
      isPol$mass <- switch(
        pola,
        pos = mapply(shewhart4hrms::get_mass, isPol$formula, isPol$adduct, MoreArgs = list(charge = 1)),
        neg = mapply(shewhart4hrms::get_mass, isPol$formula, isPol$adduct, MoreArgs = list(charge = -1))
      )
 
      isPol <- isPol[, c("name", "mass", "rt")]
      colnames(isPol) <- c("IS", "mz_calc", "rt_known")
      polRes <- merge(polRes, isPol)
      polRes$delta_mz_mDa <- abs(polRes$mz - polRes$mz_calc) * 1000 
      polRes$delta_rt_min <- abs(polRes$rt - polRes$rt_known)
      
      # compute peak width
      polRes$peak_width_min <- polRes$peak_end - polRes$peak_start
      
      polRes$time <- as.POSIXct(polRes$mtime, origin = "1970-01-01 00:00.00 UTC", tz = "Europe/Berlin")
      
      polResGraph <- switch(pola, pos = read.csv(POS_RES_GRAPH), neg = read.csv(NEG_RES_GRAPH))
      polResGraph <- polResGraph[!is.element(polResGraph$samp, polRes$samp), ]
      polResGraph <- rbind(polResGraph, polRes)
      write.csv(
        polResGraph, 
        file = switch(pola, pos = POS_RES_GRAPH, neg = NEG_RES_GRAPH), 
        row.names = F
      )
      
      cat(paste0("Completed ", pola, "\n"), file = LOGFILE, append = T)
      TRUE
    }

    if (length(new_files("pos")) > 0) {
      prog$set(detail = "Processing pos", value = 1)
      toProc <- paste(new_files("pos"), collapse = ", ")
      succPos <- process_new_files("pos")
      if (succPos) {
        showNotification(
          paste("Update successful for", toProc)
        )
      }
        
    } else {
      showNotification("No new pos files found.")
    }
    
    if (length(new_files("neg")) > 0) {
      prog$set(detail = "Processing neg", value = 2) 
      toProc <- paste(new_files("neg"), collapse = ", ")
      succ <- process_new_files("neg")
      if (succ) {
        showNotification(
          paste("Update successful for", toProc)
        )
      }
    } else {
      showNotification("No new neg files found.")
    }
    
    
    prog$set(detail = "Complete", value = 3)
    
    # only way to invalidate plots that was found to work is force change of pol
    currentPol <- isolate(input$pol)
    otherPol <- ifelse(currentPol == "pos", "neg", "pos")
    updateRadioButtons(session, "pol", selected = otherPol)
    updateRadioButtons(session, "pol", selected = currentPol)
    

    prog$close()
   
  })
  
  
  
  output$inten <- renderPlot(makeTrendPlot(dat(), fixDates(input$dates), "int_h"))
  output$intenHist <- renderPlot(makeBellCurve(dat(), fixDates(input$dates), "int_h", 
                                              input$pol))
  
  output$intenA <- renderPlot(makeTrendPlot(dat(), fixDates(input$dates), "int_a"))
  output$intenAHist <- renderPlot(makeBellCurve(dat(), fixDates(input$dates), "int_a", 
                                                 input$pol))
  
  output$mDa <- renderPlot(makeTrendPlot(dat(), fixDates(input$dates), "delta_mz_mDa"))
  output$mDaHist <- renderPlot(makeBellCurve(dat(), fixDates(input$dates), "delta_mz_mDa"))
  
  output$RT <- renderPlot(makeTrendPlot(dat(), fixDates(input$dates), "delta_rt_min"))
  output$RTHist <- renderPlot(makeBellCurve(dat(), fixDates(input$dates), "delta_rt_min"))
  
  output$width <- renderPlot(makeTrendPlot(dat(), fixDates(input$dates), "peak_width_min"))
  output$widthHist <- renderPlot(makeBellCurve(dat(), fixDates(input$dates), "peak_width_min"))
}

shinyApp(ui = ui, server = server)

