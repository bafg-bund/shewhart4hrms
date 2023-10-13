
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


#' A reference class for shewhart4hrms
#'
#' A class for processing mzXML data to extract peak information for 
#' internal standards 
#'
#' @details Initialize a Report by calling Report$new() this starts a new Report with the default settings 
#' (see settings field, below). Settings can be changed by calling the \code{changeSettings} method.
#' Load files by calling the \code{addRawFiles} method and process files with the \code{process_all}
#' method. After processing the Report can be saved for viewing in the visuallization tool (susS app).
#'
#' @field rawFiles A character vector with file paths in the order in which they are to be displayed.
#'
#' @field settings A list with settings for all processing steps.
#' The settings list contains the following fields and defaults:
#' \code{list(rttolm = 1, mztolu = 0.05, mztolu_fine = 0.005, pol = "pos", rt_res = 1.5, EIC_extraction = 0.05, baseline_noise_MS1 = 0.5, sn = 3, rtoffset = 0,
#' IS_rtoffset = 0, ISrttolm = 1, area_threshold = 1, height_threshold = 1, 
#' use_int_threshold = "area"}.
#'
#' \code{rttolm}: retention time tolerance for the suspect search. 
#' \code{mztolu}: m/z tolerance for the spectrum extraction (MS2 precursor mass).
#' \code{mztolu_fine}: m/z tolerance for precursor mass in MS1.  
#' \code{chromatography} Is the chromatographic method allowed from the SDB, use DOI of paper in 
#' which method is described.
#' \code{pol}: polarity. 
#' \code{rt_res}: Resolution of chrom. peaks (see \code{\link{ms2_search}}) in min.
#' Peaks with RT difference > \code{rt_res} are considered to be from different substances. 
#' \code{EIC_extraction}: Extraction width to produce EIC, affects integration. 
#' \code{baseline_noise_MS1}: Signals under this intensity are ignored and also not
#' included in the saved spectra. 
#' \code{sn}: signal-to-noise ratio limit for peak integration. 
#' \code{rtoffset}: Retention time offset between samples and database. 
#' \code{ISrttolm}: RT tolerance for IS peak finding and integration. 
#' \code{rtTolReinteg}: is the retention time tolerance for reintegration.
#' \code{mzTolReinteg}: is the m/z tolerance for reintegration. 
#' \code{area_threshold}: Peak area intensity threshold. 
#' \code{height_threshold}: Peak height intensity threshold.
#' \code{use_int_threshold}: can be either "area", "height", or "none"
#' \code{peaksPerPeak}: 
#' \code{mustFindChromPeak}: default FALSE, if TRUE, peaks without area are deleted
#' 
#' @field rawFilesCompl \code{data.frame} with processed files and date of processing.
#' @field peakList \code{data.frame} with all suspect search results.
#' @field currentPeakID numeric to keep track of peakIDs, do not change.
#' @field IS \code{data.frame} with list of IS to process, this has to be imported from a csv file 
#' (comma sep) using this method \code{addIS()}. The csv file should have the columns name, formula,
#' rt, adduct. Where formula is in the form e.g. for CBZ: C14 13CH12N 15NO, and adduct is in the 
#' form [M+H]+ or [M]+.
#' @field ISresults \code{data.frame} with results of IS processings.
#'
#' @import parallel
#' @import shiny
#' @import ggplot2
#' @import stringr
#' @import dplyr 
#' @return an object of Report class
#' @export Report
#' @exportClass Report
Report <- setRefClass(
  "Report",
  fields = list(
    rawFiles = "character",
    rawData = "list",
    settings = "list",
    rawFilesCompl = "data.frame",
    currentISpeakID = "numeric",
    IS = "data.frame",  # table to keep track of internal standards to be evaluated.
    ISresults = "data.frame"
  ),
  methods = list(
    initialize = function(...) {
      # Standard settings ####
      settings <<- list( 
        rttolm = 1, mztolu = 0.05,
        mztolu_fine = 0.005,
        pol = "pos", 
        rt_res = 0.333, 
        EIC_extraction = 0.02,
        EIC_extraction_range = c(0.01, 0.1),
        baseline_noise_MS1 = 0.6,
        sn = 2, rtoffset = 0, IS_rtoffset = 0, ISrttolm = 1, ISmztol = 0.005,
        rtTolReinteg = 1,  # min rt tolerance for reintegrating peaks
        mzTolReinteg = 0.005,
        mztolu_ms2 = 0.015,
        area_threshold = 1,  # peak area intensity threshold
        height_threshold = 1,  # peak height intensity threshold
        use_int_threshold = "area",
        peaksPerPeak = 10
      )
      
      rawFilesCompl <<- data.frame(path = character(), 
                                   date = character(), 
                                   stringsAsFactors = FALSE)
   
      # start counting for peakID
      currentISpeakID <<- 1L
      IS <<- data.frame(name = character(), formula = character(), rt = numeric(),
                        default = logical(), adduct = character(), 
                        stringsAsFactors = FALSE)
      ISresults <<- data.frame(samp = character(), IS = character(),
                               mz = numeric(), rt = numeric(), int_h = integer(),
                               int_a = integer(), peak_start = numeric(), 
                               peak_end = numeric(),
                               ISpeakID = integer(), 
                               eic_extraction_width = numeric(), 
                               stringsAsFactors = FALSE)
     
      callSuper(...)
    },
    addIS = function(dialog = TRUE, fileName = NULL) {
      "Include a list of internal standards in the report. Dialog indicates use of interactive file
      choosing dialog. The file should be a csv with 4 columns: 'name', 'formula', 'rt', 'adduct'."
      if (dialog)
        fileName <- rstudioapi::selectFile(path = getwd(), filter = "*.csv")
      df <- read.csv(fileName, stringsAsFactors = FALSE)
      # check that the table is ok
      stopifnot(all.equal(colnames(df), c("name", "formula", "rt", "adduct")))
      stopifnot(inherits(df$rt, "numeric"))
      
      test <- rbind(IS, df)
      
      if (any(duplicated(test$name)))
        stop("There are duplicated IS")
      
      IS <<- test
    },
    addRawFilesDir = function(dialog = TRUE, dir_path = NULL) {
      "Add all raw files in the chosen directory"
      if (dialog)
        dir_path <- rstudioapi::selectDirectory(path = "~", caption = "Select directory containing mzXML files")
      
      newFiles <- list.files(dir_path, pattern = "*\\.mzX?ML$", full.names = TRUE)
      # you are not allowed to have samples with the same name.
      testMe <- append(basename(rawFiles), basename(newFiles))
      if (anyDuplicated(testMe)) {
        warning("All file names must be unique, ignoring input")
        dupFiles <- testMe[duplicated(testMe)]
        newFiles <- newFiles[basename(newFiles) != dupFiles]
      }
      # check that files exist
      rawFiles <<- append(rawFiles, newFiles)
    },
    addRawFiles = function(dialog = TRUE, file_list = NULL) {
      "Add raw files, (mzML or mzXML) Dialog indicates use of interactive file choosing dialog.
      Does not load files to RAM, only marks file path."
      if (dialog && .Platform$OS.type == "windows") {
        newFiles <- choose.files(
          default = getwd(),
          filters = matrix(
            data = c("mzML or mzXML files", "*.mzML;*.mzXML"),
            ncol = 2
          )
        )
      } else if (dialog && .Platform$OS.type == "unix") {
        newFiles <- rstudioapi::selectFile(path = "~", filter = "mzXML file (*.mzXML)")
      } else if (!dialog) {
        newFiles <- normalizePath(file_list)
      } else {
        stop("No files selected")
      }
      # you are not allowed to have samples with the same name.
      testMe <- append(basename(rawFiles), basename(newFiles))
      if (anyDuplicated(testMe)) {
        warning("All file names must be unique, ignoring input")
        dupFiles <- testMe[duplicated(testMe)]
        newFiles <- newFiles[basename(newFiles) != dupFiles]
      }
      # check that files exist
      rawFiles <<- append(rawFiles, newFiles)
    },
    remRawFiles = function(indices) {
      "Delete files based on their indices"
      rawFilesCompl <<- rawFilesCompl[rawFilesCompl$path %notin% rawFiles[indices], ]
      ISresults <<- ISresults[!(ISresults$samp %in% basename(rawFiles[indices])), ]
      rawFiles <<- rawFiles[-indices]
    },
    moveRawFile = function(index, direction) {
      "Move a raw file in direction 'up' or 'down'"
      if (length(rawFiles) == 1) return(NULL)
      stopifnot(index %in% seq_along(rawFiles))
      stopifnot(direction %in% c("up", "down"))
      indices <- seq_along(rawFiles)
      endInd <- length(rawFiles)
      if (index == 1 && direction == "up") return(NULL)
      if (index == endInd && direction == "down") return(NULL)
      
      newIndices <- switch(
        direction,
        up = replace(indices, (index-1):index, index:(index-1)),
        down = replace(indices, index:(index+1), (index+1):index))
      
      rawFiles <<- rawFiles[newIndices]
    },
    changeSettings = function(parameter, value) {
      "Change any setting by name and then value. See settings field for more details."
      stopifnot(parameter %in% names(settings))
      # run various tests
      if (parameter == "EIC_extraction_range" && (length(value) != 2 || !is.numeric(value)))
        stop("EIC_extraction_range must be a length 2 numeric vector")
      if (parameter == "pol" && (length(value) != 1 || !(value %in% c("pos", "neg"))))
        stop("Polarity must 'pos' or 'neg'")
      if (parameter == "use_int_threshold" && (length(value) != 1 || !(value %in% c("area", "height"))))
        stop("Polarity must 'area' or 'height'")
      settings[[parameter]] <<- value
    },
    saveSettings = function(path = getwd()) {
      "Save settings to a json file. Name of file is generated from current time."
      jsonStr <- jsonlite::toJSON(settings, pretty = TRUE)
      write(jsonStr, file = file.path(
        path, paste0(format(Sys.time(), "%y%m%d-%H%M"), "_Settings.json"))
      )
    },
    loadSettings = function() {
      "Load previously settings file from current working directory. This will fail if there is more
      than one file present in the directory."
      path <- choose.files(
        default = getwd(),
        filters = matrix(
          data = c("JSON settings files", "*.json"),
          ncol = 2
        )
      )
      settings <<- jsonlite::fromJSON(path)
    },
    loadData = function(all = FALSE, indices = NULL) {
      "Load files that still need to be processed into RAM for fast access, if all = TRUE then
      all files are loaded regardless, if indices is an integer vector, only these samps will
      be loaded"
      #browser()
      to_process <- if (all) {
        rawFiles
      } else if (!is.null(indices)) {
        rawFiles[indices]
      } else {
        setdiff(rawFiles, rawFilesCompl$path)
      }
      
      to_process <- setdiff(to_process, names(rawData))  # only those which are not already loaded
      
      if (length(to_process) == 0)
        return(NULL)
      #browser()
      rawData_temp <- lapply(
        to_process, 
        function(x) try(suppressMessages(xcms::xcmsRaw(x, includeMSn = FALSE)))
      )
      names(rawData_temp) <- to_process
      if (!all(vapply(rawData_temp, inherits, what = "xcmsRaw", logical(1)))) {
        badSamp <- to_process[vapply(rawData_temp, inherits, what = "try-error", logical(1))]
        warning("Error in loading sample(s): ", paste0(basename(badSamp), collape = ", "))
        return(FALSE)
      }
        
      rawData <<- append(rawData, rawData_temp)
      
      message(paste0(basename(to_process), collape = ", "))
      return(TRUE)
    },
    
    clearData = function(indices = NULL) {
      "Remove data from RAM to clear memory, use indices of raw files"
      if (is.null(indices)) {
        rawData <<- vector("list", 0)
      } else {
        paths_to_clear <- rawFiles[indices]
        rawData[paths_to_clear] <<- NULL
      }
      gc()
    },
    getPeak = function(rawLinki, comp_mzi, comp_rti, minIndi, maxIndi, width,
                       mztoli, rttoli) {
      "Internal function to integrate peaks during processing"
      
      getPeakAtWidth <- function(widthi) {
        resul <- shewhart4hrms::peakpicking_BfG_cpp(
          i = comp_mzi - widthi / 2,
          rawData = rawLinki, mz_step = widthi,
          rt_min_scan = minIndi,
          rt_max_scan = maxIndi,
          sn = settings$sn, int_threshold = settings$baseline_noise_MS1,
          NoiseScans = 60, peakwidth_min = 4, 
          peakwidth_max = 100,  # large values chosen
          maxPeaksPerSignal = settings$peaksPerPeak,  
          precursormzTol = 5  # used to get MS2 scan, unnecessary
        )
        if (is.null(resul))
          resul <- matrix(nrow = 0, ncol = 16)
        # remove unknown extra column after switching to cpp
        resul <- resul[,-4, drop = FALSE]  
        colnames(resul) <- c("exactmass", "scantime", "peak_intens", 
                             "maxima", "scantimeleft_end", "scantimeright_end",
                             "left_end", "right_end", "noisedeviation", 
                             "peakArea", "FWHM_left", "FWHM_right", "noiselevel",
                             "i", "ms2scan")
        resul <- as.data.frame(resul)
        resul$scantime_min <- resul[, "scantime"] / 60
        # filter according to known rt from MS2 spectrum
        resul <- resul[(abs(resul$scantime_min - comp_rti) < rttoli) &
                         (abs(resul$exactmass - comp_mzi) < mztoli), ]
        if (nrow(resul) != 0)
          resul$e_width <- widthi
        resul
      }
      resu <- getPeakAtWidth(width)
      # if no peak was found, try the whole range and take the result
      # closest to the target width
      
      if (nrow(resu) == 0 && diff(settings$EIC_extraction_range) != 0) {
        ra <- settings$EIC_extraction_range
        widths <- seq(ra[1], ra[2], length.out = 10)
        # remove width already tested
        widths <- widths[-which(sapply(widths, all.equal, target = width) == "TRUE")]
        # reorder range by distance to target
        widths <- widths[order(abs(widths - width))]
        for (w in widths) {
          resu <- getPeakAtWidth(w)
          if (nrow(resu) != 0) {
            resu$e_width <- w
            break
          }
        }
      }
      resu
    },
    # Call suspect search for all (remaining) based on settings
    process_all = function() {
      "Process all currently unprocessed files in the report object. Previously recorded false
      positives (without specified sample) are deleted by default."
      # find out which are left to process
      to_process <- setdiff(rawFiles, rawFilesCompl$path)
      
      if (nrow(IS) == 0)
        stop("No internal standards loaded")
      
      # loop through each file and get IS data
      for (datFile in to_process) {
        #browser()
        # load this file
        worked <- .self$loadData(indices = which(datFile == rawFiles))
        if (!worked)
          next
        message(paste("Processing", basename(datFile)))
        stopifnot(inherits(rawData[[datFile]], "xcmsRaw"))
        rawLink <- rawData[[datFile]]
        
        # Integrate IS ####
        for (k in seq_len(nrow(IS))) {
          thisCharge <- switch(settings$pol, pos = 1, neg = -1)
          isRtSec <- (IS[k, "rt"] + settings$IS_rtoffset) * 60
          minISInd <- which.min(abs(rawLink@scantime - (isRtSec - (settings$ISrttolm * 60)/2)))
          maxISInd <- which.min(abs(rawLink@scantime - (isRtSec + (settings$ISrttolm * 60)/2)))
          
          IS_mz <- shewhart4hrms::get_mass(formula = IS[k, "formula"], charge = thisCharge,
                                         adduct = IS[k, "adduct"])
          IS_rt <- IS[k, "rt"] + settings$IS_rtoffset
          IS_res <- .self$getPeak(rawLink, IS_mz, IS_rt, minISInd, maxISInd, 
                                  settings$EIC_extraction,
                                  settings$ISmztol,
                                  settings$ISrttolm)
          
          if (nrow(IS_res) == 0) {
            message(paste(IS[k, "name"], "not found."), appendLF = FALSE)
            next
          }
          # get highest peak
          IS_res <- IS_res[which.max(IS_res$peak_intens), ]
          
          IS_res2 <- data.frame(samp = basename(datFile), IS = IS[k, "name"],
                                mz = IS_res$exactmass, rt = IS_res$scantime_min,
                                int_h = as.integer(round(IS_res$peak_intens)),
                                int_a = as.integer(round(IS_res$peakArea)),
                                peak_start = round(IS_res$scantimeleft_end / 60, 2),
                                peak_end = round(IS_res$scantimeright_end / 60, 2),
                                eic_extraction_width = IS_res$e_width,
                                ISpeakID = as.integer(currentISpeakID),
                                stringsAsFactors = FALSE)
          # export IS results 
          ISresults <<- rbind(ISresults, IS_res2)
          currentISpeakID <<- currentISpeakID + 1L
        }
        
        # update rawfiles completed list
        processed <- data.frame(path = datFile, date = date(), stringsAsFactors = FALSE)
        # update raw file complete
        rawFilesCompl <<- rbind(rawFilesCompl, processed)
        rm(rawLink)
        .self$clearData(indices = which(datFile == rawFiles))
        message("Complete")
      }
     
      
      message(sprintf("Processing completed %i files", length(to_process)))
    },


    deleteBelowIntThresh = function() {
      switch(settings$use_int_threshold,
             "area" = {
               # if no area available, use height threshold
               #browser()
               toThrow <- ifelse(is.na(peakList$int_a), 
                                 peakList$int_h < settings$height_threshold,
                                 peakList$int_a < settings$area_threshold)
               delIDs <- peakList[toThrow, "peakID"]
               .self$delPeakID(delIDs)
               message(paste(length(delIDs), 
                             "peaks were below area threshold", 
                             settings$area_threshold, "and deleted."))
             },
             "height" = {
               toThrow <- ifelse(is.na(peakList$int_h), 
                                 TRUE,
                                 peakList$int_h < settings$height_threshold)
               delIDs <- peakList[toThrow, "peakID"]
               .self$delPeakID(delIDs)
               message(paste(length(delIDs), 
                             "peaks were below height threshold", settings$height_threshold, "and deleted."))
             },
             "none" = {NULL},
             stop(
               "settings$use_int_threshold must be one of area, height or none"
             )
      )
      
    },
   
    
    clearAndSave = function(dialog = TRUE, nameReport = NULL) {
      "Clear data from RAM and save report as .report file in the current working directory. 
      Use nameReport to give a different location and different name as in 
      *.clearAndSave('D:\\exampleFolder\\example'). Note: The folder must exist beforehand.
      To read the file again, use the function shewhart4hrms::loadReport"
      
      if (dialog) {
        nameReport <- rstudioapi::selectFile("Save File as...", label = "Save", existing = F,
                                             filter = "DBscreening report file (*.report)")
        if (!grepl("\\.report$", nameReport))
          nameReport <- paste0(nameReport, ".report")
      } else if (is.null(nameReport)) {
        nameReport <- stringr::str_match(deparse(sys.call()), "^(.*)\\$clearAndSave")[,2]
        nameReport <- paste0(nameReport, ".report")
      }
      .self$clearData()
      saveRDS(.self, file = nameReport)
    },
    
    view = function() {
      "View results using shiny."
      
      require(shiny)
      require(ggplot2)
      
      nameReport <- stringr::str_match(deparse(sys.call()), "^(.*)\\$view")[,2]
      
      app <- shinyApp(
        # --UI-- ####
        ui = navbarPage(
          paste("dbscreening -", nameReport),
          # Overview ####
          tabPanel(
            "Overview",
            # Show a table of info
            verticalLayout(
              tableOutput("susSFileInfo"),
              hr(),
              DT::dataTableOutput("dataFileInfo")
            )
          ),
          # IS ####
          tabPanel(
            "IS",
            verticalLayout(
              DT::dataTableOutput("ISlist"),
              DT::dataTableOutput("ISresults")
            )
          ),
          # View IS ####
          tabPanel(
            "View IS",
            verticalLayout(
              selectInput("ISnameChosen", label = "IS name",
                          choices = unique(.self$IS$name)),
              splitLayout(plotOutput("ISarea",
                                     dblclick = "ISarea_dblclick",
                                     brush = brushOpts(
                                       id = "ISarea_brush",
                                       resetOnNew = TRUE)),
                          plotOutput("ISheight",
                                     dblclick = "ISheight_dblclick",
                                     brush = brushOpts(
                                       id = "ISheight_brush",
                                       resetOnNew = TRUE))
              ),
              splitLayout(plotOutput("ISwidth",
                                     dblclick = "ISwidth_dblclick",
                                     brush = brushOpts(
                                       id = "ISwidth_brush",
                                       resetOnNew = TRUE)),
                          plotOutput("ISrt",
                                     dblclick = "ISrt_dblclick",
                                     brush = brushOpts(
                                       id = "ISrt_brush",
                                       resetOnNew = TRUE))
              )
              
            )
          ),
         
          # Settings ####
          tabPanel(
            "Settings",
            tags$pre(textOutput("settingsOut"))
          )
          
        ),
        server = function(input, output, session) {
          
          home <- c('Home'="~")
   
          
          output$susSFileInfo <- renderTable({
            if (nrow(.self$rawFilesCompl) >= 1) {
              lastDate <- max(.self$rawFilesCompl$date)
            } else {
              lastDate <- NA
            }
            # make dataframe with info
            data.frame(
              Last_modification = lastDate,
              No._datafiles = length(.self$rawFiles),
              No._processed_datafiles = nrow(.self$rawFilesCompl),
              stringsAsFactors = FALSE
            )
          }, align = "c", spacing = "s")
          
          output$dataFileInfo <- DT::renderDataTable({
            isComplete <- .self$rawFiles %in% .self$rawFilesCompl$path
            df <- data.frame(
              Files = .self$rawFiles,
              Complete = isComplete,
              stringsAsFactors = FALSE
            )
            df$Process_date <- ifelse(
              df$Complete,
              .self$rawFilesCompl[.self$rawFilesCompl$path == df$Files, "date"],
              NA
            )
            df$Files <- stringr::str_wrap(df$Files, 20)
            df
          },
          options = list(pageLength = 25))
          
          # IS ####
          output$ISlist <- DT::renderDataTable({
            data.frame(IS_name = .self$IS$name,
                       Formula = .self$IS$formula,
                       RT = .self$IS$rt,
                       stringsAsFactors = FALSE
            )
          }, selection = "single")
          
          output$ISresults <- DT::renderDataTable({
            data.frame(
              IS = .self$ISresults$IS,
              mz = round(.self$ISresults$mz, 4),
              rt = round(.self$ISresults$rt, 2),
              height = round(.self$ISresults$int_h),
              area = round(.self$ISresults$int_a),
              width_min = round(.self$ISresults$peak_end - .self$ISresults$peak_start, 2),
              sample = stringr::str_wrap(.self$ISresults$samp, 10),
              stringsAsFactors = FALSE
            )
          }, options = list(pageLength = 25))
          
          # View IS ####
          # get some general plots about the behavior of IS, barplots over all samples, similar to check_IS
          ISnameSelected <- reactive({
            selIS <- input$ISlist_rows_selected
            if (is.null(selIS)) {
              selIS <- .self$IS[1, "name"]
            } else {
              selIS <- .self$IS[selIS, "name"]
            }
            selIS
          })
          
          
          observeEvent({
            input$ISlist_rows_selected
            input$dataFile
          },
          {
            
            updateSelectInput(session, "ISnameChosen", selected = ISnameSelected())
          })
          
          trendBreaks <- function(limits) {
            
            if (length(limits) <= 40) 
              limits else limits[seq(1, length(limits), length.out = 40)]
          }
          trendLabels <- function(breaks) {
            stringr::str_match(breaks, "^(.*)\\.mzX?ML$")[, 2]
          }
          rangesIS <- reactiveValues(x = NULL, y = NULL)
          output$ISarea <- renderPlot({
            
            toPlot <- .self$ISresults[.self$ISresults$IS == input$ISnameChosen, ]
            newPlot <- ggplot(toPlot, aes(x = samp, y = int_a)) +
              geom_bar(stat = "identity") +
              scale_x_discrete(limits = basename(.self$rawFiles), breaks = trendBreaks,
                               labels = trendLabels) +
              ylab("Intensity (area, cps)") + 
              theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
                    axis.title.x = element_blank())
            if (!is.null(rangesIS$x)) {
              from <- round(rangesIS$x[1])
              to <- round(rangesIS$x[2])
              newPlot <- newPlot +
                scale_x_discrete(limits = basename(.self$rawFiles)[from:to],
                                 breaks = trendBreaks, labels = trendLabels)
            }
            newPlot
          })
          
          observeEvent(input$ISarea_dblclick, {
            brush <- input$ISarea_brush
            if (!is.null(brush)) {
              rangesIS$x <- c(brush$xmin, brush$xmax)
              rangesIS$y <- c(brush$ymin, brush$ymax)
              
            } else {
              rangesIS$x <- NULL
              rangesIS$y <- NULL
            }
          })
          rangesISheight <- reactiveValues(x = NULL, y = NULL)
          output$ISheight <- renderPlot({
            
            toPlot <- .self$ISresults[.self$ISresults$IS == input$ISnameChosen, ]
            newPlot <- ggplot(toPlot, aes(x = samp, y = int_h)) +
              geom_bar(stat = "identity") +
              scale_x_discrete(limits = basename(.self$rawFiles), breaks = trendBreaks,
                               labels = trendLabels) +
              ylab("Intensity (height, cps)")  +
              theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
                    axis.title.x = element_blank())
            if (!is.null(rangesIS$x)) {
              from <- round(rangesIS$x[1])
              to <- round(rangesIS$x[2])
              newPlot <- newPlot +
                scale_x_discrete(limits = basename(.self$rawFiles)[from:to],
                                 breaks = trendBreaks, labels = trendLabels)
            }
            newPlot
          })
          observeEvent(input$ISheight_dblclick, {
            brush <- input$ISheight_brush
            if (!is.null(brush)) {
              rangesIS$x <- c(brush$xmin, brush$xmax)
              rangesIS$y <- c(brush$ymin, brush$ymax)
              
            } else {
              rangesIS$x <- NULL
              rangesIS$y <- NULL
            }
          })
          
          output$ISwidth <- renderPlot({
            
            toPlot <- .self$ISresults[.self$ISresults$IS == input$ISnameChosen, ]
            toPlot$width <- (toPlot$peak_end - toPlot$peak_start)
            
            newPlot <- ggplot(toPlot, aes(x = samp, y = width, group = 1)) +
              geom_line() +
              geom_point() +
              scale_x_discrete(limits = basename(.self$rawFiles), breaks = trendBreaks,
                               labels = trendLabels) +
              ylim(0, NA) +
              ylab("Peak width (min.)")  + 
              theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
                    axis.title.x = element_blank())
            if (!is.null(rangesIS$x)) {
              from <- round(rangesIS$x[1])
              to <- round(rangesIS$x[2])
              newPlot <- newPlot +
                scale_x_discrete(limits = basename(.self$rawFiles)[from:to],
                                 breaks = trendBreaks, labels = trendLabels)
            }
            newPlot
          })
          
          output$ISrt <- renderPlot({
            
            toPlot <- .self$ISresults[.self$ISresults$IS == input$ISnameChosen, ]
            newPlot <- ggplot(toPlot, aes(x = samp, y = rt, group = 1)) +
              geom_line() +
              geom_point() +
              scale_x_discrete(limits = basename(.self$rawFiles), breaks = trendBreaks,
                               labels = trendLabels) +
              ylim(0, NA) +
              ylab("RT (min., cps)") + 
              theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
                    axis.title.x = element_blank())
            if (!is.null(rangesIS$x)) {
              from <- round(rangesIS$x[1])
              to <- round(rangesIS$x[2])
              newPlot <- newPlot +
                scale_x_discrete(limits = basename(.self$rawFiles)[from:to],
                                 breaks = trendBreaks, labels = trendLabels)
            }
            newPlot
          })
          observeEvent(input$ISrt_dblclick, {
            brush <- input$ISrt_brush
            if (!is.null(brush)) {
              rangesIS$x <- c(brush$xmin, brush$xmax)
              rangesIS$y <- c(brush$ymin, brush$ymax)
              
            } else {
              rangesIS$x <- NULL
              rangesIS$y <- NULL
            }
          })

          # Settings ####
          output$settingsOut <- renderPrint({
            .self$settings
          })
          
          
        },
        
        onStart = function() {
          message(paste("Currently viewing", nameReport))
          onStop(function() {
            message("Visualization closed")
            detach("package:shinyFiles", unload = TRUE)
          })
        }
      )
      runApp(app)
    }
  )
)

