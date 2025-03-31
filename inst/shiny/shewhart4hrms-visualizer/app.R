
if (!exists("globalwd"))
  stop("Path to shewhart4hrms directory not found, use shewhart4hrms::newDirTree('path-to-shewart-directory')
       to create a new folder and viewShewhart('path-to-shewart-directory') to open the app")

startTimeValue <- Sys.time() - 1000
ui <- fixedPage(
  fixedRow(
    column(3, h2(paste(
      "shewhart4hrms", basename(globalwd)
    ))),
    column(
      6,
      sliderInput(
        inputId = "dates",
        label = "",
        timeFormat = "%F",
        value = c(startTimeValue, startTimeValue),
        min = startTimeValue,
        max = startTimeValue,
        width = "100%"
      )
    ),
    column(
      1,
      radioButtons(
        "pol",
        "",
        choices = c("pos", "neg"),
        selected = "pos",
        inline = F
      ),
      offset = 1
    ),
    column(
      1, 
      style = "margin-top: 25px;", 
      actionButton(inputId = "refresh", label = "Refresh")
    )
  ),
  fixedRow(column(
    12,
    uiOutput("tabSet")
    
  ))
)

server <- function(input, output, session) {
  
  checkShewartDir(globalwd)
  
  filePaths <- reactiveVal()
  
  observe({
    filePaths(newFilePaths(globalwd, input$pol))
  })
  
  observe({
    
    prog <- Progress$new(session, min = 0, max = 3)
    prog$set(message = "Processing new IS files", detail = "looking for new files")

    checkFileNamesHavePolShiny(filePaths(), prog)
    
    if (!isInitiated(filePaths())) {
      prog$set(detail = paste("Initiating database for ESI", isolate(input$pol)), value = 1)
      initializeResultsTable(filePaths())
    }
    
    newMeasFiles <- getNewMeasFiles(filePaths())
    if (length(newMeasFiles) > 0) {
      prog$set(detail = paste("Processing new files for ESI", isolate(input$pol)), value = 2)
      succPos <- processNewFiles(filePaths())
      if (succPos) {
        showNotification(
          paste("Update successful for", paste(newMeasFiles, collapse = ", "))
        )
      }
        
    } else {
      showNotification(paste("No new files found in ESI", isolate(input$pol)))
    }
    
    prog$close()
   
  }) %>% bindEvent(input$refresh)
  
  results <- reactive({
    input$refresh  
    requireFileExists(filePaths()$results)
    results <- getResults(filePaths())
    req(nrow(results) > 0)
    results
  })
  
  warningLevels <- reactive({
    input$refresh
    requireFileExists(filePaths()$warningLevels)
    warningLevels <- read.csv(filePaths()$warningLevels)
    subset(warningLevels, polarity == isolate(input$pol))
  })
  
  observe({
    updateSliderInput(session, "dates", min = min(results()$time), max = max(results()$time),
                      value = c(min(results()$time), max(results()$time)))
  }) %>% bindEvent(input$pol, input$refresh)
  
  plotHeightSpec <- reactive({
    numIs <- length(unique(results()$name_is))
    paste0(numIs * 200, "px")
  })
  
  output$tabSet <- renderUI({
    makeTabSetPanel(plotHeightSpec = plotHeightSpec())
  })
  
  output$intensityTimeSeries <- renderPlot(makeTimeSeries(results(), warningLevels(), input$dates, "intensity"))
  output$intensityBellCurve <- renderPlot(makeHistogram(results(), input$dates, "intensity"))
  
  output$areaTimeSeries <- renderPlot(makeTimeSeries(results(), warningLevels(), input$dates, "area"))
  output$areaBellCurve <- renderPlot(makeHistogram(results(), input$dates, "area"))
  
  output$delta_mz_mDaTimeSeries <- renderPlot(makeTimeSeries(results(), warningLevels(), input$dates, "delta_mz_mDa"))
  output$delta_mz_mDaBellCurve <- renderPlot(makeHistogram(results(), input$dates, "delta_mz_mDa"))
  
  output$delta_rt_minTimeSeries <- renderPlot(makeTimeSeries(results(), warningLevels(), input$dates, "delta_rt_min"))
  output$delta_rt_minBellCurve <- renderPlot(makeHistogram(results(), input$dates, "delta_rt_min"))
  
  output$peak_width_minTimeSeries <- renderPlot(makeTimeSeries(results(), warningLevels(), input$dates, "peak_width_min"))
  output$peak_width_minBellCurve <- renderPlot(makeHistogram(results(), input$dates, "peak_width_min"))
}

shinyApp(ui = ui, server = server)

# Copyright 2025 Bundesanstalt für Gewässerkunde
# This file is part of shewhart4hrms
