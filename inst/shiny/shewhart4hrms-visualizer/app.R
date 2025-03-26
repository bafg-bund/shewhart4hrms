



if (!exists("globalwd"))
  stop("Path to shewhart4hrms directory not found, use shewhart4hrms::newDirTree('path-to-shewart-directory')
       to create a new folder and viewShewhart('path-to-shewart-directory') to open the app")

startTimeValue <- Sys.time() - 1000
ui <- fluidPage(
  fluidRow(
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
      actionButton(inputId = "newFiles", label = "Refresh")
    )
  ),
  fluidRow(column(
    12,
    tabsetPanel(
      tabPanel(
        "Intensity_cps",
        column(8, plotOutput("inten")),
        column(4, plotOutput("intenHist"))
      ),
      tabPanel(
        "Intensity_area",
        column(8, plotOutput("intenA")),
        column(4, plotOutput("intenAHist"))
      ),
      tabPanel(
        "Mass shift mDa",
        column(8, plotOutput("mDa")),
        column(4, plotOutput("mDaHist"))
      ),
      tabPanel(
        "RT shift min", 
        column(8, plotOutput("RT")),
        column(4, plotOutput("RTHist"))
      ),
      tabPanel(
        "Width min", 
        column(8, plotOutput("width")),
        column(4, plotOutput("widthHist"))
      )
    )
  ))
)

## server ####
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
    filesHavePol <- polInFileName(filePaths())
    req(any(filesHavePol))
    
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
    
    prog$set(detail = "Complete", value = 3)

    prog$close()
   
  }) %>% bindEvent(input$newFiles)
  
  results <- reactive({
    input$newFiles
    resExists <- file.exists(filePaths()$results)
    validate(need(resExists, "No data found"))
    req(resExists)
    df <- getResults(filePaths())
    df
  })
  
  observe({
    updateSliderInput(session, "dates", min = min(results()$time), max = max(results()$time),
                      value = c(min(results()$time), max(results()$time)))
  }) %>% 
    bindEvent(input$pol, input$newFiles)
  
  output$inten <- renderPlot(makeTrendPlot(results(), input$dates, "int_h"))
  output$intenHist <- renderPlot(makeBellCurve(
    results(),
    input$dates,
    "int_h",
    filePaths()
  ))
  
  output$intenA <- renderPlot(makeTrendPlot(results(), input$dates, "int_a"))
  output$intenAHist <- renderPlot(makeBellCurve(
    results(),
    input$dates,
    "int_a",
    filePaths()
  ))
  
  output$mDa <- renderPlot(makeTrendPlot(results(), input$dates, "delta_mz_mDa"))
  output$mDaHist <- renderPlot(makeBellCurve(
    results(),
    input$dates,
    "delta_mz_mDa",
    filePaths()
  ))
  
  output$RT <- renderPlot(makeTrendPlot(results(), input$dates, "delta_rt_min"))
  output$RTHist <- renderPlot(makeBellCurve(
    results(),
    input$dates,
    "delta_rt_min",
    filePaths()
  ))
  
  output$width <- renderPlot(makeTrendPlot(results(), fixDates(input$dates), "peak_width_min"))
  output$widthHist <- renderPlot(makeBellCurve(
    results(),
    fixDates(input$dates),
    "peak_width_min",
    filePaths()
  ))
}

shinyApp(ui = ui, server = server)

# Copyright 2020-2025 Bundesanstalt für Gewässerkunde
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
