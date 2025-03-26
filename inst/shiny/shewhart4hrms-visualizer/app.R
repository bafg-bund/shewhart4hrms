



library(shiny)
library(ggplot2)

# Setup ####

filesPos <- newFilePaths(globalwd, "pos")

if (!exists("globalwd"))
  stop("Path to shewhart4hrms directory not found, use shewhart4hrms::create_dir('path-to-shewart-directory')
       to create a new folder and viewShewhart('path-to-shewart-directory') to open the app")

checkShewartDir(globalwd)
checkForPosResults(globalwd)

# Load data ####

dat1 <- getResults(filesPos)
# App ####
## ui ####
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
        value = c(min(dat1$time), max(dat1$time)),
        min = min(dat1$time),
        max = max(dat1$time),
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
        fluidRow(plotOutput("inten")),
        fluidRow(plotOutput("intenHist"))
      ),
      tabPanel(
        "Intensity_area",
        fluidRow(plotOutput("intenA")),
        fluidRow(plotOutput("intenAHist"))
      ),
      tabPanel("Mass shift mDa", fluidRow(plotOutput("mDa")), fluidRow(plotOutput("mDaHist"))),
      tabPanel("RT shift min", fluidRow(plotOutput("RT")), fluidRow(plotOutput("RTHist"))),
      tabPanel("Width min", fluidRow(plotOutput("width")), fluidRow(plotOutput("widthHist")))
    )
  ))
)

## server ####
server <- function(input, output, session) {
  
  filePaths <- reactiveVal()
  
  observe({
    filePaths(newFilePaths(globalwd, input$pol))
  })
  
  observe({
    
    prog <- Progress$new(session, min = 0, max = 2)
    prog$set(message = "Processing new IS files", detail = "looking for new files")

    checkFileNamesHavePolShiny(filePaths(), prog)
    filesHavePol <- polInFileName(filePaths())
    req(all(filesHavePol))
    
    newMeasFiles <- getNewMeasFiles(filePaths())
    if (length(newMeasFiles) > 0) {
      prog$set(detail = paste("Processing new files for ESI", isolate(input$pol)), value = 1)
      succPos <- processNewFiles(filePaths())
      if (succPos) {
        showNotification(
          paste("Update successful for", paste(newMeasFiles, collapse = ", "))
        )
      }
        
    } else {
      showNotification(paste("No new files found in ESI", isolate(input$pol)))
    }
    
    prog$set(detail = "Complete", value = 2)

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
  
  output$inten <- renderPlot(makeTrendPlot(results(), fixDates(input$dates), "int_h"))
  output$intenHist <- renderPlot(makeBellCurve(
    results(),
    fixDates(input$dates),
    "int_h",
    input$pol,
    filePaths()
  ))
  
  output$intenA <- renderPlot(makeTrendPlot(results(), fixDates(input$dates), "int_a"))
  output$intenAHist <- renderPlot(makeBellCurve(
    results(),
    fixDates(input$dates),
    "int_a",
    input$pol,
    filePaths()
  ))
  
  output$mDa <- renderPlot(makeTrendPlot(results(), fixDates(input$dates), "delta_mz_mDa"))
  output$mDaHist <- renderPlot(makeBellCurve(
    results(),
    fixDates(input$dates),
    "delta_mz_mDa",
    input$pol,
    filePaths()
  ))
  
  output$RT <- renderPlot(makeTrendPlot(results(), fixDates(input$dates), "delta_rt_min"))
  output$RTHist <- renderPlot(makeBellCurve(
    results(),
    fixDates(input$dates),
    "delta_rt_min",
    input$pol,
    filePaths()
  ))
  
  output$width <- renderPlot(makeTrendPlot(results(), fixDates(input$dates), "peak_width_min"))
  output$widthHist <- renderPlot(makeBellCurve(
    results(),
    fixDates(input$dates),
    "peak_width_min",
    input$pol,
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
