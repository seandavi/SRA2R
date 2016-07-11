library(shiny)
library(ggplot2)
# Coverage plot representing pileup depth of given SRA accession and locationon genome
ui <- shinyUI(fluidPage(
  titlePanel("Coverage Plot"),
  
  sidebarLayout(
    sidebarPanel(
      textInput("acc",
                "Accession",
                ""),
      uiOutput("ref"),
      numericInput("start",
                   "Start",
                   1,
                   min = 1),
      numericInput("stop",
                "Stop",
                0,
                min = 0),
      numericInput("minDepth",
                   "Minimum Pileup Depth",
                   0, 
                   min = 0),
      actionButton("go", "Plot")
    ),
    
    mainPanel(
      plotOutput("depthPlot")
    )
  )
))

# Define server logic required to draw a histogram
server <- shinyServer(function(input, output, clientData, session) {
  v = reactiveValues(run = FALSE)
  observeEvent(input$go, {
    v$run = input$go
  })
  output$ref <- renderUI({
    if (validAccession(input$acc)) {
      refs = getReference(input$acc)
      selectInput("ref", "Chromosome Name", sort(refs$CommonNames))
    } else {
      selectInput("ref", "Chromosome Name", list())
    }
  })
  
  output$depthPlot <- renderPlot({
    if(v$run == FALSE) return()
    isolate(if(!is.null(input$acc) &&
               !is.null(input$ref)) {
      r = getPileUp(input$acc, input$ref, input$start, input$stop, input$minDepth)
      maxDepth = max(r$PileupDepth)
      factor = 7
      col = as.factor(-floor(r$PileupDepth / maxDepth * factor))
      r = cbind(r, col)
      ggplot(data = r, aes(x = ReferencePosition, y = PileupDepth)) + geom_bar(stat = "identity", aes(fill=r$col)) + guides(fill=FALSE)
    })
  })
})

shinyApp(ui = ui, server = server)

