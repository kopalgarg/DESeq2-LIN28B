library(shiny)
library(shinyWidgets)
require(tidyverse)
require(annotables)
require(grImport)
require(BuenColors) 
require(cowplot)
require(edgeR) 
library(ggplot2)
library(magrittr)
library(shinythemes)
library(ggrepel)
library(plotly)


cell_type = c('CD34', 'BFUE', 'CFUE', 'Pro', 
              'eBaso', 'lBaso', 'Poly', 'Ortho')
data = fread('https://raw.githubusercontent.com/kopal-garg/resources/master/DE_analysis/df.tsv')

ui <- fluidPage(
  header <- shinydashboard::dashboardHeader(title = "DE Analysis: Cord- vs Peripheral-Blood" )  ,
  titlePanel('', windowTitle = 'DE analysis'),
  dropdown(
    
    selectInput(inputId = 'cell_type',
                label = 'Select Cell Type:',
                choices = cell_type),
    style = "unite", icon = icon("gear"),
    status = "danger", width = "300px",
    animate = animateOptions(
      enter = animations$fading_entrances$fadeInLeftBig,
      exit = animations$fading_exits$fadeOutRightBig
    )
  ),
  
  plotOutput(outputId = 'volcano_plot', width = 1000, height = 500)
)

server <- function(input, output, session) {
  
  plot_fig <- function(data, cell_type){
    data = data %>% filter(cell_type == !!cell_type)
    threshold_OE <- data$padj < 0.05
    data$threshold <- threshold_OE 
    data <- data[order(data$padj), ] 
    data$genelabels <- ""
    data$genelabels[1:10] <-T
  
    
  g = ggplot(data) +
      geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
      xlab("log2 fold change") + 
      ylab("-log10 adjusted p-value") + 
      theme(legend.position = "none",
            legend.title = element_blank()) +
      labs(colour = NULL) + 
      ggtitle(cell_type) +
      pretty_plot()
  
    ggplotly(g, tooltip = c("symbol"))
  }
  output$volcano_plot <- renderPlot(plot_fig(df, input$cell_type))
  
}

shinyApp(ui = ui, server = server)


