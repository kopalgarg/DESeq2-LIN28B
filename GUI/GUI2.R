library(tidyverse)
library(plotly)
library(BuenColors)
library(data.table)
require(dplyr)
require(ggplot2)
library(shiny)
library(shinyWidgets)
library(shinythemes)
library(ggrepel)
library(plotly)

cell_type = c('CD34', 'BFUE', 'CFUE', 'Pro', 
              'eBaso', 'lBaso', 'Poly', 'Ortho')
data = fread('https://raw.githubusercontent.com/kopal-garg/resources/master/DE_analysis/df.tsv')
data = data %>% distinct() %>% group_by(symbol, ensg, cell_type, baseMean,
                                        log2FoldChange) %>% mutate(padj = median(padj),
                                                                   pvalue = median(pvalue)) %>% distinct()

ui <- fluidPage(
  titlePanel(""),
  fluidRow(
    sidebarPanel(selectInput(inputId = 'cell_type',
                             label = 'Select Cell Type:',
                             choices = cell_type, selected = 'Pro')),
    column(
      width = 12,
      plotlyOutput("volcanoPlot", height = "700px"),
      downloadButton(outputId = 'download_file', 'Data'),
      downloadButton(outputId = 'download_plot', 'Plot')
    )
  )
)

server <- function(input, output) {
  
  output$volcanoPlot <- renderPlotly({
    cell_type=input$cell_type
    differentialExpressionResults <-
      data %>%
      filter(cell_type == input$cell_type) %>%
      mutate(
        minusLog10Pvalue = -log10(padj),
        tooltip = (paste(symbol,"\n", ensg, "\n",
                         "log2FoldChange: ",log2FoldChange, "\n",
                         "adjusted p-val: ", padj, "\n",
                         "p-val: ", pvalue, "\n",
                         "Base Mean", baseMean ))
      ) %>%
      mutate(group = case_when(-log10(padj) <= 3 ~ '-log10(padj) < 3',
                               abs(log2FoldChange)<= 2 ~'abs(log2FoldChange) < 2',
                               -log10(padj) >= 3 & log2FoldChange <= 2 ~ 'upregulated',
                               -log10(padj) >= 3 & log2FoldChange >= 2 ~ 'downregulated'
                               
      ))
    
    plot <- differentialExpressionResults %>%
      ggplot(aes(x = log2FoldChange,
                 y = minusLog10Pvalue,
                 text = tooltip,
                 color = group,
                 key = row.names(differentialExpressionResults))) +
      geom_point(shape = 20, alpha = 0.5) +
      scale_color_manual(values = c('blue','blue', 'red', 'green')) +
      xlab("log fold change") +
      labs(title = paste("Volcano Plot for Differential Expression (CB - PB)", "\n",
                         cell_type),
           subtitle = 'CB - PB') +
      ylab("-log10(P-value)") + 
      pretty_plot() +
      geom_vline(xintercept = -2, colour="black", linetype="dashed") +
      geom_vline(xintercept = 0, colour="gray", linetype="dashed") +
      geom_vline(xintercept = 2, colour="black", linetype="dashed") +
      geom_hline(yintercept = 3, colour="black", linetype="dashed") +
      scale_y_continuous(trans = "log1p") +
      theme(
        plot.title = element_text(hjust = 0.5, size = 10),
        legend.title = element_blank(),
        plot.subtitle = element_text(hjust = 0.3 ))
    plot %>%
      ggplotly(tooltip = "tooltip") %>%
      layout(dragmode = "select") 
  })
  
  # download dataset
  output$download_file <- downloadHandler(
    filename = function() {
      paste(input$cell_type, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(data %>% filter(cell_type==input$cell_type), file, row.names = FALSE)
    }
  )
  # download plot
  output$download_plot <- downloadHandler(
    filename = function() { paste(input$cell_type, '.png', sep='') },
    content = function(file) {
      device <- function(..., width, height) grDevices::png(..., width = width, height = height, res = 300, units = "in")
      
      # save the non-plotly plot
      cell_type=input$cell_type
      differentialExpressionResults <-
        data %>%
        filter(cell_type == input$cell_type) %>%
        mutate(
          minusLog10Pvalue = -log10(padj),
          tooltip = (paste(symbol,"\n", ensg, "\n",
                           "log2FoldChange: ",log2FoldChange, "\n",
                           "adjusted p-val: ", padj, "\n",
                           "p-val: ", pvalue, "\n",
                           "Base Mean", baseMean ))
        ) %>%
        mutate(group = case_when(-log10(padj) <= 3 ~ '-log10(padj) < 3',
                                 abs(log2FoldChange)<= 2 ~'abs(log2FoldChange) < 2',
                                 -log10(padj) >= 3 & log2FoldChange <= 2 ~ 'upregulated',
                                 -log10(padj) >= 3 & log2FoldChange >= 2 ~ 'downregulated'
                                 
        ))
      
      plot <- differentialExpressionResults %>%
        ggplot(aes(x = log2FoldChange,
                   y = minusLog10Pvalue,
                   text = tooltip,
                   color = group,
                   key = row.names(differentialExpressionResults))) +
        geom_point(shape = 20, alpha = 0.5) +
        scale_color_manual(values = c('blue','blue', 'red', 'green')) +
        xlab("log fold change") +
        labs(title = paste("Volcano Plot for Differential Expression (CB - PB)", "\n",
                           cell_type)) +
        ylab("-log10(P-value)") + 
        pretty_plot() +
        geom_vline(xintercept = -2, colour="black", linetype="dashed") +
        geom_vline(xintercept = 0, colour="gray", linetype="dashed") +
        geom_vline(xintercept = 2, colour="black", linetype="dashed") +
        geom_hline(yintercept = 3, colour="black", linetype="dashed") +
        scale_y_continuous(trans = "log1p") +
        theme(
          plot.title = element_text(hjust = 0.5, size = 10),
          legend.title = element_blank(),
          plot.subtitle = element_text(hjust = 0.3 ))
      
      ggsave(file, plot = plot, device = device)
    }
  )
  
  
}

shinyApp(ui, server, options = list(height = 600))