library(shiny)
library(tidyverse)
library(highcharter)

cell_type = c('CD34', 'BFUE', 'CFUE', 'Pro', 
              'eBaso', 'lBaso', 'Poly', 'Ortho')
data = fread('https://raw.githubusercontent.com/kopal-garg/resources/master/DE_analysis/df.tsv')


ui <- fluidPage(
  titlePanel(""),
  fluidRow(
    column(
      width = 7,
      highchartOutput("volcanoPlot", height = "500px")
    )
  )
)

server <- function(input, output) {
  
  output$volcanoPlot <- renderHighchart({
    differentialExpressionResults <-
      data %>%
      filter(cell_type == input$cell_type) %>%
      mutate(
        minusLog10Pvalue = -log10(padj),
        tooltip = (symbol)
      ) %>%
      mutate(group = case_when(padj < 0.05 & abs(log2FoldChange)< 1.5 ~ 'Significant',
                               padj > 0.05 & abs(log2FoldChange)> 1.5 ~ 'FoldChange',
                               padj < 0.05 & abs(log2FoldChange)> 1.5 ~ 'Significant & Fold Change',
                               TRUE ~ "Not Significant"
      ))
    
    plot <- highchart() %>%
      hc_chart(animation = TRUE, zoomType = "xy") %>%
      hc_xAxis(title = list(text = "log fold change")) %>%
      hc_yAxis(title = list(text = "-log10(P-value)")) %>%
      hc_colors(c('rgba(67, 67, 72, 0.6)', 'rgba(124, 181, 236, 0.6)')) %>%
      hc_legend(layout = "horizontal", align = "center", verticalAlign = "bottom") %>%
      hc_exporting(enabled = TRUE, filename = "plot.png") %>%
      hc_tooltip(
        animation = FALSE,
        formatter = JS("function() { return (this.point.tooltip) }")) %>%
      hc_plotOptions(series = list(
        cursor = "pointer",
        point = list(events = list(click = JS("function() { Shiny.onInputChange('volcanoPlot_clicked', this.id) }")))))
    
     seriesData <- differentialExpressionResults %>%
        select(id = group, x = log2FoldChange, y = minusLog10Pvalue, tooltip) %>%
        list_parse
      
      plot <- plot %>%
        hc_add_series(
          data = seriesData,
          type = "scatter",
          marker = list(symbol = "circle"),
          stickyTracking = T
        )

    plot
  })
  
}

shinyApp(ui, server, options = list(height = 575))