#C6_Shiny app
library(shiny)
library(tidyverse)
library(car)
library(stringr)
library(plotly)
library(vroom)
library(patchwork)
library(gridExtra)


colors.2 <- c("#B266FF",
              "#CC99FA",
              "darkgreen", 
              "#00FF7F", 
              "orange", 
              "brown", 
              "grey", 
              "black", 
              "#660033", 
              "#99004C", 
              "#8B4513", 
              "#D2691E", 
              'red', 
              "blue")


parseFile <- function(pathname, t = F){
  #return as a head(df) object first
  test.file <- read_csv(pathname)
  if ( ncol(test.file) == 2) {
    top <- test.file[0:9,]  #non-numeric data; top 9 rows of data file
    bottom <- test.file[10:nrow(test.file),] #rest of data file, where numeric data is present
    names(bottom) <- c("time", "data") #time in col1, data ready to be split using commands below
    time <- bottom$time[3:nrow(bottom)] 
    cols <- unlist(strsplit(as.character(bottom[1,2]), split = ",")) #phyco, chlorophyll, turb, ...
    data <- bottom[3:nrow(bottom),2] %>% #splits all data values from '1,2,3,4,5' format --> separate columns with each value in a different column/feature
      separate(data, into =cols[1:8], sep = ",", convert = TRUE)
  } else if (ncol(test.file) == 13) {
    top <- test.file[0:9,]
    bottom <- test.file[10:nrow(test.file),]
    time <- bottom[3:nrow(bottom),1] 
    cols <- as.character(bottom[1,2:9])
    data <- as.data.frame(sapply(bottom[3:nrow(bottom), 2:9], as.double)) #data is already in split up format, command coerces all values to double dtype
    names(data) <- cols
    data <- tibble(data)
  }
  
  if (t == T) 
    data[,"ts"] = time[1,1]
  
  return(data)
  
}

get_mid <- function(data) {
  depth <- pull(data, Depth)
  idxs <- which(depth == max(depth))
  return (as.integer(mean(idxs)))
}


#can you select two points at once and save them, for each cast or find cast sizes automatically
ui <- fluidPage(
  headerPanel("C6"), 
  tabsetPanel(
    tabPanel("Offsets",
        verbatimTextOutput("note_offset"),
        hr(),
        sidebarLayout(
          sidebarPanel(
            fileInput("upload", label = "Choose C6.csv file", multiple = FALSE,
                      accept = c("text/csv", 
                                 "text/comma-separated-values,text/plain", 
                                 ".csv")),
          strong("Profile 1: "), 
          numericInput("p1_min", "Start", value = 0, width = '60%'),  
          numericInput("p1_max", "End ", value = 0, width = '60%'),
          strong("Profile 2: "),
          numericInput("p2_min", "Start", value = 0, width = '60%'), 
          numericInput("p2_max", "End ", value = 0, width = '60%'),
          hr(), 
          actionButton("plot", "Plot")
          ), 
          mainPanel(
            plotlyOutput("depthGraph")
          )
        ),
        verbatimTextOutput("note_downloads"),
        dateInput("date", "Sample Date:", NA),
        downloadButton("up_download", "Download UPCAST.csv"), 
        downloadButton("down_download", "Download DOWNCAST.csv"), 
        downloadButton("plot_download", "Download Plots"), 
        verbatimTextOutput("check")
      ), 
    tabPanel("Plots", 
             plotOutput("up"), 
             hr(), 
             plotOutput("down")
             )
  )
)


server <- function(input, output, session) {
  
  #react elements
  #note: react elements req() carries over to output functions where react object is called. initialize react objects before using. (e.g: data = data())
  data <- reactive({
    req(input$upload)
    parseFile(input$upload$datapath)
  })
  
  file_name <- reactive({
    req(input$upload)
    strsplit(input$upload$name, ".csv")[[1]]
  })

  #max
  P1 <- reactive({
    req(input$upload)
    req(input$p1_min != 0)
    req(input$p1_max != 0)
    DATA <- data()[input$p1_min:input$p1_max,]
    mid_idx <- get_mid(DATA)
    DATA$profile <- NA
    DATA$cast <- NA
    DATA$cast[1:mid_idx] = "DOWN"
    DATA$cast[mid_idx:nrow(DATA)] = "UP"
    DATA$profile <- "P1"
    DATA
    
  })
  
  P2 <- reactive({
    req(input$upload)
    req(input$p2_min != 0)
    req(input$p2_max != 0)
    DATA <- data()[input$p2_min:input$p2_max,]
    mid_idx <- get_mid(DATA)
    DATA$profile <- NA
    DATA$cast <- NA
    DATA$cast[1:mid_idx] = "DOWN"
    DATA$cast[mid_idx:nrow(DATA)] = "UP"
    DATA$profile <- "P2"
    DATA
  })
  
  UPCAST <- reactive({
    P1 <- P1()
    P2 <- P2()
    UPCAST <- rbind(P1[which(P1$cast == "UP"),], P2[which(P1$cast == "UP"),]) 
    long <- UPCAST %>%
      pivot_longer(cols = names(UPCAST)[c(1:6,8)], names_to = "Variable", values_to = "Values") %>%
      mutate(color = interaction(Variable, profile, sep = "."))
    long
  })
  
  DOWNCAST <- reactive({
    P1 <- P1()
    P2 <- P2()
    DOWNCAST <- rbind(P1[which(P1$cast == "DOWN"),], P2[which(P1$cast == "DOWN"),])
    long <- DOWNCAST %>%
      pivot_longer(cols = names(DOWNCAST)[c(1:6,8)], names_to = "Variable", values_to = "Values") %>%
      mutate(color = interaction(Variable, profile, sep = "."))
    long
  })
  
  up_plot <- reactive({
    ggplot(UPCAST(), aes(y = Depth,x = Values, color = color, shape = profile)) + 
      geom_point(size = 2) + 
      facet_wrap(~Variable, scales = "free_x") + 
      theme_classic() + scale_x_continuous(position = "top") + scale_y_reverse() + 
      ylab("Depth (m)") + 
      ggtitle("UP-Cast Profiles") +
      guides(color = F) + 
      scale_color_manual(values = colors.2)
  })
  
  down_plot <- reactive({
    ggplot(DOWNCAST(), aes(y = Depth,x = Values, color = color, shape = profile)) + 
      geom_point(size = 2) + 
      facet_wrap(~Variable, scales = "free_x") + 
      theme_classic() + scale_x_continuous(position = "top") + scale_y_reverse() + 
      ylab("Depth (m)") + 
      ggtitle("DOWN-Cast Profiles") +
      guides(color = F) + 
      scale_color_manual(values = colors.2)
  })
  
  
  #output
  output$note_offset <- renderText({
    "Notes:\n- Load in File using 'Browse'\n- Profile(s) start and end at the beginning and end of each 'V'\n- Plots will be outputted under 'Plots' tab" 
    
  })
  output$note_downloads <- renderText({
    "Download data and Plots here; make sure profile indices above are filled in first\n- Set input date below to the date C6 was used"
  })
  
  output$depthGraph <- renderPlotly({
    data <- data()
    data$index <- seq(1:nrow(data))
    p <- plot_ly(data = data, x = ~index, y = ~Depth, 
                 type = 'scatter', mode = 'markers', 
                 text = ~paste("</br>Depth:", Depth, 
                               "</br>Idx:", index),
                 hoverinfo = 'text') %>% layout(
                   yaxis = list(
                     autorange = "reversed"
                   ), 
                   xaxis = list(
                     side = "top"
                   )
                 )
    p
  })
  
  observeEvent(input$plot, {
    output$up <- renderPlot({
      up_plot() 
    })
    
    output$down <- renderPlot({
      down_plot()
    })
  })
  
  output$up_download <- downloadHandler(
    filename = function() {
      paste("Processed_UPCAST_",file_name(), ".csv", sep = "" )
    }, 
    content = function(file) {
      write.csv( cbind(UPCAST()[,-ncol(UPCAST())], time = input$date), file, row.names = F)
    }
    
  )
  
  output$down_download <- downloadHandler(
    filename = function() {
      paste("Processed_DOWNCAST_",file_name(), ".csv", sep = "" )
    }, 
    content = function(file) {
      write.csv(cbind(DOWNCAST()[,-ncol(DOWNCAST( ))], time = input$date), file, row.names = F)
    }
    
  )
  
  output$plot_download <- downloadHandler(
    filename = function() {
      paste("CastPlots_",file_name(), ".png", sep = "" )
    }, 
    content = function(file) {
      
      png(file, width = 1000, height = 700, units = "px")
      print(grid.arrange(up_plot(), down_plot(), nrow = 2 ))
      dev.off()
    }
  )

}




shinyApp(ui, server)



#ggsave(file, plot = grid.arrange(up_plot(), down_plot(), nrow = 2 ), device = "png")
#make it so the total plots are custoizable in terms of which one is plotted via dropdown menu. 
#In addition should all data be downloaded at once? Or should this be UI customizable.