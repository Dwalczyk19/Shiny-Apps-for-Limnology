#C6_Shiny app
library(shiny)
library(tidyverse)
library(car)
library(stringr)
library(plotly)
library(vroom)



parseFile <- function(pathname, t = F){
  #return as a head(df) object first
  test.file <- read_csv(pathname) #assume all data files are given as csv. (That's what the C6 writes to)
  if ( ncol(test.file) == 2) {
    top <- test.file[0:9,]  #non-numeric data; top 9 rows of data file
    bottom <- test.file[10:nrow(test.file),] #rest of data file, where numeric data is present
    names(bottom) <- c("time", "data") #time in col1, data ready to be split using commands below
    time <- bottom$time[3:nrow(bottom)] 
    cols <- unlist(strsplit(as.character(bottom[1,2]), split = ",")) #phyco, chlorophyll, turb, ...
    data <- bottom[3:nrow(bottom),2] %>% #splits all data values from '1,2,3,4,5' format --> separate columns with each value in a different column/feature
      separate(data, into =cols[1:8], sep = ",", convert = TRUE)
    #names(data)[5:6] <- c("FI1", "FI2") 
    names(data) <- gsub("`", "", names(data))
    data <- data.frame(data)
  } else if (ncol(test.file) == 13) {
    top <- test.file[0:9,]
    bottom <- test.file[10:nrow(test.file),]
    time <- bottom[3:nrow(bottom),1] 
    cols <- as.character(bottom[1,2:9])
    data <- as.data.frame(sapply(bottom[3:nrow(bottom), 2:9], as.double)) #data is already in split up format, command coerces all values to double dtype
    names(data) <- cols
    #names(data)[5:6] <- c("FI1", "FI2")
    names(data) <- gsub("`", "", names(data))
    data <- data.frame(data)
  }
  
  if (t == T) 
    data[,"ts"] = time[1,1]
  
  return(data)
  
}

create_df <- function(data, p1start, p1end, p2start, p2end) { #add cast column for profiles, return a wide format as well just in case (return a vector with two elemetns)
  p1 <- data[p1start:p1end,] %>%
    mutate(profile = "P1")#check if up or down
  p2 <- data[p2start:p2end,] %>%
    mutate(profile = "P2")
  
  data.wide <- rbind(p1,p2)
  feats <- names(data.wide)[-which(names(data.wide) %in% c("Depth", "profile"))]
  data.long <- data.wide %>%
    pivot_longer(cols = all_of(feats), names_to = "Variable", values_to = "Values")
  data.long$color <- interaction(data.long$Variable, data.long$profile, sep = ".")
  return (list(data.wide, data.long))
}

make_plot <- function(data.long) {
  colors <- c("#B266FF","orange", "darkgreen", "#660033", "brown", "#8B4513","red", "#CC99FA", "#8B4513", "#00FF7F", "#99004C", "#D2691E", "black","blue")
  p <- ggplot(data.long, aes(y = Depth,x = Values, color = color, shape = profile)) +  
    geom_point(size = 2) +
    facet_wrap(~Variable, scales = "free") + #free_x or free
    theme_bw() + 
    scale_x_continuous(position = "top") + scale_y_reverse() +
    guides(color = F) + 
    labs(y = "Depth (m)", x = "Values (RFU)") + 
    theme(legend.title = element_text(size = 14), 
          strip.text = element_text(size = 14, face ="bold")) + 
    scale_color_manual(values = setNames(colors, unique(data.long$color)))
  return (p)
}


ui <- navbarPage(
  strong("C6 - Plotting Fluorescence targets"),
  id = "navbar", 
  tags$head(
    tags$style(
      HTML(".navbar-brand{ font-size: 24px;}" #.navbar specified by id, brand might be navbar specific because it's not in base CSS/HTML. different from like h1() or body()
           
           )
    )
  ),
  tabPanel("Download", 
           fluidRow(column(12, align = "left",  strong("Instructions:"),
                           br(),
                           HTML("(1) Please begin with selecting a C6 data file to be processed."),
                           br(),
                           HTML("(2) Each 'V' in the graph created below is a profile, grouped by two casts. Please select 'Start' and 'End' indices for one cast per profile."))
           ),
           br(),
           sidebarLayout(
             sidebarPanel(
               fileInput("upload", label = "C6 File", multiple = FALSE,
                         accept = c("text/csv", 
                                    "text/comma-separated-values,text/plain", 
                                    ".csv")),
               strong("Profile 1: "),
               numericInput("p1_start", "Start", value = 0, width = '60%'),  
               numericInput("p1_end", "End ", value = 0, width = '60%'),
               strong("Profile 2: "),
               numericInput("p2_start", "Start", value = 0, width = '60%'), 
               numericInput("p2_end", "End ", value = 0, width = '60%')
             ), 
             mainPanel(
               plotlyOutput("idx_check")
             )
           ), 
           br(),
           fluidRow(
             column(12, align = "center", 
                    tags$head(tags$style(".load{background-color:#5a90ed;} .load{color:white;}")), #specify download buttons
                    downloadButton("download_data", "Download Data", class = "load"),
                    downloadButton("download_plot", "Download Plot(s)", class = "load")
             )
           )
  ),
  tabPanel("Plotting",
           column(12, align = "center",
             plotOutput("plot_data", height = "100vh", width = "100vh")
           )
  )
)
server <- function(input, output, session) {
  data <- reactive({
    parseFile(input$upload$datapath)
  })
  
  data.long <- reactive({
    req(input$p1_start != 0 & input$p1_end != 0 & input$p2_start !=0 & input$p2_end != 0)
    create_df(data(), input$p1_start, input$p1_end, input$p2_start, input$p2_end)[[2]]
  })
  
  output$idx_check <- renderPlotly({
    if (!is.null(input$upload)) {
      idx <- 1:nrow(data())
      plot_ly(data = data(), x = ~idx, y = ~Depth, 
              type = "scatter", mode = "markers", text = ~paste("Index:", idx), 
              hoverinfo = "text") %>%
        layout(yaxis = list(autorange = "reversed", title = "Depth (m)"), xaxis = list(side = "top", title = "Time (s)"))
    }
  })
  
  output$plot_data <- renderPlot({
    if (!is.null(input$upload)) {
      if (input$p1_start != 0 & input$p1_end != 0 & input$p2_start !=0 & input$p2_end != 0) {
        plot <- make_plot(data.long())
        plot
      }

    }
  })
  
  output$download_plot <- downloadHandler(
    filename = "C6_plots.png", 
    content = function(file) {
      png(file, width = 700, height = 700, units = "px")
      print(make_plot(data.long()))
      dev.off()
    }
  )
  
  output$download_data <- downloadHandler(
    filename = "C6_dataLong.csv", 
    content = function(file) {
      load <- data.long() %>%
        mutate(color = NULL)
      write.csv(load, file, row.names = F)
    }
  )
  
}

shinyApp(ui, server)