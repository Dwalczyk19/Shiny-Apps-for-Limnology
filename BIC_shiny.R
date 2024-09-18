#BIC_Shiny
library(shiny)
library(tidyverse)
library(stats)
library(readxl)
library(plotly)

# '4' separate sections: 
#1 - Find offset row 
#2 - UV Check
#3 - Regression profiles (this section is added incase linear relationship is not specified by first .5/1m)
#4 - Calculated KD values

#initial parsing
parse_file <- function(pathname) {
  data <- read_csv(pathname, skip = 9, show_col_types = F) #top of df rows 1:9
  data$...16 = NULL
  data$Edz305nm <- as.numeric(data$Edz305nm)
  return(as.data.frame(data)) #[26:nrow(data),]))
}

#Lake name
get_name <- function(file) {
  space <- unlist(gregexpr(" ", file))[1]
  lake_file <- as.character(file )
  sample_date <- substr(lake_file, 1, space-1)
  lake_name <- substr(lake_file, space+1, nchar(lake_file)-4)
  return(lake_name)
}

#UV check (1)
plot_deck <- function(data_react, idx, max_idx, lake_name) {
  data <- data_react
  EDZ <- data[,c(2:5)]
  ED0 <- data[,c(10:13)]
  deck_cell <- (EDZ/ED0) * 100
  deck_cell$Depth <- data$Depth
  deck_cell_long <- deck_cell[idx:max_idx, ] %>%
    pivot_longer(cols = names(deck_cell)[1:4], values_to = "Values", names_to = "Variable") %>%
    mutate(Values = log(Values))
  #plotly
  p <- plot_ly(data = deck_cell_long, x = ~Values, y = ~Depth,
               color = ~Variable, symbol = ~Variable,
               type = "scatter", mode = "markers") %>%
      layout(
        xaxis = list(side = 'top', title = "Log(% of deck cell)"),   # Move x-axis to the top
        yaxis = list(autorange = 'reversed', title = "Depth"),  #Flip y-axis
        title = list(text = "Cell Deck % - UV Check",x =0,  y = 0.95) , 
        margin = list(t = 100)
      )
  
  return(p)
}

#UV Check (2)
plot_subsurface <- function(data_react, idx, max_idx, lake_name) {
  data <- data_react
  EDZ <- data[,c(2:5)]
  ED0 <- data[,c(10:13)]
  deck_cell <- (EDZ/ED0) * 100
  deck_cell$Depth <- data$Depth
  subsurface_means <- colMeans( deck_cell[(idx-28):(idx-24),] )
  albedo <- 100 - subsurface_means
  subsurface <- data.frame(matrix(ncol = 4, nrow = nrow(deck_cell)))
  names(subsurface) <- c("Edz305nm", "Edz320nm", "Edz380nm", "EdzPAR")
  subsurface$Edz305nm <- ( deck_cell$Edz305nm / subsurface_means[1] ) * 100
  subsurface$Edz320nm <- ( deck_cell$Edz320nm / subsurface_means[2] ) * 100
  subsurface$Edz380nm <- ( deck_cell$Edz380nm / subsurface_means[3] ) * 100
  subsurface$EdzPAR <- ( deck_cell$EdzPAR / subsurface_means[4] ) * 100
  subsurface$Depth <- data$Depth
  subsurface_long <- subsurface[idx:max_idx,] %>%
    pivot_longer(cols = names(subsurface)[1:4], values_to = "Values", names_to = "Variable")
  p <- plot_ly(data = subsurface_long, x = ~log(Values), y = ~Depth,
               color = ~Variable, symbol = ~Variable,
               type = "scatter", mode = "markers") %>%
    layout(
      xaxis = list(side = 'top', title = "Log(% of Subsurface)"),   # Move x-axis to the top
      yaxis = list(autorange = 'reversed', title = "Depth"),  #Flip y-axis
      title = list(text = "Subsurface % - UV Check",x =0,  y = 0.95) , 
      margin = list(t = 100)
    )
  
  return (p)
  
}

R2_base <- function(data) {
  max_idx <- c(0.5, 1)
  lr_data.5 <- data[which(data$Depth <= .5),]
  lr_data1 <- data[which(data$Depth <= 1),]
  
  lr_.5 <-  summary(lm(log(lr_data.5[,1])~lr_data.5[,2]))
  lr_1 <- summary(lm(log(lr_data1[,1])~lr_data1[,2]))
  
  if  (lr_.5[["r.squared"]] > lr_1[["r.squared"]] ) {
    return ( c(round(lr_.5[["r.squared"]],3), -round(lr_.5[["coefficients"]][2],3), .5) ) #KD_arr[3] = .5, for display
  } else {
    return ( c(round(lr_1[["r.squared"]],3), -round(lr_1[["coefficients"]][2],3), 1)) #KD_arr[3] = 1, for display
  }
}

R2 <- function(data, min_x, max_x){ #use min and max depth value as domain
  if (max_x == -1) { #PAR
    lr_data <- data
  } else {
    lr_data <- data[which( (data$Depth >= min_x) & (data$Depth <= max_x) ), ]
  }
  
  lr <- summary( lm(log(lr_data[,1])~lr_data[,2]) )
  return( c(round(lr[["r.squared"]],3), -round(lr[["coefficients"]][2],3)) )
}

ui <- fluidPage(
  headerPanel("BIC - Calculating 1%, 10% and normal extinction coefficients (Kd)"), 
  tabsetPanel(
    tabPanel("Offset Row",
       sidebarLayout(
         sidebarPanel(width = 10,
           fileInput("upload", label = "BIC File", multiple = FALSE,
                     accept = c("text/csv", 
                                "text/comma-separated-values,text/plain", 
                                ".csv")), 
           numericInput("offset_idx", "Pick the row index where EDZ (r) and Depth (b) change systematically:", value = 0, width = "60%")
         ), 
         mainPanel(
           verbatimTextOutput("text_check"),
           plotlyOutput("depth_EDZ_graph")
         )
       )
    ), 
    tabPanel("UV Check",
             verbatimTextOutput("Note_UV"),
             hr(),
             plotlyOutput("UV_deck"),
             hr(),
             plotlyOutput("UV_subsurface"),
             ),
    
    tabPanel("Regression Profiles",
             verbatimTextOutput("Note_Regress"),
             hr(), 
             textOutput("check"),
             splitLayout(
               cellWidths = c("50%","50%"), plotOutput("p1", brush = "p1_brush"), plotOutput("p2", brush = "p2_brush")
             ), 
             splitLayout(
               cellwidths = c("50%", "50%"), verbatimTextOutput("p1_R2"), verbatimTextOutput("p2_R2")
             ),
             splitLayout(
               cellWidths = c("50%","50%"), plotOutput("p3", brush = "p3_brush"), plotOutput("p4", brush = "p4_brush")
             ), 
             splitLayout(
               cellwidths = c("50%", "50%"), verbatimTextOutput("p3_R2"), verbatimTextOutput("p4_R2")
             )

    ),
    tabPanel("Kd Calculated", 
             verbatimTextOutput("KD_check"),
             downloadButton("downloadData", "Download"))
  )
)

server <- function(input, output, session) {
  
  #set INITIAL CONDITION BEFORE ANY INPUT HAS EVEN BEEN INPUTTED. Helps to avoid any potential problems
  #--------REACTIVE--------
  #used for labeling
  lake_name <- reactive({
    req(input$upload)
    get_name(input$upload$name)
  })
  #full original data + EDZ + ED0
  input_data <- reactive({
    req(input$upload)
    parse_file(input$upload$datapath)
  })
  #corrected data
  data <- reactive({
    req(input$upload)
    req(input$offset_idx)
    #browser() #acts as a breakpoint 
    data <- input_data()
    data$Depth <- data$Depth - data$Depth[input$offset_idx] #corrected depth
    data
  })
  #irradiance @ depth 'z' = microwats/cm^2/nm, represents the power per unit area of electromagnetic radiation 
  EDZ <- reactive({
    req(input$upload)
    req(input$offset_idx)
    data()[,2:5]
  })

  imx <- reactive({ 
    req(input$upload)
    req(input$offset_idx)
    which(data()[,"Depth"] == max(data()[,"Depth"], na.rm = T))
  })
  


  
  #--------OUTPUT--------
  #Depth graph to get offset_idx
  output$depth_EDZ_graph <- renderPlotly({
    data_long <- input_data()[,c(3,7)] %>%
      mutate(Depth_plot = Depth, idx = seq(nrow(input_data()))) %>%
      pivot_longer(cols = c(Edz320nm, Depth), names_to = "Variable", values_to = "Values") 
    
    plot_ly(data_long, x = ~idx, y = ~Values, 
            type = "scatter", mode = "markers",
            color = ~Variable, 
            colors = c('#636EFA', '#EF553B'))
    
  })
  
  #printing lake name
  output$text_check <- renderText({
    req(input$upload)
    lake_name()
  })
  
  output$Note_UV <- renderText({
    "Note: Smaller wavelengths will trail off shallower than PAR."
    })
  
  #Plot cell_deck % UV check
  output$UV_deck <- renderPlotly({ 
    req(input$upload)
    req(input$offset_idx)
    plot_deck(data(), input$offset_idx, imx(), lake_name())
  })
  
  #Plot subsurface % UV check
  output$UV_subsurface <- renderPlotly({ #renderPlotly
    req(input$upload)
    req(input$offset_idx)
    plot_subsurface(data(), input$offset_idx, imx(), lake_name())
  })
  
  output$Note_Regress <- renderText({
    "Highlight the points where the highest linear relationship is.\nOtherwise, base linear relationships are bound to .5 and 1m for 305, 320 & 380nm"
  })
  
  #-----REGRESSION-----
  #Regress Profiles
  colors <- list("Edz305nm" = "#621F87", 
                 "Edz320nm" = "#2E1F87", 
                 "Edz380nm" = "blue", 
                 "EdzPAR" = "red")
  
  
  #305nm
  KD305 <- reactiveVal()
  output$p1 <- renderPlot({
    req(input$upload)
    req(input$offset_idx)
    arr <- data()[input$offset_idx:imx(),c(2,7)] 
    ggplot(arr, aes(x = Depth, y = log(Edz305nm))) + 
      geom_point(color = colors$Edz305nm) + 
      labs(y = "Edz305nm", title = "Edz305nm")
  }) 
  
  output$p1_R2 <- renderPrint({
    req(input$upload)
    req(input$offset_idx)
    arr <- data()[input$offset_idx:imx(),c(2,7)] #as nm increases, change 2-->3, then -->4, etc..
    check <- input$p1_brush
    if (is.null(check)) {
      KD_arr <- R2_base(arr)
      KD305(KD_arr[2])
      paste("R2:",KD_arr[1], ", Kd:", KD_arr[2], ", Depth-Range: [ 0 ,",KD_arr[3],"]")
    } else {
      x_min <- check$xmin
      x_max <- check$xmax
      KD_arr <- R2(arr, x_min, x_max)
      KD305(KD_arr[2])
      paste("R2:",KD_arr[1], ", Kd:", KD_arr[2],", Depth-Range (m): [",round(x_min,3),",",round(x_max,3),"]")
    }

  })
  
  #320nm
  KD320 <- reactiveVal()
  output$p2 <- renderPlot({
    req(input$upload)
    req(input$offset_idx)
    arr <- data()[input$offset_idx:imx(),c(3,7)] 
    ggplot(arr, aes(x = Depth, y = log(Edz320nm))) + 
      geom_point(color = colors$Edz320nm) + 
      labs(y = "Edz320nm", title = "Edz320nm")
  })
  
  output$p2_R2 <- renderPrint({
    req(input$upload)
    req(input$offset_idx)
    arr <- data()[input$offset_idx:imx(),c(3,7)] 
    check <- input$p2_brush
    if (is.null(check)) {
      KD_arr <- R2_base(arr)
      KD320(KD_arr[2])
      paste("R2:",KD_arr[1], ", Kd:", KD_arr[2], ", Depth-Range: [ 0 ,",KD_arr[3],"]")
    } else {
      x_min <- check$xmin
      x_max <- check$xmax
      KD_arr <- R2(arr, x_min, x_max)
      KD320(KD_arr[2])
      paste("R2:",KD_arr[1], ", Kd:", KD_arr[2],", Depth-Range (m): [",round(x_min,3),",",round(x_max,3),"]")
    }
    

  })
  
  #380nm
  KD380 <- reactiveVal()
  output$p3 <- renderPlot({
    req(input$upload)
    req(input$offset_idx)
    arr <- data()[input$offset_idx:imx(),c(4,7)] 
    ggplot(arr, aes(x = Depth, y = log(Edz380nm))) + 
      geom_point(color = colors$Edz380nm) + 
      labs(y = "Edz380nm", title = "Edz380nm")
  })
  
  output$p3_R2 <- renderPrint({
    req(input$upload)
    req(input$offset_idx)
    arr <- data()[input$offset_idx:imx(),c(4,7)] 
    check <- input$p3_brush
    if (is.null(check)) {
      KD_arr <- R2_base(arr)
      KD380(KD_arr[2])
      paste("R2:",KD_arr[1], ", Kd:", KD_arr[2], ", Depth-Range: [ 0 ,",KD_arr[3],"]")
    } else {
      x_min <- check$xmin
      x_max <- check$xmax
      KD_arr <- R2(arr, x_min, x_max)
      KD380(KD_arr[2])
      paste("R2:",KD_arr[1], ", Kd:", KD_arr[2],", Depth-Range (m): [",round(x_min,3),",",round(x_max,3),"]")
    }
    
    
  })
  
  #PAR - base calculation uses whole depth profile for PAR
  KDPAR <- reactiveVal()
  output$p4 <- renderPlot({
    req(input$upload)
    req(input$offset_idx)
    arr <- data()[input$offset_idx:imx(),c(5,7)] 
    ggplot(arr, aes(x = Depth, y = log(EdzPAR))) + 
      geom_point(color = colors$EdzPAR) + 
      labs(y = "EdzPAR", title = "EdzPAR")
  })
  
  output$p4_R2 <- renderPrint({
    req(input$upload)
    req(input$offset_idx)
    arr <- data()[input$offset_idx:imx(),c(5,7)] 
    check <- input$p4_brush
    if (is.null(check)) {
      KD_arr <- R2(arr, 0, -1)
      KDPAR(KD_arr[2])
      paste("R2:",KD_arr[1], ", Kd:", KD_arr[2], ", Depth-Range: [ 0 ,",round(arr[nrow(arr), 'Depth'],3),"]")
    } else {
      x_min <- check$xmin
      x_max <- check$xmax
      KD_arr <- R2(arr, x_min, x_max)
      KDPAR(KD_arr[2])
      paste("R2:",KD_arr[1], ", Kd:", KD_arr[2],", Depth-Range (m): [",round(x_min,3),",",round(x_max,3),"]")
    }
    
    
  })
  return_data <- reactive({
    KD <- c(KD305(), KD320(), KD380(), KDPAR())
    data.frame(wavelength = c("305", "320", "380", "PAR"), Kd = KD, 
               Kd10 = log(10) / KD, KD1 = log(100) / KD)
              
  })
  
  output$KD_check <- renderPrint({
    return_data()
  })
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("KD_",lake_name(), "_",Sys.Date(), ".csv", sep = "" )
    }, 
    content = function(file) {
      write.csv(return_data(), file, row.names = F)
    }
    
  )
                            
  
}

shinyApp(ui, server)







#https://mastering-shiny.org/action-graphics.html
#https://plotly.com/r/subplots/
#https://stackoverflow.com/questions/34384907/how-can-put-multiple-plots-side-by-side-in-shiny-r
#also run on IDEA cluster under folder ~/ShinyApps/BIC_profile/app.R
#create new objects for each graph. Problem with obtaining brushed points is that the index has to correspond

#to the data being plotted. The data being plotted is in long format so create long dataframes
#for each regression wavelength
#