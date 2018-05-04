#30/4/18
library(shiny)
#source("global.R")
markers <- mget(load("data/all_markers.RData")) 
markers_labels <- mget(load("data/all_markers_selin_labels.RData"))
# version of markers "../..Hydra/sjessa/single_cell/2018-02_mouse_dev/objects/all_markers.RData"
references <- mget(load("data/references.RData"))
markers_top50 <- mget(load("data/all_markers_selin_labelstop50.RData"))


#### -------------- UI ------------------- ####
ui <- fluidPage(
   
   # Application title
   titlePanel("Systematic Cell Type Identification"),
   
   # Sidebar with input select for input data vs reference
   sidebarLayout(
      sidebarPanel(
        
        helpText("Here you can select two datasets from this repository to comare marker genes identified
                 in different datasets or upload a new dataset to compare to any dataset in this repository"),
        
        ## Sample_1 Selection
        selectInput(inputId = "s1",
        label = h3("Select a dataset for Sample_1"),
        choices = list("Zhang_2015", "Lein_2007", "Hochgerner_2018(early)", "Hochgerner_2018(mature)", 
                       "Usoskin_2015", "Daneman_2010", "Habib_2016(a)", "Habib_2016(b)", "Chen_2017", 
                       "Gokce_2018", "Cahoy_2008", "Zeisal_2015", "Haber_2017",  "Shekhar_2016", 
                       "Zhong_2018", "Newman_2015", "Zeisal_2018", "Darmanis_2015",  "Macosko_2015", 
                       
                       #our markers version: 2018-02_mouse_dev
                       "ct_e12", "ct_e15", "ct_p0", "po_e12", "po_e15", "po_p0",
                       "ct_e12_top50", "ct_e15_top50", "ct_p0_top50", "po_e12_top50", "po_e15_top50", "po_p0_top50"
                       ), 
          selected = "Zhang_2015"),
        
        tags$hr(), ##horizontal line
      
        ## Sample_2 selection
        selectInput(inputId = "s2",
                  label = h3("Select a dataset for Sample_2"),
                  choices = list("Zhang_2015", "Lein_2007", "Hochgerner_2018(early)", "Hochgerner_2018(mature)", 
                                 "Usoskin_2015", "Daneman_2010", "Habib_2016(a)", "Habib_2016(b)", "Chen_2017", 
                                 "Gokce_2018", "Cahoy_2008", "Zeisal_2015", "Haber_2017",  "Shekhar_2016", 
                                 "Zhong_2018", "Newman_2015", "Zeisal_2018", "Darmanis_2015",  "Macosko_2015", 
                                 
                                 #our markers version: 2018-02_mouse_dev
                                 "ct_e12", "ct_e15", "ct_p0", "po_e12", "po_e15", "po_p0",
                                 "ct_e12_top50", "ct_e15_top50", "ct_p0_top50", "po_e12_top50", "po_e15_top50", "po_p0_top50"), 
                  selected = "ct_e12"),
        
        ##or upload a sample
        # code snippet from https://shiny.rstudio.com/articles/upload.html
        
        fileInput(inputId = "file1", 
                  label = h3("Or upload a dataset for Sample_2"),
                  accept = c("text/csv", "text/comma-separated-values,text/plain",
                             ".csv")),
        helpText("Files can be .txt or .csv. Must have columns titled 'cluster' (with cluster numbers or 
                 labels in the column), 'external_gene_name', and 'ensembl_gene_id' (with lists of genes)"),
        radioButtons(inputId = "sep", label = "Delimiter",
                     choices = c(Comma = ",",
                                 Semicolon = ";",
                                 Tab = "\t"),
                     selected = ","),
        ##
        tags$hr(),
        
        
        ## Choose which will be the denominator for heatmap proportions
        radioButtons(inputId = "denominator",
                           label = "Choose denominator",
                           choices = list("Sample_1", "Sample_2"),
                           selected = "Sample_1"),
        
        ## Option to show labels on our datasets
        radioButtons(inputId = "label_option",
                     label = "Add labels",
                     choices = list("Labels", "No_labels"),
                     selected = "Labels"),
        
        ## Use gene symbols or ensemble IDs
        radioButtons(inputId = "sym_ens",
                     label = "Use Gene Symbols or Ensembl IDs?",
                     choices = list("Gene_Symbols", "Ensemble_IDs"),
                     selected = "Gene_Symbols"),
        
        tags$hr(),
        
        ###these are the clusters that are compared
        helpText("The first selection will determine the barplot, the second selection is intersected to display the shared genes."),
        htmlOutput("selectUI"),
        htmlOutput("selectC2")
        
   ),
          # Show a heatplot of shared genes
      mainPanel(
        fluidRow(plotOutput("heatPlot")),
        #add a line and white space between plots
        tags$hr(),
        br(), br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(), 
        fluidRow(plotOutput("barPlot")),
        tags$hr(),
        textOutput("text1"),
        tags$hr(),
        textOutput("refs1"),
        textOutput("refs2")
      ))
      
   ### sidebar for cluster analysis
   
   )




#### -------------- SERVER ------------------- ####
server <- function(input, output) {
  source('global.R')
  
  ## Draw the heatmaps
  output$heatPlot <- renderPlot({
     
     #Show labels option
     if (input$label_option=="Labels"){
      sample_1 <- #reactive({
        switch(input$s1,
                      "Zhang_2015" = zhang_2015, "Lein_2007" = lein_2007, "Hochgerner_2018(early)" = hg_early_2018, 
                      "Hochgerner_2018(mature)" = hg_mature_2018, "Usoskin_2015" = usoskin_2015, "Daneman_2010" = daneman_2010,
                      "Habib_2016(a)" = habib_2016_1, "Habib_2016(b)" = habib_2016_2, "Chen_2017" = chen_2017, 
                      "Gokce_2018" = gokce_2018, "Cahoy_2008" = cahoy_2008, "Zeisal_2015" = zeisal_2015, 
                      "Haber_2017" = haber_2017, "Shekhar_2016" = shekhar_2016, "Zhong_2018" = zhong_2018, 
                      "Newman_2015" = newman_2015, "Zeisal_2018" = zeisal_2018, "Darmanis_2015" = darmanis_2018, 
                      "Macosko_2015" = macosko_2015,
                      #Our datasets with Selin's labels
                      "ct_e12" = markers_ct_e12_labels, "ct_e15" = markers_ct_e15_labels, "ct_p0" = markers_ct_p0_labels,
                      "po_e12" = markers_po_e12_labels, "po_e15" = markers_po_e15_labels, "po_p0" = markers_po_p0_labels,
               "ct_e12_top50" = markers_ct_e12_labels_top50, "ct_e15_top50" = markers_ct_e15_labels_top50, 
               "ct_p0_top50" = markers_ct_p0_labels_top50, "po_e12_top50" = markers_po_e12_labels_top50, 
               "po_e15_top50" = markers_po_e15_labels_top50, "po_p0_top50" = markers_po_p0_labels_top50
                    )
      
      if (is.null(input$file1)==TRUE){
        sample_2 <- switch(input$s2,
                         "Zhang_2015" = zhang_2015, "Lein_2007" = lein_2007, "Hochgerner_2018(early)" = hg_early_2018, 
                         "Hochgerner_2018(mature)" = hg_mature_2018, "Usoskin_2015" = usoskin_2015, "Daneman_2010" = daneman_2010,
                         "Habib_2016(a)" = habib_2016_1, "Habib_2016(b)" = habib_2016_2, "Chen_2017" = chen_2017, 
                         "Gokce_2018" = gokce_2018, "Cahoy_2008" = cahoy_2008, "Zeisal_2015" = zeisal_2015, 
                         "Haber_2017" = haber_2017, "Shekhar_2016" = shekhar_2016, "Zhong_2018" = zhong_2018, 
                         "Newman_2015" = newman_2015, "Zeisal_2018" = zeisal_2018, "Darmanis_2015" = darmanis_2018, 
                         "Macosko_2015" = macosko_2015,
                         #Our datasets with Selin's labels
                         "ct_e12" = markers_ct_e12_labels, "ct_e15" = markers_ct_e15_labels, "ct_p0" = markers_ct_p0_labels,
                         "po_e12" = markers_po_e12_labels, "po_e15" = markers_po_e15_labels, "po_p0" = markers_po_p0_labels,
               "ct_e12_top50" = markers_ct_e12_labels_top50, "ct_e15_top50" = markers_ct_e15_labels_top50, "ct_p0_top50" = markers_ct_p0_labels_top50, 
               "po_e12_top50" = markers_po_e12_labels_top50, "po_e15_top50" = markers_po_e15_labels_top50, "po_p0_top50" = markers_po_p0_labels_top50
      )
      }
      else if (is.null(input$file1)==FALSE){
        sample_2 <- read.csv(input$file1$datapath)
      }
      
      
     }
     
     else if (input$label_option=="No_labels"){
       sample_1 <- switch(input$s1,
                          "Zhang_2015" = zhang_2015,  "Lein_2007" = lein_2007,  "Hochgerner_2018(early)" = hg_early_2018, 
                          "Hochgerner_2018(mature)" = hg_mature_2018, "Usoskin_2015" = usoskin_2015, "Daneman_2010" = daneman_2010,
                          "Habib_2016(a)" = habib_2016_1, "Habib_2016(b)" = habib_2016_2, "Chen_2017" = chen_2017, 
                          "Gokce_2018" = gokce_2018, "Cahoy_2008" = cahoy_2008, "Zeisal_2015" = zeisal_2015, 
                          "Haber_2017" = haber_2017, "Shekhar_2016" = shekhar_2016, "Zhong_2018" = zhong_2018, 
                          "Newman_2015" = newman_2015, "Zeisal_2018" = zeisal_2018, "Darmanis_2015" = darmanis_2018, 
                          "Macosko_2015" = macosko_2015,
                          #Our datasets without labels
                          "ct_e12" = markers_ct_e12, "ct_e15" = markers_ct_e15, "ct_p0" = markers_ct_p0, "po_e12" = markers_po_e12,
                          "po_e15" = markers_po_e15, "po_p0" = markers_po_p0, "ct_e12_top50" = markers_ct_e12_labels_top50, 
                          "ct_e15_top50" = markers_ct_e15_labels_top50, "ct_p0_top50" = markers_ct_p0_labels_top50, "po_e12_top50" = markers_po_e12_labels_top50, 
                          "po_e15_top50" = markers_po_e15_labels_top50, "po_p0_top50" = markers_po_p0_labels_top50)
       
       if (is.null(input$file1)==TRUE){
        sample_2 <- switch(input$s2,
                          "Zhang_2015" = zhang_2015, "Lein_2007" = lein_2007, "Hochgerner_2018(early)" = hg_early_2018, 
                          "Hochgerner_2018(mature)" = hg_mature_2018, "Usoskin_2015" = usoskin_2015, "Daneman_2010" = daneman_2010,
                          "Habib_2016(a)" = habib_2016_1, "Habib_2016(b)" = habib_2016_2, "Chen_2017" = chen_2017, "Gokce_2018" = gokce_2018, "Cahoy_2008" = cahoy_2008, 
                          "Zeisal_2015" = zeisal_2015, "Haber_2017" = haber_2017, 
                          "Shekhar_2016" = shekhar_2016, "Zhong_2018" = zhong_2018, "Newman_2015" = newman_2015, 
                          "Zeisal_2018" = zeisal_2018, "Darmanis_2015" = darmanis_2018, "Macosko_2015" = macosko_2015,
                          #Our datasets without labels
                          "ct_e12" = markers_ct_e12, "ct_e15" = markers_ct_e15, "ct_p0" = markers_ct_p0, "po_e12" = markers_po_e12,
                          "po_e15" = markers_po_e15,  "po_p0" = markers_po_p0,
                          "ct_e12_top50" = markers_ct_e12_labels_top50, "ct_e15_top50" = markers_ct_e15_labels_top50, "ct_p0_top50" = markers_ct_p0_labels_top50, 
                          "po_e12_top50" = markers_po_e12_labels_top50, "po_e15_top50" = markers_po_e15_labels_top50, "po_p0_top50" = markers_po_p0_labels_top50)
       }
       else if (is.null(input$file1)==FALSE){
        sample_2 <- read.csv(input$file1$datapath)
       }
       }
     
    ## Radio buttons determine the denominator and which gene IDs wil be used
    if (input$sym_ens=="Gene_Symbols"){
      markers_overlap <- computeMarkersOverlap(sample_1 = sample_1,
                                               sample_2 = sample_2,
                                               denom = input$denominator)
    }
    else if (input$sym_ens=="Ensemble_IDs") {
      markers_overlap <- computeMarkersOverlap_ensid(sample_1 = sample_1,
                                                     sample_2 = sample_2,
                                                     denom = input$denominator)
    }
    
    ## return the names of the sample chosen
    sample_name_s1 <- switch(input$s1,
                              "Zhang_2015" = "Zhang_2015", "Lein_2007" = "Lein_2007", 
                              "Hochgerner_2018(early)" = "Hochgerner_2018(early)", "Hochgerner_2018(mature)" = "Hochgerner_2018(mature)", 
                              "Usoskin_2015" = "Usoskin_2015", "Daneman_2010" = "Daneman_2010", "Habib_2016(a)" = "Habib_2016(a)",
                              "Habib_2016(b)" = "Habib_2016(b)", "Chen_2017" = "Chen_2017", "Gokce_2018" = "Gokce_2018",
                              "Cahoy_2008" = "Cahoy_2008",  "Zeisal_2015" = "Zeisal_2015",  "Haber_2017" = "Haber_2017", 
                              "Shekhar_2016" = "Shekhar_2016", "Zhong_2018" = "Zhong_2018", "Newman_2015" = "Newman_2015",
                              "Zeisal_2018" = "Zeisal_2018",  "Darmanis_2015" = "Darmanis_2015",  "Macosko_2015" = "Macosko_2015", 
                              
                              #our markers version: 2018-02_mouse_dev
                              "ct_e12" = "ct_e12", "ct_e15" = "ct_e15", "ct_p0" = "ct_p0", "po_e12" = "po_e12", "po_e15" = "po_e15", "po_p0" = "po_p0",
                              "ct_e12_top50" ="ct_e12_top50", "ct_e15_top50" = "ct_e15_top50", "ct_p0_top50" = "ct_p0_top50", 
                              "po_e12_top50" = "po_e12_top50", "po_e15_top50" = "po_e15_top50", "po_p0_top50" = "po_p0_top50") 
    
    if (is.null(input$file1)==TRUE){
      sample_name_s2 <- switch(input$s2,
                             "Zhang_2015" = "Zhang_2015", "Lein_2007" = "Lein_2007", 
                             "Hochgerner_2018(early)" = "Hochgerner_2018(early)", "Hochgerner_2018(mature)" = "Hochgerner_2018(mature)", 
                             "Usoskin_2015" = "Usoskin_2015", "Daneman_2010" = "Daneman_2010", "Habib_2016(a)" = "Habib_2016(a)",
                             "Habib_2016(b)" = "Habib_2016(b)", "Chen_2017" = "Chen_2017", "Gokce_2018" = "Gokce_2018",
                             "Cahoy_2008" = "Cahoy_2008",  "Zeisal_2015" = "Zeisal_2015",  "Haber_2017" = "Haber_2017", 
                             "Shekhar_2016" = "Shekhar_2016", "Zhong_2018" = "Zhong_2018", "Newman_2015" = "Newman_2015",
                             "Zeisal_2018" = "Zeisal_2018",  "Darmanis_2015" = "Darmanis_2015",  "Macosko_2015" = "Macosko_2015", 
                             
                             #our markers version: 2018-02_mouse_dev
                             "ct_e12" = "ct_e12", "ct_e15" = "ct_e15", "ct_p0" = "ct_p0", "po_e12" = "po_e12", "po_e15" = "po_e15", "po_p0" = "po_p0",
                             "ct_e12_top50" ="ct_e12_top50", "ct_e15_top50" = "ct_e15_top50", "ct_p0_top50" = "ct_p0_top50", 
                             "po_e12_top50" = "po_e12_top50", "po_e15_top50" = "po_e15_top50", "po_p0_top50" = "po_p0_top50") 
    }
    else if (is.null(input$file1)==FALSE){
      sample_name_s2 <- "User Input"
    }
    
    
    
    p <- heatmapMarkersOverlap(markers_overlap, sample_1 = sample_1, sample_2 = sample_2)
    final_p <- p + labs(x=sample_name_s2, y =sample_name_s1)
    final_p
     }, height=700, units="px") # end of heatmap arguments
   

     
     
   ## Draw the barplots
  
  ##But first, render the UI for Cluster 1 selection
  output$selectUI <- renderUI({
    #Show labels option
    if (input$label_option=="Labels"){
      sample_1 <- switch(input$s1,
                         "Zhang_2015" = zhang_2015, "Lein_2007" = lein_2007, "Hochgerner_2018(early)" = hg_early_2018, 
                         "Hochgerner_2018(mature)" = hg_mature_2018, "Usoskin_2015" = usoskin_2015, "Daneman_2010" = daneman_2010,
                         "Habib_2016(a)" = habib_2016_1, "Habib_2016(b)" = habib_2016_2, "Chen_2017" = chen_2017, 
                         "Gokce_2018" = gokce_2018, "Cahoy_2008" = cahoy_2008, "Zeisal_2015" = zeisal_2015, 
                         "Haber_2017" = haber_2017, "Shekhar_2016" = shekhar_2016, "Zhong_2018" = zhong_2018, 
                         "Newman_2015" = newman_2015, "Zeisal_2018" = zeisal_2018, "Darmanis_2015" = darmanis_2018, 
                         "Macosko_2015" = macosko_2015,
                         #Our datasets with Selin's labels
                         "ct_e12" = markers_ct_e12_labels, "ct_e15" = markers_ct_e15_labels, "ct_p0" = markers_ct_p0_labels,
                         "po_e12" = markers_po_e12_labels, "po_e15" = markers_po_e15_labels, "po_p0" = markers_po_p0_labels,
                         "ct_e12_top50" = markers_ct_e12_labels_top50, "ct_e15_top50" = markers_ct_e15_labels_top50, 
                         "ct_p0_top50" = markers_ct_p0_labels_top50, "po_e12_top50" = markers_po_e12_labels_top50, 
                         "po_e15_top50" = markers_po_e15_labels_top50, "po_p0_top50" = markers_po_p0_labels_top50
      )
      if (is.null(input$file1)==TRUE){
        sample_2 <- switch(input$s2,
                         "Zhang_2015" = zhang_2015, "Lein_2007" = lein_2007, "Hochgerner_2018(early)" = hg_early_2018, 
                         "Hochgerner_2018(mature)" = hg_mature_2018, "Usoskin_2015" = usoskin_2015, "Daneman_2010" = daneman_2010,
                         "Habib_2016(a)" = habib_2016_1, "Habib_2016(b)" = habib_2016_2, "Chen_2017" = chen_2017, 
                         "Gokce_2018" = gokce_2018, "Cahoy_2008" = cahoy_2008, "Zeisal_2015" = zeisal_2015, 
                         "Haber_2017" = haber_2017, "Shekhar_2016" = shekhar_2016, "Zhong_2018" = zhong_2018, 
                         "Newman_2015" = newman_2015, "Zeisal_2018" = zeisal_2018, "Darmanis_2015" = darmanis_2018, 
                         "Macosko_2015" = macosko_2015,
                         #Our datasets with Selin's labels
                         "ct_e12" = markers_ct_e12_labels, "ct_e15" = markers_ct_e15_labels, "ct_p0" = markers_ct_p0_labels,
                         "po_e12" = markers_po_e12_labels, "po_e15" = markers_po_e15_labels, "po_p0" = markers_po_p0_labels,
                         "ct_e12_top50" = markers_ct_e12_labels_top50, "ct_e15_top50" = markers_ct_e15_labels_top50, "ct_p0_top50" = markers_ct_p0_labels_top50, 
                         "po_e12_top50" = markers_po_e12_labels_top50, "po_e15_top50" = markers_po_e15_labels_top50, "po_p0_top50" = markers_po_p0_labels_top50
      )
      }
      else if (is.null(input$file1)==FALSE){
        sample_2 <- read.csv(input$file1$datapath)
      }
    }
    
    else if (input$label_option=="No_labels"){
      sample_1 <- switch(input$s1,
                         "Zhang_2015" = zhang_2015, "Lein_2007" = lein_2007, "Hochgerner_2018(early)" = hg_early_2018, 
                         "Hochgerner_2018(mature)" = hg_mature_2018, "Usoskin_2015" = usoskin_2015,
                         "Daneman_2010" = daneman_2010,  "Habib_2016(a)" = habib_2016_1, "Habib_2016(b)" = habib_2016_2, 
                         "Chen_2017" = chen_2017, "Gokce_2018" = gokce_2018, "Cahoy_2008" = cahoy_2008,
                         "Zeisal_2015" = zeisal_2015, "Haber_2017" = haber_2017, "Shekhar_2016" = shekhar_2016, 
                         "Zhong_2018" = zhong_2018, "Newman_2015" = newman_2015, "Zeisal_2018" = zeisal_2018, 
                         "Darmanis_2015" = darmanis_2018, "Macosko_2015" = macosko_2015,
                         #Our datasets without labels
                         "ct_e12" = markers_ct_e12, "ct_e15" = markers_ct_e15, "ct_p0" = markers_ct_p0,
                         "po_e12" = markers_po_e12, "po_e15" = markers_po_e15, "po_p0" = markers_po_p0,
                         "ct_e12_top50" = markers_ct_e12_labels_top50, "ct_e15_top50" = markers_ct_e15_labels_top50, 
                         "ct_p0_top50" = markers_ct_p0_labels_top50, "po_e12_top50" = markers_po_e12_labels_top50, 
                         "po_e15_top50" = markers_po_e15_labels_top50, "po_p0_top50" = markers_po_p0_labels_top50
      )
      if (is.null(input$file1)==TRUE){
        sample_2 <- switch(input$s2,
                         "Zhang_2015" = zhang_2015, "Lein_2007" = lein_2007, "Hochgerner_2018(early)" = hg_early_2018, 
                         "Hochgerner_2018(mature)" = hg_mature_2018, "Usoskin_2015" = usoskin_2015,
                         "Daneman_2010" = daneman_2010,  "Habib_2016(a)" = habib_2016_1, "Habib_2016(b)" = habib_2016_2, 
                         "Chen_2017" = chen_2017, "Gokce_2018" = gokce_2018, "Cahoy_2008" = cahoy_2008,
                         "Zeisal_2015" = zeisal_2015, "Haber_2017" = haber_2017, "Shekhar_2016" = shekhar_2016, 
                         "Zhong_2018" = zhong_2018, "Newman_2015" = newman_2015, "Zeisal_2018" = zeisal_2018, 
                         "Darmanis_2015" = darmanis_2018, "Macosko_2015" = macosko_2015,
                         #Our datasets without labels
                         "ct_e12" = markers_ct_e12, "ct_e15" = markers_ct_e15, "ct_p0" = markers_ct_p0, "po_e12" = markers_po_e12,
                         "po_e15" = markers_po_e15, "po_p0" = markers_po_p0, 
                         "ct_e12_top50" = markers_ct_e12_labels_top50, "ct_e15_top50" = markers_ct_e15_labels_top50, 
                         "ct_p0_top50" = markers_ct_p0_labels_top50, "po_e12_top50" = markers_po_e12_labels_top50, 
                         "po_e15_top50" = markers_po_e15_labels_top50, "po_p0_top50" = markers_po_p0_labels_top50
      )
      }
      else if (is.null(input$file1)==FALSE){
        sample_2 <- read.csv(input$file1$datapath)
      }
    }
    
    cluster_a_options <- switch(input$denominator, #cluster a can be from sample 1 or sample 2
                              "Sample_1" = names(split(sample_1, f = sample_1$cluster)),
                              "Sample_2" = names(split(sample_2, f = sample_2$cluster)))
    selectInput(inputId = "cluster_name",label = "Choose a cluster", choices = cluster_a_options)
  })
  
  
  ### Then render the UI for a cluster 2 selection 
  output$selectC2 <- renderUI({#Show labels option
    if (input$label_option=="Labels"){
      sample_1 <- switch(input$s1,
                         "Zhang_2015" = zhang_2015, "Lein_2007" = lein_2007, "Hochgerner_2018(early)" = hg_early_2018, 
                         "Hochgerner_2018(mature)" = hg_mature_2018, "Usoskin_2015" = usoskin_2015, "Daneman_2010" = daneman_2010,
                         "Habib_2016(a)" = habib_2016_1, "Habib_2016(b)" = habib_2016_2, "Chen_2017" = chen_2017, 
                         "Gokce_2018" = gokce_2018, "Cahoy_2008" = cahoy_2008, "Zeisal_2015" = zeisal_2015, 
                         "Haber_2017" = haber_2017, "Shekhar_2016" = shekhar_2016, "Zhong_2018" = zhong_2018, 
                         "Newman_2015" = newman_2015, "Zeisal_2018" = zeisal_2018, "Darmanis_2015" = darmanis_2018, 
                         "Macosko_2015" = macosko_2015,
                         #Our datasets with Selin's labels
                         "ct_e12" = markers_ct_e12_labels, "ct_e15" = markers_ct_e15_labels, "ct_p0" = markers_ct_p0_labels,
                         "po_e12" = markers_po_e12_labels, "po_e15" = markers_po_e15_labels, "po_p0" = markers_po_p0_labels,
                         "ct_e12_top50" = markers_ct_e12_labels_top50, "ct_e15_top50" = markers_ct_e15_labels_top50, 
                         "ct_p0_top50" = markers_ct_p0_labels_top50, "po_e12_top50" = markers_po_e12_labels_top50, 
                         "po_e15_top50" = markers_po_e15_labels_top50, "po_p0_top50" = markers_po_p0_labels_top50
      )
      if (is.null(input$file1)==TRUE){
        sample_2 <- switch(input$s2,
                         "Zhang_2015" = zhang_2015, "Lein_2007" = lein_2007, "Hochgerner_2018(early)" = hg_early_2018, 
                         "Hochgerner_2018(mature)" = hg_mature_2018, "Usoskin_2015" = usoskin_2015, "Daneman_2010" = daneman_2010,
                         "Habib_2016(a)" = habib_2016_1, "Habib_2016(b)" = habib_2016_2, "Chen_2017" = chen_2017, 
                         "Gokce_2018" = gokce_2018, "Cahoy_2008" = cahoy_2008, "Zeisal_2015" = zeisal_2015, 
                         "Haber_2017" = haber_2017, "Shekhar_2016" = shekhar_2016, "Zhong_2018" = zhong_2018, 
                         "Newman_2015" = newman_2015, "Zeisal_2018" = zeisal_2018, "Darmanis_2015" = darmanis_2018, 
                         "Macosko_2015" = macosko_2015,
                         #Our datasets with Selin's labels
                         "ct_e12" = markers_ct_e12_labels, "ct_e15" = markers_ct_e15_labels, "ct_p0" = markers_ct_p0_labels,
                         "po_e12" = markers_po_e12_labels, "po_e15" = markers_po_e15_labels, "po_p0" = markers_po_p0_labels,
                         "ct_e12_top50" = markers_ct_e12_labels_top50, "ct_e15_top50" = markers_ct_e15_labels_top50, "ct_p0_top50" = markers_ct_p0_labels_top50, 
                         "po_e12_top50" = markers_po_e12_labels_top50, "po_e15_top50" = markers_po_e15_labels_top50, "po_p0_top50" = markers_po_p0_labels_top50
      )
      }
      else if  (is.null(input$file1)==FALSE){
      sample_2 <- read.csv(input$file1$datapath)
      }
    }
    
    else if (input$label_option=="No_labels"){
      sample_1 <- switch(input$s1,
                         "Zhang_2015" = zhang_2015, "Lein_2007" = lein_2007, "Hochgerner_2018(early)" = hg_early_2018, 
                         "Hochgerner_2018(mature)" = hg_mature_2018, "Usoskin_2015" = usoskin_2015,
                         "Daneman_2010" = daneman_2010,  "Habib_2016(a)" = habib_2016_1, "Habib_2016(b)" = habib_2016_2, 
                         "Chen_2017" = chen_2017, "Gokce_2018" = gokce_2018, "Cahoy_2008" = cahoy_2008,
                         "Zeisal_2015" = zeisal_2015, "Haber_2017" = haber_2017, "Shekhar_2016" = shekhar_2016, 
                         "Zhong_2018" = zhong_2018, "Newman_2015" = newman_2015, "Zeisal_2018" = zeisal_2018, 
                         "Darmanis_2015" = darmanis_2018, "Macosko_2015" = macosko_2015,
                         #Our datasets without labels
                         "ct_e12" = markers_ct_e12, "ct_e15" = markers_ct_e15, "ct_p0" = markers_ct_p0,
                         "po_e12" = markers_po_e12, "po_e15" = markers_po_e15, "po_p0" = markers_po_p0,
                         "ct_e12_top50" = markers_ct_e12_labels_top50, "ct_e15_top50" = markers_ct_e15_labels_top50, 
                         "ct_p0_top50" = markers_ct_p0_labels_top50, "po_e12_top50" = markers_po_e12_labels_top50, 
                         "po_e15_top50" = markers_po_e15_labels_top50, "po_p0_top50" = markers_po_p0_labels_top50
      )
      if (is.null(input$file1)==TRUE){
        sample_2 <- switch(input$s2,
                         "Zhang_2015" = zhang_2015, "Lein_2007" = lein_2007, "Hochgerner_2018(early)" = hg_early_2018, 
                         "Hochgerner_2018(mature)" = hg_mature_2018, "Usoskin_2015" = usoskin_2015,
                         "Daneman_2010" = daneman_2010,  "Habib_2016(a)" = habib_2016_1, "Habib_2016(b)" = habib_2016_2, 
                         "Chen_2017" = chen_2017, "Gokce_2018" = gokce_2018, "Cahoy_2008" = cahoy_2008,
                         "Zeisal_2015" = zeisal_2015, "Haber_2017" = haber_2017, "Shekhar_2016" = shekhar_2016, 
                         "Zhong_2018" = zhong_2018, "Newman_2015" = newman_2015, "Zeisal_2018" = zeisal_2018, 
                         "Darmanis_2015" = darmanis_2018, "Macosko_2015" = macosko_2015,
                         #Our datasets without labels
                         "ct_e12" = markers_ct_e12, "ct_e15" = markers_ct_e15, "ct_p0" = markers_ct_p0, "po_e12" = markers_po_e12,
                         "po_e15" = markers_po_e15, "po_p0" = markers_po_p0, 
                         "ct_e12_top50" = markers_ct_e12_labels_top50, "ct_e15_top50" = markers_ct_e15_labels_top50, 
                         "ct_p0_top50" = markers_ct_p0_labels_top50, "po_e12_top50" = markers_po_e12_labels_top50, 
                         "po_e15_top50" = markers_po_e15_labels_top50, "po_p0_top50" = markers_po_p0_labels_top50
      )
      }
      if (is.null(input$file1)==FALSE){
        sample_2 <- read.csv(input$file1$datapath)
      }
    }
    cluster_b_options <- switch(input$denominator,
                                 "Sample_1" = names(split(sample_2, f = sample_2$cluster)),
                                 "Sample_2" = names(split(sample_1, f = sample_1$cluster)))
    
    selectInput(inputId = "cluster_2_name",label = "Choose a cluster", choices = cluster_b_options)
  })
  
  
  
  ### Now actually draw the barplot
   output$barPlot <- renderPlot({
   
   #Show labels option
     if (input$label_option=="Labels"){
       sample_1 <- switch(input$s1,
                "Zhang_2015" = zhang_2015, "Lein_2007" = lein_2007, "Hochgerner_2018(early)" = hg_early_2018, 
                "Hochgerner_2018(mature)" = hg_mature_2018, "Usoskin_2015" = usoskin_2015, "Daneman_2010" = daneman_2010,
                "Habib_2016(a)" = habib_2016_1, "Habib_2016(b)" = habib_2016_2, "Chen_2017" = chen_2017, 
                "Gokce_2018" = gokce_2018, "Cahoy_2008" = cahoy_2008, "Zeisal_2015" = zeisal_2015, 
                "Haber_2017" = haber_2017, "Shekhar_2016" = shekhar_2016, "Zhong_2018" = zhong_2018, 
                "Newman_2015" = newman_2015, "Zeisal_2018" = zeisal_2018, "Darmanis_2015" = darmanis_2018, 
                "Macosko_2015" = macosko_2015,
                #Our datasets with Selin's labels
                "ct_e12" = markers_ct_e12_labels, "ct_e15" = markers_ct_e15_labels, "ct_p0" = markers_ct_p0_labels,
                "po_e12" = markers_po_e12_labels, "po_e15" = markers_po_e15_labels, "po_p0" = markers_po_p0_labels,
                "ct_e12_top50" = markers_ct_e12_labels_top50, "ct_e15_top50" = markers_ct_e15_labels_top50, 
                "ct_p0_top50" = markers_ct_p0_labels_top50, "po_e12_top50" = markers_po_e12_labels_top50, 
                "po_e15_top50" = markers_po_e15_labels_top50, "po_p0_top50" = markers_po_p0_labels_top50
         )
       if (is.null(input$file1)==TRUE){
        sample_2 <- switch(input$s2,
                "Zhang_2015" = zhang_2015, "Lein_2007" = lein_2007, "Hochgerner_2018(early)" = hg_early_2018, 
                "Hochgerner_2018(mature)" = hg_mature_2018, "Usoskin_2015" = usoskin_2015, "Daneman_2010" = daneman_2010,
                "Habib_2016(a)" = habib_2016_1, "Habib_2016(b)" = habib_2016_2, "Chen_2017" = chen_2017, 
                "Gokce_2018" = gokce_2018, "Cahoy_2008" = cahoy_2008, "Zeisal_2015" = zeisal_2015, 
                "Haber_2017" = haber_2017, "Shekhar_2016" = shekhar_2016, "Zhong_2018" = zhong_2018, 
                "Newman_2015" = newman_2015, "Zeisal_2018" = zeisal_2018, "Darmanis_2015" = darmanis_2018, 
                "Macosko_2015" = macosko_2015,
                #Our datasets with Selin's labels
                "ct_e12" = markers_ct_e12_labels, "ct_e15" = markers_ct_e15_labels, "ct_p0" = markers_ct_p0_labels,
                "po_e12" = markers_po_e12_labels, "po_e15" = markers_po_e15_labels, "po_p0" = markers_po_p0_labels,
                "ct_e12_top50" = markers_ct_e12_labels_top50, "ct_e15_top50" = markers_ct_e15_labels_top50, "ct_p0_top50" = markers_ct_p0_labels_top50, 
                "po_e12_top50" = markers_po_e12_labels_top50, "po_e15_top50" = markers_po_e15_labels_top50, "po_p0_top50" = markers_po_p0_labels_top50
         )
       }
       else if (is.null(input$file1)==FALSE){
         sample_2 <- read.csv(input$file1$datapath)
       }
     }
   
   else if (input$label_option=="No_labels"){
     sample_1 <- switch(input$s1,
                        "Zhang_2015" = zhang_2015, "Lein_2007" = lein_2007, "Hochgerner_2018(early)" = hg_early_2018, 
                        "Hochgerner_2018(mature)" = hg_mature_2018, "Usoskin_2015" = usoskin_2015,
                        "Daneman_2010" = daneman_2010,  "Habib_2016(a)" = habib_2016_1, "Habib_2016(b)" = habib_2016_2, 
                        "Chen_2017" = chen_2017, "Gokce_2018" = gokce_2018, "Cahoy_2008" = cahoy_2008,
                        "Zeisal_2015" = zeisal_2015, "Haber_2017" = haber_2017, "Shekhar_2016" = shekhar_2016, 
                        "Zhong_2018" = zhong_2018, "Newman_2015" = newman_2015, "Zeisal_2018" = zeisal_2018, 
                        "Darmanis_2015" = darmanis_2018, "Macosko_2015" = macosko_2015,
                        #Our datasets without labels
                        "ct_e12" = markers_ct_e12, "ct_e15" = markers_ct_e15, "ct_p0" = markers_ct_p0,
                        "po_e12" = markers_po_e12, "po_e15" = markers_po_e15, "po_p0" = markers_po_p0,
                        "ct_e12_top50" = markers_ct_e12_labels_top50, "ct_e15_top50" = markers_ct_e15_labels_top50, 
                        "ct_p0_top50" = markers_ct_p0_labels_top50, "po_e12_top50" = markers_po_e12_labels_top50, 
                        "po_e15_top50" = markers_po_e15_labels_top50, "po_p0_top50" = markers_po_p0_labels_top50
       )
     if (is.null(input$file1)==TRUE){
      sample_2 <- switch(input$s2,
              "Zhang_2015" = zhang_2015, "Lein_2007" = lein_2007, "Hochgerner_2018(early)" = hg_early_2018, 
              "Hochgerner_2018(mature)" = hg_mature_2018, "Usoskin_2015" = usoskin_2015,
              "Daneman_2010" = daneman_2010,  "Habib_2016(a)" = habib_2016_1, "Habib_2016(b)" = habib_2016_2, 
              "Chen_2017" = chen_2017, "Gokce_2018" = gokce_2018, "Cahoy_2008" = cahoy_2008,
              "Zeisal_2015" = zeisal_2015, "Haber_2017" = haber_2017, "Shekhar_2016" = shekhar_2016, 
              "Zhong_2018" = zhong_2018, "Newman_2015" = newman_2015, "Zeisal_2018" = zeisal_2018, 
              "Darmanis_2015" = darmanis_2018, "Macosko_2015" = macosko_2015,
              #Our datasets without labels
              "ct_e12" = markers_ct_e12, "ct_e15" = markers_ct_e15, "ct_p0" = markers_ct_p0, "po_e12" = markers_po_e12,
              "po_e15" = markers_po_e15, "po_p0" = markers_po_p0, 
              "ct_e12_top50" = markers_ct_e12_labels_top50, "ct_e15_top50" = markers_ct_e15_labels_top50, 
              "ct_p0_top50" = markers_ct_p0_labels_top50, "po_e12_top50" = markers_po_e12_labels_top50, 
              "po_e15_top50" = markers_po_e15_labels_top50, "po_p0_top50" = markers_po_p0_labels_top50
       )
     }
     else if (is.null(input$file1)==FALSE){
       sample_2 <- read.csv(input$file1$datapath)
     }
   }
     

     p <- cluster_barplot(sample_1 = sample_1, sample_2 = sample_2, 
                     denom = input$denominator, clust_select =input$cluster_name) 
     show(p)
     #decide the inputselect option for clust_select
   }) 
   
   ##display the genes
   output$text1 <- renderText({
     if (input$label_option=="Labels"){
       sample_1 <- switch(input$s1,
                          "Zhang_2015" = zhang_2015, "Lein_2007" = lein_2007, "Hochgerner_2018(early)" = hg_early_2018, 
                          "Hochgerner_2018(mature)" = hg_mature_2018, "Usoskin_2015" = usoskin_2015, "Daneman_2010" = daneman_2010,
                          "Habib_2016(a)" = habib_2016_1, "Habib_2016(b)" = habib_2016_2, "Chen_2017" = chen_2017, 
                          "Gokce_2018" = gokce_2018, "Cahoy_2008" = cahoy_2008, "Zeisal_2015" = zeisal_2015, 
                          "Haber_2017" = haber_2017, "Shekhar_2016" = shekhar_2016, "Zhong_2018" = zhong_2018, 
                          "Newman_2015" = newman_2015, "Zeisal_2018" = zeisal_2018, "Darmanis_2015" = darmanis_2018, 
                          "Macosko_2015" = macosko_2015,
                          #Our datasets with Selin's labels
                          "ct_e12" = markers_ct_e12_labels, "ct_e15" = markers_ct_e15_labels, "ct_p0" = markers_ct_p0_labels,
                          "po_e12" = markers_po_e12_labels, "po_e15" = markers_po_e15_labels, "po_p0" = markers_po_p0_labels,
                          "ct_e12_top50" = markers_ct_e12_labels_top50, "ct_e15_top50" = markers_ct_e15_labels_top50, 
                          "ct_p0_top50" = markers_ct_p0_labels_top50, "po_e12_top50" = markers_po_e12_labels_top50, 
                          "po_e15_top50" = markers_po_e15_labels_top50, "po_p0_top50" = markers_po_p0_labels_top50
       )
       if (is.null(input$file1)==TRUE){
        sample_2 <- switch(input$s2,
                          "Zhang_2015" = zhang_2015, "Lein_2007" = lein_2007, "Hochgerner_2018(early)" = hg_early_2018, 
                          "Hochgerner_2018(mature)" = hg_mature_2018, "Usoskin_2015" = usoskin_2015, "Daneman_2010" = daneman_2010,
                          "Habib_2016(a)" = habib_2016_1, "Habib_2016(b)" = habib_2016_2, "Chen_2017" = chen_2017, 
                          "Gokce_2018" = gokce_2018, "Cahoy_2008" = cahoy_2008, "Zeisal_2015" = zeisal_2015, 
                          "Haber_2017" = haber_2017, "Shekhar_2016" = shekhar_2016, "Zhong_2018" = zhong_2018, 
                          "Newman_2015" = newman_2015, "Zeisal_2018" = zeisal_2018, "Darmanis_2015" = darmanis_2018, 
                          "Macosko_2015" = macosko_2015,
                          #Our datasets with Selin's labels
                          "ct_e12" = markers_ct_e12_labels, "ct_e15" = markers_ct_e15_labels, "ct_p0" = markers_ct_p0_labels,
                          "po_e12" = markers_po_e12_labels, "po_e15" = markers_po_e15_labels, "po_p0" = markers_po_p0_labels,
                          "ct_e12_top50" = markers_ct_e12_labels_top50, "ct_e15_top50" = markers_ct_e15_labels_top50, "ct_p0_top50" = markers_ct_p0_labels_top50, 
                          "po_e12_top50" = markers_po_e12_labels_top50, "po_e15_top50" = markers_po_e15_labels_top50, "po_p0_top50" = markers_po_p0_labels_top50
       )
       }
       else if (is.null(input$file1)==FALSE){
         sample_2 <- read.csv(input$file1$datapath)
       }
     }
     
     else if (input$label_option=="No_labels"){
       sample_1 <- switch(input$s1,
                          "Zhang_2015" = zhang_2015, "Lein_2007" = lein_2007, "Hochgerner_2018(early)" = hg_early_2018, 
                          "Hochgerner_2018(mature)" = hg_mature_2018, "Usoskin_2015" = usoskin_2015,
                          "Daneman_2010" = daneman_2010,  "Habib_2016(a)" = habib_2016_1, "Habib_2016(b)" = habib_2016_2, 
                          "Chen_2017" = chen_2017, "Gokce_2018" = gokce_2018, "Cahoy_2008" = cahoy_2008,
                          "Zeisal_2015" = zeisal_2015, "Haber_2017" = haber_2017, "Shekhar_2016" = shekhar_2016, 
                          "Zhong_2018" = zhong_2018, "Newman_2015" = newman_2015, "Zeisal_2018" = zeisal_2018, 
                          "Darmanis_2015" = darmanis_2018, "Macosko_2015" = macosko_2015,
                          #Our datasets without labels
                          "ct_e12" = markers_ct_e12, "ct_e15" = markers_ct_e15, "ct_p0" = markers_ct_p0,
                          "po_e12" = markers_po_e12, "po_e15" = markers_po_e15, "po_p0" = markers_po_p0,
                          "ct_e12_top50" = markers_ct_e12_labels_top50, "ct_e15_top50" = markers_ct_e15_labels_top50, 
                          "ct_p0_top50" = markers_ct_p0_labels_top50, "po_e12_top50" = markers_po_e12_labels_top50, 
                          "po_e15_top50" = markers_po_e15_labels_top50, "po_p0_top50" = markers_po_p0_labels_top50
       )
       if (is.null(input$file1)==TRUE){
        sample_2 <- switch(input$s2,
                          "Zhang_2015" = zhang_2015, "Lein_2007" = lein_2007, "Hochgerner_2018(early)" = hg_early_2018, 
                          "Hochgerner_2018(mature)" = hg_mature_2018, "Usoskin_2015" = usoskin_2015,
                          "Daneman_2010" = daneman_2010,  "Habib_2016(a)" = habib_2016_1, "Habib_2016(b)" = habib_2016_2, 
                          "Chen_2017" = chen_2017, "Gokce_2018" = gokce_2018, "Cahoy_2008" = cahoy_2008,
                          "Zeisal_2015" = zeisal_2015, "Haber_2017" = haber_2017, "Shekhar_2016" = shekhar_2016, 
                          "Zhong_2018" = zhong_2018, "Newman_2015" = newman_2015, "Zeisal_2018" = zeisal_2018, 
                          "Darmanis_2015" = darmanis_2018, "Macosko_2015" = macosko_2015,
                          #Our datasets without labels
                          "ct_e12" = markers_ct_e12, "ct_e15" = markers_ct_e15, "ct_p0" = markers_ct_p0, "po_e12" = markers_po_e12,
                          "po_e15" = markers_po_e15, "po_p0" = markers_po_p0, 
                          "ct_e12_top50" = markers_ct_e12_labels_top50, "ct_e15_top50" = markers_ct_e15_labels_top50, 
                          "ct_p0_top50" = markers_ct_p0_labels_top50, "po_e12_top50" = markers_po_e12_labels_top50, 
                          "po_e15_top50" = markers_po_e15_labels_top50, "po_p0_top50" = markers_po_p0_labels_top50
       )
       }
       else if (is.null(input$file1)==FALSE){
         sample_2 <- read.csv(input$file1$datapath)
       }
     }
     
     cluster_options_c2 <- switch(input$denominator,
                                  "Sample_1" = names(split(sample_2, f = sample_2$cluster)),
                                  "Sample_2" = names(split(sample_1, f = sample_1$cluster)))
     
     
     shared_genes <- marker_return(sample_1 = sample_1, sample_2 = sample_2, 
                                   cluster_s1 = input$cluster_name, cluster_s2 = input$cluster_2_name,
                                   genelist_type = input$sym_ens)
     paste0(shared_genes[1], "\n")
     #paste0(shared_genes[2])
     #paste0(shared_genes[3])
     
   })
      
   
   
   ##Write out all of the reference information
   output$refs1 <- renderText({
     ## Attach reference labels to the two samples
       sample_1 <- switch(input$s1,
                          "Zhang_2015" = "https://doi.org/10.1523/JNEUROSCI.4506-14.2015",
                          "Lein_2007" = "https://doi.org/10.1038/nature05453", 
                          "Hochgerner_2018(early)" = "https://doi.org/10.1038/s41593-017-0056-2", 
                          "Hochgerner_2018(mature)" = "https://doi.org/10.1038/s41593-017-0056-2", 
                          "Usoskin_2015" = "https://doi.org/10.1038/nn.3881", 
                          "Daneman_2010" = "https://doi.org/10.1371/journal.pone.0013741",
                          "Habib_2016(a)" = "https://doi.org/10.1126/science.aad7038", 
                          "Habib_2016(b)" = "https://doi.org/10.1126/science.aad7038", 
                          "Chen_2017" = "https://doi.org/10.1016/j.celrep.2017.03.004", 
                          "Gokce_2018" = "https://doi.org/10.1016/j.celrep.2016.06.059", 
                          "Cahoy_2008" = "https://doi.org/10.1523/JNEUROSCI.4178-07.2008", 
                          "Zeisal_2015" = "https://doi.org/10.1126/science.aaa1934", 
                          "Haber_2017" = "https://doi.org/10.1038/nature24489", 
                          "Shekhar_2016" = "https://doi.org/10.1016/j.cell.2016.07.054.", 
                          "Zhong_2018" = "https://doi.org/10.1038/nature24489", 
                          "Newman_2015" = "https://doi.org/10.1038/nature25980", 
                          "Zeisal_2018" = "https://doi.org/10.1101/294918", 
                          "Darmanis_2015" = "https://doi.org/10.1073/pnas.1507125112", 
                          "Macosko_2015" = "https://doi.org/10.1016/j.cell.2015.05.002",
                          #Our datasets with Selin's labels
                          "ct_e12" = "Hydra/sjessa/single_cell/2018-02_mouse_dev/objects/all_markers.RData", 
                          "ct_e15" = "Hydra/sjessa/single_cell/2018-02_mouse_dev/objects/all_markers.RData", 
                          "ct_p0" = "Hydra/sjessa/single_cell/2018-02_mouse_dev/objects/all_markers.RData",
                          "po_e12" = "Hydra/sjessa/single_cell/2018-02_mouse_dev/objects/all_markers.RData", 
                          "po_e15" = "Hydra/sjessa/single_cell/2018-02_mouse_dev/objects/all_markers.RData", 
                          "po_p0" = "Hydra/sjessa/single_cell/2018-02_mouse_dev/objects/all_markers.RData",
                          "ct_e12_top50" = "Hydra/sjessa/single_cell/2018-02_mouse_dev/objects/all_markers.RData",
                          "ct_e15_top50" = "Hydra/sjessa/single_cell/2018-02_mouse_dev/objects/all_markers.RData", 
                          "ct_p0_top50" = "Hydra/sjessa/single_cell/2018-02_mouse_dev/objects/all_markers.RData",
                          "po_e12_top50" = "Hydra/sjessa/single_cell/2018-02_mouse_dev/objects/all_markers.RData", 
                          "po_e15_top50" = "Hydra/sjessa/single_cell/2018-02_mouse_dev/objects/all_markers.RData", 
                          "po_p0_top50" = "Hydra/sjessa/single_cell/2018-02_mouse_dev/objects/all_markers.RData"
       )
       
     
     paste0( "Sample_1 source: ", sample_1)
   })   
   output$refs2 <- renderText({
     if (is.null(input$file1)==TRUE){
      sample_2 <- switch(input$s2,
                        "Zhang_2015" = "https://doi.org/10.1523/JNEUROSCI.4506-14.2015",
                        "Lein_2007" = "https://doi.org/10.1038/nature05453", 
                        "Hochgerner_2018(early)" = "https://doi.org/10.1038/s41593-017-0056-2", 
                        "Hochgerner_2018(mature)" = "https://doi.org/10.1038/s41593-017-0056-2", 
                        "Usoskin_2015" = "https://doi.org/10.1038/nn.3881", 
                        "Daneman_2010" = "https://doi.org/10.1371/journal.pone.0013741",
                        "Habib_2016(a)" = "https://doi.org/10.1126/science.aad7038", 
                        "Habib_2016(b)" = "https://doi.org/10.1126/science.aad7038", 
                        "Chen_2017" = "https://doi.org/10.1016/j.celrep.2017.03.004", 
                        "Gokce_2018" = "https://doi.org/10.1016/j.celrep.2016.06.059", 
                        "Cahoy_2008" = "https://doi.org/10.1523/JNEUROSCI.4178-07.2008", 
                        "Zeisal_2015" = "https://doi.org/10.1126/science.aaa1934", 
                        "Haber_2017" = "https://doi.org/10.1038/nature24489", 
                        "Shekhar_2016" = "https://doi.org/10.1016/j.cell.2016.07.054.", 
                        "Zhong_2018" = "https://doi.org/10.1038/nature24489", 
                        "Newman_2015" = "https://doi.org/10.1038/nature25980", 
                        "Zeisal_2018" = "https://doi.org/10.1101/294918", 
                        "Darmanis_2015" = "https://doi.org/10.1073/pnas.1507125112", 
                        "Macosko_2015" = "https://doi.org/10.1016/j.cell.2015.05.002",
                        #Our datasets with Selin's labels
                        "ct_e12" = "Hydra/sjessa/single_cell/2018-02_mouse_dev/objects/all_markers.RData", 
                        "ct_e15" = "Hydra/sjessa/single_cell/2018-02_mouse_dev/objects/all_markers.RData", 
                        "ct_p0" = "Hydra/sjessa/single_cell/2018-02_mouse_dev/objects/all_markers.RData",
                        "po_e12" = "Hydra/sjessa/single_cell/2018-02_mouse_dev/objects/all_markers.RData", 
                        "po_e15" = "Hydra/sjessa/single_cell/2018-02_mouse_dev/objects/all_markers.RData", 
                        "po_p0" = "Hydra/sjessa/single_cell/2018-02_mouse_dev/objects/all_markers.RData",
                        "ct_e12_top50" = "Hydra/sjessa/single_cell/2018-02_mouse_dev/objects/all_markers.RData",
                        "ct_e15_top50" = "Hydra/sjessa/single_cell/2018-02_mouse_dev/objects/all_markers.RData", 
                        "ct_p0_top50" = "Hydra/sjessa/single_cell/2018-02_mouse_dev/objects/all_markers.RData",
                        "po_e12_top50" = "Hydra/sjessa/single_cell/2018-02_mouse_dev/objects/all_markers.RData", 
                        "po_e15_top50" = "Hydra/sjessa/single_cell/2018-02_mouse_dev/objects/all_markers.RData", 
                        "po_p0_top50" = "Hydra/sjessa/single_cell/2018-02_mouse_dev/objects/all_markers.RData"
     )
     }
     else if (is.null(input$file1)==FALSE){
       sample_2 <- "User defined"
     }
     paste0("Sample_2 source: ", sample_2)
   }) 
     

   
   
   }## end of server


# Run the application 
shinyApp(ui = ui, server = server)



# paper <- reactive({
#  switch(input$s1,
#                "Zhang_2015" = a("10.1523/JNEUROSCI.4506-14.2015", href = "https://dx.doi.org/10.1523%2FJNEUROSCI.4506-14.2015"),
#               "Lein_2007" = a("doi:10.1038/nature05453", href = "https://www.nature.com/articles/nature05453"))#, 
#"Hochgerner_2018(early)" = hg_early_2018, 
#"Hochgerner_2018(mature)" = hg_mature_2018, 
#"Usoskin_2015" = usoskin_2015,
#"Daneman_2010" = daneman_2010,
#"Habib_2016(a)" = habib_2016_1, 
#"Habib_2016(b)" = habib_2016_2, 
#"Chen_2017" = chen_2017, 
#"Gokce_2018" = gokce_2018, 
#"Cahoy_2008" = cahoy_2008,
#"Zeisal_2015" = zeisal_2015, 
#"Haber_2017" = haber_2017, 
#"Shekhar_2016" = shekhar_2016, 
#"Zhong_2018" = zhong_2018, 
#"Newman_2015" = newman_2015, 
#"Zeisal_2018" = zeisal_2015, 
#"Darmanis_2015" = darmanis_2018, 
#"Macosko_2015" = macosko_2015,)
# output$DOI <- renderUI({
#   taglist("URL link", paper)
#  })



