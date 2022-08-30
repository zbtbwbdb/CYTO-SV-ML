#
#    # Main panel for displaying outputs ----
ui<- dashboardPage(

     ################################################################################################################### 
     
     #1. head information
     dashboardHeader(title = "CYTO-SV-ML pipeline for prediction of somatic cytogenetic results through genomic structural variant classification", titleWidth = 1200),
     
     ###################################################################################################################
     
     #2. set up input widget    
     dashboardSidebar(
     # action button for ploting    
                actionButton("Make_simulated_SV",'Run_CYTO_SV_ML'),
     # Input: Slider for SV type                    
           selectInput(inputId = "sv_type",
                       label = "SV type:",
                       choices = c("DEL", "DUP", "INV","TRS"),
                       selected = "DEL"),                     
     # Input: Slider for SV chromosome at 1st breakpoint
           selectInput(inputId = "sv_chr",
                       label = "SV chromosome:",
                       choices = paste('chr',c(1:22,'X','Y'),sep=''),
                       selected = "chr1"),        
     # Input: Slider for SV chromosome at 2nd breakpoint
           selectInput(inputId = "sv_chr2",
                       label = "SV 2nd chromosome (TRS only):",
                       choices = c(paste('chr',c(1:22,'X','Y'),sep=''), 'ns'),
                       selected = "chr1"),                     
     # Input: Slider for SV chromosome location                      
           selectInput(inputId = "cytoband",
                       label = "SV chrom arm:",
                       choices = c("p", "q","ns"),
                       selected = "ns"),                                                                            
     # Input: Slider for SV feature for axis X                     
           selectInput(inputId = "feature_x",
                       label = "Pairplot X-axis:",
                       choices = c("sv_chr", "sv_chr2", "sv_type","sv_bp_st_ci_range","sv_bp_end_ci_range","sv_read_ratio","sv_read_diff","sv_bp_st_cc1","sv_bp_end_cc1","prediction_TA","prediction_TG","prediction_TS"),
                       selected = "sv_read_ratio"),                  
     # Input: Slider for SV feature at axis Y                      
           selectInput(inputId = "feature_y",
                       label = "Pairplot Y-axis:",
                       choices = c("sv_chr", "sv_chr2", "sv_type","sv_bp_st_ci_range","sv_bp_end_ci_range","sv_read_ratio","sv_read_diff","sv_bp_st_cc1","sv_bp_end_cc1","prediction_TA","prediction_TG","prediction_TS"),
                       selected = "sv_bp_st_ci_range"),                  
     # Input: Slider for SV feature at axis Z                      
#           selectInput(inputId = "feature_z",
#                       label = "SV feature at axis Z:",
#                       choices = c("sv_chr", "sv_chr2", "sv_type","sv_bp_st_ci_range","sv_bp_end_ci_range","sv_read_ratio","sv_read_diff","sv_bp_st_cc1","sv_bp_end_cc1"),
#                       selected = "sv_type"),
     sidebarMenu(
          id = "tabs",
            menuItem("Abbreviation definition", icon = icon("info"), tabName = "Abbreviation definition"),           
            menuItem(print("bp: break point")),            
            menuItem(print("st: 1st breakpoint")),
            menuItem(print("ci: confidence interval")),   
            menuItem(print("cc: sequence complexity")),  
            menuItem(print("TA: True systematic artifact SV")),  
            menuItem(print("TG: True germline SV")), 
            menuItem(print("TS: True somatic cytogenetic SV")),    
            menuItem(print("UC: Unclassified SV"))                                                           
            )           
#           menuItem("Contact", icon = icon("phone"), tabName = "contact")
#           )
      ), 
                       
     ###################################################################################################################
     
     #3. layout for main page of plots and tables                       
     dashboardBody(
#            useShinyjs(),     
#            fluidPage(
#            fluidRow(id = "mainContent",
#                 column(12, h1("Main Content"))
#                 ),
#                 hidden(fluidRow(id = "contact", h1("Contact Info")))
#            ),
              #3.1 logo
              includeScript('www/addlogo.js'),
              includeCSS('www/custom.css'),
#              tags$style(type="text/css", "#timestamp {white-space: pre-wrap; word-break: break-word;}"),
              tags$style(type="text/css", "#timestamp {white-space: wrap; word-break: break-word;}", "#TxtOut {white-space: normal;}", "#value{ height: 200px; font-family: monospace;}"),
              #3.2 pop message
              #tags$head(tags$script(HTML('Shiny.addCustomMessageHandler("jsCode",function(message) {eval(message.value);});'))),

              #3.3 main layout for plots and tables
              tabsetPanel(
                         id = "tabselected",              
                        fluidRow(                 
                        column(width = 6, 
                              box(
                              title = "SV summary", width = NULL,height=280, solidHeader = TRUE,status = "primary",
                              wellPanel(htmlOutput("text"),height=280),
                              style = "height:280px; overflow-y: scroll;overflow-x: scroll;"                        
                              ),  
                              box(
                              title = "SV feature pair-plot", width = NULL, solidHeader = TRUE,status = "warning",
                              plotlyOutput("scatter2", width = 600,height=325),    
                              style = "height:280px; overflow-y: scroll;overflow-x: scroll;"                                 
                              )),                               
                        column(width = 6,
                              box(
                              title = "SV prediction plot", width = NULL, solidHeader = TRUE,status = "primary",
                              plotlyOutput("scatter1", width = 620,height=620),
                              style = "height:280px; overflow-y: scroll;overflow-x: scroll;"                                 
                              ))
                              ),  
                        fluidRow(                                           
                        box(
                        title = "SV data summary", width = 12,    solidHeader = TRUE,  status = "primary",                          
                                                lineupOutput("lineup1"))
                              )
                        )  
))
