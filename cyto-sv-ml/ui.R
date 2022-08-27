#
#    # Main panel for displaying outputs ----
ui<- dashboardPage(

     ################################################################################################################### 
     
     #1. head information
     dashboardHeader(title = "The Demonstration of CYTO-SV-ML pipeline for cytogenetic somatic SV classification", titleWidth = 800),
     
     ###################################################################################################################
     
     #2. set up input widget    
     dashboardSidebar(
     # action button for ploting    
                actionButton("Make_simulated_SV",'Run_simulation'),
     # Input: Slider for SV type                    
           selectInput(inputId = "sv_type",
                       label = "SV type:",
                       choices = c("DEL", "DUP", "INV","BND"),
                       selected = "DEL"),                     
     # Input: Slider for SV chromosome at 1st breakpoint
           selectInput(inputId = "sv_chr",
                       label = "SV chrom at 1st breakpoint:",
                       choices = paste('chr',c(1:22,'X','Y'),sep=''),
                       selected = "chr1"),        
     # Input: Slider for SV chromosome at 2nd breakpoint
           selectInput(inputId = "sv_chr2",
                       label = "SV chrom at 2nd breakpoint:",
                       choices = paste('chr',c(1:22,'X','Y'),sep=''),
                       selected = "chr1"),                     
     # Input: Slider for SV chromosome location                      
           selectInput(inputId = "cytoband",
                       label = "SV chrom arm:",
                       choices = c("p", "q","ns"),
                       selected = "ns"),                                                                            
     # Input: Slider for SV feature for axis X                     
           selectInput(inputId = "feature_x",
                       label = "Pairplot X-axis:",
                       choices = c("sv_chr", "sv_chr2", "sv_type","sv_bp_st_ci_range","sv_bp_end_ci_range","sv_read_ratio","sv_read_diff","sv_bp_st_cc1","sv_bp_end_cc1"),
                       selected = "sv_read_ratio"),                  
     # Input: Slider for SV feature at axis Y                      
           selectInput(inputId = "feature_y",
                       label = "Pairplot Y-axis:",
                       choices = c("sv_chr", "sv_chr2", "sv_type","sv_bp_st_ci_range","sv_bp_end_ci_range","sv_read_ratio","sv_read_diff","sv_bp_st_cc1","sv_bp_end_cc1"),
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
              tags$head(tags$script(HTML('Shiny.addCustomMessageHandler("jsCode",function(message) {eval(message.value);});'))),

              #3.3 main layout for plots and tables
              tabsetPanel(
                         id = "tabselected",              
                        fluidRow(                 
                        column(width = 6, 
                              box(
                              title = "SV class prediction", width = NULL, solidHeader = TRUE,status = "primary",
                              wellPanel(htmlOutput("text"))                         
                              ),  
                              box(
                              title = "SV class prediction", width = NULL, solidHeader = TRUE,status = "warning",
                              plotlyOutput("scatter2", width = 600,height=315)                              
                              )),                               
                        column(width = 6,
                              box(
                              title = "SV class prediction", width = NULL, solidHeader = TRUE,status = "primary",
                              plotlyOutput("scatter1", width = 600,height=595)
                              ))
                              ),  
                        fluidRow(                                           
                        box(
                        title = "SV feature summary", width = 12,    solidHeader = TRUE,  status = "primary",                          
                                                lineupOutput("lineup1"))
                              )
                        )  
))
