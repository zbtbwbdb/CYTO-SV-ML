#
#    # Main panel for displaying outputs ----
ui<- dashboardPage(

     ################################################################################################################### 
     
     #1. head information
     dashboardHeader(title = "CYTO-SV-ML pipeline for cytogenetic somatic SV classification", titleWidth = 600),
     
     ###################################################################################################################
     
     #2. set up input widget    
     dashboardSidebar(
     # action button for ploting    
                actionButton("Make_simulated_SV",'Plot_SV_data'),
     # Input: Slider for SV type                    
           selectInput(inputId = "sv_type",
                       label = "SV type:",
                       choices = c("DEL", "DUP", "INV","BND"),
                       selected = "DEL"),                     
     # Input: Slider for SV chromosome at 1st breakpoint
           selectInput(inputId = "sv_chr",
                       label = "SV chromosome at 1st breakpoint:",
                       choices = c(1:22,'X',"Y"),
                       selected = "1"),        
     # Input: Slider for SV chromosome at 2nd breakpoint
           selectInput(inputId = "sv_chr2",
                       label = "SV chromosome at 2nd breakpoint:",
                       choices = c(1:22,'X',"Y"),
                       selected = "1"),                     
     # Input: Slider for SV chromosome location                      
           selectInput(inputId = "cytoband",
                       label = "SV chromosome location:",
                       choices = c("p", "q","ns"),
                       selected = "ns"),                                                                            
     # Input: Slider for SV feature for axis X                     
           selectInput(inputId = "feature_x",
                       label = "SV feature at axis X:",
                       choices = c("sv_chr", "sv_chr2", "sv_type","sv_st_ci_range","sv_end_ci_range","sv_read_ratio","sv_read_diff"),
                       selected = "sv_chr"),                  
     # Input: Slider for SV feature at axis Y                      
           selectInput(inputId = "feature_y",
                       label = "SV feature at axis Y:",
                       choices = c("sv_chr", "sv_chr2", "sv_type","sv_st_ci_range","sv_end_ci_range","sv_read_ratio","sv_read_diff"),
                       selected = "sv_chr2"),                  
     # Input: Slider for SV feature at axis Z                      
           selectInput(inputId = "feature_z",
                       label = "SV feature at axis Z:",
                       choices = c("sv_chr", "sv_chr2", "sv_type","sv_st_ci_range","sv_end_ci_range","sv_read_ratio","sv_read_diff"),
                       selected = "sv_type"), 
                       
     ###################################################################################################################
     
     #3. layout for main page of plots and tables                       
     dashboardBody(
     
              #3.1 logo
              includeScript('www/addlogo.js'),
              includeCSS('www/custom.css'),
              tags$style(type="text/css", "#timestamp {white-space: pre-wrap;}"),
              
              #3.2 pop message
              tags$head(tags$script(HTML('Shiny.addCustomMessageHandler("jsCode",function(message) {eval(message.value);});'))),
              
              #3.3 main layout for plots and tables
              tabsetPanel(
                         id = "tabselected",
                         tabPanel(
                                  'Structural variant classification summary', value = 1,
                                  fluidRow(
                                    column(width = 4,DT::dataTableOutput("disttable1"))),
                                  fluidRow(
                                    column(width = 8,plotlyOutput("distPlot1", width = 800, height = 400)))
                                  )
                         )
                 )  
)