#
#    # Main panel for displaying outputs ----
ui<- dashboardPage(

     ################################################################################################################### 
     
     #1. head information
     dashboardHeader(title = "One Year Event Free Survival Prediction Calculator", titleWidth = 600),
     
     ###################################################################################################################
     
     #2. set up input widget    
     dashboardSidebar(
     # action button for ploting    
                actionButton("makePlot",'Plot Data'),
     # Input: Slider for Recipient age group
           selectInput(inputId = "agegp",
                       label = "Recipient age group:",
                       choices = c("18-29", "30-39","40-49", "50-59","60+"),
                       selected = "18-29"),        
     # Input: Slider for Recipient gender 
           selectInput(inputId = "sex",
                       label = "Recipient gender:",
                       choices = c("Male", "Female"),
                       selected = "Male"),                     
     # Input: Slider for race group                      
           selectInput(inputId = "Ethnicity",
                       label = "race group:",
                       choices = c("Caucasian", "Hispanic", "African-American","Asian/Pacific_Islander", "Native_American","Other_Race","Unknown_Race"),
                       selected = "Caucasian"),                                                       
     # Input: Slider for Karnofsky score                    
           selectInput(inputId = "kpslt90",
                       label = "Karnofsky score:",
                       choices = c("<90", "90+","unknown_score"),
                       selected = "<90"),                          
     # Input: Slider for disease type at HSCT                     
           selectInput(inputId = "disease",
                       label = "disease type at HSCT:",
                       choices = c("Acute_myeloid_leukemia", "Acute_lymphoblastic_leukemia", "Chronic_myeloid_leukemia","Myelodysplastic_syndrome"),
                       selected = "Acute_myeloid_leukemia"),                  
     # Input: Slider for Disease status prior to HSCT                      
           selectInput(inputId = "disstat",
                       label = "Disease status prior to HSCT:",
                       choices = c("Early", "Intermediate", "Advanced"),
                       selected = "Early"),                  
     # Input: Slider for Time from diagnosis to HSCT                      
           selectInput(inputId = "intdxtxgp",
                       label = "Time from diagnosis to HSCT:",
                       choices = c("<6_months", "6-12_months", "12+_months", "Missing_date"),
                       selected = "<6_months"), 
     # Input: Slider for Conditioning intensity                      
           selectInput(inputId = "condint",
                       label = "Conditioning intensity:",
                       choices = c("MAC", "RIC"),
                       selected = "MAC"), 
     # Input: Slider for Total body irradiation                     
           selectInput(inputId = "tbi",
                       label = "Total body irradiation:",
                       choices = c("Yes", "No","Missing"),
                       selected = "Yes"), 
     # Input: Slider for Graft type                     
           selectInput(inputId = "graftype",
                       label = "Graft type:",
                       choices = c("Bone_Marrow", "Peripheral_Blood"),
                       selected = "Bone_Marrow"),                   
     # Input: Slider for Donor age in decades                      
           textInput(inputId ="dnrage", label = "Donor age in decades:", value = 20),                   
     # Input: Slider for Tacrolimus (FK506) based GVHD prophylaxis                    
           selectInput(inputId = "tac",
                       label = "Tacrolimus (FK506) based GVHD prophylaxis:",
                       choices = c("Yes", "No"),
                       selected = "Yes"), 
     # Input: Slider for Cyclosporine A (CSA) based GVHD prophylaxis                     
           selectInput(inputId = "csa",
                       label = "Cyclosporine A (CSA) based GVHD prophylaxis:",
                       choices = c("Yes", "No"),
                       selected = "No"), 
     # Input: Slider for Non-FK506/Non-CSA based GVHD prophylaxis                     
           selectInput(inputId = "gvhdo",
                       label = "Non-FK506/Non-CSA based GVHD prophylaxis:",
                       choices = c("Yes", "No"),
                       selected = "No"), 
     # Input: Slider for In vivo T-cell depletio                     
           selectInput(inputId = "invivo_tcd",
                       label = "In vivo T-cell depletion:",
                       choices = c("Yes", "No"),
                       selected = "Yes"), 
     # Input: Slider for plot subgroup labelYear of HSCTputId = "yeargp",
           selectInput(inputId = "yeargp",
                       label = "Year of HSCT:",
                       choices = c("1999-2002", "2003-2006","2007-2011","2012-2014","2015-2016"),
                       selected = "2012-2014"), 
     # Input: Slider for TRecipient Cytomegalovirus                      
           selectInput(inputId = "rcmv",
                       label = "TRecipient Cytomegalovirus:",
                       choices = c("Positive", "Negative","Unknown"),
                       selected = "Yes")), 
                       
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
                                  'Survival statistics summary of prediction', value = 1,
                                  fluidRow(
                                    column(width = 4,DT::dataTableOutput("disttable1")),  column(width = 8,plotlyOutput("distPlot1", width = 800, height = 400)),
                                  fluidRow(
                                    column(width = 12,DT::dataTableOutput("disttable2")))
                                          )
                                  )
                         )
                 )  
)