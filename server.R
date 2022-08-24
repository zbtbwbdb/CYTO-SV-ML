server<-function(input, output) 
    {
    ################################################################################################################### 
    
    # 1. load environment with the data below:
    # label_hash: dictionary for transformation of clinical variable code ( from App widget input or encode3.merge.txt
    # all_donor_db: data frame with patient ID and age information for all the donor available in our database
    # optimaldonor: data frame with optimal donor of BART3 model
    # recipient_db: data frame from encode3.merge.txt simplified with only donors of real transplant (one donor for each recipient)
    # post: pre-compiled BART3 model
     
#     data("EFS_data3")
     load("EFS_data3.RData")
     removeModal()
     
    ################################################################################################################### 
    
    #2. assign the input values from App widget
    checkagegp <- eventReactive(input$makePlot, {input$agegp})
    checksex <- eventReactive(input$makePlot, {input$sex})
    checkEthnicity <- eventReactive(input$makePlot, {input$Ethnicity})
    checkkpslt90 <- eventReactive(input$makePlot, {input$kpslt90})
    checkdisease <- eventReactive(input$makePlot, {input$disease})
    checkdisstat <- eventReactive(input$makePlot, {input$disstat})
    checkintdxtxgp <- eventReactive(input$makePlot, {input$intdxtxgp})
    checkcondint <- eventReactive(input$makePlot, {input$condint})
    checktbi <- eventReactive(input$makePlot, {input$tbi})
    checkgraftype <- eventReactive(input$makePlot, {input$graftype})
    checkdnrage <- eventReactive(input$makePlot, {input$dnrage})
    checktac <- eventReactive(input$makePlot, {input$tac})
    checkcsa <- eventReactive(input$makePlot, {input$csa})      
    checkgvhdo <- eventReactive(input$makePlot, {input$gvhdo})
    checkinvivo_tcd <- eventReactive(input$makePlot, {input$invivo_tcd})                    
    checkyeargp <- eventReactive(input$makePlot, {input$yeargp})
    checkrcmv <- eventReactive(input$makePlot, {input$rcmv})  
    
    ###################################################################################################################   
    
    #3. wrap prediction function for post: pre-compiled BART3 model
      bart3_pred<-function(x)
           {         
           #3.1 predict the clustered subset with input sample 
           x_pred = predict(post, x, mc.cores=1)
           #3.2 re-ogranize prediction data into data frame    
           x_p.mean<-as.data.frame(x_pred$prob.test.mean,Stringasfactor=F)
           x_p.lower <- as.data.frame(x_pred$prob.test.lower,Stringasfactor=F)
           x_p.upper <- as.data.frame(x_pred$prob.test.upper,Stringasfactor=F)   
           input_p<-rbind(t(x_p.mean),t(x_p.lower),t(x_p.upper))
           #3.3 return data frame
           return (input_p)
           }
      
    ###################################################################################################################  
    #4. generate the input data subset for recipient list table 
    input_subset_tmp<- reactive({
          #4.1 load assigned input values
          agegp<-checkagegp()
          sex<-checksex()
          Ethnicity<-checkEthnicity()
          kpslt90<-checkkpslt90() 
          disease<-checkdisease() 
          disstat<-checkdisstat()
          intdxtxgp<-checkintdxtxgp() 
          condint<-checkcondint()
          tbi<-checktbi()
          graftype<-checkgraftype() 
          dnrage<-checkdnrage()
          tac<-checktac() 
          csa<-checkcsa()      
          gvhdo<-checkgvhdo()
          invivo_tcd<-checkinvivo_tcd()                  
          yeargp<-checkyeargp()
          rcmv<-checkrcmv() 
          
          #4.2 transform input values into data frame with the same format as optimaldonor (ready as BART3 model input)
          input_dp<-as.data.frame(t(c(label_hash[[agegp]],as.numeric(sapply(label_hash[[Ethnicity]], function(x) str_split(x, ",")[[1]],USE.NAMES=FALSE)),label_hash[[sex]],as.numeric(sapply(label_hash[[kpslt90]], function(x) str_split(x, ",")[[1]],USE.NAMES=FALSE)),as.numeric(sapply(label_hash[[disease]], function(x) str_split(x, ",")[[1]],USE.NAMES=FALSE)),label_hash[[disstat]],as.numeric(sapply(label_hash[[intdxtxgp]], function(x) str_split(x, ",")[[1]],USE.NAMES=FALSE)),label_hash[[condint]],as.numeric(sapply(label_hash[paste("tbi_",tbi,sep="")], function(x) str_split(x, ",")[[1]],USE.NAMES=FALSE)),label_hash[[tac]],label_hash[[csa]],label_hash[[gvhdo]],label_hash[[invivo_tcd]],label_hash[[graftype]],label_hash[[yeargp]],as.numeric(sapply(label_hash[paste("rcmv_",rcmv,sep="")], function(x) str_split(x, ",")[[1]],USE.NAMES=FALSE)),as.numeric(dnrage)/10)),Stringasfactor=F)
          colnames(input_dp)<-colnames(optimaldonor)
          
          #4.3 set up recipient database and selected variables
          tmp_db<-recipient_db
        	race<-colnames(input_dp)[which(input_dp[1,2:8]==1)+1]
        	disease<-colnames(input_dp)[which(input_dp[1,13:16]==1)+12]
        	disstat<-input_dp[1,17]
         
          #4.4 randomize 50 recipient from database based on disease type and status, and race 
          all_data_sel<-tmp_db[tmp_db[,11]==label_hash[[race]] & tmp_db[,8]==label_hash[[disease]] & tmp_db[,9] %in% sapply(label_hash[[as.character(disstat)]], function(x) str_split(x, ",")[[1]],USE.NAMES=FALSE),]
          # additional resampling when not enough samples in the database
          if (dim(all_data_sel)[1]<50)
                { 
                n<-50-dim(all_data_sel)[1]       
                all_data_tmp2<-tmp_db[tmp_db[,11]==1 & tmp_db[,8]==label_hash[[disease]] & tmp_db[,9] %in% sapply(label_hash[[as.character(disstat)]], function(x) str_split(x, ",")[[1]],USE.NAMES=FALSE),]
                all_data_tmp3<-all_data_tmp2[sample(seq(1:dim(all_data_tmp2)[1]),n),]
                all_data_tmp<-rbind(all_data_tmp,all_data_tmp3)  
                }  
          else
                {
                  all_data_tmp<-all_data_sel[sample(seq(1:dim(all_data_sel)[1]),50),]
                }
              
          #4.5 re-orgnize the data frame and return output
### The changes of line 102-105 for adding donor id column in the table          
          all_data_tmp<-all_data_tmp[,c(1,8,9,11,12,13,7,6,4)]
          colnames(all_data_tmp)<-c('crid_set', 'disease_set','disstat_set','racegp_set','hisp_set','sex_set','agegp_set','dnrage_set','nmdp_did')
          all_data_tmp[,10:42]<-input_dp[rep(1,50),]
          all_data_tmp[,42]<-all_data_tmp['dnrage_set']/10  
          all_data_tmp
          })
    
    ###################################################################################################################  
     
    #5. generate recipient list table with the input data subset from 4.
    output$disttable1 <- renderDataTable( 
          input_subset_tmp()[,1:4], extensions = 'Buttons', 
          options = list(dom = 'Bfrtip',
                         buttons = c('copy', 'csv', 'excel', 'pdf', 'print')),
          rownames= FALSE
          )
                                          
    ###################################################################################################################   
    
    #6. detect selected sample for visualiztion of BART3 prediction plot and pop up messege & generate table for significantly optimal donor available
    observe({
          #6.1 detect the selected row with recipient ID information
          req(input$disttable1_rows_selected)
          selRow <- input_subset_tmp()[input$disttable1_rows_selected,] 
          id<<-selRow[[1]]
          
          
          #6.2 pull out recipient clinical information and donor age information          
          r_var<<-input_subset_tmp()[input_subset_tmp()[,1]==id,]
### The changes of line 132-133, 138-140 for adding donor id information          
          r_var2<-as.numeric(all_donor_db[all_donor_db[,1]==id,2])/10
          r_var3<-all_donor_db[all_donor_db[,1]==id,2]        
          #6.3 re-organize the data frame of recipient & donors and re-assign the donor age information
          v_t<-length(r_var2)+1
          dfs0<-r_var[rep(seq_len(nrow(r_var)), each = v_t), ]
          dfs0[2:v_t,'dnrage']<-r_var2
          dfs0[2:v_t,'nmdp_did']<-r_var3
          dfs0_tmp<-dfs0[1:min(dim(dfs0)[1],50),10:42]
          for (i in 1:dim(dfs0_tmp)[2])
                {
                  dfs0_tmp[,i]<-as.numeric(dfs0_tmp[,i])
                }

          #6.4 BART3 prediction for data frame of selected recipient & donors
          dfs0_pred = bart3_pred(dfs0_tmp)
          
          #6.5. transform prediction data for plot
          output$distPlot1 <- renderPlotly({  
          
          #6.5.1 re-organize the BART3 prediction, label the transplant donor and order by risk probability
          dfs0_pred_tmp<-as.data.frame(t(dfs0_pred),stringsAsFactors=F)
          colnames(dfs0_pred_tmp)<-c('mean','lower','upper')
###  The changes of line 155 for index switch due to additional donor id column           
          dfs1<-cbind(dfs0[1:min(dim(dfs0)[1],50),10:42],dfs0_pred_tmp)
          dfs1$group=c('patient_with_transplant_donor',rep('patient_with_non-transplant_donor',min(dim(dfs0)[1],50)-1)) 
          dfs2<<-dfs1[order(dfs1$mean),]
          dfs2$casenum=1:nrow(dfs2)  
              
          #6.5.2 bart3 donor comparison ( the diff value will be checked later with the cutoff for pop up messege & generate table for significantly optimal donor available )
          pred_diff<<- dfs0_pred[1,1]-min(dfs0_pred[1,])

          #6.5.3 extract clinical variables for hover information of disease / dissat / dnrage on the plot    
          pred_mean<-dfs2[,'mean']
          pred_upper<-dfs2[,'upper']
          pred_lower<-dfs2[,'lower']
          disease<-colnames(which(dfs2[1,13:16]==1)+12)
          disstat<-dfs2[,'disstat']
          race<-colnames(which(dfs2[1,2:8]==1)+1)
          dnrage<-dfs2[,'dnrage']

          ###################################################################################################################            
          #6.5.4 Plot posterior mean of selected recipient & donors from BART3 model    
###   The changes of line 175 for correcting the bug of confident interval by mislabelling dfs2 as dfs1 for ribbons data source          
          p2<- dfs2 %>% plot_ly(type = 'scatter',mode = 'markers',x=~dfs2$casenum, y=~mean,color=~factor(group),colors = c( "blue", "red"),hovertemplate = paste("prediction_mean:",pred_mean,"<br>upper_CI:",pred_upper,"<br>lower_CI:",pred_lower,"<br>recipient disease:",disease,"<br>recipient disease stage:",disstat,"<br>recipient race:",race,"<br>donor age:", dnrage)) %>% add_ribbons(data=dfs2, ymin =~lower , ymax = ~upper, line = list(color = 'rgba(7, 164, 181, 0.05)'), fillcolor = 'rgba(7, 164, 181, 0.2)', name = '95% CI')
          p2 <- p2  %>% layout(yaxis = list(title = "One year death/GVHD risk probability", range = c(0,1)))  %>% layout(xaxis = list(title = "Sample No.", range = c(0,51)))
          p2      })
          ###################################################################################################################     
           
          #6.5.5 pop up messege & generate table for significantly optimal donor available
          if(as.numeric(pred_diff) >= 0.05)
                  {
###    The changes of line 185-189 for the pop up messege of detecting significant better donor            
                  # pop up messege 
                  showModal(modalDialog(
                            title = "optimal donor detected",
                            "please double check whether any donor is significantly better than the donor selected for transplant",
                            easyClose = TRUE,
                            footer = NULL))
                  # generate table
###    The changes of line 192 for selectively showing donor information in the table                   
                  output$disttable2 <- renderDataTable(dfs2[,33:38], extensions = 'Buttons', 
                                                         options = list(dom = 'Bfrtip',
                                                                        buttons = c('copy', 'csv', 'excel', 'pdf', 'print'), scrollX = TRUE),
                                                         rownames= FALSE                                              
                                                         )   
                  }
                  else
                  {     
###    The changes of line 201 for selectively showing donor information in the table                       
                    output$disttable2 <- renderDataTable(dfs2[,33:38], extensions = 'Buttons', 
                                                         options = list(dom = 'Bfrtip',
                                                                        buttons = c('copy', 'csv', 'excel', 'pdf', 'print'), scrollX = TRUE),
                                                         rownames= FALSE                                              
                    )            
                  }  
###    The changes of line 208 for cleanup of public variable pred_diff ( still need fix the initiation of this variable )             
          rm(pred_diff)
    })       
}

