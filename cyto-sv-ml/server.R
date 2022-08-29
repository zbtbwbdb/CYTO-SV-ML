server<-function(input, output) 
    {
    ################################################################################################################### 
    
    # 1. load environment with the data below:
    # label_hash: dictionary for transformation of clinical variable code ( from App widget input or encode3.merge.txt
    # all_donor_db: data frame with patient ID and age information for all the donor available in our database
    # optimaldonor: data frame with optimal donor of BART3 model
    # recipient_db: data frame from encode3.merge.txt simplified with only donors of real transplant (one donor for each recipient)
    # post: pre-compiled BART3 model

     #data("data/cyto_sv_ml")
     load("data/cyto_sv_ml.RData")
     removeModal()
     
    ################################################################################################################### 
    
    #2. assign the input values from App widget
    check_sv_chr <- eventReactive(input$Make_simulated_SV, {input$sv_chr})
    check_sv_chr2 <- eventReactive(input$Make_simulated_SV, {input$sv_chr2})
    check_sv_type <- eventReactive(input$Make_simulated_SV, {input$sv_type})   
    check_cytoband <- eventReactive(input$Make_simulated_SV, {input$cytoband})
    check_feature_x <- eventReactive(input$Run_simulation, {input$feature_x})
    check_feature_y <- eventReactive(input$Run_simulation, {input$feature_y})
#    check_feature_z <- eventReactive(input$Run_simulation, {input$feature_z})
 
    ###################################################################################################################  
    #3. generate the input data subset for recipient list table 
    input_subset_tmp<- reactive({
          #3.1 load assigned input values
            sv_chr<-check_sv_chr()
            sv_chr2<-check_sv_chr2() 
#            sv_chr<-paste('chr',sv_chr,sep='')               
#            sv_chr2<-paste('chr',sv_chr2,sep='')           
            sv_type<-check_sv_type()
            sv_chr_arm<-check_cytoband()    
            #3.2 generate the simulated examples from random resampling
#            bp_min=genome_coor[genome_coor$chr==sv_chr & genome_coor$chr_arm==sv_chr_arm,3]
#            bp_max=genome_coor[genome_coor$chr==sv_chr & genome_coor$chr_arm==sv_chr_arm,4]           
            cytoband_cen1=genome_coor[genome_coor$chr==sv_chr,2]
            cytoband_cen2=genome_coor[genome_coor$chr==sv_chr,3]            
      	    cytoband_end=genome_coor[genome_coor$chr==sv_chr,4]
      	    if (sv_chr_arm=='p')
      	          {                         
                        cytoband_range=c(0,cytoband_cen1)
      	    } else if (sv_chr_arm=='q') {
                        cytoband_range=c(cytoband_cen2,cytoband_end)                                       
                  } else {
                        cytoband_range=c(0,cytoband_end)
                  }
      	    print(cytoband_range)  
#            if (sv_type=='BND')
#                  {
#                  benchmark_sv=read.table("data/cnv_benchmark_all.vcf",sep="\t",header=T)
#                  sv_simu_database=read.table("data/cnv_simu.vcf",sep="\t",header=T)
#                  input_sv=sv_simu_database[sv_simu_database$Chrom==sv_chr & sv_simu_database$chr2==sv_chr2 & sv_simu_database$bp_st %in% cytoband_range,]
#                  input_sv=sv_simu_database[sv_simu_database$Chrom==sv_chr & sv_simu_database$chr2==sv_chr2,]
#            } else {
#                  benchmark_sv=read.table("data/trs_benchmark_all.vcf",sep="\t",header=T)
#                  sv_simu_database=read.table("data/trs_simu.vcf",sep="\t",header=T)
#                  input_sv=sv_simu_database[sv_simu_database$sv_chr==sv_chr & sv_simu_database$sv_type==sv_type & sv_simu_database$bp_st %in% cytoband_range,]
#                  input_sv=sv_simu_database[sv_simu_database$sv_chr==sv_chr & sv_simu_database$sv_type==sv_type,]
#                  }
#            input_sv=input_sv[sample(1:dim(input_sv)[1],1),]      
#            write.table(input_sv,"data/input_sv.vcf",sep="\t",col.names=T,row.names=F)
#            #3.3 model prediction by system call python script   
#            pred_sv=system2("python","data/cyto-sv-ml.py","input_sv.vcf",sv_type, stdout = TRUE, stderr = TRUE)
#            input_sv['prediction']=int(pred_sv)
            if ( sv_type=='BND')
                  {
              all_input_sv=all_sv[all_sv$label=='UC' & all_sv$sv_chr==sv_chr & all_sv$sv_chr2==sv_chr2 & all_sv$sv_type==sv_type & ((all_sv$sv_bp_st>=cytoband_range[1] & all_sv$sv_bp_st<=cytoband_range[2]) | (all_sv$sv_bp_end>=cytoband_range[1] & all_sv$sv_bp_end<=cytoband_range[2])),]
            } else {
              all_input_sv=all_sv[all_sv$label=='UC' & all_sv$sv_chr==sv_chr & all_sv$sv_type==sv_type & ((all_sv$sv_bp_st>=cytoband_range[1] & all_sv$sv_bp_st<=cytoband_range[2]) | (all_sv$sv_bp_end>=cytoband_range[1] & all_sv$sv_bp_end<=cytoband_range[2])),]
                  }
            print(dim(all_input_sv))  
            if(dim(all_input_sv)[1]==0)
                  {
                    ###    the pop up messege of detecting none sv            
                    # pop up messege 
                    showModal(modalDialog(
                    title = "warning",
                    "No such SVs in current WGS database found and will change Chrom arm or Chrom2 to non-specific ",
                    easyClose = TRUE,
                    footer = NULL))
                    all_input_sv1=all_sv[all_sv$label=='UC' & all_sv$sv_chr==sv_chr & all_sv$sv_chr2==sv_chr2 & all_sv$sv_type==sv_type,]
                    all_input_sv=all_input_sv1
                    if(dim(all_input_sv1)[1]==0) {
                          all_input_sv1=all_sv[all_sv$label=='UC' & all_sv$sv_chr==sv_chr & all_sv$sv_type==sv_type & ((all_sv$sv_bp_st>=cytoband_range[1] & all_sv$sv_bp_st<=cytoband_range[2]) | (all_sv$sv_bp_end>=cytoband_range[1] & all_sv$sv_bp_end<=cytoband_range[2])),]
                          all_input_sv=all_input_sv1
                          if(dim(all_input_sv1)[1]==0) {
                                all_input_sv1=all_sv[all_sv$label=='UC' & all_sv$sv_chr==sv_chr & all_sv$sv_type==sv_type,]  
                                all_input_sv=all_input_sv1
                                }
                          }
                  }
            print(dim(all_input_sv)) 
            input_sv=all_input_sv[sample(1:dim(all_input_sv)[1],1),]          
            all_benchmark_sv=all_sv[all_sv$label %in% c('TA','TG','TS'), ]
            print(dim(all_benchmark_sv)) 
            benchmark_sv=all_benchmark_sv[all_benchmark_sv$sv_chr==sv_chr & all_benchmark_sv$sv_chr2==sv_chr2 & all_benchmark_sv$sv_type==sv_type & ((all_benchmark_sv$sv_bp_st>=cytoband_range[1] & all_benchmark_sv$sv_bp_st<=cytoband_range[2]) | (all_benchmark_sv$sv_bp_end>=cytoband_range[1] & all_benchmark_sv$sv_bp_end<=cytoband_range[2])),]          
            print(dim(benchmark_sv)) 
            if (dim(benchmark_sv)[1]>100)
                  {
                  benchmark_sv=benchmark_sv[sample(1:dim(benchmark_sv)[1],100),]       
            }  else {
                  benchmark_sv0=all_benchmark_sv[all_benchmark_sv$sv_chr==sv_chr & all_benchmark_sv$sv_type==sv_type & ((all_benchmark_sv$sv_bp_st>=cytoband_range[1] & all_benchmark_sv$sv_bp_st<=cytoband_range[2]) | (all_benchmark_sv$sv_bp_end>=cytoband_range[1] & all_benchmark_sv$sv_bp_end<=cytoband_range[2])),]          
                  if (dim(benchmark_sv0)[1]<=100)
                        {                  
                        benchmark_sv0=all_benchmark_sv[all_benchmark_sv$sv_chr==sv_chr & ((all_benchmark_sv$sv_bp_st>=cytoband_range[1] & all_benchmark_sv$sv_bp_st<=cytoband_range[2]) | (all_benchmark_sv$sv_bp_end>=cytoband_range[1] & all_benchmark_sv$sv_bp_end<=cytoband_range[2])),]                           
                        if (dim(benchmark_sv)[1]>100)
                              {                        
                              benchmark_sv=rbind(benchmark_sv,benchmark_sv0[sample(1:dim(benchmark_sv0)[1],100),])   
                              }                        
                  } else {
                        benchmark_sv=rbind(benchmark_sv,benchmark_sv0[sample(1:dim(benchmark_sv0)[1],100),])     
                        }   
            } 
                        
            #3.5 return data frame
            return (list(input_sv,benchmark_sv,all_benchmark_sv))      
            })

    ###################################################################################################################  
     
    #4. generate summary text for input sv data
      output$text <- renderUI({
            input_sv=input_subset_tmp()[[1]]      
            input_sv_type=ifelse(input_sv$sv_type=='BND','Translocation', ifelse(input_sv$sv_type=='DEL', "Deletion", ifelse(input_sv$sv_type=='DUP', 'Duplication',ifelse(input_sv$sv_type=='INV', 'Inversion'))))
            input_sv_class=ifelse(input_sv$predict_label=='TA','Systematic artifact SV', ifelse(input_sv$predict_label=='TG', "True germline SV", ifelse(input_sv$predict_label=='TS', 'True somatic cytogenetic SV')))
            HTML(paste(paste0("This SV is <b>", input_sv_type, "</b> with 1st break point at chromsome <b>", input_sv$sv_chr, "</b> position <b>", input_sv$sv_bp_st, "</b> (CI=",input_sv$sv_bp_st_ci0,",",input_sv$sv_bp_st_ci1,") and 2nd break point at chromsome <b>", input_sv$sv_chr2, "</b> position <b>", input_sv$sv_bp_end,"</b> (CI=",input_sv$sv_bp_st_ci0,",",input_sv$sv_bp_st_ci1,"), with the coverage of <b>",input_sv$sv_read_r, "</b> reference reads and <b>",input_sv$sv_read_a, "</b> alternated reads."),paste0("The classification by CYTO-SV-ML pipeline is <b>",input_sv_class,"</b> with prediction probility -- System artifact SV: <b>", round(input_sv$prediction_TA,2),"</b>, True germline SV: <b>", round(input_sv$prediction_TG,2),"</b>, True cytogenetic somatic SV: <b>", round(input_sv$prediction_TS,2),"</b>."), paste0(""), paste0("(<u>The pair-plot and table at below are based on 100 benchmark SVs close to  genomic region of the simulated SV. The 3D dot plot at the right is based on simulated SV and benchmark SVs.</u>)"),sep="<br/>"))
            })                       
                                     
    ###################################################################################################################  
    #5. generate summary plot for input and benchmark sv data     
      output$scatter1 <-  renderPlotly({
            input_sv=input_subset_tmp()[[1]]        
            all_benchmark_sv= input_subset_tmp()[[3]] 
            input_sv$predict_max<-input_sv$predict_max/100
            all_sel_sv=rbind(input_sv,all_benchmark_sv)               
            all_sel_sv$alphag <- ifelse(all_sel_sv$label %in% c('TS','TA','TS'),0.2,1)
            p2<- all_sel_sv %>% plot_ly(type = 'scatter3d',mode = 'markers',x=~-log(prediction_TA,10), y=~-log(prediction_TG,10),z=~-log(prediction_TS,10),size=~-log(predict_max,10)*10, color=~factor(label),  colors =c(rgb(169,169,169, alpha = 0.01 * 255, maxColorValue = 255),rgb(127,255,0, alpha = 0.01 * 255, maxColorValue = 255), rgb(30,144,255, alpha = 0.01 * 255, maxColorValue = 255),rgb(255,0,0, alpha = 0.6 * 255, maxColorValue = 255)) ,marker=list(opacity=all_sel_sv$alphag),hovertemplate = paste("sv_type:",all_sel_sv$sv_type,"<br>sv_chr:",all_sel_sv$sv_chr,"<br>sv_chr2:",all_sel_sv$sv_chr2,"<br>sv_bp_st:",all_sel_sv$sv_bp_st,"<br>sv_bp_end:",all_sel_sv$sv_bp_end, "<br>TA_SV_prob:", round(all_sel_sv$prediction_TA,2),"<br>TG_SV_prob:", round(all_sel_sv$prediction_TG,2),"<br>TS_SV_prob:", round(all_sel_sv$prediction_TS,2),"<br>sv_database:",all_sel_sv$sv_database),opacity=0.36)#as.numeric(all_sel_sv$prediction_TS))#ifelse(all_sel_sv$label %in% c('TS','UC'),0,1)) 
            # p2 <-p2 %>%  add_trace(input_sv, x=~-log(prediction_TA,10), y=~-log(prediction_TG,10),z=~-log(prediction_TS,10), type = "scatter3d", mode = "markers",marker_symbol = 'star', marker_size = 15, color = I("red"), inherit = FALSE)
            p2 <- p2 %>% layout(scene = list(xaxis = list(title ="prediction_TA (-log)"), yaxis = list(title = "prediction_TG (-log)"),zaxis = list(title = "prediction_TS (-log)")))
            p2
    
      }) 
    
    ###################################################################################################################              
    #6. generate pair-plot for sv data
    output$scatter2 <- renderPlotly({
      req(input$feature_x,input$feature_y)          
      input_sv=input_subset_tmp()[[1]]
      benchmark_sv= input_subset_tmp()[[2]]              
      sel_sv=rbind(input_sv,benchmark_sv)   
      rownames(sel_sv)<-c('Unclassified_SV', paste('Benchmark SV_',1:dim(benchmark_sv)[1],sep=''))
      col_index<<-colnames(input_sv)
      # re-arrange data.frame for feature X Y Z 
      x=which(colnames(sel_sv)==input$feature_x)
      y=which(colnames(sel_sv)==input$feature_y)            
      x_y=which(1:dim(sel_sv)[2] %nin% c(x,y))
      all_sv_sel<-sel_sv[,c(x,y,x_y)]            
      
      #6.3 Plot prediction           
      p3<- all_sv_sel %>% plot_ly(type = 'scatter',mode = 'markers',x=~all_sv_sel[,1], y=~all_sv_sel[,2],size=~-log(all_sv_sel$predict_max,10)*10, color=~factor(all_sv_sel$label),  colors = c(rgb(169,169,169, alpha = 0.01 * 255, maxColorValue = 255),rgb(127,255,0, alpha = 0.01 * 255, maxColorValue = 255), rgb(30,144,255, alpha = 0.01 * 255, maxColorValue = 255),rgb(255,0,0, alpha = 0.6 * 255, maxColorValue = 255)), hovertemplate = paste("sv_type:",all_sv_sel$sv_type,"<br>sv_chr:",all_sv_sel$sv_chr,"<br>sv_chr2:",all_sv_sel$sv_chr2,"<br>sv_bp_st:",all_sv_sel$sv_bp_st,"<br>sv_bp_end:",all_sv_sel$sv_bp_end,"<br>sv_read_ref:",all_sv_sel$sv_read_r,"<br>sv_read_alt:",all_sv_sel$sv_read_a,"<br>sv_database:",all_sv_sel$sv_database),opacity=0.6)
      p3 <- p3 %>% layout(xaxis = list(title = as.character(input$feature_x)), yaxis = list(title = as.character(input$feature_y)))
      p3   
            })

    ###################################################################################################################              
    #7. generate summary table for sv data
       output$lineup1 <- renderLineup({
            req(input$feature_x,input$feature_y)          
            input_sv=input_subset_tmp()[[1]]
            benchmark_sv= input_subset_tmp()[[2]]              
            sel_sv=rbind(input_sv,benchmark_sv)   
            rownames(sel_sv)<-c('Unclassified_SV', paste('Benchmark SV_',1:dim(benchmark_sv)[1],sep=''))
            col_index<<-colnames(input_sv)
            # re-arrange data.frame for feature X Y Z 
            x=which(colnames(sel_sv)==input$feature_x)
            y=which(colnames(sel_sv)==input$feature_y)            
            x_y=which(1:dim(sel_sv)[2] %nin% c(x,y))
            all_sv_sel<-sel_sv[,c(x,y,x_y)]               
            all_sv_sel2<-all_sv_sel[,match(col_index,colnames(all_sv_sel))]
            all_sv_sel2<-all_sv_sel2[,colnames(all_sv_sel2) %nin% c('sv_chr','sv_chr2')]            
            lineup(all_sv_sel2, width = "100%")
             })
             
    }
