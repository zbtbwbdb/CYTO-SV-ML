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
    check_sv_chr <- eventReactive(input$Make_simulated_SV, {input$sv_chr})
    check_sv_chr2 <- eventReactive(input$Make_simulated_SV, {input$sv_chr2})
    check_sv_type <- eventReactive(input$Make_simulated_SV, {input$sv_type})   
    check_cytoband <- eventReactive(input$Make_simulated_SV, {input$cytoband})
    check_feature_x <- eventReactive(input$Plot_SV_data, {input$feature_x})
    check_feature_y <- eventReactive(input$Plot_SV_data, {input$feature_y})
    check_feature_z <- eventReactive(input$Plot_SV_data, {input$feature_z})
 
    ###################################################################################################################  
    #4. generate the input data subset for recipient list table 
    input_subset_tmp<- reactive({
          #4.1 load assigned input values
            sv_chr<-check_sv_chr()
            sv_type<-check_sv_type()
            sv_chr_arm<-check_cytoband()

            #3.1 generate the simulated examples from random resampling
            bp_min=genome_coor[genome_coor$chr==sv_chr & genome_coor$chr_arm=sv_chr_arm,3]
            bp_max=genome_coor[genome_coor$chr==sv_chr & genome_coor$chr_arm=sv_chr_arm,4]
            if (sv_type=='BND')
                  {
                  benchmark_sv=read.table("cnv_benchmark_all.vcf",sep="\t",header=T)
                  input_sv=sv_simu_database[sv_simu_database$Chrom==sv_chr & sv_simu_database$chr2==sv_chr2 & sv_simu_database$bp_st<=bp_max & sv_simu_database$bp_st>=bp_min,]
            } else {
                  benchmark_sv=read.table("trs_benchmark_all.vcf",sep="\t",header=T)
                  input_sv=sv_simu_database[sv_simu_database$sv_chr==sv_chr & sv_simu_database$sv_type==sv_type & sv_simu_database$500.inv_sv_all_bp_st<=bp_max & sv_simu_database$500.inv_sv_all_bp_st>=bp_min,]
                  }
            write.table(input_sv,"input_sv.vcf",sep="\t",col.names=T,row.names=F)
            #3.2 model prediction by system call python script   
            pred_sv=system("python","cyto-sv-ml.py","input_sv.vcf", stdout = TRUE, stderr = TRUE)
            input_sv['prediction']=int(pred_sv)
            #3.3 return data frame
            return (list(input_sv,benchmark_sv))     
            })
    
    ###################################################################################################################  
     
    #5. generate recipient list table with the input data subset from 4.
    output$disttable1 <- renderDataTable( 
          input_subset_tmp()[,c()], extensions = 'Buttons', 
          options = list(dom = 'Bfrtip',
                         buttons = c('copy', 'csv', 'excel', 'pdf', 'print')),
          rownames= FALSE
          )
                                          
    ###################################################################################################################   
    
    #6. detect selected sample for visualiztion of BART3 prediction plot and pop up messege & generate table for significantly optimal donor available
   output$distPlot1 <- renderPlot({
            feature_x <- check_feature_x() 
            feature_y <- check_feature_y() 
            feature_z <- check_feature_z()            
            input_sv=input_subset_tmp()[[1]]
            benchmark_sv= input_subset_tmp()[[2]]  
            input_sv_cr= input_sv[,colnames(input_sv) %in% colnames(benchmark_sv)]  
            all_sv=rbind(input_sv_cr,benchmark_sv)                                                             
          ###################################################################################################################            
          #6.5.4 Plot posterior mean of selected recipient & donors from BART3 model    
###   The changes of line 175 for correcting the bug of confident interval by mislabelling dfs2 as dfs1 for ribbons data source          
          p2<- all_sv %>% plot_ly(type = 'scatter',mode = 'markers',x=~all_sv[,1], y=~all_sv[,2],z=all_sv[,3],symbol=~factor(label),color=~factor(prediction),size=~factor(sv_class), colors = c( "blue", "red"),hovertemplate = paste("prediction_mean:",pred_mean,"<br>upper_CI:",pred_upper,"<br>lower_CI:",pred_lower,"<br>recipient disease:",disease,"<br>recipient disease stage:",disstat,"<br>recipient race:",race,"<br>donor age:", dnrage)) %>% add_ribbons(data=dfs2, ymin =~lower , ymax = ~upper, line = list(color = 'rgba(7, 164, 181, 0.05)'), fillcolor = 'rgba(7, 164, 181, 0.2)', name = '95% CI')
          #p2 <- p2  %>% layout(yaxis = list(title = "One year death/GVHD risk probability", range = c(0,1)))  %>% layout(xaxis = list(title = "Sample No.", range = c(0,51)))
          p2      })
    }
