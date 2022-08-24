library(shiny)
#library(package)
library(ggforce)
library(ggplot2)
library(gridExtra)
library(DALEX)
library(hash)
library(DT)
library(cluster)
library(remotes)
library(stringr)
#install_github("rsparapa/bnptools/BART3")
library(BART3)
library(shinydashboard)
library(plotly)
#
#options(encoding = 'UTF-8')
#options(timeout= 1200)

# #setwd("C:/Users/tzhang/Desktop")
# data<-read.table("test.cr.cr.cr.cr.cr.cr.cr.cr_0.01",sep="\t",header=F,stringsAsFactors = F)
# data[,8]<-1
# colnames(data)<-c('Variant_type','Genomad_AF_log','VAF','COSMIC_variant','MDS_gene','Variant_source','samples', 'All')

# save(data, file='./data/MDSdata2.RData')


# data preparation
# data_mds<-read.table("./data/test.all_g.mds.sim",sep="\t",header=F,stringsAsFactors = F)
# data_mds_coding<-data_mds[data_mds[,2]!='noncoding',]
# data<-read.table("./data/test.cr.cr.cr.cr.cr.cr.cr.cr_0.01",sep="\t",header=F,stringsAsFactors = F)
# MDS_gene<-read.table("./data/MDS_gene_conservative",sep="\t",header=F,stringsAsFactors = F)
# data[,8]<-1
# colnames(data)<-c('Variant_type','Genomad_AF_log','VAF','COSMIC_variant','MDS_gene','Variant_source','samples', 'All')
# 
# data_mds_coding[,9]<-1
# colnames(data_mds_coding)<-c('MDS_prognostic_gene','Variant_type','Genomad_AF_log','VAF','COSMIC_variant','MDS_gene','Variant_source','samples', 'All')
# 
# save(data, data_mds_coding, MDS_gene, file='./data/MDSdata2.RData')

