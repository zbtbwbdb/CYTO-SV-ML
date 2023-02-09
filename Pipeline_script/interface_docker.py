library(optparse)
library(cluster)
library(reghelper)
library(bios2mds)
library(Hmisc)
library(stringr)
library(GOSemSim)
library(org.Hs.eg.db)
library('biomaRt')
################################################################################################################################################################################

option_list<-list(make_option(c("-i", "--input_file"),type="character",help="input file", default=NA,metavar="filename"), make_option(c("-o", "--output_file"),type="character",help="output file", default=NA,metavar="filename"))
opt_parser<-OptionParser(option_list=option_list)
opt=parse_args(opt_parser)
input<-as.character(opt$input_file)
opc<-as.character(opt$output_file)

# load the optimal model

# load the training data tranformer

# read in all WGS SV data

# transform all the original SV data

# predict all the transformed SV data

# attach the prediction label to all the original SV data

# save all the SV data

