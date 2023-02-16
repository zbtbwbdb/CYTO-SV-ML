library(optparse)

################################################################################################################################################################################

option_list<-list(make_option(c("-i", "--input_file"),type="character",help="input file", default=NA,metavar="filename"), make_option(c("-o", "--output_file"),type="character",help="output file", default=NA,metavar="filename"), make_option(c("-c", "--cyto_band_file"),type="character",help="cyto band file", default=NA,metavar="filename"))
opt_parser<-OptionParser(option_list=option_list)
opt=parse_args(opt_parser)
input<-as.character(opt$input_file)
output<-as.character(opt$output_file)
cyto_band_file<-as.character(opt$cyto_band_file)

# read in trs and nontrs data
all_sv=read.table(paste(input,'.all_pred',sep=''),sep='\t',header=T)
genome_coor=read.table(cyto_band_file,sep='\t',header=T)

# save the SV data
save(all_sv,genome_coor, file = output)
