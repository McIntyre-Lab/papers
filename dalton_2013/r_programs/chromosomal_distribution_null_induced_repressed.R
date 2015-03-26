#===============================================================================
# 
#   DESCRIPTION: This script creates a barchart of showing the distribution of
#   induced or repressed genes accross the different chromosomes.
# 
#       OPTIONS: ---
#  REQUIREMENTS: ---
#        AUTHOR: Justin Fear (JMF), jfear@ufl.edu
#  ORGANIZATION: University of Florida
#       CREATED: 08/27/2012 12:59:04 PM EDT
#      REVISION:  ---
#===============================================================================

library(RSvgDevice)
PROJ <- "/home/jfear/mclab/arbeitman_fru_network";
REPORTS <- paste(PROJ,"reports_external",sep="/");
INPUT <- paste(REPORTS,"chrom_enrichment_tables_null_v2.csv",sep="/")

mydata <- read.csv(INPUT,header=TRUE);

####
ind <- subset(mydata, mydata$flag_sig == "1" & mydata$flag_chrom == "1" & name == "induced")
small.ind <- ind[5:6]
row.names(small.ind) <- ind[[2]]

svg(paste(REPORTS,"/chrom_dist/null_induced_v2.svg",sep=""))
barplot(t(small.ind),beside=TRUE,legend.text=c("Observed","Expected"),main="Null Induced",ylim=c(0,200),args.legend=list(x="topright",inset=c(0,-0.1)))
dev.off()

####
rep <- subset(mydata, mydata$flag_sig == "1" & mydata$flag_chrom == "1" & name == "repressed")
small.rep <- rep[5:6]
row.names(small.rep) <- rep[[2]]

svg(paste(REPORTS,"/chrom_dist/null_repressed_v2.svg",sep=""))
barplot(t(small.rep),beside=TRUE,legend.text=c("Observed","Expected"),main="Null Repressed",ylim=c(0,200),args.legend=list(x="topright",inset=c(0,-0.1)))
dev.off()
