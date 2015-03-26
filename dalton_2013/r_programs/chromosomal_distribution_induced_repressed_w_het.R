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
INPUT <- paste(REPORTS,"chrom_enrichment_tables_w_het.csv",sep="/")

mydata <- read.csv(INPUT,header=TRUE);

####
female.ind <- subset(mydata, mydata$sex == "female" & mydata$flag_sig == "1" & mydata$flag_chrom == "1"& name == "induced")
female <- female.ind[6:7]
row.names(female) <- female.ind[[3]]

svg(paste(REPORTS,"/chrom_dist/female_induced_w_het.svg",sep=""))
barplot(t(female),beside=TRUE,legend.text=c("Observed","Expected"),main="Female Induced",ylim=c(0,100),args.legend=list(x="topright",inset=c(0,-0.1)))
dev.off()

####
male.ind <- subset(mydata, mydata$sex == "male" & mydata$flag_sig == "1" & mydata$flag_chrom == "1" & name == "induced")
male <- male.ind[6:7]
row.names(male) <- male.ind[[3]]

svg(paste(REPORTS,"/chrom_dist/male_induced_w_het.svg",sep=""))
barplot(t(male),beside=TRUE,legend.text=c("Observed","Expected"),main="Male Induced",ylim=c(0,300),args.legend=list(x="topright",inset=c(0,-0.1)))
dev.off()

####
female.rep <- subset(mydata, mydata$sex == "female" & mydata$flag_sig == "1" & mydata$flag_chrom == "1"& name == "repressed")
female <- female.rep[6:7]
row.names(female) <- female.rep[[3]]

svg(paste(REPORTS,"/chrom_dist/female_repressed_w_het.svg",sep=""))
barplot(t(female),beside=TRUE,legend.text=c("Observed","Expected"),main="Female Repressed",ylim=c(0,130),args.legend=list(x="topright",inset=c(0,-0.1)))
dev.off()

####
male.rep <- subset(mydata, mydata$sex == "male" & mydata$flag_sig == "1" & mydata$flag_chrom == "1" & name == "repressed")
male <- male.rep[6:7]
row.names(male) <- male.rep[[3]]

svg(paste(REPORTS,"/chrom_dist/male_repressed_w_het.svg",sep=""))
barplot(t(male),beside=TRUE,legend.text=c("Observed","Expected"),main="Male Repressed",ylim=c(0,150),args.legend=list(x="topright",inset=c(0,-0.1)))
dev.off()
