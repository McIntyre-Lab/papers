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
INPUT <- paste(REPORTS,"chrom_enrichment_tables_by_fru.csv",sep="/")

mydata <- read.csv(INPUT,header=TRUE);

#### 
female.a.ind <- subset(mydata, mydata$sex == "female" & mydata$flag_sig == "1" & mydata$flag_chrom == "1"& (name == "a_ind" ))
female.a <- female.a.ind[6:7]
row.names(female.a) <- female.a.ind[[3]]

svg(paste(REPORTS,"/chrom_dist/female_fru_a_ind.svg",sep=""))
barplot(t(female.a),beside=TRUE,legend.text=c("Observed","Expected"),main="Female Fru A Induced",ylim=c(0,100),args.legend=list(x="topright",inset=c(0,-0.1)))
dev.off()


female.b.ind <- subset(mydata, mydata$sex == "female" & mydata$flag_sig == "1" & mydata$flag_chrom == "1"& (name == "b_ind" ))
female.b <- female.b.ind[6:7]
row.names(female.b) <- female.b.ind[[3]]

svg(paste(REPORTS,"/chrom_dist/female_fru_b_ind.svg",sep=""))
barplot(t(female.b),beside=TRUE,legend.text=c("Observed","Expected"),main="Female Fru B Induced",ylim=c(0,100),args.legend=list(x="topright",inset=c(0,-0.1)))
dev.off()


female.c.ind <- subset(mydata, mydata$sex == "female" & mydata$flag_sig == "1" & mydata$flag_chrom == "1"& (name == "c_ind" ))
female.c <- female.c.ind[6:7]
row.names(female.c) <- female.c.ind[[3]]

svg(paste(REPORTS,"/chrom_dist/female_fru_c_ind.svg",sep=""))
barplot(t(female.c),beside=TRUE,legend.text=c("Observed","Expected"),main="Female Fru C Induced",ylim=c(0,100),args.legend=list(x="topright",inset=c(0,-0.1)))
dev.off()

#### 
male.a.ind <- subset(mydata, mydata$sex == "male" & mydata$flag_sig == "1" & mydata$flag_chrom == "1"& (name == "a_ind" ))
male.a <- male.a.ind[6:7]
row.names(male.a) <- male.a.ind[[3]]

svg(paste(REPORTS,"/chrom_dist/male_fru_a_ind.svg",sep=""))
barplot(t(male.a),beside=TRUE,legend.text=c("Observed","Expected"),main="Male Fru A Induced",ylim=c(0,300),args.legend=list(x="topright",inset=c(0,-0.1)))
dev.off()


male.b.ind <- subset(mydata, mydata$sex == "male" & mydata$flag_sig == "1" & mydata$flag_chrom == "1"& (name == "b_ind" ))
male.b <- male.b.ind[6:7]
row.names(male.b) <- male.b.ind[[3]]

svg(paste(REPORTS,"/chrom_dist/male_fru_b_ind.svg",sep=""))
barplot(t(male.b),beside=TRUE,legend.text=c("Observed","Expected"),main="Male Fru B Induced",ylim=c(0,300),args.legend=list(x="topright",inset=c(0,-0.1)))
dev.off()


male.c.ind <- subset(mydata, mydata$sex == "male" & mydata$flag_sig == "1" & mydata$flag_chrom == "1"& (name == "c_ind" ))
male.c <- male.c.ind[6:7]
row.names(male.c) <- male.c.ind[[3]]

svg(paste(REPORTS,"/chrom_dist/male_fru_c_ind.svg",sep=""))
barplot(t(male.c),beside=TRUE,legend.text=c("Observed","Expected"),main="Male Fru C Induced",ylim=c(0,300),args.legend=list(x="topright",inset=c(0,-0.1)))
dev.off()

#### 
female.a.rep <- subset(mydata, mydata$sex == "female" & mydata$flag_sig == "1" & mydata$flag_chrom == "1"& (name == "a_rep" ))
female.a <- female.a.rep[6:7]
row.names(female.a) <- female.a.rep[[3]]

svg(paste(REPORTS,"/chrom_dist/female_fru_a_rep.svg",sep=""))
barplot(t(female.a),beside=TRUE,legend.text=c("Observed","Expected"),main="Female Fru A Repressed",ylim=c(0,100),args.legend=list(x="topright",inset=c(0,-0.1)))
dev.off()


female.b.rep <- subset(mydata, mydata$sex == "female" & mydata$flag_sig == "1" & mydata$flag_chrom == "1"& (name == "b_rep" ))
female.b <- female.b.rep[6:7]
row.names(female.b) <- female.b.rep[[3]]

svg(paste(REPORTS,"/chrom_dist/female_fru_b_rep.svg",sep=""))
barplot(t(female.b),beside=TRUE,legend.text=c("Observed","Expected"),main="Female Fru B Repressed",ylim=c(0,100),args.legend=list(x="topright",inset=c(0,-0.1)))
dev.off()


female.c.rep <- subset(mydata, mydata$sex == "female" & mydata$flag_sig == "1" & mydata$flag_chrom == "1"& (name == "c_rep" ))
female.c <- female.c.rep[6:7]
row.names(female.c) <- female.c.rep[[3]]

svg(paste(REPORTS,"/chrom_dist/female_fru_c_rep.svg",sep=""))
barplot(t(female.c),beside=TRUE,legend.text=c("Observed","Expected"),main="Female Fru C Repressed",ylim=c(0,100),args.legend=list(x="topright",inset=c(0,-0.1)))
dev.off()

#### 
male.a.rep <- subset(mydata, mydata$sex == "male" & mydata$flag_sig == "1" & mydata$flag_chrom == "1"& (name == "a_rep" ))
male.a <- male.a.rep[6:7]
row.names(male.a) <- male.a.rep[[3]]

svg(paste(REPORTS,"/chrom_dist/male_fru_a_rep.svg",sep=""))
barplot(t(male.a),beside=TRUE,legend.text=c("Observed","Expected"),main="Male Fru A Repressed",ylim=c(0,100),args.legend=list(x="topright",inset=c(0,-0.1)))
dev.off()


male.b.rep <- subset(mydata, mydata$sex == "male" & mydata$flag_sig == "1" & mydata$flag_chrom == "1"& (name == "b_rep" ))
male.b <- male.b.rep[6:7]
row.names(male.b) <- male.b.rep[[3]]

svg(paste(REPORTS,"/chrom_dist/male_fru_b_rep.svg",sep=""))
barplot(t(male.b),beside=TRUE,legend.text=c("Observed","Expected"),main="Male Fru B Repressed",ylim=c(0,100),args.legend=list(x="topright",inset=c(0,-0.1)))
dev.off()


male.c.rep <- subset(mydata, mydata$sex == "male" & mydata$flag_sig == "1" & mydata$flag_chrom == "1"& (name == "c_rep" ))
male.c <- male.c.rep[6:7]
row.names(male.c) <- male.c.rep[[3]]

svg(paste(REPORTS,"/chrom_dist/male_fru_c_rep.svg",sep=""))
barplot(t(male.c),beside=TRUE,legend.text=c("Observed","Expected"),main="Male Fru C Repressed",ylim=c(0,100),args.legend=list(x="topright",inset=c(0,-0.1)))
dev.off()
