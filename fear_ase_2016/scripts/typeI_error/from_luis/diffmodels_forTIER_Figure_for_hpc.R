#Program used to analyze the simulated data sets 
#generated by R program simulatingdata.R
#explained in Synthetic Data Set Plan_LLN08072014.docx
#For two scenarios 
#"NOAI_YESBIAS"
#"YESAI_YESBIAS"

#HPC=				0		#If this flag is 1 it used the path in the HPC
#path=				ifelse(HPC==1,"","/Users/luis/projects/lauren/")
##pathfigs=		paste(path,"paper/figs/",sep="")
#pathprograms=	paste(path,"rprograms/",sep="")
#pathdata=		paste(path,"data/forTIERfigure/",sep="")
#pathresults=	paste(path,"simulationresults/forTIERfigure/",sep="")
#
#source(paste(pathprograms,"AI_poissongamma_functions.R",sep=""))#Subrutines for Poissong Gamma
#
#
#scenario=	"NOAI"
#
#
#
#
#q_TRUE=		0.35	#c(0.35,0.375,0.4,0.425,0.45,0.46,0.47,0.48,0.49,0.50,0.51,0.52,0.53,0.54,0.55,0.575,0.6,0.625,0.65)
#					#Run this program 	19 times! 
#					#Once qith each value of q_TRUE above
#fortitleq_true=paste("q_true_0d",q_TRUE*1000,sep="")
#
#
#fileinname<-paste("simulateddata_forTIER_Figure",scenario,"Poisson",fortitleq_true,sep="")
#
#filein=		paste(pathdata,fileinname,".csv",sep="")
#fileout=	paste(pathresults,"results_",fileinname,".csv",sep="")
			
# Get command line args with filenames for processing
args <- commandArgs(TRUE)

# Import Function
#source("/scratch/lfs/mcintyre/luis/projects/SNP/rprograms/Subroutines_model2_experimental.R")#Subrutines in Rita s paper
source("/scratch/lfs/mcintyre/luis/projects/lauren/rprograms/AI_poissongamma_functions.R")#Subrutines for Poissong Gamma
			
# File paths
con= file(args[1],"r")
fileout= args[2]
q_TRUE = as.numeric(args[3])
			

headers_out=paste("FUSION_ID",
"S_RNA","TotalRNA","PER_S_RNA",
"estimate_binomial_test","q025_binomial_test","q975_binomial_test","binomial_test_pvalue","flag_AI_binomial_test",
"q_true",
"mean_q_true","q025_q_true","q975_q_true","Bayesianpvalue_q_true","flag_AI_q_true",
"mean_q_ahalf",	"q025_q_ahalf","q975_q_ahalf","Bayesianpvalue_q_ahalf","flag_AI_q_ahalf",
"q_sim_below_01per",
"mean_q_sim_below_01per","q025_q_sim_below_01per","q975_q_sim_below_01per","Bayesianpvalue_q_sim_below_01per","flag_AI_q_sim_below_01per",
"q_sim_above_01per",
"mean_q_sim_above_01per","q025_q_sim_above_01per","q975_q_sim_above_01per","Bayesianpvalue_q_sim_above_01per","flag_AI_q_sim_above_01per",
"q_sim_below_02per",
"mean_q_sim_below_02per","q025_q_sim_below_02per","q975_q_sim_below_02per","Bayesianpvalue_q_sim_below_02per","flag_AI_q_sim_below_02per",
"q_sim_above_02per",
"mean_q_sim_above_02per","q025_q_sim_above_02per","q975_q_sim_above_02per","Bayesianpvalue_q_sim_above_02per","flag_AI_q_sim_above_02per",
"q_sim_below_05per",
"mean_q_sim_below_05per","q025_q_sim_below_05per","q975_q_sim_below_05per","Bayesianpvalue_q_sim_below_05per","flag_AI_q_sim_below_05per",
"q_sim_above_05per",
"mean_q_sim_above_05per","q025_q_sim_above_05per","q975_q_sim_above_05per","Bayesianpvalue_q_sim_above_05per","flag_AI_q_sim_above_05per",
"q_sim_below_10per",
"mean_q_sim_below_10per","q025_q_sim_below_10per","q975_q_sim_below_10per","Bayesianpvalue_q_sim_below_10per","flag_AI_q_sim_below_10per",
"q_sim_above_10per",
"mean_q_sim_above_10per","q025_q_sim_above_10per","q975_q_sim_above_10per","Bayesianpvalue_q_sim_above_10per","flag_AI_q_sim_above_10per",
sep=",")


#con=		file(filein,"r")
percentages=c(1,2,5,10)
qsims=q_TRUE*(1+c(-1,1)*rep(percentages,each=2)/100)
q_sims_fornames=paste(c("below","above"),rep(percentages,each=2),"per",sep="")

cat(headers_out,file=fileout,append=FALSE,sep="\n")


newline<-readLines(con,n=1) #moving the pointer away from the headers
headers_in=strsplit(newline,split=",")[[1]]



means=qs025=qs975=flag_AI=rep(NA,11)
names(qs025)=names(qs975)=names(means)=names(flag_AI)=c("q_true",
"q_ahalf","Confidence",q_sims_fornames)
Bayesianpvalue=rep(NA,10)
names(Bayesianpvalue)=c("q_true","q_ahalf",q_sims_fornames)	
mydata=rep(NA,length(headers_in))
names(mydata)=headers_in
m=1

#newline<-readLines(con,n=8)
while(length(newline) != 0 ){
newline<-readLines(con,n=1)	
if (length(newline) == 0){break}
mydata<-as.vector(strsplit(newline,split=",")[[1]])
names(mydata)=headers_in
	
print(paste("------------Analyzing fusion--------------",m));m=m+1


x	<-	as.numeric(mydata[c("RNA_Mel_all_Rep1","RNA_Mel_all_Rep2","RNA_Mel_all_Rep3")])
x2	<-	as.numeric(mydata[c("RNA_Sim_all_Rep1","RNA_Sim_all_Rep2","RNA_Sim_all_Rep3")])
q_true=as.numeric(mydata["q_true"])


#Poisson gamma with q_true###############################
tem_PG=gibbs_poissongamma(nsim=10000,nburnin=100,lag=5,
x=as.numeric(x),y=as.numeric(x2),both=c(0,0,0),
a_mu=	.1,
b_mu=	.1,
a_alpha=50,
b_alpha=50,
a_beta=	1/2,
b_beta=	1/2,
q_=		q_true,abundance=FALSE)
pps=		tem_PG$alphas/(1+tem_PG$alphas)

qs025["q_true"]=	quantile(pps,c(.025),na.rm=TRUE)
qs975["q_true"]=	quantile(pps,c(.975),na.rm=TRUE)
means["q_true"]=	mean(pps)
Bayesianpvalue["q_true"]=2*min(mean(pps<1/2),mean(pps>1/2))
flag_AI["q_true"]=ifelse(Bayesianpvalue["q_true"]<0.05,1,0)

#Poisson Gammanaive i.e. q =1/2:#############################
tem_PGnaive=gibbs_poissongamma(nsim=10000,nburnin=100,lag=5,
x=as.numeric(x),y=as.numeric(x2),both=c(0,0,0),
a_mu=.1,
b_mu=.1,
a_alpha=50,
b_alpha=50,
a_beta=1/2,
b_beta=1/2,
q_=0.5,abundance=FALSE)
pps=tem_PGnaive$alphas/(1+tem_PGnaive$alphas)

quants=quantile(pps,c(.025,0.975),na.rm=TRUE)
qs025["q_ahalf"]=	quants[1]
qs975["q_ahalf"]=	quants[2]
means["q_ahalf"]=	mean(pps)
Bayesianpvalue["q_ahalf"]=2*min(mean(pps<1/2),mean(pps>1/2))
flag_AI["q_ahalf"]=ifelse(Bayesianpvalue["q_ahalf"]<0.05,1,0)

#############################################
#Poisson gamma with different q_sims
forout=c()
for(ii in 1:length(qsims)){
tem_PG=gibbs_poissongamma(nsim=10000,nburnin=100,lag=5,
x=as.numeric(x),y=as.numeric(x2),both=c(0,0,0),
a_mu=	.1,
b_mu=	.1,
a_alpha=50,
b_alpha=50,
a_beta=	1/2,
b_beta=	1/2,
q_=		qsims[ii],abundance=FALSE)
pps=		tem_PG$alphas/(1+tem_PG$alphas)

qlabel=	q_sims_fornames[ii]
quants=quantile(pps,c(.025,0.975),na.rm=TRUE)
qs025[qlabel]=quants[1]
qs975[qlabel]=	quants[2]
means[qlabel]=	mean(pps)
Bayesianpvalue[qlabel]=2*min(mean(pps<1/2),mean(pps>1/2))
flag_AI[qlabel]=ifelse(Bayesianpvalue[qlabel]<0.05,1,0)
forout=c(forout,qsims[ii],means[qlabel],qs025[qlabel],qs975[qlabel],Bayesianpvalue[qlabel],flag_AI[qlabel])
}
#################################################################
#Confidence interval
n_i=x+x2
theta_star=			1-sum(x)/sum(n_i)
binomial_test=	binom.test(x=sum(n_i-x), n=sum(n_i), p=0.5, conf.level=0.95)


binomial_test_pvalue=binomial_test$p.value

#lims=				theta_star+c(-1,1)*1.965*sqrt(theta_star*(1-theta_star)/sum(n_i))
qs025["Confidence"]=binomial_test$conf.int[1]
qs975["Confidence"]=binomial_test$conf.int[2]
means["Confidence"]=binomial_test$estimate
flag_AI["Confidence"]=ifelse(binomial_test_pvalue<= 0.05,1,0)



#################################Saving results
SNPout=paste(mydata["FUSION_ID"],
paste(
round(
c(sum(x2),sum(n_i),round(sum(x2)/sum(n_i),3),
means["Confidence"],qs025["Confidence"],qs975["Confidence"],binomial_test_pvalue,flag_AI["Confidence"],
q_true,means["q_true"],qs025["q_true"],qs975["q_true"],Bayesianpvalue["q_true"],flag_AI["q_true"],
means["q_ahalf"],qs025["q_ahalf"],qs975["q_ahalf"],Bayesianpvalue["q_ahalf"],flag_AI["q_ahalf"],
forout)
,3),
collapse=",")
,sep=",")

print(SNPout)
cat(SNPout,file=fileout,append=TRUE,sep="\n")


}
close(con)

