# Get command line args with file names for processing
args <- commandArgs(TRUE)

# MBE PG functions
source("/scratch/lfs/mcintyre/cegs_ase_paper/scripts/r_lib/Subroutines_model2_experimental.R")
source("/scratch/lfs/mcintyre/cegs_ase_paper/scripts/r_lib/AI_poissongamma_functions.R")

# Prepare output file
fileout = args[2]
headers_out = "line,mating_status,fusion_id,q4_mean_theta,q4_q025,q4_q975,q5_mean_theta,q5_q025,q5_q975,q6_mean_theta,q6_q025,q6_q975"
cat(headers_out,file=fileout,append=FALSE,sep="\n")

# Make Connection to input
con = file(args[1],"r")
newline<-readLines(con,n=1) # Go to header line
headers_in=strsplit(newline,split=",")[[1]]

mydata=rep(NA,length(headers_in))
names(mydata)=headers_in
m=1

while(length(newline) != 0 ){
    # Move past fusions that are not flag_analyze
    flaganalyze=0
    while(flaganalyze==0){
        newline<-readLines(con,n=1) 
        if(length(newline)==0){break}
        mydata<-as.vector(strsplit(newline,split=",")[[1]])
        names(mydata)=headers_in
        flaganalyze=as.numeric(mydata["flag_analyze"])
    }
    if(length(newline)==0){break}

    print(paste("------------Analyzing",mydata['line'],mydata['mating_status'],mydata['fusion_id'], "--------------"));
    m=m+1
    X_RNA <- as.numeric(mydata[c("LINE_TOTAL_1","LINE_TOTAL_2","LINE_TOTAL_3","LINE_TOTAL_4","LINE_TOTAL_5","LINE_TOTAL_6")])
    X_RNA <- X_RNA[!is.na(X_RNA)]
    Y_RNA <- as.numeric(mydata[c("TESTER_TOTAL_1","TESTER_TOTAL_2","TESTER_TOTAL_3","TESTER_TOTAL_4","TESTER_TOTAL_5","TESTER_TOTAL_6")])
    Y_RNA <- Y_RNA[!is.na(Y_RNA)]
    n_i <- X_RNA + Y_RNA

    # Create storage vectors for Quantiles and posterior means for theta (proportion of male reads) under different models
    means=q_025=q_975=rep(NA,3)
    names(means)=names(q_025)=names(q_975)=c("PG_4","PG_5","PG_6")

    #Poisson Gamma models q = .4
    tem_PG=gibbs_poissongamma(nsim=1000,nburnin=1000,lag=10,x=X_RNA,y=Y_RNA,both=c(0,0,0),a_mu=1/2,b_mu=1/2,a_alpha=1/2,b_alpha=1/2,a_beta=1/2,b_beta=1/2,
    q_=.4,#<-------As before but with q=q_sim 
    abundance=FALSE)

    thetas=tem_PG$alphas/(1+tem_PG$alphas)
    q_025["PG_4"]= quantile(thetas,c(.025),na.rm=TRUE)
    q_975["PG_4"]= quantile(thetas,c(.975),na.rm=TRUE)
    means["PG_4"]= mean(thetas)

    #Poisson Gamma models q = .5
    tem_PG=gibbs_poissongamma(nsim=1000,nburnin=1000,lag=10,x=X_RNA,y=Y_RNA,both=c(0,0,0),a_mu=1/2,b_mu=1/2,a_alpha=1/2,b_alpha=1/2,a_beta=1/2,b_beta=1/2,
    q_=.5,#<-------As before but with q=q_sim 
    abundance=FALSE)

    thetas=tem_PG$alphas/(1+tem_PG$alphas)
    q_025["PG_5"]= quantile(thetas,c(.025),na.rm=TRUE)
    q_975["PG_5"]= quantile(thetas,c(.975),na.rm=TRUE)
    means["PG_5"]= mean(thetas)

    #Poisson Gamma models q = .6
    tem_PG=gibbs_poissongamma(nsim=1000,nburnin=1000,lag=10,x=X_RNA,y=Y_RNA,both=c(0,0,0),a_mu=1/2,b_mu=1/2,a_alpha=1/2,b_alpha=1/2,a_beta=1/2,b_beta=1/2,
    q_=.6,#<-------As before but with q=q_sim 
    abundance=FALSE)

    thetas=tem_PG$alphas/(1+tem_PG$alphas)
    q_025["PG_6"]= quantile(thetas,c(.025),na.rm=TRUE)
    q_975["PG_6"]= quantile(thetas,c(.975),na.rm=TRUE)
    means["PG_6"]= mean(thetas)

    # Create output and write to table
    SNPout = paste(mydata["line"],mydata["mating_status"],mydata["fusion_id"],paste(round(c(means["PG_4"],q_025["PG_4"],q_975["PG_4"]),3),collapse=","),paste(round(c(means["PG_5"],q_025["PG_5"],q_975["PG_5"]),3),collapse=","),paste(round(c(means["PG_6"],q_025["PG_6"],q_975["PG_6"]),3),collapse=","),sep=",")
    cat(SNPout,file=fileout,append=TRUE,sep="\n")
}
