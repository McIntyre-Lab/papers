# Get command line args with file names for processing
args <- commandArgs(TRUE)

# MBE PG functions
source("/scratch/lfs/mcintyre/cegs_ase_paper/scripts/r_lib/Subroutines_model2_experimental.R")
source("/scratch/lfs/mcintyre/cegs_ase_paper/scripts/r_lib/AI_poissongamma_functions.R")

# Prepare output file
fileout = args[2]
headers_out = "line,mating_status,fusion_id,qsim_mean_theta,qsim_q025,qsim_q975,Bayesianpvalue_qsim,flag_AI_qsim"
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
    qsim <- as.numeric(mydata["qsim_tester"])
    n_i <- X_RNA + Y_RNA

    #Poisson Gamma with q random and obtained from the DNA
    #Poisson Gamma models q = .4
    PGSIM <- gibbs_poissongamma(nsim=1000,
                                nburnin=1000,
                                lag=10,
                                x=X_RNA,
                                y=Y_RNA,
                                both=c(0,0,0),
                                a_mu=1/2,
                                b_mu=1/2,
                                a_alpha=1/2,
                                b_alpha=1/2,
                                a_beta=1/2,
                                b_beta=1/2,
                                q_=qsim,
                                abundance=FALSE
                               )

    pps <- PGSIM$alphas/(1+PGSIM$alphas)
    qs025 <- quantile(pps,c(.025),na.rm=TRUE)
    qs975 <- quantile(pps,c(.975),na.rm=TRUE)
    means <- mean(pps)
    Bayesianpvalue <- 2*min(mean(pps<1/2),mean(pps>1/2))
    flag_AI <- ifelse(Bayesianpvalue<0.05,1,0)

    # Create output and write to table
    SNPout = paste(mydata["line"],mydata["mating_status"],mydata["fusion_id"],paste(round(c(means,qs025,qs975,Bayesianpvalue),3),collapse=","),flag_AI,sep=",")
    cat(SNPout,file=fileout,append=TRUE,sep="\n")
}
