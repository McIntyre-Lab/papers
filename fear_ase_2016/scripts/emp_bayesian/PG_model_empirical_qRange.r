# Get command line args with file names for processing
args <- commandArgs(TRUE)

# MBE PG functions
source("/scratch/lfs/mcintyre/cegs_ase_paper/scripts/r_lib/Subroutines_model2_experimental.R")
source("/scratch/lfs/mcintyre/cegs_ase_paper/scripts/r_lib/AI_poissongamma_functions.R")

nameFunc = function(q, pre='', post=''){
    # Function that takes a vector and merges on pre and post to every member
    # of that vector and returns a vector of these strings.

    # Convert to string and remove decimal point
    qStr = sub('\\.', '', toString(q))

    # Concatenate on pre and post
    combine = paste(pre, qStr, post, sep='')
    return(combine)
}

# Generate range of q's
#qRange = seq(.01, 0.99, by=0.01)
qRange = c(0.4, 0.5, 0.6)

# Prepare names for header
n1 = sapply(qRange, nameFunc, pre='q', post='_mean_theta')
n2 = sapply(qRange, nameFunc, pre='q', post='_q025')
n3 = sapply(qRange, nameFunc, pre='q', post='_q975')
nVector = as.vector(t(cbind(n1, n2, n3)))

# Prepare output file
fileout = args[2]
headers_out = paste("line,mating_status,fusion_id,", paste(nVector, collapse=','), sep='')
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
    means=q_025=q_975=rep(NA,length(qRange))
    n4 = sapply(qRange, nameFunc, pre='PG_')
    names(means)=names(q_025)=names(q_975)=n4

    for(i in 1:length(qRange)){
        # Iterate over values of q
        q = qRange[i]
        n = n4[i]

        #Poisson Gamma models q
        tem_PG = gibbs_poissongamma(nsim=1000, nburnin=1000, lag=10, x=X_RNA, y=Y_RNA, both=c(0,0,0), a_mu=1/2, b_mu=1/2, a_alpha=1/2, b_alpha=1/2, a_beta=1/2, b_beta=1/2, 
                                  q_=q, abundance=FALSE)

        thetas = tem_PG$alphas/(1+tem_PG$alphas)
        q_025[n] = quantile(thetas, c(.025), na.rm=TRUE)
        q_975[n] = quantile(thetas, c(.975), na.rm=TRUE)
        means[n] = mean(thetas)
    }

    # Combine into an output table
    outTable = rbind(q_025, means, q_975)

    # Create output and write to table
    SNPout = paste(mydata["line"], mydata["mating_status"], mydata["fusion_id"], paste(round(outTable, 3), collapse=","), sep=",")
    cat(SNPout, file=fileout, append=TRUE, sep="\n")
}
