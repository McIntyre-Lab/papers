#Program used to analyze the simulated data sets 
#generated by R program simulatingdata.R
#explained in Synthetic Data Set Plan_LLN08072014.docx
#For two scenarios 
#"NOAI_YESBIAS"
#"YESAI_YESBIAS"
#scenario=    "NOAI"

# Get command line args with filenames for processing
# args[1] = input
# args[2] = output
# args[3] = qTrue
args <- commandArgs(TRUE)

# Import Function
source("/scratch/lfs/mcintyre/cegs_ase_paper/scripts/r_lib/Subroutines_model2_experimental.R")#Subrutines from MBE paper
source("/scratch/lfs/mcintyre/cegs_ase_paper/scripts/r_lib/AI_poissongamma_functions.R")#Subrutines for Poissong Gamma

# Input File 
con <- file(args[1],"r")
newline <- readLines(con,n=1) #moving the pointer away from the headers
headers_in <- strsplit(newline,split=",")[[1]]
mydata <- rep(NA,length(headers_in))
names(mydata) <- headers_in
m <- 1

# Output File
fileout <- args[2]
headers_out <- paste("FUSION_ID","S_RNA","TotalRNA","PER_S_RNA",
                     "estimate_binomial_test","q025_binomial_test","q975_binomial_test","binomial_test_pvalue","flag_AI_binomial_test",
                     "q_true","mean_q_true","q025_q_true","q975_q_true","Bayesianpvalue_q_true","flag_AI_q_true",
                     "mean_q4","q025_q4","q975_q4","Bayesianpvalue_q4","flag_AI_q4",
                     "mean_q5","q025_q5","q975_q5","Bayesianpvalue_q5","flag_AI_q5",
                     "mean_q6","q025_q6","q975_q6","Bayesianpvalue_q6","flag_AI_q6",
                     "q_sim_below_01per","mean_q_sim_below_01per","q025_q_sim_below_01per","q975_q_sim_below_01per","Bayesianpvalue_q_sim_below_01per","flag_AI_q_sim_below_01per",
                     "q_sim_above_01per","mean_q_sim_above_01per","q025_q_sim_above_01per","q975_q_sim_above_01per","Bayesianpvalue_q_sim_above_01per","flag_AI_q_sim_above_01per",
                     "q_sim_below_02per","mean_q_sim_below_02per","q025_q_sim_below_02per","q975_q_sim_below_02per","Bayesianpvalue_q_sim_below_02per","flag_AI_q_sim_below_02per",
                     "q_sim_above_02per","mean_q_sim_above_02per","q025_q_sim_above_02per","q975_q_sim_above_02per","Bayesianpvalue_q_sim_above_02per","flag_AI_q_sim_above_02per",
                     "q_sim_below_05per","mean_q_sim_below_05per","q025_q_sim_below_05per","q975_q_sim_below_05per","Bayesianpvalue_q_sim_below_05per","flag_AI_q_sim_below_05per",
                     "q_sim_above_05per","mean_q_sim_above_05per","q025_q_sim_above_05per","q975_q_sim_above_05per","Bayesianpvalue_q_sim_above_05per","flag_AI_q_sim_above_05per",
                     "q_sim_below_10per","mean_q_sim_below_10per","q025_q_sim_below_10per","q975_q_sim_below_10per","Bayesianpvalue_q_sim_below_10per","flag_AI_q_sim_below_10per",
                     "q_sim_above_10per","mean_q_sim_above_10per","q025_q_sim_above_10per","q975_q_sim_above_10per","Bayesianpvalue_q_sim_above_10per","flag_AI_q_sim_above_10per",
                     "q_sim_below_20per","mean_q_sim_below_20per","q025_q_sim_below_20per","q975_q_sim_below_20per","Bayesianpvalue_q_sim_below_20per","flag_AI_q_sim_below_20per",
                     "q_sim_above_20per","mean_q_sim_above_20per","q025_q_sim_above_20per","q975_q_sim_above_20per","Bayesianpvalue_q_sim_above_20per","flag_AI_q_sim_above_20per",
                     sep=",")

cat(headers_out,file=fileout,append=FALSE,sep="\n")

# Get Bias parameter
#q_TRUE <- 0.35    #c(0.35,0.375,0.4,0.425,0.45,0.46,0.47,0.48,0.49,0.50,0.51,0.52,0.53,0.54,0.55,0.575,0.6,0.625,0.65)
q_TRUE <- as.numeric(args[3])
            
# Get misspecification values
percentages <- c(1,2,5,10,20)
qsims <- q_TRUE*(1+c(-1,1)*rep(percentages,each=2)/100)
q_sims_fornames <- paste(c("below","above"),rep(percentages,each=2),"per",sep="")

# Set defaults
means <- qs025 <- qs975 <- flag_AI <- rep(NA,15)
names(qs025) <- names(qs975) <- names(means) <- names(flag_AI) <- c("q_true","q4","q5","q6","Confidence",q_sims_fornames)
Bayesianpvalue <- rep(NA,14)
names(Bayesianpvalue) <- c("q_true","q4","q5","q6",q_sims_fornames)    

while(length(newline) != 0 ){
    newline <- readLines(con,n=1)    
    if (length(newline) == 0){break}
    mydata <- as.vector(strsplit(newline,split=",")[[1]])
    names(mydata) <- headers_in
    print(paste("------------Analyzing fusion--------------",m));m=m+1
    x <- as.numeric(mydata[c("line_RNA_Rep1","line_RNA_Rep2","line_RNA_Rep3")])
    x2 <- as.numeric(mydata[c("tester_RNA_Rep1","tester_RNA_Rep2","tester_RNA_Rep3")])
    q_true <- as.numeric(mydata["qTrue"])

    # Make sure there are reads for both
    if(sum(x) > 0 && sum(x2) > 0){
        # Confidence interval
        n_i <- x+x2
        theta_star <- 1-sum(x)/sum(n_i)
        binomial_test <- binom.test(x=sum(n_i-x), n=sum(n_i), p=0.5, conf.level=0.95)
        binomial_test_pvalue <- binomial_test$p.value
        qs025["Confidence"] <- binomial_test$conf.int[1]
        qs975["Confidence"] <- binomial_test$conf.int[2]
        means["Confidence"] <- binomial_test$estimate
        flag_AI["Confidence"] <- ifelse(binomial_test_pvalue<= 0.05,1,0)

        # Poisson gamma with q_true
        tem_PG <- gibbs_poissongamma(nsim=10000,nburnin=100,lag=5,
                                     x=as.numeric(x),y=as.numeric(x2),both=c(0,0,0),
                                     a_mu=.1,
                                     b_mu=.1,
                                     a_alpha=50,
                                     b_alpha=50,
                                     a_beta=1/2,
                                     b_beta=1/2,
                                     q_=q_true,abundance=FALSE
                                    )

        pps <- tem_PG$alphas/(1+tem_PG$alphas)
        qs025["q_true"] <- quantile(pps,c(.025),na.rm=TRUE)
        qs975["q_true"] <- quantile(pps,c(.975),na.rm=TRUE)
        means["q_true"] <- mean(pps)
        Bayesianpvalue["q_true"] <- 2*min(mean(pps<1/2),mean(pps>1/2))
        flag_AI["q_true"] <- ifelse(Bayesianpvalue["q_true"]<0.05,1,0)

        #Poisson Gamma naive i.e. q = 0.4
        tem_PGnaive <- gibbs_poissongamma(nsim=10000,nburnin=100,lag=5,
                                          x=as.numeric(x),y=as.numeric(x2),both=c(0,0,0),
                                          a_mu=.1,
                                          b_mu=.1,
                                          a_alpha=50,
                                          b_alpha=50,
                                          a_beta=1/2,
                                          b_beta=1/2,
                                          q_=0.4,abundance=FALSE
                                         )

        pps <- tem_PGnaive$alphas/(1+tem_PGnaive$alphas)
        quants <- quantile(pps,c(.025,0.975),na.rm=TRUE)
        qs025["q4"] <- quants[1]
        qs975["q4"] <- quants[2]
        means["q4"] <- mean(pps)
        Bayesianpvalue["q4"] <- 2*min(mean(pps<1/2),mean(pps>1/2))
        flag_AI["q4"] <- ifelse(Bayesianpvalue["q4"]<0.05,1,0)

        #Poisson Gamma naive i.e. q = 1/2
        tem_PGnaive <- gibbs_poissongamma(nsim=10000,nburnin=100,lag=5,
                                          x=as.numeric(x),y=as.numeric(x2),both=c(0,0,0),
                                          a_mu=.1,
                                          b_mu=.1,
                                          a_alpha=50,
                                          b_alpha=50,
                                          a_beta=1/2,
                                          b_beta=1/2,
                                          q_=0.5,abundance=FALSE
                                         )

        pps <- tem_PGnaive$alphas/(1+tem_PGnaive$alphas)
        quants <- quantile(pps,c(.025,0.975),na.rm=TRUE)
        qs025["q5"] <- quants[1]
        qs975["q5"] <- quants[2]
        means["q5"] <- mean(pps)
        Bayesianpvalue["q5"] <- 2*min(mean(pps<1/2),mean(pps>1/2))
        flag_AI["q5"] <- ifelse(Bayesianpvalue["q5"]<0.05,1,0)

        #Poisson Gamma naive i.e. q = 0.6
        tem_PGnaive <- gibbs_poissongamma(nsim=10000,nburnin=100,lag=5,
                                          x=as.numeric(x),y=as.numeric(x2),both=c(0,0,0),
                                          a_mu=.1,
                                          b_mu=.1,
                                          a_alpha=50,
                                          b_alpha=50,
                                          a_beta=1/2,
                                          b_beta=1/2,
                                          q_=0.6,abundance=FALSE
                                         )

        pps <- tem_PGnaive$alphas/(1+tem_PGnaive$alphas)
        quants <- quantile(pps,c(.025,0.975),na.rm=TRUE)
        qs025["q6"] <- quants[1]
        qs975["q6"] <- quants[2]
        means["q6"] <- mean(pps)
        Bayesianpvalue["q6"] <- 2*min(mean(pps<1/2),mean(pps>1/2))
        flag_AI["q6"] <- ifelse(Bayesianpvalue["q6"]<0.05,1,0)

        #Poisson gamma with different q_sims
        forout=c()
        for(ii in 1:length(qsims)){
            tem_PG <- gibbs_poissongamma(nsim=10000,nburnin=100,lag=5,
                                         x=as.numeric(x),y=as.numeric(x2),both=c(0,0,0),
                                         a_mu=.1,
                                         b_mu=.1,
                                         a_alpha=50,
                                         b_alpha=50,
                                         a_beta=1/2,
                                         b_beta=1/2,
                                         q_=qsims[ii],abundance=FALSE
                                        )

            pps <- tem_PG$alphas/(1+tem_PG$alphas)
            qlabel <- q_sims_fornames[ii]
            quants <- quantile(pps,c(.025,0.975),na.rm=TRUE)
            qs025[qlabel] <- quants[1]
            qs975[qlabel] <- quants[2]
            means[qlabel] <- mean(pps)
            Bayesianpvalue[qlabel] <- 2*min(mean(pps<1/2),mean(pps>1/2))
            flag_AI[qlabel] <- ifelse(Bayesianpvalue[qlabel]<0.05,1,0)
            forout <- c(forout,qsims[ii],means[qlabel],qs025[qlabel],qs975[qlabel],Bayesianpvalue[qlabel],flag_AI[qlabel])
        }

        #Saving results
        SNPout=paste(mydata["fusion_id"],paste(round(c(sum(x2),sum(n_i),round(sum(x2)/sum(n_i),3),
                                                       means["Confidence"],qs025["Confidence"],qs975["Confidence"],
                                                       binomial_test_pvalue,flag_AI["Confidence"],q_true,means["q_true"],
                                                       qs025["q_true"],qs975["q_true"],Bayesianpvalue["q_true"],flag_AI["q_true"],
                                                       means["q4"],qs025["q4"],qs975["q4"],Bayesianpvalue["q4"],flag_AI["q4"],
                                                       means["q5"],qs025["q5"],qs975["q5"],Bayesianpvalue["q5"],flag_AI["q5"],
                                                       means["q6"],qs025["q6"],qs975["q6"],Bayesianpvalue["q6"],flag_AI["q6"],
                                                       forout),3),collapse=","),sep=",")
    }
    else {
        SNPout=paste(mydata["fusion_id"], paste(rep(NA, 8), collapse=','), 
                     mydata["qTrue"], paste(rep(NA, 80), collapse=','),sep=",")
    }
    print(SNPout)
    cat(SNPout,file=fileout,append=TRUE,sep="\n")
}
close(con)
