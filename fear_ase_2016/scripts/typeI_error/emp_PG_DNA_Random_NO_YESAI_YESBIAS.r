# Get command line args with filenames for processing
args <- commandArgs(TRUE)

# Import Function
source("/scratch/lfs/mcintyre/cegs_ase_paper/scripts/r_lib/Subroutines_model2_experimental.R")#Subrutines in Rita s paper
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
headers_out <- paste("FUSION_ID", "S_RNA","TotalRNA","PER_S_RNA","q_simulation",
                  "mean_PG_q_true",  "q025_PG_q_true","q975_PG_q_true","Bayesianpvalue_PG_q_true","flag_AI_PG_q_true",
                  "mean_PG_q4", "q025_PG_q4","q975_PG_q4","Bayesianpvalue_PG_q4","flag_AI_PG_q4",
                  "mean_PG_q5", "q025_PG_q5","q975_PG_q5","Bayesianpvalue_PG_q5","flag_AI_PG_q5",
                  "mean_PG_q6", "q025_PG_q6","q975_PG_q6","Bayesianpvalue_PG_q6","flag_AI_PG_q6",
                  "estimate_binomial_test","q025_binomial_test","q975_binomial_test","binomial_test_pvalue","flag_AI_binomial_test", sep=",")
cat(headers_out,file=fileout,append=FALSE,sep="\n")

# Create output holders
means = qs025 = qs975 = flag_AI = rep(NA,5)
names(qs025) = names(qs975) = names(means) = names(flag_AI) = c("PG_q_true","PG_q4","PG_q5","PG_q6","Confidence")
Bayesianpvalue <- rep(NA,4)
names(Bayesianpvalue) <- c("PG_q_true","PG_q5","PG_q5","PG_q6") 

# Iterate over rest of the file
while(length(newline) != 0 ){
    newline <- readLines(con,n=1) 
    if (length(newline) == 0){break}
    mydata <- as.vector(strsplit(newline,split=",")[[1]])
    names(mydata) <- headers_in
       
    print(paste("------------Analyzing fusion--------------",m));m=m+1
    x <- as.numeric(mydata[c("RNA_Mel_all_Rep1","RNA_Mel_all_Rep2","RNA_Mel_all_Rep3")])
    x2 <- as.numeric(mydata[c("RNA_Sim_all_Rep1","RNA_Sim_all_Rep2","RNA_Sim_all_Rep3")])
    n_i <- x+x2
    q_true <- as.numeric(mydata["q_true"])

    # Binomail Test
    binomial_test <- binom.test(x=sum(n_i-x), n=sum(n_i), p=0.5, conf.level=0.95)
    binomial_test_pvalue <- binomial_test$p.value
    qs025["Confidence"] <- binomial_test$conf.int[1]
    qs975["Confidence"] <- binomial_test$conf.int[2]
    means["Confidence"] <- binomial_test$estimate
    flag_AI["Confidence"] <- ifelse(binomial_test_pvalue<= 0.05,1,0)

    # Poisson gamma with q=0.4
    tem_PGnaive <- gibbs_poissongamma(nsim=10000,nburnin=100,lag=5,
                                   x=as.numeric(x),y=as.numeric(x2),both=c(0,0,0),
                                   a_mu=.1,
                                   b_mu=.1,
                                   a_alpha=50,
                                   b_alpha=50,
                                   a_beta=1/2,
                                   b_beta=1/2,
                                   q_=0.4,abundance=FALSE)

    pps <- tem_PGnaive$alphas/(1+tem_PGnaive$alphas)

    quants <- quantile(pps,c(.025,0.975),na.rm=TRUE)
    qs025["PG_q4"] <- quants[1]
    qs975["PG_q4"] <- quants[2]
    means["PG_q4"] <- mean(pps)
    Bayesianpvalue["PG_q4"] <- 2*min(mean(pps<1/2),mean(pps>1/2))
    flag_AI["PG_q4"] <- ifelse(Bayesianpvalue["PG_q4"]<0.05,1,0)

    # Poisson gamma with q=1/2
    tem_PGnaive <- gibbs_poissongamma(nsim=10000,nburnin=100,lag=5,
                                   x=as.numeric(x),y=as.numeric(x2),both=c(0,0,0),
                                   a_mu=.1,
                                   b_mu=.1,
                                   a_alpha=50,
                                   b_alpha=50,
                                   a_beta=1/2,
                                   b_beta=1/2,
                                   q_=0.5,abundance=FALSE)

    pps <- tem_PGnaive$alphas/(1+tem_PGnaive$alphas)

    quants <- quantile(pps,c(.025,0.975),na.rm=TRUE)
    qs025["PG_q5"] <- quants[1]
    qs975["PG_q5"] <- quants[2]
    means["PG_q5"] <- mean(pps)
    Bayesianpvalue["PG_q5"] <- 2*min(mean(pps<1/2),mean(pps>1/2))
    flag_AI["PG_q5"] <- ifelse(Bayesianpvalue["PG_q5"]<0.05,1,0)

    # Poisson gamma with q=0.6
    tem_PGnaive <- gibbs_poissongamma(nsim=10000,nburnin=100,lag=5,
                                   x=as.numeric(x),y=as.numeric(x2),both=c(0,0,0),
                                   a_mu=.1,
                                   b_mu=.1,
                                   a_alpha=50,
                                   b_alpha=50,
                                   a_beta=1/2,
                                   b_beta=1/2,
                                   q_=0.6,abundance=FALSE)

    pps <- tem_PGnaive$alphas/(1+tem_PGnaive$alphas)

    quants <- quantile(pps,c(.025,0.975),na.rm=TRUE)
    qs025["PG_q6"] <- quants[1]
    qs975["PG_q6"] <- quants[2]
    means["PG_q6"] <- mean(pps)
    Bayesianpvalue["PG_q6"] <- 2*min(mean(pps<1/2),mean(pps>1/2))
    flag_AI["PG_q6"] <- ifelse(Bayesianpvalue["PG_q6"]<0.05,1,0)

    # model 2 Poi5amma with q_true
    tem_PG <- gibbs_poissongamma(nsim=10000,nburnin=100,lag=5,
                              x=as.numeric(x),y=as.numeric(x2),both=c(0,0,0),
                              a_mu= .1,
                              b_mu= .1,
                              a_alpha=50,
                              b_alpha=50,
                              a_beta= 1/2,
                              b_beta= 1/2,
                              q_=  q_true,abundance=FALSE)
    pps <-   tem_PG$alphas/(1+tem_PG$alphas)

    qs025["PG_q_true"] <-  quantile(pps,c(.025),na.rm=TRUE)
    qs975["PG_q_true"] <-  quantile(pps,c(.975),na.rm=TRUE)
    means["PG_q_true"] <-  mean(pps)
    Bayesianpvalue["PG_q_true"] <- 2*min(mean(pps<1/2),mean(pps>1/2))
    flag_AI["PG_q_true"] <- ifelse(Bayesianpvalue["PG_q_true"]<0.05,1,0)

    #################################################################
    #Confidence interval
    theta_star <-    1-sum(x)/sum(n_i)

    # Create Output
    SNPout <- paste(mydata["FUSION_ID"],
                    paste(round(c(sum(x2),sum(n_i),round(sum(x2)/sum(n_i),3),q_true,
                    means["PG_q_true"],qs025["PG_q_true"],qs975["PG_q_true"],Bayesianpvalue["PG_q_true"],flag_AI["PG_q_true"],
                    means["PG_q4"],qs025["PG_q4"],qs975["PG_q4"],Bayesianpvalue["PG_q4"],flag_AI["PG_q4"],
                    means["PG_q5"],qs025["PG_q5"],qs975["PG_q5"],Bayesianpvalue["PG_q5"],flag_AI["PG_q5"],
                    means["PG_q6"],qs025["PG_q6"],qs975["PG_q6"],Bayesianpvalue["PG_q6"],flag_AI["PG_q6"],
                    means["Confidence"],qs025["Confidence"],qs975["Confidence"],binomial_test_pvalue,flag_AI["Confidence"]
                    ),3),collapse=","),sep=",")

    print(SNPout)
    cat(SNPout,file=fileout,append=TRUE,sep="\n")
}
close(con)
