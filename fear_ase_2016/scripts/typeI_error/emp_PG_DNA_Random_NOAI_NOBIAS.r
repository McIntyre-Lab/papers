# Get command line args with filenames for processing
# args[1] - input file
# args[2] - output file

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
headers_out <- paste("FUSION_ID","S_RNA","TotalRNA","PROP_S_RNA","S_DNA","TotalDNA","PROP_S_DNA","q_true",
                     "estimate_binomial_test","q025_binomial_test","q975_binomial_test","binomial_test_pvalue","flag_AI_binomial_test",
                     "mean_NB_p_random_DNA","q025_NB_p_random_DNA","q975_NB_p_random_DNA","Bayesianpvalue_NB_p_random_DNA","flag_AI_NB_p_random_DNA",
                     "mean_PG_q_random_DNA","q025_PG_p_random_DNA","q975_PG_q_random_DNA","Bayesianpvalue_PG_q_random_DNA","flag_AI_PG_q_random_DNA",
                     "mean_PG_q4","q025_PG_q4","q975_PG_q4","Bayesianpvalue_PG_q4","flag_AI_PG_q4",
                     "mean_PG_q5","q025_PG_q5","q975_PG_q5","Bayesianpvalue_PG_q5","flag_AI_PG_q5",
                     "mean_PG_q6","q025_PG_q6","q975_PG_q6","Bayesianpvalue_PG_q6","flag_AI_PG_q6",
                     sep=",")

cat(headers_out,file=fileout,append=FALSE,sep="\n")

# Create output holders
ps <- rep(NA,2)

flag_AI = Bayesianpvalue = means = qs025 = qs975 = rep(NA,6)
names(flag_AI) = names(Bayesianpvalue) = names(means) = names(qs025) = names(qs975) = c("PG_q4","PG_q5","PG_q5","Confidence","PG_RandomDNA","NB")

while(length(newline) != 0 ){
    newline<-readLines(con,n=1)    
    if (length(newline) == 0){break}

    mydata <- as.vector(strsplit(newline,split=",")[[1]])
    names(mydata) <- headers_in
            
    print(paste("------------Analyzing fusion--------------",m));
    m <- m+1
    x <- as.numeric(mydata[c("line_RNA_Rep1","line_RNA_Rep2","line_RNA_Rep3")])
    x2 <- as.numeric(mydata[c("tester_RNA_Rep1","tester_RNA_Rep2","tester_RNA_Rep3")])
    n_i <- x+x2
    y <- as.numeric(mydata[c("line_DNA_Rep1","line_DNA_Rep2","line_DNA_Rep3")])
    y2 <- as.numeric(mydata[c("tester_DNA_Rep1","tester_DNA_Rep2","tester_DNA_Rep3")])
    m_i <- y+y2
    nRNA <- x+x2
    nDNA <- y+y2

    if(sum(x) > 0 and sum(x2) > 0){
        # Binomail Test
        binomial_test <- binom.test(x=sum(n_i-x), n=sum(n_i), p=0.5, conf.level=0.95)
        binomial_test_pvalue <- binomial_test$p.value
        qs025["Confidence"] <- binomial_test$conf.int[1]
        qs975["Confidence"] <- binomial_test$conf.int[2]
        means["Confidence"] <- binomial_test$estimate
        flag_AI["Confidence"] <- ifelse(binomial_test_pvalue<= 0.05,1,0)
        ps["p_sim"] <- 1/2

        # Negative Binomial
        NB <- gibbs_sampler(x=x,
                            n_i=n_i,
                            y=y,
                            m_i=m_i,
                            ##Hyperparameters
                            tt= sum(as.numeric(nRNA)), #prior variance of theta is p(1-p)/(t+1)
                            v=1, #prior variance of p is 1/(4(v+1))
                            sigma_MH=.04,  #Standard deviacion in the proposed distribution in the Metropolis Hasting algorithm when simulating p
                            sigma_theta=0.04,#Standard deviations of the sigma_thetas in the MH. I uses these only when t is large Parameers of the gamma distribution of lambda and lambda delta
                            a_lambda=0.5,
                            b_lambda=0.5,
                            a_delta=0.5,
                            b_delta=0.5,
                            ##Gibbs Sampler parameters
                            burn_in=1000,
                            storing_every=10,
                            m_simulations=1000,
                            figures=0
                           )

        qs025["NB"] <- quantile(NB$theta,c(.025),na.rm=TRUE)[1]
        qs975["NB"] <- quantile(NB$theta,c(.975),na.rm=TRUE)
        means["NB"] <- mean(NB$theta)
        Bayesianpvalue["NB"] <- 2*min(mean(NB$theta<1/2),mean(NB$theta>1/2))
        flag_AI["NB"] <- ifelse(Bayesianpvalue["NB"]<0.05,1,0)

        #Poisson Gamma with q random and obtained from the DNA
        PGDNA <- gibbs_poissongamma_q_random(nsim=1000,nburnin=100,lag=10,
                                                  x=as.numeric(x),y=as.numeric(x2),
                                                  both=c(0,0,0),
                                                  x_DNA=as.numeric(y2),y_DNA=as.numeric(y),
                                                  a_mu=0.1,
                                                  b_mu=0.1,
                                                  a_alpha=50,
                                                  b_alpha=50,
                                                  a_beta=1/2,
                                                  b_beta=1/2,
                                                  abundance=FALSE, 
                                                  sigma_MH=0.1, #SD of the normal dist to simulate from p in the MH step
                                                  v=1, #parameter of the beta prior distribution for p 
                                                  qs=NB$p
                                                 )
        pps <- PGDNA$alphas/(1+PGDNA$alphas)
        qs025["PG_RandomDNA"] <- quantile(pps,c(.025),na.rm=TRUE)
        qs975["PG_RandomDNA"] <- quantile(pps,c(.975),na.rm=TRUE)
        means["PG_RandomDNA"] <- mean(pps)
        Bayesianpvalue["PG_RandomDNA"] <- 2*min(mean(pps<1/2),mean(pps>1/2))
        flag_AI["PG_RandomDNA"] <- ifelse(Bayesianpvalue["PG_RandomDNA"]<0.05,1,0)

        # Poisson gamma with q=0.4
        PG4 <- gibbs_poissongamma(nsim=10000,nburnin=100,lag=10,
                                     x=as.numeric(x),y=as.numeric(x2),both=c(0,0,0),
                                     a_mu=0.1,
                                     b_mu=0.1,
                                     a_alpha=50,
                                     b_alpha=50,
                                     a_beta=1/2,
                                     b_beta=1/2,
                                     q_=0.4,        
                                     abundance=FALSE
                                    )

        pps <- PG4$alphas/(1+PG4$alphas)
        qs025["PG_q4"] <- quantile(pps,c(.025),na.rm=TRUE)
        qs975["PG_q4"] <- quantile(pps,c(.975),na.rm=TRUE)
        means["PG_q4"] <- mean(pps)
        Bayesianpvalue["PG_q4"] <- 2*min(mean(pps<1/2),mean(pps>1/2))
        flag_AI["PG_q4"] <- ifelse(Bayesianpvalue["PG_q4"]<0.05,1,0)

        # Poisson gamma with q=1/2
        PG5 <- gibbs_poissongamma(nsim=10000,nburnin=100,lag=10,
                                     x=as.numeric(x),y=as.numeric(x2),both=c(0,0,0),
                                     a_mu=0.1,
                                     b_mu=0.1,
                                     a_alpha=50,
                                     b_alpha=50,
                                     a_beta=1/2,
                                     b_beta=1/2,
                                     q_=1/2,        
                                     abundance=FALSE
                                    )

        pps <- PG5$alphas/(1+PG5$alphas)
        qs025["PG_q5"] <- quantile(pps,c(.025),na.rm=TRUE)
        qs975["PG_q5"] <- quantile(pps,c(.975),na.rm=TRUE)
        means["PG_q5"] <- mean(pps)
        Bayesianpvalue["PG_q5"] <- 2*min(mean(pps<1/2),mean(pps>1/2))
        flag_AI["PG_q5"] <- ifelse(Bayesianpvalue["PG_q5"]<0.05,1,0)

        # Poisson gamma with q=0.6
        PG6 <- gibbs_poissongamma(nsim=10000,nburnin=100,lag=10,
                                  x=as.numeric(x),y=as.numeric(x2),both=c(0,0,0),
                                  a_mu=0.1,
                                  b_mu=0.1,
                                  a_alpha=50,
                                  b_alpha=50,
                                  a_beta=1/2,
                                  b_beta=1/2,
                                  q_=0.6,        
                                  abundance=FALSE
                                 )

        pps <- PG6$alphas/(1+PG6$alphas)
        qs025["PG_q6"] <- quantile(pps,c(.025),na.rm=TRUE)
        qs975["PG_q6"] <- quantile(pps,c(.975),na.rm=TRUE)
        means["PG_q6"] <- mean(pps)
        Bayesianpvalue["PG_q6"] <- 2*min(mean(pps<1/2),mean(pps>1/2))
        flag_AI["PG_q6"] <- ifelse(Bayesianpvalue["PG_q6"]<0.05,1,0)

        # Create Output
        SNPout <- paste(mydata["FUSION_ID"],
                        paste(round(c(sum(x2),sum(n_i),round(sum(x2)/sum(n_i),3),sum(y2),sum(m_i),round(sum(y2)/sum(m_i),3),ps["p_sim"],
                        means["Confidence"],qs025["Confidence"],qs975["Confidence"],binomial_test_pvalue,flag_AI["Confidence"],
                        means["NB"],qs025["NB"],qs975["NB"],Bayesianpvalue["NB"],flag_AI["NB"],
                        means["PG_RandomDNA"],qs025["PG_RandomDNA"],qs975["PG_RandomDNA"],Bayesianpvalue["PG_RandomDNA"],flag_AI["PG_RandomDNA"],
                        means["PG_q4"],qs025["PG_q4"],qs975["PG_q4"],Bayesianpvalue["PG_q4"],flag_AI["PG_q4"],
                        means["PG_q5"],qs025["PG_q5"],qs975["PG_q5"],Bayesianpvalue["PG_q5"],flag_AI["PG_q5"],
                        means["PG_q6"],qs025["PG_q6"],qs975["PG_q6"],Bayesianpvalue["PG_q6"],flag_AI["PG_q6"]
                       ),3),collapse=","),sep=",")
        print(SNPout)
        cat(SNPout,file=fileout,append=TRUE,sep="\n")
    }
}
