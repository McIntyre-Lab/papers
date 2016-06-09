#Simulating dataset according to 
#Synthetic Data Set Plan_RMG07302014_v2 proposed by Rita August 2014
#This is a version of the program simulatingdata
#but without using deseq to compute the sample means of the real dataset
#I updaed this because when I tried to run simulatingdata.R for Justin in Nov 17, 2014
#The version of DESeq in my computer had been updated and functions on the original program
#were obsolete
#Since at the end, we sample from Poisson we do not require DESeq to estimate the 
#overdispersion parameters
 
setwd('/home/jfear/mclab/cegs_ase_paper/pipeline_output/typeI_error/output')
fname <- '/home/jfear/mclab/cegs_ase_paper/pipeline_output/typeI_error/input/r101_RNA_sim_DNA_cnts.csv'

####################################################################
##"OBTAINING" THE SIMULATION TRUE PARAMETERS FROM THE REAL DATASET
set.seed(0)
n.genes <- 10^4        #Number of simulated exons
mydatatem <- read.csv(fname)
mydata <- mydatatem[sort(sample(1:nrow(mydatatem),n.genes)),]
rm(mydatatem)

countsTable <- mydata[c( "line_RNA_Rep1","line_RNA_Rep2","line_RNA_Rep3","tester_RNA_Rep1","tester_RNA_Rep2","tester_RNA_Rep3" )]
meanM <- apply(countsTable[,c("line_RNA_Rep1","line_RNA_Rep2","line_RNA_Rep3")],1,mean)
meanS <- apply(countsTable[,c("tester_RNA_Rep1","tester_RNA_Rep2","tester_RNA_Rep3")],1,mean)
summeans <- meanM+meanS
##########################################################################

# Initialize output
out <- data.frame(FUSION_ID=mydata$fusion_id)                    

##Simulating data according to the different values of bias
qs <- c(0.35, 0.375, 0.4, 0.425, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.575, 0.6, 0.625, 0.65)

for(j in 1:length(qs)){
    q_true <- qs[j]            #Choose this value to select the bias

    R=1
    out$q_true=c(rep(q_true,10000))
    b <- 1/((R-1)*(1-out$q_true)+1)
    a <- R*b

    out$truemeanM <- a*summeans
    out$truemeanS <- b*summeans

    out$tester_RNA_Rep3 <- out$tester_RNA_Rep2 <- out$tester_RNA_Rep1 <- out$line_RNA_Rep3 <- out$line_RNA_Rep2 <- out$line_RNA_Rep1 <- rep(NA,nrow(out))
    names(out)

    for(i in 1:nrow(out)){
        out[i,c("line_RNA_Rep1","line_RNA_Rep2","line_RNA_Rep3")] <- rpois(3,(1-out$q_true)[i]*out$truemeanM[i])
        out[i,c("tester_RNA_Rep1","tester_RNA_Rep2","tester_RNA_Rep3")] <- rpois(3,out$q_true[i]*out$truemeanS[i])
    }

    par(mfrow=c(2,1))
    sumM <- rowSums(out[,c("line_RNA_Rep1","line_RNA_Rep2","line_RNA_Rep3")])
    sumS <- rowSums(out[,c("tester_RNA_Rep1","tester_RNA_Rep2","tester_RNA_Rep3")])
    plot((sumM/(1-out$q_true))/(sumS/(out$q_true)),ylab="Mel/Sim after bias correction");abline(h=R,col="red",lwd=3)
    plot(log10(summeans),log10((sumM/3)+(sumS/3)),ylab="sum of simulated sample",xlab="Real data set sum of sample means");abline(a=0,b=1,col="red",lwd=3)

    write.csv(out,file=paste("./r101_simulateddata_forTIER_FigureNOAIPoisson_q_true_0d",q_true*1000,".csv",sep=""),quote=FALSE,row.names = FALSE)
}
