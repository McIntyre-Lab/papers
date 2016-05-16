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
fname <- '/home/jfear/mclab/cegs_ase_paper/pipeline_output/typeI_error/input/onebigdataset_4luispaper.csv'

####################################################################
##"OBTAINING" THE SIMULATION TRUE PARAMETERS FROM THE REAL DATASET
set.seed(0)
n.genes <- 10^4        #Number of simulated exons
mydatatem <- read.csv(fname)
mydatatem <- mydatatem[which(mydatatem$in_final_MBE==1),]
mydata <- mydatatem[sort(sample(1:nrow(mydatatem),n.genes)),]
rm(mydatatem)

countsTable <- mydata[c( "RNA_Mel_all_Rep1","RNA_Mel_all_Rep2","RNA_Mel_all_Rep3","RNA_Sim_all_Rep1","RNA_Sim_all_Rep2","RNA_Sim_all_Rep3" )]
vonds <- c( "M", "M", "M", "S", "S", "S" ) 

meanM <- apply(countsTable[,c("RNA_Mel_all_Rep1","RNA_Mel_all_Rep2","RNA_Mel_all_Rep3")],1,mean)
meanS <- apply(countsTable[,c("RNA_Sim_all_Rep1","RNA_Sim_all_Rep2","RNA_Sim_all_Rep3")],1,mean)
##########################################################################
##Simulating data according to the different scenarios
scenario <-     "NOAI_NOBIAS"        #Choose one of these three scenarios
                #"NOAI_YESBIAS"
                #"YESAI_YESBIAS"
                
out <- data.frame(FUSION_ID=mydata$FUSION_ID)                    

summeans <- meanM+meanS
q_truefirst5000 <- 0.45            #Choose this value to select the bias

if(scenario=="NOAI_NOBIAS"){R=1;out$q_true=1/2}
if(scenario=="NOAI_YESBIAS"){R=1;out$q_true=c(rep(q_truefirst5000,5000),rep(1-q_truefirst5000,5000))}
if(scenario=="YESAI_YESBIAS"){R=1.5;out$q_true=c(rep(q_truefirst5000,5000),rep(1-q_truefirst5000,5000))}
fortitle <- ifelse(R==1,R,paste(floor(R),"d",(R-floor(R))*10,sep=""))

out$trueoverdispersion <- 0#(overdispersionM+overdispersionS)/2);
b <- 1/((R-1)*(1-out$q_true)+1)
a <- R*b

out$truemeanM <- a*summeans
out$truemeanS <- b*summeans


out$RNA_Sim_all_Rep3 <- out$RNA_Sim_all_Rep2 <- out$RNA_Sim_all_Rep1 <- out$RNA_Mel_all_Rep3 <- out$RNA_Mel_all_Rep2 <- out$RNA_Mel_all_Rep1 <- rep(NA,nrow(out))
names(out)

for(i in 1:nrow(out)){
    out[i,c("RNA_Mel_all_Rep1","RNA_Mel_all_Rep2","RNA_Mel_all_Rep3")] <- rpois(3,(1-out$q_true)[i]*out$truemeanM[i])

    out[i,c("RNA_Sim_all_Rep1","RNA_Sim_all_Rep2","RNA_Sim_all_Rep3")] <- rpois(3,out$q_true[i]*out$truemeanS[i])
}
par(mfrow=c(2,1))
sumM <- rowSums(out[,c("RNA_Mel_all_Rep1","RNA_Mel_all_Rep2","RNA_Mel_all_Rep3")])
sumS <- rowSums(out[,c("RNA_Sim_all_Rep1","RNA_Sim_all_Rep2","RNA_Sim_all_Rep3")])
plot((sumM/(1-out$q_true))/(sumS/(out$q_true)),ylab="Mel/Sim after bias correction");abline(h=R,col="red",lwd=3)
plot(log10(summeans),log10((sumM/3)+(sumS/3)),ylab="sum of simulated sample",xlab="Real data set sum of sample means");abline(a=0,b=1,col="red",lwd=3)

if(scenario != "NOAI_NOBIAS"){
    write.csv(out,file=paste("./simulated_",scenario,"Requal",fortitle,"Poisson_q_truefirst5000_0d",q_truefirst5000*100,".csv",sep=""),quote=FALSE,row.names = FALSE)
} else{
    # If scenario == "NOAI_NOBIAS" we also generate DNA counts
    countsTable_DNA=mydata[c("DNA_Mel_all_Rep1","DNA_Mel_all_Rep2","DNA_Mel_all_Rep3","DNA_Sim_all_Rep1","DNA_Sim_all_Rep2","DNA_Sim_all_Rep3" )]

    # Luis was originally using these sum statments, but I think it should have
    # been a mean statment so I went a head and changed it.
    #sumM_DNA <- rowSums(countsTable_DNA[,c("DNA_Mel_all_Rep1","DNA_Mel_all_Rep2","DNA_Mel_all_Rep3")])
    #sumS_DNA <- rowSums(countsTable_DNA[,c("DNA_Sim_all_Rep1","DNA_Sim_all_Rep2","DNA_Sim_all_Rep3")])

    sumM_DNA <- apply(countsTable_DNA[,c("DNA_Mel_all_Rep1","DNA_Mel_all_Rep2","DNA_Mel_all_Rep3")],1, mean)
    sumS_DNA <- apply(countsTable_DNA[,c("DNA_Sim_all_Rep1","DNA_Sim_all_Rep2","DNA_Sim_all_Rep3")],1,mean)

    summeans_DNA <- sumM_DNA+sumS_DNA
    out_DNA <- data.frame(trueoverdispersion_DNA=rep(0,nrow(out)))#(overdispersionM_DNA+overdispersionS_DNA)/2)
    
    # Simulating DNA with no systematic bias (q_true=1/2) and no AI (R=1) 
    q_true_DNA <- 1/2
    out_DNA$q_true_DNA <- q_true_DNA
    R_DNA <- 1
    b_DNA <- 1/((R_DNA-1)*q_true_DNA+1)
    a_DNA <- R_DNA*b_DNA

    out_DNA$truemeanM_DNA <- a_DNA*summeans_DNA
    out_DNA$truemeanS_DNA <- b_DNA*summeans_DNA

    out_DNA$DNA_Sim_all_Rep3 <- out_DNA$DNA_Sim_all_Rep2 <- out_DNA$DNA_Sim_all_Rep1 <- out_DNA$DNA_Mel_all_Rep3 <- out_DNA$DNA_Mel_all_Rep2 <- out_DNA$DNA_Mel_all_Rep1 <- rep(NA,nrow(out))
    head(out_DNA)

    for(i in 1:nrow(out)){
        out_DNA[i,c("DNA_Mel_all_Rep1","DNA_Mel_all_Rep2","DNA_Mel_all_Rep3")] <- rpois(3,(1-out_DNA$q_true_DNA[i])*out_DNA$truemeanM_DNA[i])
        
        out_DNA[i,c("DNA_Sim_all_Rep1","DNA_Sim_all_Rep2","DNA_Sim_all_Rep3")] <- rpois(3,out_DNA$q_true_DNA[i]*out_DNA$truemeanS_DNA[i])
    }
        
    par(mfrow=c(2,1))
    sumM_DNA <- rowSums(out_DNA[,c("DNA_Mel_all_Rep1","DNA_Mel_all_Rep2","DNA_Mel_all_Rep3")])
    sumS_DNA <- rowSums(out_DNA[,c("DNA_Sim_all_Rep1","DNA_Sim_all_Rep2","DNA_Sim_all_Rep3")])
    plot((sumM_DNA/out_DNA$q_true_DNA)/(sumS_DNA/(1-out_DNA$q_true_DNA)),ylab="DNA Mel/Sim after bias correction");abline(h=R,col="red",lwd=3)
    plot(log10(summeans_DNA),log10((sumM_DNA/3)+(sumS_DNA/3)),ylab="DNA sum of simulated sample",xlab="Real data set sum of sample means");abline(a=0,b=1,col="red",lwd=3)
    head(cbind(out,out_DNA))
    write.csv(cbind(out,out_DNA),file=paste("./simulated_",scenario,"Requal",R,"Poisson_q_truefirst5000_0d",q_truefirst5000*100,".csv",sep=""),quote=FALSE,row.names = FALSE)
}
