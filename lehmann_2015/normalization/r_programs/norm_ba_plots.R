library(ggplot2)
library(gridExtra)

# Findout where MCLAB is mounted at
mclab = Sys.getenv("MCLAB")

# Create a storage list to hold the plots
plots = list()

# For the mated dataset figure out all of the pairwise combinations
mydata <- read.csv('/home/jfear/tmp/mmethod.csv', header=TRUE)
mynames <- names(mydata)[3:dim(mydata)[2]]
mycomb = combn(mynames,2)

for (i in 1:dim(mycomb)[2]){
    comp = mycomb[,i]
    mydata$mean = rowMeans(cbind(mydata[,comp[1]],mydata[,comp[2]]))
    mydata$diff = mydata[,comp[1]] - mydata[,comp[2]]
    plots[[length(plots) + 1 ]] <- qplot(mean,diff,data=mydata,main=paste("Mated ",comp[1], " vs ", comp[2], sep=""))
    
}


# For the virgin dataset figure out all of the pairwise combinations
mydata <- read.csv('/home/jfear/tmp/mmethod.csv', header=TRUE)
mynames <- names(mydata)[3:dim(mydata)[2]]
mycomb = combn(mynames,2)

for (i in 1:dim(mycomb)[2]){
    comp = mycomb[,i]
    mydata$mean = rowMeans(cbind(mydata[,comp[1]],mydata[,comp[2]]))
    mydata$diff = mydata[,comp[1]] - mydata[,comp[2]]
    plots[[length(plots) + 1 ]] <- qplot(mean,diff,data=mydata,main=paste("Virgin ",comp[1], " vs ", comp[2], sep=""))
    
}


png(paste(mclab,"/cegs_sergey/r_plots/apn_M_ba_plots.png", sep=""))
grid.arrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], ncol=2, main=textGrob("&myline", gp=gpar(fontsize=20)))
dev.off()
