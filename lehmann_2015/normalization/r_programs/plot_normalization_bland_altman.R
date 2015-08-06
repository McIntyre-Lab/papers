library(ggplot2)

mclab <- Sys.getenv('MCLAB')
mydata <- read.csv('/home/jfear/tmp.csv')

line <- mydata$line[1]
mating_status <- mydata$mating_status[1]
mytype <- mydata$name[1]
mynames <- names(mydata)[5:dim(mydata)[2]]
mycomb <- combn(mynames,2)

for (i in 1:dim(mycomb)[2]){
    comp <- mycomb[,i]
    mydata$mean <- rowMeans(cbind(mydata[,comp[1]],mydata[,comp[2]]))
    mydata$diff <- mydata[,comp[1]] - mydata[,comp[2]]
    png(paste(mclab,'/cegs_sergey/reports/line_normalization/ba_plots/bland_altman_',mytype,'_',line,'_',mating_status,'_',comp[1],'_',comp[2],'.png',sep=""))
    p <- ggplot(mydata,aes(mean,diff)) + geom_point() + ggtitle(paste("Bland-Altman Plots",mytype,"\n",line,mating_status,comp[1],comp[2]))
    print(p)
    dev.off()
}
