library(ggplot2)

mclab <- Sys.getenv('MCLAB')
mydata <- read.csv('/home/jfear/tmp.csv')

mating_status <- mydata$mating_status[1]
mytype <- mydata$name[1]
mynames <- names(mydata)[4:dim(mydata)[2]-1]
mycomb <- combn(mynames,2)

for (i in 1:dim(mycomb)[2]){
# The above will plot all 4k+ graphs, below will plot just r101
#for (i in 1:92){
    comp <- mycomb[,i]
    mydata$mean <- rowMeans(cbind(mydata[,comp[1]],mydata[,comp[2]]))
    mydata$diff <- mydata[,comp[1]] - mydata[,comp[2]]
    png(paste(mclab,'/cegs_sergey/reports/line_normalization/line_ba_plots/bland_altman_',mytype,'_',mating_status,'_',comp[1],'_',comp[2],'.png',sep=""))
    p <- ggplot(mydata,aes(mean,diff)) + geom_point() + ggtitle(paste("Bland-Altman Plots",mytype,"\n",mating_status,comp[1],comp[2])) + ylim(-6,6)
    print(p)
    dev.off()
}
