# Get command line arguments
#args[1] = input CSV of concentration and RPKM
#args[2] = output PDF plots

args <- commandArgs(TRUE)
cname <- args[1]
pname <- args[2]

data=read.csv(cname)
attach(data)

pdf(pname)
par(mfrow=c(3,3))
for(i in unique(sample)){
  plot(data[sample==i,]$log_Mix1_adj_conc,data[sample==i, ]$log_rpkm,main=i,xlab="log_Mix1_adj_conc",ylab="log(RPKM)",cex.main=.8, xlim=c(0,30), ylim=c(0,30))
}
dev.off()

