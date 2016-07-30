# Plotting the smoothed curves of pi, theta, and D on one plot per chromosome
plot(x=chrX$BIN_START, y=chrX$PI, type="n", main="X" ,  col="#003399",ylim=c(0,0.02), xlim=c(0,3e+07))
with(chrX, lines(loess.smooth(chrX$BIN_START, chrX$PI, span=0.09), col = "blue"))


#PI
pi10 <- read.table("/ufgi/SHARE/Mcintyre_Lab/ethanol/Sim_Pop_Gen/output/nuc_div/pi_10kb_genome_filter_nolab.windowed.pi", header=TRUE)

#Subset the 10kb data
chr2L_p <- subset(pi10, pi10$CHROM=="2L")
chr2R_p <- subset(pi10, pi10$CHROM=="2R")
chr3L_p <- subset(pi10, pi10$CHROM=="3L")
chr3R_p <- subset(pi10, pi10$CHROM=="3R")
chr4_p <- subset(pi10, pi10$CHROM=="4")
chrX_p <- subset(pi10, pi10$CHROM=="X")

#THETA
thet<- read.table('/ufgi/SHARE/Mcintyre_Lab/ethanol/Sim_Pop_Gen/output/snp_density/snpden_10kb_genome_filter_nolab_theta_4-20.tsv', header=FALSE)

chr2L_t <- subset(thet, thet$V1=='2L')
chr2R_t <- subset(thet, thet$V1=='2R')
chr3L_t <- subset(thet, thet$V1=='3L')
chr3R_t <- subset(thet, thet$V1=='3R')
chr4_t <- subset(thet, thet$V1=='4')
chrX_t <- subset(thet, thet$V1=='X')

#TAJIMA'S D
td10kb <- read.table('/ufgi/SHARE/Mcintyre_Lab/ethanol/Sim_Pop_Gen/output/tajima_d/filter_nolab_tajd_10kb.Tajima.D', header=TRUE)

# Plot the 10kb window TsD
chr2L_d <- subset(td10kb, td10kb$CHROM=='2L')
chr2R_d <- subset(td10kb, td10kb$CHROM=='2R')
chr3L_d <- subset(td10kb, td10kb$CHROM=='3L')
chr3R_d <- subset(td10kb, td10kb$CHROM=='3R')
chr4_d <- subset(td10kb, td10kb$CHROM=='4')
chrX_d <- subset(td10kb, td10kb$CHROM=='X')



colfunc <- colorRampPalette(c("#C90072", "white"))
colfunc(10)
plot(rep(1,10),col=colfunc(10),pch=19,cex=3)

dev.off


par(mfrow=c(2,3))
  

#ChrX
plot(x=chrX_p$BIN_START, y=chrX_p$PI, type="n", main="X" ,  col="#003399",ylim=c(0,0.024), xlim=c(0,3e+07), xlab=NA, ylab="PI and Theta")
  #axis(side=2, labels=TRUE)
  #axis(side=1, labels=TRUE)
  #axis(side=3, labels=FALSE)
  with(chrX_p, lines(loess.smooth(chrX_p$BIN_START, chrX_p$PI, span=0.005), lwd=4, col="#0061A1"))
  with(chrX_t, lines(loess.smooth(chrX_t$V2, chrX_t$V5, span=0.005), lwd=4, col="#AA8E39"))
par(new=T)  
plot(x=chrX_d$BIN_START, y=chrX_d$TajimaD ,  type="n",  col="#C90072", xlab=NA, ylab=NA, axes=F,   ylim=c(-4,4),  xlim=c(0,3e+07))
  with(chrX_d, lines(loess.smooth(chrX_d$BIN_START, chrX_d$TajimaD, span=0.005),  lwd=4,  col="#C90072"))
  axis(side=4)


# Chr2L


plot(x=chr2L_p$BIN_START, y=chr2L_p$PI, type="n", main="2L" ,  col="#003399",ylim=c(0,0.024), xlim=c(0,3e+07),xlab=NA, ylab="PI and Theta")
  #axis(side=2, labels=TRUE)
  #axis(side=1, labels=TRUE)
  with(chr2L_p, lines(loess.smooth(chr2L_p$BIN_START, chr2L_p$PI, span=0.005), lwd=4, col="#0061A1"))
  with(chr2L_t, lines(loess.smooth(chr2L_t$V2, chr2L_t$V5, span=0.005), lwd=4, col="#AA8E39"))
par(new=T)
plot(x=chr2L_d$BIN_START, y=chr2L_d$TajimaD ,  type="n",  col="#C90072", xlab=NA, ylab=NA, axes=F,   ylim=c(-4,4),  xlim=c(0,3e+07))
  with(chr2L_d, lines(loess.smooth(chr2L_d$BIN_START, chr2L_d$TajimaD, span=0.005), ylab="D", lwd=4,  col="#C90072"))
  axis(side=4, labels=TRUE)

   
 
  # Chr2R
plot(x=chr2R_p$BIN_START, y=chr2R_p$PI, type="n", main="2R" ,  col="#003399",ylim=c(0,0.024), xlim=c(0,3e+07), xlab=NA, ylab="PI and Theta")
  #axis(side=2, labels=TRUE)
  #axis(side=1, labels=TRUE)
  with(chr2R_p, lines(loess.smooth(chr2R_p$BIN_START, chr2R_p$PI, span=0.005), lwd=4, col="#0061A1"))
  with(chr2R_t, lines(loess.smooth(chr2R_t$V2, chr2R_t$V5, span=0.005), lwd=4, col="#AA8E39"))
par(new=T)
plot(x=chr2R_d$BIN_START, y=chr2R_d$TajimaD ,  type="n",  col="#C90072", xlab=NA, ylab=NA, axes=F,   ylim=c(-4,4),  xlim=c(0,3e+07))
  with(chr2R_d, lines(loess.smooth(chr2R_d$BIN_START, chr2R_d$TajimaD, span=0.005),  lwd=4,  col="#C90072"))
  axis(side=4)

  

# Chr3L
plot(x=chr3L_p$BIN_START, y=chr3L_p$PI, type="n", main="3L" ,  col="#003399",ylim=c(0,0.024),  xlim=c(0,3e+07), xlab=NA, ylab="PI and Theta")
  #axis(side=2, labels=TRUE)
  #axis(side=1, labels=TRUE)
  with(chr3L_p, lines(loess.smooth(chr3L_p$BIN_START, chr3L_p$PI, span=0.005), lwd=4, col="#0061A1"))
  with(chr3L_t, lines(loess.smooth(chr3L_t$V2, chr3L_t$V5, span=0.005), lwd=4, col="#AA8E39"))
par(new=T)
plot(x=chr3L_d$BIN_START, y=chr3L_d$TajimaD ,  type="n",  col="#C90072", xlab=NA, ylab=NA, axes=F,   ylim=c(-4,4),  xlim=c(0,3e+07))
  with(chr3L_d, lines(loess.smooth(chr3L_d$BIN_START, chr3L_d$TajimaD, span=0.005),  lwd=4,  col="#C90072"))
  axis(side=4)

  
# Chr3R
plot(x=chr3R_p$BIN_START, y=chr3R_p$PI, type="n", main="3R" ,  col="#003399",ylim=c(0,0.024), xlim=c(0,3e+07), xlab=NA, ylab="PI and Theta")
  #axis(side=2, labels=TRUE)
  #axis(side=1, labels=TRUE)
  with(chr3R_p, lines(loess.smooth(chr3R_p$BIN_START, chr3R_p$PI, span=0.005), lwd=4, col="#0061A1"))
  with(chr3R_t, lines(loess.smooth(chr3R_t$V2, chr3R_t$V5, span=0.005), lwd=4, col="#AA8E39"))
par(new=T)  
plot(x=chr3R_d$BIN_START, y=chr3R_d$TajimaD ,  type="n",  col="#C90072", xlab=NA, ylab=NA, axes=F,   ylim=c(-4,4),  xlim=c(0,3e+07))
  with(chr3R_d, lines(loess.smooth(chr3R_d$BIN_START, chr3R_d$TajimaD, span=0.005),  lwd=4,  col="#C90072"))
  axis(side=4)

  

  
# Chr4
plot(x=chr4_p$BIN_START, y=chr4_p$PI, type="n", main="4" ,  col="#003399",ylim=c(0,0.010), xlim=c(0,6e+06), xlab=NA, ylab="PI and Theta")
  #axis(side=2, labels=TRUE)
  #axis(side=1, labels=TRUE)
  with(chr4_p, lines(loess.smooth(chr4_p$BIN_START, chr4_p$PI, span=0.07), lwd=4, col="#0061A1"))
  with(chr4_t, lines(loess.smooth(chr4_t$V2, chr4_t$V5, span=0.07), lwd=4, col="#AA8E39"))
par(new=T)
plot(x=chr4_d$BIN_START, y=chr4_d$TajimaD ,  type="n",  col="#C90072", xlab=NA, ylab=NA, axes=F,   ylim=c(-4,4),  xlim=c(0,6e+06))
  with(chr4_d, lines(loess.smooth(chr4_d$BIN_START, chr4_d$TajimaD, span=0.07),  lwd=4,  col="#C90072"))
  axis(side=4)

