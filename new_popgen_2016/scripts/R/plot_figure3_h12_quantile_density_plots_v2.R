# H12 quantile plots


ld_h_2l <- read.table('/home/fnew/ufgi_share/SHARE/McIntyre_Lab/ethanol/Sim_Pop_Gen/output/LD_vs_h12/chr2L_h12_and_ld.txt', header=TRUE)
ld_h_2r <- read.table('/home/fnew/ufgi_share/SHARE/McIntyre_Lab/ethanol/Sim_Pop_Gen/output/LD_vs_h12/chr2R_h12_and_ld.txt', header=TRUE)
ld_h_3l <- read.table('/home/fnew/ufgi_share/SHARE/McIntyre_Lab/ethanol/Sim_Pop_Gen/output/LD_vs_h12/chr3L_h12_and_ld.txt', header=TRUE)
ld_h_3r <- read.table('/home/fnew/ufgi_share/SHARE/McIntyre_Lab/ethanol/Sim_Pop_Gen/output/LD_vs_h12/chr3R_h12_and_ld.txt', header=TRUE)
ld_h_4 <- read.table('/home/fnew/ufgi_share/SHARE/McIntyre_Lab/ethanol/Sim_Pop_Gen/output/LD_vs_h12/chr4_h12_and_ld.txt', header=TRUE)
ld_h_x <- read.table('/home/fnew/ufgi_share/SHARE/McIntyre_Lab/ethanol/Sim_Pop_Gen/output/LD_vs_h12/chrX_h12_and_ld.txt', header=TRUE)

## Subset the groups
## Chrom2L
ld_h_2l_top <- ld_h_2l[ which(ld_h_2l$flag_top_25 ==1),]
ld_h_2l_bot <- ld_h_2l[ which(ld_h_2l$flag_bot_25 ==1),]
ld_h_2l_mid <- ld_h_2l[ which(ld_h_2l$flag_mid_50 ==1),]
## Chrom 2R
ld_h_2r_top <- ld_h_2r[ which(ld_h_2r$flag_top_25 ==1),]
ld_h_2r_bot <- ld_h_2r[ which(ld_h_2r$flag_bot_25 ==1),]
ld_h_2r_mid <- ld_h_2r[ which(ld_h_2r$flag_mid_50 ==1),]
## Chrom 3L
ld_h_3l_top <- ld_h_3l[ which(ld_h_3l$flag_top_25 ==1),]
ld_h_3l_bot <- ld_h_3l[ which(ld_h_3l$flag_bot_25 ==1),]
ld_h_3l_mid <- ld_h_3l[ which(ld_h_3l$flag_mid_50 ==1),]
## Chrom 3R
ld_h_3r_top <- ld_h_3r[ which(ld_h_3r$flag_top_25 ==1),]
ld_h_3r_bot <- ld_h_3r[ which(ld_h_3r$flag_bot_25 ==1),]
ld_h_3r_mid <- ld_h_3r[ which(ld_h_3r$flag_mid_50 ==1),]
## Chrom X
ld_h_x_top <- ld_h_x[ which(ld_h_x$flag_top_25 ==1),]
ld_h_x_bot <- ld_h_x[ which(ld_h_x$flag_bot_25 ==1),]
ld_h_x_mid <- ld_h_x[ which(ld_h_x$flag_mid_50 ==1),]
## Chrom 4
ld_h_4_top <- ld_h_4[ which(ld_h_4$flag_top_25 ==1),]
ld_h_4_bot <- ld_h_4[ which(ld_h_4$flag_bot_25 ==1),]
ld_h_4_mid <- ld_h_4[ which(ld_h_4$flag_mid_50 ==1),]

par(mfrow=c(1,1))
quartz(title="PCoA",12,6)
par(mfrow=c(1,2),oma=c(5,0,0,0),xpd=NA)

plot(1:3,4:6,main="plot 1")

plot(1:3,4:6,main="plot 2")
legend(-0.5,3.5,ncol=3,c("0-1 km","1-5 km","outside barrier"), 
       fill=c("green","orange","#C90072"), title="Fetch")

## PLOT
par(mfrow=c(2,6), xpd=FALSE)
#H2/H1
plot(density(ld_h_x_top$H2H1, bw=0.04), ylim=c(0,5.5),lwd=c(2.5), xlim=c(0,1.0), main="X\nH2/H1")
  lines(density(ld_h_x_mid$H2H1), lwd=c(2.5), col="#1F6CC2")
  lines(density(ld_h_x_bot$H2H1), lwd=c(2.5), col="#C90072")
  #  legend("topright", c("TOP 25%","MIDDLE 50%","BOTTOM 25%"), cex=c(.40),inset=c(0.01),  lty=c(1,1), lwd=c(2.5,2.5),col=c("black","#1F6CC2","#C90072"))

plot(density(ld_h_2l_top$H2H1, bw=0.04), ylim=c(0,5.5),lwd=c(2.5), xlim=c(0,1.0), main="2L\nH2/H1")
  lines(density(ld_h_2l_mid$H2H1), lwd=c(2.5), col="#1F6CC2")
  lines(density(ld_h_2l_bot$H2H1), lwd=c(2.5), col="#C90072")
#  legend("topleft", c("TOP 25%","MIDDLE 50%","BOTTOM 25%"), cex=c(.40), inset=c(0.01), lty=c(1,1), lwd=c(2.5,2.5),col=c("black","#1F6CC2","#C90072"))

plot(density(ld_h_2r_top$H2H1, bw=0.04), ylim=c(0,5.5),lwd=c(2.5), xlim=c(0,1.0), main="2R\nH2/H1")
  lines(density(ld_h_2r_mid$H2H1), lwd=c(2.5), col="#1F6CC2")
  lines(density(ld_h_2r_bot$H2H1), lwd=c(2.5), col="#C90072")
#  legend("topleft", c("TOP 25%","MIDDLE 50%","BOTTOM 25%"), cex=c(.40), inset=c(0.01),  lty=c(1,1), lwd=c(2.5,2.5),col=c("black","#1F6CC2","#C90072"))

plot(density(ld_h_3l_top$H2H1, bw=0.04), ylim=c(0,5.5),lwd=c(2.5), xlim=c(0,1.0), main="3L\nH2/H1")
  lines(density(ld_h_3l_mid$H2H1), lwd=c(2.5), col="#1F6CC2")
  lines(density(ld_h_3l_bot$H2H1), lwd=c(2.5), col="#C90072")
#  legend("topleft", c("TOP 25%","MIDDLE 50%","BOTTOM 25%"), cex=c(.40), inset=c(0.01), lty=c(1,1), lwd=c(2.5,2.5),col=c("black","#1F6CC2","#C90072"))

plot(density(ld_h_3r_top$H2H1, bw=0.04), ylim=c(0,5.5),lwd=c(2.5), xlim=c(0,1.0), main="3R\nH2/H1")
  lines(density(ld_h_3r_mid$H2H1), lwd=c(2.5), col="#1F6CC2")
  lines(density(ld_h_3r_bot$H2H1), lwd=c(2.5), col="#C90072")
#  legend("topleft", c("TOP 25%","MIDDLE 50%","BOTTOM 25%"), cex=c(.40), inset=c(0.01), lty=c(1,1), lwd=c(2.5,2.5),col=c("black","#1F6CC2","#C90072"))


plot(density(ld_h_4_top$H2H1, bw=0.04), ylim=c(0,5.5),lwd=c(2.5), xlim=c(0,1.0), main="4\nH2/H1")
  lines(density(ld_h_4_mid$H2H1), lwd=c(2.5), col="#1F6CC2")
  lines(density(ld_h_4_bot$H2H1), lwd=c(2.5), col="#C90072")
#  legend("topright", c("TOP 25%","MIDDLE 50%","BOTTOM 25%"), cex=c(.40),inset=c(0.01),  lty=c(1,1), lwd=c(2.5,2.5),col=c("black","#1F6CC2","#C90072"))

#LD
  
plot(density(ld_h_x_top$LD, bw=0.008), lwd=c(2.5), ylim=c(0,25), xlim=c(0,0.7), main="X\nLD")
  lines(density(ld_h_x_mid$LD), lwd=c(2.5), col="#1F6CC2")
  lines(density(ld_h_x_bot$LD), lwd=c(2.5), col="#C90072")
  #  legend("topright", c("Top 25%","Middle 50%","Bottom 25%"), lty=c(1,1),cex=c(.50), lwd=c(2.5,2.5),col=c("black","#1F6CC2","#C90072"))
  
plot(density(ld_h_2l_top$LD, bw=0.008), lwd=c(2.5), ylim=c(0,25), xlim=c(0,0.7), main="2L\nLD")
  lines(density(ld_h_2l_mid$LD), lwd=c(2.5), col="#1F6CC2")
  lines(density(ld_h_2l_bot$LD), lwd=c(2.5), col="#C90072")
#  legend("topright", c("Top 25%","Middle 50%","Bottom 25%"), lty=c(1,1),cex=c(.50), lwd=c(2.5,2.5),col=c("black","#1F6CC2","#C90072"))

plot(density(ld_h_2r_top$LD, bw=0.008), lwd=c(2.5), ylim=c(0,25), xlim=c(0,0.7), main="2R\nLD")
  lines(density(ld_h_2r_mid$LD), lwd=c(2.5), col="#1F6CC2")
  lines(density(ld_h_2r_bot$LD), lwd=c(2.5), col="#C90072")
#  legend("topright", c("Top 25%","Middle 50%","Bottom 25%"), lty=c(1,1),cex=c(.50), lwd=c(2.5,2.5),col=c("black","#1F6CC2","#C90072"))

plot(density(ld_h_3l_top$LD, bw=0.008), lwd=c(2.5), ylim=c(0,25), xlim=c(0,0.7), main="3L\nLD")
  lines(density(ld_h_3l_mid$LD), lwd=c(2.5), col="#1F6CC2")
  lines(density(ld_h_3l_bot$LD), lwd=c(2.5), col="#C90072")
#  legend("topright", c("Top 25%","Middle 50%","Bottom 25%"), lty=c(1,1),cex=c(.50), lwd=c(2.5,2.5),col=c("black","#1F6CC2","#C90072"))

plot(density(ld_h_3r_top$LD, bw=0.008), lwd=c(2.5), ylim=c(0,25), xlim=c(0,0.7), main="3R\nLD")
  lines(density(ld_h_3r_mid$LD), lwd=c(2.5), col="#1F6CC2")
  lines(density(ld_h_3r_bot$LD), lwd=c(2.5), col="#C90072")
#  legend("topright", c("Top 25%","Middle 50%","Bottom 25%"), lty=c(1,1),cex=c(.50), lwd=c(2.5,2.5),col=c("black","#1F6CC2","#C90072"))

plot(density(ld_h_4_top$LD, bw=0.008), lwd=c(2.5), ylim=c(0,25), xlim=c(0,0.7), main="4\nLD")
  lines(density(ld_h_4_mid$LD), lwd=c(2.5), col="#1F6CC2")
  lines(density(ld_h_4_bot$LD), lwd=c(2.5), col="#C90072")
#  legend("topright", c("Top 25%","Middle 50%","Bottom 25%"), lty=c(1,1),cex=c(.50), lwd=c(2.5,2.5),col=c("black","#1F6CC2","#C90072"))
