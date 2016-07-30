## Figure 2 is LD vs Tajima's D
## Panel A will show drop in LD for the chromosomes
## Panel B will show LD vs D (100kb, 10kb, and 1kb windows), including the regression lines for negative and positive data

# LD pairwise comparisons across the genome
ld2l<- read.table('/ufgi$/SHARE/McIntyre_Lab/ethanol/Sim_Pop_Gen/output/LD/plink_LD_output/plink_2L_LD.ld',  header=TRUE)
ld2r<- read.table('/ufgi$/SHARE/McIntyre_Lab/ethanol/Sim_Pop_Gen/output/LD/plink_LD_output/plink_2R_LD.ld',  header=TRUE)
ld3l<- read.table('/ufgi$/SHARE/McIntyre_Lab/ethanol/Sim_Pop_Gen/output/LD/plink_LD_output/plink_3L_LD.ld',  header=TRUE)
ld3r<- read.table('/ufgi$/SHARE/McIntyre_Lab/ethanol/Sim_Pop_Gen/output/LD/plink_LD_output/plink_3R_LD.ld',  header=TRUE)
ld4<- read.table('/ufgi$/SHARE/McIntyre_Lab/ethanol/Sim_Pop_Gen/output/LD/plink_LD_output/plink_4_LD.ld',  header=TRUE)
ldx<- read.table('/ufgi$/SHARE/McIntyre_Lab/ethanol/Sim_Pop_Gen/output/LD/plink_LD_output/plink_X_LD.ld',  header=TRUE)


ld2l$dist <- (ld2l$BP_B - ld2l$BP_A)
ld2r$dist <- (ld2r$BP_B - ld2r$BP_A)
ld3l$dist <- (ld3l$BP_B - ld3l$BP_A)
ld3r$dist <- (ld3r$BP_B - ld3r$BP_A)
ld4$dist <- (ld4$BP_B - ld4$BP_A)
ldx$dist <- (ldx$BP_B - ldx$BP_A)

## Part 1: plot_figure2_ld_decay_rank_means_plots.R
##-----------------------------------------------------------------
## Part 2: Plot LD vs D and add regression lines for x<0 and x>0
  

# LD vs D in windows
chr4 <- read.table("/ufgi$/SHARE/McIntyre_Lab/ethanol/Sim_Pop_Gen/output/LD/ld_and_D/chr4_ld_and_td.txt", header=TRUE)
chr2l <- read.table("/ufgi$/SHARE/McIntyre_Lab/ethanol/Sim_Pop_Gen/output/LD/ld_and_D/chr2l_ld_and_td.txt", header=TRUE)
chrX <- read.table("/ufgi$/SHARE/McIntyre_Lab/ethanol/Sim_Pop_Gen/output/LD/ld_and_D/chrX_ld_and_td.txt", header=TRUE)
chr3L <- read.table("/ufgi$/SHARE/McIntyre_Lab/ethanol/Sim_Pop_Gen/output/LD/ld_and_D/chr3L_ld_and_td.txt", header=TRUE)
chr2R <- read.table("/ufgi$/SHARE/McIntyre_Lab/ethanol/Sim_Pop_Gen/output/LD/ld_and_D/chr2R_ld_and_td.txt", header=TRUE)
chr3R <- read.table("/ufgi$/SHARE/McIntyre_Lab/ethanol/Sim_Pop_Gen/output/LD/ld_and_D/chr3R_ld_and_td.txt", header=TRUE)



chr2L_pos <- chr2l[ which(chr2l$TajimaD > 0),]
chr2L_neg <- chr2l[ which(chr2l$TajimaD < 0),]
chr2R_pos <- chr2R[ which(chr2R$TajimaD > 0),]
chr2R_neg <- chr2R[ which(chr2R$TajimaD < 0),]
chr3L_pos <- chr3L[ which(chr3L$TajimaD > 0),]
chr3L_neg <- chr3L[ which(chr3L$TajimaD < 0),]
chr3R_pos <- chr3R[ which(chr3R$TajimaD > 0),]
chr3R_neg <- chr3R[ which(chr3R$TajimaD < 0),]
chr4_pos <- chr4[ which(chr4$TajimaD > 0),]
chr4_neg <- chr4[ which(chr4$TajimaD < 0),]
chrX_pos <- chrX[ which(chrX$TajimaD > 0),]
chrX_neg <- chrX[ which(chrX$TajimaD < 0),]

dev.off
par(mfrow=c(1,1))
par(xpd=FALSE)
par(mfrow=c(2,3))
par(mfrow=c(1,3))

#X
pos <- lm(chrX_pos$R2~ chrX_pos$TajimaD)
neg <- lm(chrX_neg$R2 ~ chrX_neg$TajimaD)
summary(lm(chrX_pos$R2 ~ chrX_pos$TajimaD))$r.squared #0.02488322
summary(lm(chrX_pos$R2 ~ chrX_pos$TajimaD))$adj.r.squared #0.0242463
summary(lm(chrX_neg$R2 ~ chrX_neg$TajimaD))$r.squared #0.3984849
summary(lm(chrX_neg$R2 ~ chrX_neg$TajimaD))$adj.r.squared #0.397371


plot(x=chrX$TajimaD, y=chrX$R2, xlab="D", ylab="LD", pch=20, xlim=c(-4,4), ylim=c(0,0.6), main="X")
abline(pos, lwd=3, lty=2,  col="#C90072", xlim=c(0,4))
abline(neg, lwd=3, lty=2,  col="#1F6CC2", xlim=c(-4,0))
text(3,0.59, expression(R^{2} ~  "= 0.0248"), col="#C90072")
text(3,0.54, expression(R^{2} ~ "= 0.3973"),  col="#1F6CC2")


#2L
pos <- lm(chr2L_pos$ld~ chr2L_pos$TajimaD)
neg <- lm(chr2L_neg$ld ~ chr2L_neg$TajimaD)
 summary(lm(chr2L_pos$ld ~ chr2L_pos$TajimaD))$r.squared #0.1176369
 summary(lm(chr2L_pos$ld ~ chr2L_pos$TajimaD))$adj.r.squared #0.1172018
 summary(lm(chr2L_neg$ld ~ chr2L_neg$TajimaD))$r.squared #0.3099885
 summary(lm(chr2L_neg$ld ~ chr2L_neg$TajimaD))$adj.r.squared #0.3067027


plot(x=chr2l$TajimaD, y=chr2l$ld, xlab="D", ylab="LD", pch=20, xlim=c(-4,4), ylim=c(0,0.6), main="2L")
  abline(pos, lwd=3, lty=2, col="#C90072", xlim=c(0,4))
  abline(neg, lwd=3, lty=2,  col="#1F6CC2", xlim=c(-4,0))
  text(3,0.59, expression(R^{2} ~  "= 0.1176"), col="#C90072")
  text(3,0.54, expression(R^{2} ~ "= 0.3099"),  col="#1F6CC2")

#2R
pos <- lm(chr2R_pos$R2~ chr2R_pos$TajimaD)
neg <- lm(chr2R_neg$R2 ~ chr2R_neg$TajimaD)
 summary(lm(chr2R_pos$R2 ~ chr2R_pos$TajimaD))$r.squared #0.1264926
 summary(lm(chr2R_pos$R2 ~ chr2R_pos$TajimaD))$adj.r.squared #0.1260171
 summary(lm(chr2R_neg$R2 ~ chr2R_neg$TajimaD))$r.squared #0.1379368
 summary(lm(chr2R_neg$R2 ~ chr2R_neg$TajimaD))$adj.r.squared #0.13373`6


plot(x=chr2R$TajimaD, y=chr2R$R2, xlab="D", ylab="LD", pch=20, xlim=c(-4,4), ylim=c(0,0.6), main="2R")
  abline(pos, lwd=3, lty=2,  col="#C90072", xlim=c(0,4))
  abline(neg, lwd=3, lty=2,  col="#1F6CC2", xlim=c(-4,0))
  text(3,0.59, expression(R^{2} ~  "= 0.1264"), col="#C90072")
  text(3,0.54, expression(R^{2} ~ "= 0.1379"),  col="#1F6CC2")

#3L
pos <- lm(chr3L_pos$R2~ chr3L_pos$TajimaD)
neg <- lm(chr3L_neg$R2 ~ chr3L_neg$TajimaD)
  summary(lm(chr3L_pos$R2 ~ chr3L_pos$TajimaD))$r.squared #0.177305
  summary(lm(chr3L_pos$R2 ~ chr3L_pos$TajimaD))$adj.r.squared #0.1769223
  summary(lm(chr3L_neg$R2 ~ chr3L_neg$TajimaD))$r.squared #0.1836232
  summary(lm(chr3L_neg$R2 ~ chr3L_neg$TajimaD))$adj.r.squared #0.1787926
  

plot(x=chr3L$TajimaD, y=chr3L$R2, xlab="D", ylab="LD", pch=20,  xlim=c(-4,4), ylim=c(0,0.6), main="3L")
  abline(pos, lwd=3, lty=2,  col="#C90072", xlim=c(0,4))
  abline(neg, lwd=3, lty=2,  col="#1F6CC2", xlim=c(-4,0))
  text(3,0.59, expression(R^{2} ~  "= 0.1773"), col="#C90072")
  text(3,0.54, expression(R^{2} ~ "= 0.1836"),  col="#1F6CC2")
  
  
  #3R
pos <- lm(chr3R_pos$R2~ chr3R_pos$TajimaD)
neg <- lm(chr3R_neg$R2 ~ chr3R_neg$TajimaD)
  summary(lm(chr3R_pos$R2 ~ chr3R_pos$TajimaD))$r.squared #0.1515024
  summary(lm(chr3R_pos$R2 ~ chr3R_pos$TajimaD))$adj.r.squared #0.1511713
  summary(lm(chr3R_neg$R2 ~ chr3R_neg$TajimaD))$r.squared #0.4758257
  summary(lm(chr3R_neg$R2 ~ chr3R_neg$TajimaD))$adj.r.squared #0.4721601
  

plot(x=chr3R$TajimaD, y=chr3R$R2, xlab="D", ylab="LD", pch=20, xlim=c(-4,4), ylim=c(0,0.6), main="3R")
  abline(pos, lwd=3, lty=2,  col="#C90072", xlim=c(0,4))
  abline(neg, lwd=3, lty=2,  col="#1F6CC2", xlim=c(-4,0))
  text(3,0.59, expression(R^{2} ~  "= 0.1515"), col="#C90072")
  text(3,0.54, expression(R^{2} ~ "= 0.4758"),  col="#1F6CC2")
 
  

  
#4 -- There are no positive values of D for Chromosome 4!
neg <- lm(chr4_neg$ld ~ chr4_neg$TajimaD)
  summary(lm(chr4_neg$ld ~ chr4_neg$TajimaD))$r.squared #0.3561973
  summary(lm(chr4_neg$ld ~ chr4_neg$TajimaD))$adj.r.squared #0.3497593
  

plot(x=chr4$TajimaD, y=chr4$ld, xlab="D", ylab="LD", pch=20, xlim=c(-4,4), ylim=c(0,0.6), main="4")
  abline(neg, lwd=3, lty=2,  col="#1F6CC2", xlim=c(-4,0))
  text(3,0.54, expression(R^{2} ~ "= 0.3562"),  col="#1F6CC2")  
  

  
  