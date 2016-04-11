# Attempting to plot LD...

ld<- read.table('/home/fnew/dsim/ld/plink.ld',  header=TRUE)

ggplot(ld) +
  geom_line(aes(x=BP_B - BP_A, y = R2))
plot(x=ld$BP_B - ld$BP_A, y=ld$R2, type="l")


ld2l<- read.table('/home/fnew/dsim/ld/plink_2L_LD.ld',  header=TRUE)
ld2r<- read.table('/home/fnew/dsim/ld/plink_2R_LD.ld',  header=TRUE)
ld3l<- read.table('/home/fnew/dsim/ld/plink_3L_LD.ld',  header=TRUE)
ld3r<- read.table('/home/fnew/dsim/ld/plink_3R_LD.ld',  header=TRUE)
ld4<- read.table('/home/fnew/dsim/ld/plink_4_LD.ld',  header=TRUE)
ldx<- read.table('/home/fnew/dsim/ld/plink_X_LD.ld',  header=TRUE)
  

  ld2l$dist <- (ld2l$BP_B - ld2l$BP_A)
  ld2r$dist <- (ld2r$BP_B - ld2r$BP_A)
  ld3l$dist <- (ld3l$BP_B - ld3l$BP_A)
  ld3r$dist <- (ld3r$BP_B - ld3r$BP_A)
  ld4$dist <- (ld4$BP_B - ld4$BP_A)
  ldx$dist <- (ldx$BP_B - ldx$BP_A)
  
plot(x=ld4$dist, y=ld4$R2, pch=20, xlab="Distance (bp)", ylab="R2", main="LD Decay Chromosome 4")
plot(x=ld2l$dist, y=ld2l$R2, pch=20, xlab="Distance (bp)", ylab="R2", main="LD Decay Chromosome 2L")
plot(x=ld2r$dist, y=ld2r$R2, pch=20, xlab="Distance (bp)", ylab="R2", main="LD Decay Chromosome 2R")
plot(x=ld3l$dist, y=ld3l$R2, pch=20, xlab="Distance (bp)", ylab="R2", main="LD Decay Chromosome 3L")
plot(x=ld3r$dist, y=ld3r$R2, pch=20, xlab="Distance (bp)", ylab="R2", main="LD Decay Chromosome 3R")
plot(x=ldx$dist, y=ldx$R2, pch=20, xlab="Distance (bp)", ylab="R2", main="LD Decay Chromosome X")

#xlim=2500 
plot(x=ld4$dist, y=ld4$R2, pch=20, xlim=c(0,2500), xlab="Distance (bp)", ylab="R2", main="LD Decay Chromosome 4")
plot(x=ld2l$dist, y=ld2l$R2, pch=20, xlim=c(0,2500), xlab="Distance (bp)", ylab="R2", main="LD Decay Chromosome 2L")
plot(x=ld2r$dist, y=ld2r$R2, pch=20, xlim=c(0,2500), xlab="Distance (bp)", ylab="R2", main="LD Decay Chromosome 2R")
plot(x=ld3l$dist, y=ld3l$R2, pch=20, xlim=c(0,2500), xlab="Distance (bp)", ylab="R2", main="LD Decay Chromosome 3L")
plot(x=ld3r$dist, y=ld3r$R2, pch=20, xlim=c(0,2500), xlab="Distance (bp)", ylab="R2", main="LD Decay Chromosome 3R")
plot(x=ldx$dist, y=ldx$R2, pch=20, xlim=c(0,2500), xlab="Distance (bp)", ylab="R2", main="LD Decay Chromosome X")


#xlim=max x 
plot(x=ld4$dist, y=ld4$R2, pch=20, xlim=c(0,150000), xlab="Distance (bp)", ylab="R2", main="LD Decay Chromosome 4")
plot(x=ld2l$dist, y=ld2l$R2, pch=20, xlim=c(0,150000), xlab="Distance (bp)", ylab="R2", main="LD Decay Chromosome 2L")
plot(x=ld2r$dist, y=ld2r$R2, pch=20, xlim=c(0,150000), xlab="Distance (bp)", ylab="R2", main="LD Decay Chromosome 2R")
plot(x=ld3l$dist, y=ld3l$R2, pch=20, xlim=c(0,150000), xlab="Distance (bp)", ylab="R2", main="LD Decay Chromosome 3L")
plot(x=ld3r$dist, y=ld3r$R2, pch=20, xlim=c(0,150000), xlab="Distance (bp)", ylab="R2", main="LD Decay Chromosome 3R")
plot(x=ldx$dist, y=ldx$R2, pch=20, xlim=c(0,150000), xlab="Distance (bp)", ylab="R2", main="LD Decay Chromosome X")




write.table(ld2l, file="/home/fnew/dsim/ld/chr2L_LD_dist.txt",
   append=FALSE, quote = FALSE, sep="\t")