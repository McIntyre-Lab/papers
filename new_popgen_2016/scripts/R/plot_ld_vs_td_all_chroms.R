# plot LD and TD 

chr4 <- read.table("/home/fnew/dsim/vcf_split_td_ld/chr4_ld_and_td.txt", header=TRUE)

plot(x=chr4$TajimaD, y=chr4$ld, xlab="Tajima's D (10kb windows)", ylab="LD (10kb windows)", main="Tajima's D vs LD on Chromosome 4")
abline(plt1)


plt <- lm(chr4$ld ~ chr4$TajimaD)
abline(plt)
summary(lm(chr4$ld ~ chr4$TajimaD))
legend("topright", bty="n", legend=paste("R2 = 0.3498"))

## Chrom 2L
chr2l <- read.table("/home/fnew/dsim/vcf_split_td_ld/chr2l_ld_and_td.txt", header=TRUE)

plot(x=chr2l$TajimaD, y=chr2l$ld, xlab="Tajima's D (10kb windows)", ylab="LD (10kb windows)", main="Tajima's D vs LD on Chromosome 2L")
abline(plt1)


plt <- lm(chr2l$ld ~ chr2l$TajimaD)
abline(plt)
summary(lm(chr2l$ld ~ chr2l$TajimaD))
legend("topright", bty="n", legend=paste("R2 = 0.02066"))


ggplot(chrX, aes(x=chrX$TajimaD, y=chrX$R2)) +
  geom_point(shape=1) +    
  geom_smooth(method=lm)

## Chrom X
chrX <- read.table("/home/fnew/dsim/vcf_split_td_ld/chrX_ld_and_td.txt", header=TRUE)

plot(x=chrX$TajimaD, y=chrX$R2, xlab="Tajima's D (10kb windows)", ylab="LD (10kb windows)", main="Tajima's D vs LD on Chromosome X")
abline(plt1)

negX <- chrX[ which(chrX$TajimaD<0),]
posX <- chrX[ which(chrX$TajimaD>0),]

plot(x=posX$TajimaD, y=posX$R2,  xlab="Tajima's D (10kb windows)",ylab="LD (10kb windows)", main="Tajima's D vs LD on Chromosome X")



plt <- lm(chrX$R2 ~ chrX$TajimaD)
abline(plt)
summary(lm(chrX$R2 ~ chrX$TajimaD))
legend("topright", bty="n", legend=paste("R2 = 0.1259"))




## Chrom 3L
chr3L <- read.table("/home/fnew/dsim/vcf_split_td_ld/chr3L_ld_and_td.txt", header=TRUE)

plot(x=chr3L$TajimaD, y=chr3L$R2, xlab="Tajima's D (10kb windows)", ylab="LD (10kb windows)", main="Tajima's D vs LD on Chromosome 3L")
abline(plt1)


plt <- lm(chr3L$R2 ~ chr3L$TajimaD)
abline(plt)
summary(lm(chr3L$R2 ~ chr3L$TajimaD))
legend("topright", bty="n", legend=paste("R2 = 0.1259"))


## Chrom 3R
chr3R <- read.table("/home/fnew/dsim/vcf_split_td_ld/chr3R_ld_and_td.txt", header=TRUE)

plot(x=chr3R$TajimaD, y=chr3R$R2, xlab="Tajima's D (10kb windows)", ylab="LD (10kb windows)", main="Tajima's D vs LD on Chromosome 3R")
abline(plt1)


plt <- lm(chr3R$R2 ~ chr3R$TajimaD)
abline(plt)
summary(lm(chr3R$R2 ~ chr3R$TajimaD))
legend("topright", bty="n", legend=paste("R2 = 0.02953"))


## Chrom 2R
chr2R <- read.table("/home/fnew/dsim/vcf_split_td_ld/chr2R_ld_and_td.txt", header=TRUE)

plot(x=chr2R$TajimaD, y=chr2R$R2, xlab="Tajima's D (10kb windows)", ylab="LD (10kb windows)", main="Tajima's D vs LD on Chromosome 2R")
abline(plt1)


plt <- lm(chr2R$R2 ~ chr2R$TajimaD)
abline(plt)
summary(lm(chr2R$R2 ~ chr2R$TajimaD))
legend("topright", bty="n", legend=paste("R2 = 0.009303"))



## 3D surface plots 

scatterplot3d(chrX$bin, chrX$TajimaD, chrX$R2, main="Tajima's D vs LD by Position ChromX")

scatterplot3d(chr3L$bin, chr3L$TajimaD, chr3L$R2, main="Tajima's D vs LD by Position Chrom3L")

scatterplot3d(chr3R$bin, chr3R$TajimaD, chr3R$R2, main="Tajima's D vs LD by Position Chrom3R")

scatterplot3d(chr2l$bin, chr2l$TajimaD, chr2l$ld, main="Tajima's D vs LD by Position Chrom2L")

scatterplot3d(chr2R$bin, chr2R$TajimaD, chr2R$R2, main="Tajima's D vs LD by Position Chrom2R")

scatterplot3d(chr4$bin, chr4$TajimaD, chr4$R2, main="Tajima's D vs LD by Position Chrom4")

## LD vs position

plot(x=chr4$bin, y=chr4$ld, xlab="Position (10kb windows)", ylab="LD (10kb windows)", main="LD on Chromosome 4")
plot(x=chr2l$bin, y=chr2l$ld, xlab="Position (10kb windows)", ylab="LD (10kb windows)", main="LD on Chromosome 2L")
plot(x=chr2R$bin, y=chr2R$R2, xlab="Position (10kb windows)", ylab="LD (10kb windows)", main="LD on Chromosome 2R")
plot(x=chr3L$bin, y=chr3L$R2, xlab="Position (10kb windows)", ylab="LD (10kb windows)", main="LD on Chromosome 3L")
plot(x=chr3R$bin, y=chr3R$R2, xlab="Position (10kb windows)", ylab="LD (10kb windows)", main="LD on Chromosome 3R")
plot(x=chrX$bin, y=chrX$R2, xlab="Position (10kb windows)", ylab="LD (10kb windows)", main="LD on Chromosome X")



