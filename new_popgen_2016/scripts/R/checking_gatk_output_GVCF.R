## Checking the raw X vcf with lower PHREDD cut off ##
## Will plot het, snpdensity, pi, tajima's d, ld, ....   ##


#### HETEROZYGOSITY ####
het<-read.table("/home/fnew/dsim/checking_gatk_output/sim_gvcf_biall_het.het", header=TRUE)
attach(het)
het2<- het[order(F),]

plt <- plot(het2$F, ylab="F", main="Coefficient of inbreeding", ylim=c(-12,2))

abline(h=0.96, col="red")
#abline(h=0.88, col="blue")
#abline(h=0.67, col="orange")

text(x=25, y=1.5, "F=0.96, 15 gen inbreeding", cex=0.8)

ten<-read.table("/home/fnew/dsim/checking_gatk_output/10permiss.het", header=TRUE)
nomiss<-read.table("/home/fnew/dsim/checking_gatk_output/nomiss.het", header=TRUE)

attach(ten)
ten2<-ten[order(F),]
attach(nomiss)
nomiss2<-nomiss[order(F),]

plt<- plot(ten2$F, main="10 percent missing")
abline(h=0.96, col="red")
abline(h=0.88, col="blue")
abline(h=0.67, col="orange")

plt2<-plot(nomiss2$F, main="no missing data")
abline(h=0.96, col="red")
abline(h=0.88, col="blue")
abline(h=0.67, col="orange")


#### SNP DENSITY ####
xsnp <- read.table("/home/fnew/dsim/checking_gatk_output/10permiss_gvcf_10kb.snpden", header=TRUE)

plot(x=xsnp$BIN_START, y=xsnp$SNP_COUNT, type="h", main="SNP density on the X\n10kb windows",
     ylab="SNP Count", xlab="Position")

dengen <- density(xsnp$SNP_COUNT, adjust=0.5)
plot(dengen, main="SNP density across the X \n 10kb windows")
points(x=386.73, y=0, pch=17, col="red", cex=2) #mean
points(x=395.5, y=0, pch="*", col="blue", cex=4) #median

avg <- mean(xsnp$SNP_COUNT)
med <- median(xsnp$SNP_COUNT)

#### PI ####
pi <- read.table("/home/fnew/dsim/checking_gatk_output/10permiss_gvcf_10kb.windowed.pi", header=TRUE)

plot(x=pi$BIN_START, y=pi$PI, type="l", main="Nucleotide Diversity on X" , xlab="Position in 10kb windows", ylab="pi", ylim=c(0,0.015))


#### TAJIMA'S D ####
d <- read.table("/home/fnew/dsim/checking_gatk_output/10permiss_gvcf_10kb.Tajima.D", header=TRUE)

plt <- ggplot(data=d, aes(x=d$BIN_START, y=d$TajimaD)) +
  geom_bar(stat='identity') +
  ggtitle("Tajima's D along Chromosome X in bins of 10kb") +
  xlab("Position along Chromosome X") + ylab("Tajima's D in bins of 10kb")
plt

plt <- ggplot(data=d, aes(x=d$BIN_START, y=d$TajimaD)) +
  geom_bar(stat='identity') +
  ggtitle("Tajima's D along Chromosome X in bins of 10kb") +
  xlab("Position along Chromosome X") + ylab("Tajima's D in bins of 10kb")
plt

plot(x=d$BIN_START, y=d$TajimaD, type="h", ylim=c(-3,4), main="Tajima's D ChrX", ylab="D", xlab="Position in 10kb windows")
abline(h=-1.765, col="red")
abline(h=2.095, col="red")



#### LD vs TsD ####

## Chrom X
chrX <- read.table("/home/fnew/dsim/checking_gatk_output/chrX_ld_and_td.txt", header=TRUE)

plot(x=chrX$TajimaD, y=chrX$R2, xlab="Tajima's D (10kb windows)", ylab="LD (10kb windows)", main="Tajima's D vs LD on Chromosome X\nChecking")
abline(plt1)

negX <- chrX[ which(chrX$TajimaD<0),]
posX <- chrX[ which(chrX$TajimaD>0),]

plot(x=posX$TajimaD, y=posX$R2,  xlab="Tajima's D (10kb windows)",ylab="LD (10kb windows)", main="Tajima's D vs LD on Chromosome X")



plt <- lm(chrX$R2 ~ chrX$TajimaD)
abline(plt)
summary(lm(chrX$R2 ~ chrX$TajimaD))
legend("topright", bty="n", legend=paste("R2 = 0.1259"))














