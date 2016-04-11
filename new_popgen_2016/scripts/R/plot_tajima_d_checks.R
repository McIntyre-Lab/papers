tdtest <- read.table('/home/fnew/dsim/tajima_d/tajima_d_10kb_remove_83lines_10permiss_nolab.Tajima.D', header=TRUE)
tdtest100 <- read.table('/home/fnew/dsim/tajima_d/tajima_d_100kb_remove_83lines_10permiss_nolab.Tajima.D', header=TRUE)


#Test with the 82 het lines removed
plt <- ggplot(data=tdtest, aes(x=tdtest$BIN_START, y=tdtest$TajimaD, fill=tdtest$CHROM)) +
  geom_bar(stat='identity')+ scale_fill_discrete()+
  ggtitle("Tajima's D Whole Genome in bins of 100KB") +
  xlab("Position along the genome") + ylab("Tajima's D in bins of 100kb") +
  facet_grid(. ~CHROM) 
plt

dentest <- density(tdtest$TajimaD, adjust=0.4)
plot(dentest, main="Tajima's D whole genome \n 100kb window")
abline(v=-1.765, col="red")
abline(v=2.095, col="red")

# Plot the test data window TsD
chr_2L <- subset(tdtest, tdtest$CHROM=='2L')
chr_2R <- subset(tdtest, tdtest$CHROM=='2R')
chr_3L <- subset(tdtest, tdtest$CHROM=='3L')
chr_3R <- subset(tdtest, tdtest$CHROM=='3R')
chr_4 <- subset(tdtest, tdtest$CHROM=='4')
chr_X <- subset(tdtest, tdtest$CHROM=='X')

par(mfrow=c(3,2))
plot(x=chr_2L$BIN_START, y=chr_2L$TajimaD, type="h", main="Tajima's D Chr2L")
abline(h=-1.783, col="red")
abline(h=2.071, col="red")
plot(x=chr_2R$BIN_START, y=chr_2R$TajimaD, type="h",  main="Tajima's D Chr2R")
abline(h=-1.783, col="red")
abline(h=2.071, col="red")
plot(x=chr_3L$BIN_START, y=chr_3L$TajimaD, type="h",  main="Tajima's D Chr3L")
abline(h=-1.783, col="red")
abline(h=2.071, col="red")
plot(x=chr_3R$BIN_START, y=chr_3R$TajimaD, type="h",  main="Tajima's D Chr3R")
abline(h=-1.783, col="red")
abline(h=2.071, col="red")
plot(x=chr_4$BIN_START, y=chr_4$TajimaD, type="h",  main="Tajima's D Chr4")
abline(h=-1.783, col="red")
abline(h=2.071, col="red")
plot(x=chr_X$BIN_START, y=chr_X$TajimaD, type="h",  main="Tajima's D ChrX")
abline(h=-1.783, col="red")
abline(h=2.071, col="red")





plt <- ggplot(data=chr_4, aes(x=chr_4$BIN_START, y=chr_4$TajimaD)) +
  geom_bar(stat='identity') +
  ggtitle("Tajima's D along Chromosome 4 in bins of 1kb") +
  xlab("Position along Chromosome 4") + ylab("Tajima's D in bins of 1kb")
plt
plt <- ggplot(data=chr_X, aes(x=chr_X$BIN_START, y=chr_X$TajimaD)) +
  geom_bar(stat='identity') +
  ggtitle("Tajima's D along Chromosome X in bins of 1kb") +
  xlab("Position along Chromosome X") + ylab("Tajima's D in bins of 1kb")
plt
plt <- ggplot(data=chr_3L, aes(x=chr_3L$BIN_START, y=chr_3L$TajimaD)) +
  geom_bar(stat='identity') +
  ggtitle("Tajima's D along Chromosome 3L in bins of 1kb") +
  xlab("Position along Chromosome 3L") + ylab("Tajima's D in bins of 1kb")
plt
plt <- ggplot(data=chr_3R, aes(x=chr_3R$BIN_START, y=chr_3R$TajimaD)) +
  geom_bar(stat='identity') +
  ggtitle("Tajima's D along Chromosome 3R in bins of 1kb") +
  xlab("Position along Chromosome 3R") + ylab("Tajima's D in bins of 1kb")
plt

plt <- ggplot(data=chr_2L, aes(x=chr_2L$BIN_START, y=chr_2L$TajimaD)) +
  geom_bar(stat='identity') +
  ggtitle("Tajima's D along Chromosome 2L in bins of 1kb") +
  xlab("Position along Chromosome 2L") + ylab("Tajima's D in bins of 1kb")
plt
plt <- ggplot(data=chr_2R, aes(x=chr_2R$BIN_START, y=chr_2R$TajimaD)) +
  geom_bar(stat='identity') +
  ggtitle("Tajima's D along Chromosome 2R in bins of 1kb") +
  xlab("Position along Chromosome 2R") + ylab("Tajima's D in bins of 1kb")
plt




##Plot density curves for TsD
den2L <- density(chr_2L$TajimaD, adjust=0.4)
den2R <- density(chr_2R$TajimaD, adjust=0.4)
den3L <- density(chr_3L$TajimaD, adjust=0.4)
den3R <- density(chr_3R$TajimaD, adjust=0.4)
den4 <- density(chr_4$TajimaD, adjust=0.2) #1kb window, not 10kb
denX <- density(chr_X$TajimaD, adjust=0.4)

par(mfrow=c(3,2))
plot(den2L, main="Tajima's D on Chr2L \n 10kb window")
abline(v=-1.765, col="red")
abline(v=2.095, col="red")
plot(den2R, main="Tajima's D on Chr2R \n 10kb window")
abline(v=-1.765, col="red")
abline(v=2.095, col="red")
plot(den3L, main="Tajima's D on Chr3L \n 10kb window")
abline(v=-1.765, col="red")
abline(v=2.095, col="red")
plot(den3R, main="Tajima's D on Chr3R \n 10kb window")
abline(v=-1.765, col="red")
abline(v=2.095, col="red")
plot(den4, main="Tajima's D on Chr4 \n 10kb window")
abline(v=-1.765, col="red")
abline(v=2.095, col="red")
plot(denX, main="Tajima's D on ChrX \n 10kb window")
abline(v=-1.765, col="red")
abline(v=2.095, col="red")

## Test the results from removing the six Dmel lines
tdmel <- read.table('/home/fnew/dsim/tajima_d/check_contam_tsd_100kb.Tajima.D', header=TRUE)
tdmel10 <- read.table('/home/fnew/dsim/tajima_d/check_contam_tsd_10kb.Tajima.D', header=TRUE)

chr_2L <- subset(tdmel, tdmel$CHROM=='2L')
chr_2R <- subset(tdmel, tdmel$CHROM=='2R')
chr_3L <- subset(tdmel, tdmel$CHROM=='3L')
chr_3R <- subset(tdmel, tdmel$CHROM=='3R')
chr_4 <- subset(tdmel, tdmel$CHROM=='4')
chr_X <- subset(tdmel, tdmel$CHROM=='X')

par(mfrow=c(3,2))
plot(x=chr_2L$BIN_START, y=chr_2L$TajimaD, type="h", ylim=c(-3,3), main="Tajima's D Chr2L NO MEL")
abline(h=-1.783, col="red")
abline(h=2.071, col="red")
plot(x=chr_2R$BIN_START, y=chr_2R$TajimaD, type="h",ylim=c(-3,3),  main="Tajima's D Chr2R NO MEL")
abline(h=-1.783, col="red")
abline(h=2.071, col="red")
plot(x=chr_3L$BIN_START, y=chr_3L$TajimaD, type="h",ylim=c(-3,3),  main="Tajima's D Chr3L NO MEL")
abline(h=-1.783, col="red")
abline(h=2.071, col="red")
plot(x=chr_3R$BIN_START, y=chr_3R$TajimaD, type="h", ylim=c(-3,3), main="Tajima's D Chr3R NO MEL")
abline(h=-1.783, col="red")
abline(h=2.071, col="red")
plot(x=chr_4$BIN_START, y=chr_4$TajimaD, type="h", ylim=c(-3,3), main="Tajima's D Chr4 NO MEL")
abline(h=-1.783, col="red")
abline(h=2.071, col="red")
plot(x=chr_X$BIN_START, y=chr_X$TajimaD, type="h",  ylim=c(-3,3), main="Tajima's D ChrX NO MEL")
abline(h=-1.783, col="red")
abline(h=2.071, col="red")

chr_2L <- subset(tdmel10, tdmel10$CHROM=='2L')
chr_2R <- subset(tdmel10, tdmel10$CHROM=='2R')
chr_3L <- subset(tdmel10, tdmel10$CHROM=='3L')
chr_3R <- subset(tdmel10, tdmel10$CHROM=='3R')
chr_4 <- subset(tdmel10, tdmel10$CHROM=='4')
chr_X <- subset(tdmel10, tdmel10$CHROM=='X')


par(mfrow=c(3,2))
plot(x=chr_2L$BIN_START, y=chr_2L$TajimaD, type="h", ylim=c(-3,4), main="Tajima's D Chr2L NO MEL")
abline(h=-1.783, col="red")
abline(h=2.071, col="red")
plot(x=chr_2R$BIN_START, y=chr_2R$TajimaD, type="h",ylim=c(-3,4),  main="Tajima's D Chr2R NO MEL")
abline(h=-1.783, col="red")
abline(h=2.071, col="red")
plot(x=chr_3L$BIN_START, y=chr_3L$TajimaD, type="h",ylim=c(-3,4),  main="Tajima's D Chr3L NO MEL")
abline(h=-1.783, col="red")
abline(h=2.071, col="red")
plot(x=chr_3R$BIN_START, y=chr_3R$TajimaD, type="h", ylim=c(-3,4), main="Tajima's D Chr3R NO MEL")
abline(h=-1.783, col="red")
abline(h=2.071, col="red")
plot(x=chr_4$BIN_START, y=chr_4$TajimaD, type="h", ylim=c(-3,4), main="Tajima's D Chr4 NO MEL")
abline(h=-1.783, col="red")
abline(h=2.071, col="red")
plot(x=chr_X$BIN_START, y=chr_X$TajimaD, type="h",  ylim=c(-3,4), main="Tajima's D ChrX NO MEL")
abline(h=-1.783, col="red")
abline(h=2.071, col="red")