# Plot distributions of Tajima's D
# TsD calculated in windows of 50 and 100 for comparison...100 looks better....

td50 <- read.table('/home/fnew/mclab/ethanol/Sim_Pop_Gen/output/tajima_d/dsim_tajimad_50.Tajima.D', header=TRUE)
td100 <- read.table('/home/fnew/ufgi_share/SHARE/McIntyre_Lab/ethanol/Sim_Pop_Gen/output/tajima_d/filtered_vcf_TajimaD_100kb_window.Tajima.D', header=TRUE)
td10kb <- read.table('/home/fnew/dsim/tajima_d/filter_nolab_tajd_10kb.Tajima.D', header=TRUE)
td100kb <- read.table('/home/fnew/dsim/tajima_d/filter_nolab_tajd_100kb.Tajima.D', header=TRUE)
td1kb <- read.table('/home/fnew/dsim/filter_nolab_tajd_1kb.Tajima.D', header=TRUE)


#Plotting the whole genome with the 100kb window TsD
plt <- ggplot(data=td100kb, aes(x=td100kb$BIN_START, y=td100kb$TajimaD, fill=td100$CHROM)) +
  geom_bar(stat='identity')+ scale_fill_discrete()+
  ggtitle("Tajima's D Whole Genome in bins of 100KB") +
  xlab("Position along the genome") + ylab("Tajima's D in bins of 100kb") +
  facet_grid(. ~CHROM) 
plt

dengen <- density(td100kb$TajimaD, adjust=0.4)
plot(dengen, main="Tajima's D whole genome \n 100kb window")


# Plot the 1kb window TsD
chr_2L <- subset(td1kb, td1kb$CHROM=='2L')
chr_2R <- subset(td1kb, td1kb$CHROM=='2R')
chr_3L <- subset(td1kb, td1kb$CHROM=='3L')
chr_3R <- subset(td1kb, td1kb$CHROM=='3R')
chr_4 <- subset(td1kb, td1kb$CHROM=='4')
chr_X <- subset(td1kb, td1kb$CHROM=='X')


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



##Plot the 10kb window TsD
chr_2L <- subset(td100kb, td100kb$CHROM=='2L')
chr_2R <- subset(td100kb, td100kb$CHROM=='2R')
chr_3L <- subset(td100kb, td100kb$CHROM=='3L')
chr_3R <- subset(td100kb, td100kb$CHROM=='3R')
chr_4 <- subset(td100kb, td100kb$CHROM=='4')
chr_X <- subset(td100kb, td100kb$CHROM=='X')


chr_2L <- subset(td10kb, td10kb$CHROM=='2L')
chr_2R <- subset(td10kb, td10kb$CHROM=='2R')
chr_3L <- subset(td10kb, td10kb$CHROM=='3L')
chr_3R <- subset(td10kb, td10kb$CHROM=='3R')
chr_4 <- subset(td10kb, td10kb$CHROM=='4')
chr_X <- subset(td10kb, td10kb$CHROM=='X')

par(mfrow=c(3,2))
plot(x=chr_2L$BIN_START, y=chr_2L$TajimaD, type="h", ylim=c(-3,4), main="Tajima's D Chr2L")
abline(h=-1.765, col="red")
abline(h=2.095, col="red")
plot(x=chr_2R$BIN_START, y=chr_2R$TajimaD, type="h", ylim=c(-3,4), main="Tajima's D Chr2R")
abline(h=-1.765, col="red")
abline(h=2.095, col="red")
plot(x=chr_3L$BIN_START, y=chr_3L$TajimaD, type="h", ylim=c(-3,4), main="Tajima's D Chr3L")
abline(h=-1.765, col="red")
abline(h=2.095, col="red")
plot(x=chr_3R$BIN_START, y=chr_3R$TajimaD, type="h",ylim=c(-3,4),  main="Tajima's D Chr3R")
abline(h=-1.765, col="red")
abline(h=2.095, col="red")
plot(x=chr_4$BIN_START, y=chr_4$TajimaD, type="h", ylim=c(-3,4), main="Tajima's D Chr4")
abline(h=-1.765, col="red")
abline(h=2.095, col="red")
plot(x=chr_X$BIN_START, y=chr_X$TajimaD, type="h", ylim=c(-3,4), main="Tajima's D ChrX")
abline(h=-1.765, col="red")
abline(h=2.095, col="red")



plt <- ggplot(data=chr_4, aes(x=chr_4$BIN_START, y=chr_4$TajimaD)) +
  geom_bar(stat='identity') +
  ggtitle("Tajima's D along Chromosome 4 in bins of 10kb") +
  xlab("Position along Chromosome 4") + ylab("Tajima's D in bins of 10kb")
plt
plt <- ggplot(data=chr_X, aes(x=chr_X$BIN_START, y=chr_X$TajimaD)) +
  geom_bar(stat='identity') +
  ggtitle("Tajima's D along Chromosome X in bins of 10kb") +
  xlab("Position along Chromosome X") + ylab("Tajima's D in bins of 10kb")
plt
plt <- ggplot(data=chr_3L, aes(x=chr_3L$BIN_START, y=chr_3L$TajimaD)) +
  geom_bar(stat='identity') +
  ggtitle("Tajima's D along Chromosome 3L in bins of 10kb") +
  xlab("Position along Chromosome 3L") + ylab("Tajima's D in bins of 10kb")
plt
plt <- ggplot(data=chr_3R, aes(x=chr_3R$BIN_START, y=chr_3R$TajimaD)) +
  geom_bar(stat='identity') +
  ggtitle("Tajima's D along Chromosome 3R in bins of 10kb") +
  xlab("Position along Chromosome 3R") + ylab("Tajima's D in bins of 10kb")
plt

plt <- ggplot(data=chr_2L, aes(x=chr_2L$BIN_START, y=chr_2L$TajimaD)) +
  geom_bar(stat='identity') +
  ggtitle("Tajima's D along Chromosome 2L in bins of 10kb") +
  xlab("Position along Chromosome 2L") + ylab("Tajima's D in bins of 10kb")
plt
plt <- ggplot(data=chr_2R, aes(x=chr_2R$BIN_START, y=chr_2R$TajimaD)) +
  geom_bar(stat='identity') +
  ggtitle("Tajima's D along Chromosome 2R in bins of 10kb") +
  xlab("Position along Chromosome 2R") + ylab("Tajima's D in bins of 10kb")
plt


##Plot density curves for TsD
den2L <- density(chr_2L$TajimaD, adjust=0.4)
den2R <- density(chr_2R$TajimaD, adjust=0.4)
den3L <- density(chr_3L$TajimaD, adjust=0.4)
den3R <- density(chr_3R$TajimaD, adjust=0.4)
den4 <- density(chr_4$TajimaD, adjust=0.2) #1kb window, not 10kb
denX <- density(chr_X$TajimaD, adjust=0.4)

plot(den2L, main="Tajima's D on Chr2L \n 10kb window")
plot(den2R, main="Tajima's D on Chr2R \n 10kb window")
plot(den3L, main="Tajima's D on Chr3L \n 10kb window")
plot(den3R, main="Tajima's D on Chr3R \n 10kb window")
plot(den4, main="Tajima's D on Chr4 \n 1kb window")
plot(denX, main="Tajima's D on ChrX \n 10kb window")