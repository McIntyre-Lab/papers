#Plot SNP densities from 1kb, 10kb, and 100kb windows

snp1 <- read.table('/home/fnew/dsim/snp_den/filter_snp_den_1kb.snpden', header=TRUE)
snp10 <- read.table('/home/fnew/dsim/snp_den/filter_snp_den_10kb.snpden', header=TRUE)
snp100 <- read.table('/home/fnew/dsim/snp_den/filter_snp_den_100kb.snpden', header=TRUE)


#Plot the whole genome with the 100kb snp density data
plot(x=snp100$BIN_START, y=snp100$SNP_COUNT, type="h")
dengen <- density(snp100$VARIANTS.KB, adjust=0.5)
plot(dengen, main="SNP density across the genome \n 100kb windows")


plt <- ggplot(data=snp100, aes(x=snp100$BIN_START, y=snp100$VARIANTS.KB, fill=snp100$CHROM)) +
  geom_bar(stat='identity')+ scale_fill_discrete()+
  ggtitle("SNP density Whole Genome in bins of 100KB") +
  xlab("Position along the genome") + ylab("SNP Density in bins of 100kb") +
  facet_grid(. ~CHROM) 
plt

#Make boxplots?
par(mfrow=c(1,1))
plot(x=snp100$CHROM, y=snp100$SNP_COUNT, type="b", main="SNPs in 100kb windows")

#Subset the chroms
chr2L <- subset(snp100, snp100$CHROM=="2L")
chr2R <- subset(snp100, snp100$CHROM=="2R")
chr3L <- subset(snp100, snp100$CHROM=="3L")
chr3R <- subset(snp100, snp100$CHROM=="3R")
chr4 <- subset(snp100, snp100$CHROM=="4")
chrX <- subset(snp100, snp100$CHROM=="X")

plot(x=chr2L$BIN_START, y=chr2L$VARIANTS.KB, type='h')
plot(x=chr2L$BIN_START, y=chr2L$SNP_COUNT, type='h')
plot(x=chr2R$BIN_START, y=chr2R$SNP_COUNT, type='h')
plot(x=chr3L$BIN_START, y=chr3L$SNP_COUNT, type='h')
plot(x=chr3R$BIN_START, y=chr3R$SNP_COUNT, type='h')
plot(x=chr4$BIN_START, y=chr4$SNP_COUNT, type='h')
plot(x=chrX$BIN_START, y=chrX$SNP_COUNT, type='h')



# Density curves for the chromosomes
den2L <- density(chr2L$SNP_COUNT, adjust=0.5)
den2R <- density(chr2R$SNP_COUNT, adjust=0.5)
den3L <- density(chr3L$SNP_COUNT, adjust=0.5)
den3R <- density(chr3R$SNP_COUNT, adjust=0.5)
den4 <- density(chr4$SNP_COUNT, adjust=0.5) #10kb window, not 100kb
denX <- density(chrX$SNP_COUNT, adjust=0.5)

par(mfrow=c(3,2))
par(mfrow=c(1,1))
par(mfrow=c(6,1))
plot(den2L, main="Number of SNPs within 100kb Windows \n Chr2L", xlab="Number of SNPs", xlim=c(0,7000))
points(x=5009, y=0, pch=17, col="red", cex=2) #mean
points(x=5519.5, y=0, pch="*", col="blue", cex=4) #median

plot(den2R, main="Number of SNPs within 100kb Windows \n Chr2R", xlab="Number of SNPs", xlim=c(0,7000))
points(x=4847, y=0, pch=17, col="red", cex=2) #mean
points(x=5278, y=0, pch="*", col="blue", cex=4) #median

plot(den3L, main="Number of SNPs within 100kb Windows \n Chr3L", xlab="Number of SNPs", xlim=c(0,7000))
points(x=5087.55, y=0, pch=17, col="red", cex=2) #mean
points(x=5433, y=0, pch="*", col="blue", cex=4) #median

plot(den3R, main="Number of SNPs within 100kb Windows \n Chr3R", xlab="Number of SNPs", xlim=c(0,7000))
points(x=5460.4, y=0, pch=17, col="red", cex=2) #mean
points(x=5491, y=0, pch="*", col="blue", cex=4) #median

plot(den4, main="Number of SNPs within 100kb Windows \n Chr4", xlab="Number of SNPs")
points(x=2916.6, y=0, pch=17, col="red", cex=2) #mean
points(x=3119, y=0, pch="*", col="blue", cex=4) #median

plot(denX, main="Number of SNPs within 100kb Windows \n ChrX", xlab="Number of SNPs", xlim=c(0,7000))
points(x=4859.78, y=0, pch=17, col="red", cex=2) #mean
points(x=4937, y=0, pch="*", col="blue", cex=4) #median


#get the average and median for each chrom
avg2L <- mean(chr2L$SNP_COUNT)
med2L <- median(chr2L$SNP_COUNT)

avg2R <- mean(chr2R$SNP_COUNT)
med2R <- median(chr2R$SNP_COUNT)

avg3L <- mean(chr3L$SNP_COUNT)
med3L <- median(chr3L$SNP_COUNT)

avg3R <- mean(chr3R$SNP_COUNT)
med3R <- median(chr3R$SNP_COUNT)

avg4 <- mean(chr4$SNP_COUNT)
med4 <- median(chr4$SNP_COUNT)

avgX <- mean(chrX$SNP_COUNT)
medX <- median(chrX$SNP_COUNT)

#Subset the chroms for 10kb data
chr2L <- subset(snp10, snp10$CHROM=="2L")
chr2R <- subset(snp10, snp10$CHROM=="2R")
chr3L <- subset(snp10, snp10$CHROM=="3L")
chr3R <- subset(snp10, snp10$CHROM=="3R")
chr4 <- subset(snp10, snp10$CHROM=="4")
chrX <- subset(snp10, snp10$CHROM=="X")

plot(x=chr2L$BIN_START, y=chr2L$SNP_COUNT, type='h')
plot(x=chr2R$BIN_START, y=chr2R$SNP_COUNT, type='h')
plot(x=chr3L$BIN_START, y=chr3L$SNP_COUNT, type='l')
plot(x=chr3R$BIN_START, y=chr3R$SNP_COUNT, type='h')
plot(x=chr4$BIN_START, y=chr4$SNP_COUNT, type='h')
plot(x=chrX$BIN_START, y=chrX$SNP_COUNT, type='h')

#Just do chr4 with the 1kb data
chr4 <- subset(snp1, snp1$CHROM=="4")
plot(x=chr4$BIN_START, y=chr4$VARIANTS.KB, type='h')



##Plot the features
exon <- read.table('/home/fnew/dsim/snp_den/filtered_vcf_exons_snpden_10kb.snpden', header=TRUE)
intron <-read.table('/home/fnew/dsim/snp_den/filtered_vcf_introns_snpden_10kn.snpden', header=TRUE)
gene <-read.table('/home/fnew/dsim/snp_den/filtered_vcf_genes_snpden_10kb.snpden', header=TRUE)
trans <-read.table('/home/fnew/dsim/snp_den/filtered_vcf_transcript_snpden_10kb.snpden', header=TRUE)
utr3 <-read.table('/home/fnew/dsim/snp_den/filtered_vcf_3utr_snpden_10kb.snpden', header=TRUE)
utr5 <-read.table('/home/fnew/dsim/snp_den/filtered_vcf_5utr_snpden_10kb.snpden', header=TRUE)


par(mfrow=c(3,2))

plot(x=gene$CHROM, y=gene$SNP_COUNT, type="b", main="SNPs per 10kb within Genes", ylim=c(0,900))
plot(x=trans$CHROM, y=trans$SNP_COUNT, type="b", main="SNPs per 10kb within Transcripts", ylim=c(0,900))
plot(x=exon$CHROM, y=exon$SNP_COUNT, type="b", main="SNPs per 10kb within Exons", ylim=c(0,900))
plot(x=intron$CHROM, y=intron$SNP_COUNT, type="b", main="SNPs per 10kb within Introns", ylim=c(0,900))
plot(x=utr3$CHROM, y=utr3$SNP_COUNT, type="b", main="SNPs per 10kb within 3'UTRs", ylim=c(0,100))
plot(x=utr5$CHROM, y=utr5$SNP_COUNT, type="b", main="SNPs per 10kb within 5'UTRS", ylim=c(0,100))

dene <- density(exon$SNP_COUNT,adjust=0.3)
deng <- density(gene$SNP_COUNT, adjust=0.3)
dent <- density(trans$SNP_COUNT, adjust=0.3)
deni <- density(intron$SNP_COUNT, adjust=0.3)
den3 <- density(utr3$SNP_COUNT)
den5 <- density(utr5$SNP_COUNT)

plot(deng, main="Number of SNPs within 10KB\nGenes", xlab="Number of SNPs per 10KB", ylim=c(0,0.017), xlim=c(0,800))
points(x=5009, y=-1, pch=17, col="red")
plot(dent,  main="Number of SNPs within 10KB\nTranscripts", xlab="Number of SNPs per 10KB", ylim=c(0,0.017), xlim=c(0,800))
plot(dene,  main="Number of SNPs within 10KB\nExons", xlab="Number of SNPs per 10KB", ylim=c(0,0.017), xlim=c(0,800))
plot(deni,  main="Number of SNPs within 10KB\nIntrons", xlab="Number of SNPs per 10KB", ylim=c(0,0.017), xlim=c(0,800))
plot(den3,  main="Number of SNPs within 10KB\n3'UTRs", xlab="Number of SNPs per 10KB", ylim=c(0,0.06), xlim=c(0,800))
plot(den5,  main="Number of SNPs within 10KB\n5'UTRs", xlab="Number of SNPs per 10KB", ylim=c(0,0.06), xlim=c(0,800))
