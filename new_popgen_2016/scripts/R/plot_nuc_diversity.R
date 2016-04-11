#plot pi in 1kb, 10kb, and 100kb windows

pi1 <- read.table("/Users/felicianew/Documents/filt_nolab_pi_1kb.windowed.pi", header=TRUE)
pi10 <- read.table("/Users/felicianew/Documents/filt_nolab_pi_10kb.windowed.pi", header=TRUE)
pi100 <- read.table("/Users/felicianew/Documents/filt_nolab_pi_100kb.windowed.pi", header=TRUE)



#Whole genome
plot(x=pi100$BIN_START, y=pi100$PI, type="l")

plt <- ggplot(data=pi100, aes(x=pi100$BIN_START, y=pi100$PI, fill=pi100$CHROM)) +
  geom_bar(stat='identity')+ scale_fill_discrete()+
  ggtitle("Nucleotide diversity Whole Genome") +
  xlab("Position along the genome") + ylab("Nucleotide Diversity \n 100kb window") +
  facet_grid(. ~CHROM) 
plt


#Subset the 10kb data
chr2L <- subset(pi10, pi10$CHROM=="2L")
chr2R <- subset(pi10, pi10$CHROM=="2R")
chr3L <- subset(pi10, pi10$CHROM=="3L")
chr3R <- subset(pi10, pi10$CHROM=="3R")
chr4 <- subset(pi10, pi10$CHROM=="4")
chrX <- subset(pi10, pi10$CHROM=="X")

par(mfrow=c(3,2))
plot(x=chr2L$BIN_START, y=chr2L$PI, type="l", main="Nucleotide Diversity on 2L" , xlab="Position in 10kb windows", ylab="pi", ylim=c(0,0.02))
plot(x=chr2R$BIN_START, y=chr2R$PI, type="l", main="Nucleotide Diversity on 2R" , xlab="Position in 10kb windows", ylab="pi", ylim=c(0,0.02))
plot(x=chr3L$BIN_START, y=chr3L$PI, type="l", main="Nucleotide Diversity on 3L" , xlab="Position in 10kb windows", ylab="pi", ylim=c(0,0.02))
plot(x=chr3R$BIN_START, y=chr3R$PI, type="l", main="Nucleotide Diversity on 3R" , xlab="Position in 10kb windows", ylab="pi", ylim=c(0,0.02))
plot(x=chr4$BIN_START, y=chr4$PI, type="l", main="Nucleotide Diversity on 4" , xlab="Position in 10kb windows", ylab="pi")
plot(x=chrX$BIN_START, y=chrX$PI, type="l", main="Nucleotide Diversity on X" , xlab="Position in 10kb windows", ylab="pi", ylim=c(0,0.02))


# subset the 100kb data and plot
chr2L <- subset(pi100, pi100$CHROM=="2L")
chr2R <- subset(pi100, pi100$CHROM=="2R")
chr3L <- subset(pi100, pi100$CHROM=="3L")
chr3R <- subset(pi100, pi100$CHROM=="3R")
chr4 <- subset(pi100, pi100$CHROM=="4")
chrX <- subset(pi100, pi100$CHROM=="X")

par(mfrow=c(3,2))
plot(x=chr2L$BIN_START, y=chr2L$PI, type="l", main="Nucleotide Diversity on 2L" , xlab="Position in 100kb windows", ylab="pi", ylim=c(0,0.013))
plot(x=chr2R$BIN_START, y=chr2R$PI, type="l", main="Nucleotide Diversity on 2R" , xlab="Position in 100kb windows", ylab="pi", ylim=c(0,0.013))
plot(x=chr3L$BIN_START, y=chr3L$PI, type="l", main="Nucleotide Diversity on 3L" , xlab="Position in 100kb windows", ylab="pi", ylim=c(0,0.013))
plot(x=chr3R$BIN_START, y=chr3R$PI, type="l", main="Nucleotide Diversity on 3R" , xlab="Position in 100kb windows", ylab="pi", ylim=c(0,0.013))
plot(x=chr4$BIN_START, y=chr4$PI, type="l", main="Nucleotide Diversity on 4" , xlab="Position in 100kb windows", ylab="pi")
plot(x=chrX$BIN_START, y=chrX$PI, type="l", main="Nucleotide Diversity on X" , xlab="Position in 100kb windows", ylab="pi", ylim=c(0,0.013))


# check chrom 4 in the smallest windows
chr4 <- subset(pi1, pi1$CHROM=="4")
plot(x=chr4$BIN_START, y=chr4$PI, type="l", main="Nucleotide Diversity on 4" , xlab="Position in 1kb windows", ylab="pi")


## By features
piG <- read.table("/home/fnew/dsim/pi/filtered_vcf_genes_pi_10kb.windowed.pi", header=TRUE)
piT <- read.table("/home/fnew/dsim/pi/filtered_vcf_transcript_pi_10kb.windowed.pi", header=TRUE)
piI <- read.table("/home/fnew/dsim/pi/filtered_vcf_introns_pi_10kn.windowed.pi", header=TRUE)
piE <- read.table("/home/fnew/dsim/pi/filtered_vcf_exons_pi_10kb.windowed.pi", header=TRUE)
pi3 <- read.table("/home/fnew/dsim/pi/filtered_vcf_3utr_pi_10kb.windowed.pi", header=TRUE)
pi5 <- read.table("/home/fnew/dsim/pi/filtered_vcf_5utr_pi_10kb.windowed.pi", header=TRUE)


denpiG <-density(piG$PI, adjust=0.5)
denpiT <- density(piT$PI, adjust=0.5)
denpiI<-density(piI$PI, adjust=0.5)
denpiE <- density(piE$PI, adjust=0.5)
denpi3 <-density(pi3$PI, adjust=0.5)
denpi5 <- density(pi5$PI, adjust=0.5)

par(mfrow=c(3,2))
plot(denpiG, main="Nucleotide diversity \n Genes", xlab="Nucleotide Diversity", xlim=c(0,0.02) )
plot(denpiT, main="Nucleotide diversity \n Transcripts", xlim=c(0,0.02))
plot(denpiE, main="Nucleotide diversity \n Exons",xlab="Nucleotide Diversity", ylim=c(0,350), xlim=c(0,0.02))
plot(denpiI, main="Nucleotide diversity \n Introns", xlab="Nucleotide Diversity",ylim=c(0,350), xlim=c(0,0.02))
plot(denpi3, main="Nucleotide diversity \n 3'UTRs", xlab="Nucleotide Diversity",ylim=c(0,1800))
plot(denpi5, main="Nucleotide diversity \n 5'UTRs",xlab="Nucleotide Diversity", ylim=c(0,1800))
