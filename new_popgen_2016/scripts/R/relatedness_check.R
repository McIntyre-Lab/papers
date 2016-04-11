# Relatedness

rel <- read.table("/home/fnew/ufgi_share/SHARE/McIntyre_Lab/ethanol/Sim_Pop_Gen/output/filter_10permiss_nolab_nomel_relat.relatedness", header=TRUE)

rel2 <- rel[order(rel$RELATEDNESS_AJK),]
summary(rel2)

den_rel <- density(rel2$RELATEDNESS_AJK)

plot(den_rel)


#h12
h2l <- read.table("/home/fnew/ufgi_share/SHARE/McIntyre_Lab/ethanol/Sim_Pop_Gen/output/H12/chr2L_H12_output.txt", header=FALSE)
summary(h2l$V9)

h2r <- read.table("/home/fnew/ufgi_share/SHARE/McIntyre_Lab/ethanol/Sim_Pop_Gen/output/H12/chr2R_H12_output.txt", header=FALSE)
summary(h2r$V9)
h3l <- read.table("/home/fnew/ufgi_share/SHARE/McIntyre_Lab/ethanol/Sim_Pop_Gen/output/H12/chr3L_H12_output.txt", header=FALSE)
summary(h3l$V9)
h3r <- read.table("/home/fnew/ufgi_share/SHARE/McIntyre_Lab/ethanol/Sim_Pop_Gen/output/H12/chr3R_H12_output.txt", header=FALSE)
summary(h3r$V9)
hx <- read.table("/home/fnew/ufgi_share/SHARE/McIntyre_Lab/ethanol/Sim_Pop_Gen/output/H12/chrX_H12_output.txt", header=FALSE)
summary(hx$V9)
h4 <- read.table("/home/fnew/ufgi_share/SHARE/McIntyre_Lab/ethanol/Sim_Pop_Gen/output/H12/chr4_H12_output.txt", header=FALSE)
summary(h4$V9)

all <- read.table("/home/fnew/ufgi_share/SHARE/McIntyre_Lab/ethanol/Sim_Pop_Gen/output/H12/all_chromosomes_H12_output.txt", header=FALSE)
summary(all$V9)
