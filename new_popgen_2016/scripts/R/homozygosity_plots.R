#This script plots the coefficient of inbreeding. 
#We are interested in the different ways to filter based on individuals retained and level of homozygosity.

#Import the het file
het <- read.csv("/home/fnew/dsim/dsim_het.het", header=TRUE)
het2 <- read.csv("/home/fnew/dsim/resid_het_by_indiv.het", header=TRUE)


#I want the plots in one window
par(mfrow=c(1,2))


plt1 <- plot(het2$F, main="Coefficient of inbreeding\nNumber of individuals retained", ylab="F", col="#0066CC", pch=1, ylim=c(0,1))


#abline(h=0.53189) # -25 lines
#abline(h=0.62590) # -50 lines 
#abline(h=0.71039) # -100 lines
#abline(h=0.83105) # -150 lines, too many, not going to graph
#text(x=175, y=0.56, labels="25 (53%)")
#text(x=175, y=0.65, labels="50 (63%)")
#text(x=173, y=0.74, labels="100 (71%)")


#I want to do individuals retained now
abline(h=0.72296) #75 retained
abline(h=0.68622) #100 retained
abline(h=0.64630) #125 retained
abline(h=0.57344) #150 retained

text(x=175, y=0.59, labels="150 (57%)",cex=0.7)
text(x=175, y=0.66, labels="125 (65%)",cex=0.7)
text(x=175, y=0.70, labels="100 (69%)",cex=0.7)
text(x=175, y=0.74, labels="75 (72%)",cex=0.7)



plt2 <- plot(het2$F, main="Coefficient of inbreeding\nLevel of homozygosity", ylab="F", pch=1, col="#0066CC", ylim=c(0,1))

abline(h=0.50)
abline(h=0.60)
abline(h=0.70)
abline(h=0.80)

text(x=178, y=0.53, labels="50%(165)",cex=0.7)
text(x=178, y=0.63, labels="60%(140)",cex=0.7)
text(x=178, y=0.73, labels="70% (95)",cex=0.7)
text(x=177, y=0.83, labels="80% (45)",cex=0.7)


