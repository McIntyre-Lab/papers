
#project directory staph\seuong-ho\output

match<-read.csv("match.csv")
no_match<-read.csv("no_match.csv")

den_match<-density(match[,2])

den_nomatch<-density(no_match[,2])

snp<-read.csv("for_plotting_snps.csv")

library(vioplot)

x1 <- snp$sum_no_match[snp$group=="CC"]
x2 <- snp$sum_no_match[snp$group=="Pair"]
x3 <- snp$sum_no_match[snp$group=="PFGE"]
x4 <- snp$sum_no_match[snp$group=="ST"]

vioplot(x1, x2, x3,x4, names=c("CC", "Pair", "PFGE", "ST"),
        col="gold")

title("Violin Plots of Miles Per Gallon")


snp_pair<-snps[,(snps[,1]==1)]

den_pair<-density(snps)

plot(den_match$x,den_match$y,type="l", xlim=c(0,2000))

plot(den_nomatch$x,den_nomatch$y,type="l", ylim=c(0,0.010))
lines(den_match$x,den_match$y,type="l")