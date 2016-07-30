## Proc rank done in SAS to rank distances. Then means taken of the LD for each 'bin'
## Plotting these ~smoothed curves on one graph
## This is part 1 to the script "plot_figure2_for_manuscript_2.R

rank_2l<- read.table("/home/fnew/ufgi_share/SHARE/McIntyre_Lab/ethanol/Sim_Pop_Gen/output/LD/chr2L_LD_rank_means_1000.txt", header=TRUE)
rank_2r<- read.table("/home/fnew/ufgi_share/SHARE/McIntyre_Lab/ethanol/Sim_Pop_Gen/output/LD/chr2R_LD_rank_means_1000.txt", header=TRUE)
rank_3l<- read.table("/home/fnew/ufgi_share/SHARE/McIntyre_Lab/ethanol/Sim_Pop_Gen/output/LD/chr3L_LD_rank_means_1000.txt", header=TRUE)
rank_3r<- read.table("/home/fnew/ufgi_share/SHARE/McIntyre_Lab/ethanol/Sim_Pop_Gen/output/LD/chr3R_LD_rank_means_1000.txt", header=TRUE)
rank_4<- read.table("/home/fnew/ufgi_share/SHARE/McIntyre_Lab/ethanol/Sim_Pop_Gen/output/LD/chr4_LD_rank_means_100.txt", header=TRUE)
rank_x<- read.table("/home/fnew/ufgi_share/SHARE/McIntyre_Lab/ethanol/Sim_Pop_Gen/output/LD/chrX_LD_rank_means_1000.txt", header=TRUE)


par(mar=c(5.1,4.1,4.1,7.1), xpd=FALSE)
plot(x=rank_2l$median, rank_2l$mean, type="l", lwd=3, col="purple", xlab="BP", ylab="LD", main="LD Decay", ylim=c(0,0.35), xlim=c(0,1200))
  with(rank_2r, lines(rank_2r$median, rank_2r$mean, type="l", lwd=3, col="dark orange"))
  with(rank_3l, lines(rank_3l$median, rank_3l$mean, type="l", lwd=3, col="dark green"))
  with(rank_3r, lines(rank_3r$median, rank_3r$mean, type="l", lwd=3, col="dark blue"))
  with(rank_x, lines(rank_x$median, rank_x$mean, type="l", lwd=3, col="black"))
  with(rank_4, lines(rank_4$median, rank_4$mean, type="l", lwd=3, col="red"))
  legend("topright",  fill=c("purple","dark orange","dark green","dark blue","black", "red"), legend=c("2L","2R","3L","3R","X","4"))
  
  
## Check what rank 100 looks like
rank_2l<- read.table("/home/fnew/ufgi_share/SHARE/McIntyre_Lab/ethanol/Sim_Pop_Gen/output/LD/chr2L_LD_rank_means_100.txt", header=TRUE)
rank_2r<- read.table("/home/fnew/ufgi_share/SHARE/McIntyre_Lab/ethanol/Sim_Pop_Gen/output/LD/chr2R_LD_rank_means_100.txt", header=TRUE)
rank_3l<- read.table("/home/fnew/ufgi_share/SHARE/McIntyre_Lab/ethanol/Sim_Pop_Gen/output/LD/chr3L_LD_rank_means_100.txt", header=TRUE)
rank_3r<- read.table("/home/fnew/ufgi_share/SHARE/McIntyre_Lab/ethanol/Sim_Pop_Gen/output/LD/chr3R_LD_rank_means_100.txt", header=TRUE)
rank_4<- read.table("/home/fnew/ufgi_share/SHARE/McIntyre_Lab/ethanol/Sim_Pop_Gen/output/LD/chr4_LD_rank_means_100.txt", header=TRUE)
rank_x<- read.table("/home/fnew/ufgi_share/SHARE/McIntyre_Lab/ethanol/Sim_Pop_Gen/output/LD/chrX_LD_rank_means_100.txt", header=TRUE)
  
  
  par(mar=c(5.1,4.1,4.1,7.1), xpd=FALSE)
  plot(x=rank_2l$median, y=rank_2l$mean, type="l", lwd=3, col="purple", xlab="BP", ylab="LD", ylim=c(0,0.6), xlim=c(0,750))
  with(rank_2r, lines(rank_2r$median, rank_2r$mean, type="l", lwd=3, col="dark orange"))
  with(rank_3l, lines(rank_3l$median, rank_3l$mean, type="l", lwd=3, col="dark green"))
  with(rank_3r, lines(rank_3r$median, rank_3r$mean, type="l", lwd=3, col="dark blue"))
  with(rank_x, lines(rank_x$median, rank_x$mean, type="l", lwd=3, col="black"))
  with(rank_4, lines(rank_4$median, rank_4$mean, type="l", lwd=3, col="red"))
  legend("topright",  fill=c("purple","dark orange","dark green","dark blue","black", "red"), legend=c("2L","2R","3L","3R","X","4"))
  
  
  
## Try smoothing the 1000 bins lines and plotting
  
par(mar=c(5.1,4.1,4.1,7.1), xpd=FALSE)
plot(loess.smooth(x=rank_2l$median, rank_2l$mean, span=0.5), type="l", lwd=3, col="purple", xlab="BP", ylab="LD", main="LD Decay", ylim=c(0,0.35), xlim=c(0,1200))
  with(rank_2r, lines(loess.smooth(rank_2r$median, rank_2r$mean, span=0.5), type="l", lwd=3, col="dark orange"))
  with(rank_3l, lines(loess.smooth(rank_3l$median, rank_3l$mean, span=0.5), type="l", lwd=3, col="dark green"))
  with(rank_3r, lines(loess.smooth(rank_3r$median, rank_3r$mean, span=0.5), type="l", lwd=3, col="dark blue"))
  with(rank_x, lines(loess.smooth(rank_x$median, rank_x$mean, span=0.5), type="l", lwd=3, col="black"))
  with(rank_4, lines(loess.smooth(rank_4$median, rank_4$mean, span=0.5), type="l", lwd=3, col="red"))
  legend("topright",  fill=c("purple","dark orange","dark green","dark blue","black", "red"), legend=c("2L","2R","3L","3R","X","4"))
  
  