MCLAB <- Sys.getenv('MCLAB')

library(ggplot2)

infile <- paste(MCLAB,"/arbeitman/arbeitman_fru_network/reports_internal/motif_analysis/motif_distribution.csv",sep="")
indata <- read.csv(infile,header=TRUE)

mindfile <- paste(MCLAB,"/arbeitman/arbeitman_fru_network/reports_internal/motif_analysis/male_ind_motif_distribution.csv",sep="")
minddata <- read.csv(mindfile,header=TRUE)

mrepfile <- paste(MCLAB,"/arbeitman/arbeitman_fru_network/reports_internal/motif_analysis/male_rep_motif_distribution.csv",sep="")
mrepdata <- read.csv(mrepfile,header=TRUE)

findfile <- paste(MCLAB,"/arbeitman/arbeitman_fru_network/reports_internal/motif_analysis/female_ind_motif_distribution.csv",sep="")
finddata <- read.csv(findfile,header=TRUE)

frepfile <- paste(MCLAB,"/arbeitman/arbeitman_fru_network/reports_internal/motif_analysis/female_rep_motif_distribution.csv",sep="")
frepdata <- read.csv(frepfile,header=TRUE)

nindfile <- paste(MCLAB,"/arbeitman/arbeitman_fru_network/reports_internal/motif_analysis/null_male_induced_motif_distribution.csv",sep="")
ninddata <- read.csv(nindfile,header=TRUE)

nrepfile <- paste(MCLAB,"/arbeitman/arbeitman_fru_network/reports_internal/motif_analysis/null_male_repressed_motif_distribution.csv",sep="")
nrepdata <- read.csv(nrepfile,header=TRUE)


dis <- ggplot(data=indata,aes(x=motif_count,y=num_genes,fill=motif,alpha=0.2)) + geom_ribbon(aes(ymin=0,ymax=num_genes)) + ggtitle("Distirbution of FruM Motifs") + scale_alpha(guide=F) + theme(plot.background=element_blank(),panel.background=element_blank(),panel.grid=element_blank()) + scale_x_continuous(limits=c(0,20),breaks=seq(0,20,2))

mind   <- ggplot(data=minddata,aes(x=motif_count,y=num_genes,fill=motif,alpha=0.2)) + geom_ribbon(aes(ymin=0,ymax=num_genes)) + ggtitle("Distirbution of FruM Motifs\nMales Induced") + scale_alpha(guide=F) + theme(plot.background=element_blank(),panel.background=element_blank(),panel.grid=element_blank()) + scale_x_continuous(limits=c(0,20),breaks=seq(0,20,2))
mrep   <- ggplot(data=mrepdata,aes(x=motif_count,y=num_genes,fill=motif,alpha=0.2)) + geom_ribbon(aes(ymin=0,ymax=num_genes)) + ggtitle("Distirbution of FruM Motifs\nMales Repressed") + scale_alpha(guide=F) + theme(plot.background=element_blank(),panel.background=element_blank(),panel.grid=element_blank()) + scale_x_continuous(limits=c(0,20),breaks=seq(0,20,2))

find   <- ggplot(data=finddata,aes(x=motif_count,y=num_genes,fill=motif,alpha=0.2)) + geom_ribbon(aes(ymin=0,ymax=num_genes)) + ggtitle("Distirbution of FruM Motifs\nFemales Induced") + scale_alpha(guide=F) + theme(plot.background=element_blank(),panel.background=element_blank(),panel.grid=element_blank()) + scale_x_continuous(limits=c(0,20),breaks=seq(0,20,2))
frep   <- ggplot(data=frepdata,aes(x=motif_count,y=num_genes,fill=motif,alpha=0.2)) + geom_ribbon(aes(ymin=0,ymax=num_genes)) + ggtitle("Distirbution of FruM Motifs\nFemales Repressed") + scale_alpha(guide=F) + theme(plot.background=element_blank(),panel.background=element_blank(),panel.grid=element_blank()) + scale_x_continuous(limits=c(0,20),breaks=seq(0,20,2))

nind   <- ggplot(data=ninddata,aes(x=motif_count,y=num_genes,fill=motif,alpha=0.2)) + geom_ribbon(aes(ymin=0,ymax=num_genes)) + ggtitle("Distirbution of FruM Motifs\nMale Null Induced") + scale_alpha(guide=F) + theme(plot.background=element_blank(),panel.background=element_blank(),panel.grid=element_blank()) + scale_x_continuous(limits=c(0,20),breaks=seq(0,20,2))
nrep   <- ggplot(data=nrepdata,aes(x=motif_count,y=num_genes,fill=motif,alpha=0.2)) + geom_ribbon(aes(ymin=0,ymax=num_genes)) + ggtitle("Distirbution of FruM Motifs\nMale Null Repressed") + scale_alpha(guide=F) + theme(plot.background=element_blank(),panel.background=element_blank(),panel.grid=element_blank()) + scale_x_continuous(limits=c(0,20),breaks=seq(0,20,2))


ggsave(plot=dis,filename=paste(MCLAB,"/arbeitman/arbeitman_fru_network/reports_internal/motif_analysis/motif_distribution.svg",sep=""),width=8,height=4)
ggsave(plot=mind,filename=paste(MCLAB,"/arbeitman/arbeitman_fru_network/reports_internal/motif_analysis/male_ind_motif_distribution.svg",sep=""),width=8,height=4)
ggsave(plot=mrep,filename=paste(MCLAB,"/arbeitman/arbeitman_fru_network/reports_internal/motif_analysis/male_rep_motif_distribution.svg",sep=""),width=8,height=4)
ggsave(plot=find,filename=paste(MCLAB,"/arbeitman/arbeitman_fru_network/reports_internal/motif_analysis/female_ind_motif_distribution.svg",sep=""),width=8,height=4)
ggsave(plot=frep,filename=paste(MCLAB,"/arbeitman/arbeitman_fru_network/reports_internal/motif_analysis/female_rep_motif_distribution.svg",sep=""),width=8,height=4)
ggsave(plot=nind,filename=paste(MCLAB,"/arbeitman/arbeitman_fru_network/reports_internal/motif_analysis/null_male_ind_motif_distribution_v2.svg",sep=""),width=8,height=4)
ggsave(plot=nrep,filename=paste(MCLAB,"/arbeitman/arbeitman_fru_network/reports_internal/motif_analysis/null_male_rep_motif_distribution_v2.svg",sep=""),width=8,height=4)
