# Load Library
library(Vennerable)

# Where is mclab
mclab <- Sys.getenv("MCLAB")

# data munging
ctrl.F <- read.csv('/home/jfear/mclab/arbeitman/arbeitman_dsx/deg_output/control_F_bias_ctrl.tab', header=FALSE)
ctrl.M <- read.csv('/home/jfear/mclab/arbeitman/arbeitman_dsx/deg_output/control_M_bias_ctrl.tab', header=FALSE)
ctrl <- unique(rbind(ctrl.F,ctrl.M))
names(ctrl.F) <- 'Female Bias'
names(ctrl.M) <- 'Male Bias'
names(ctrl) <- 'Control'

dsxD <- read.csv('/home/jfear/mclab/arbeitman/arbeitman_dsx/deg_output/dsxDF_dsxD_bias.tab', header=FALSE)
dsxDwt <- read.csv('/home/jfear/mclab/arbeitman/arbeitman_dsx/deg_output/dsxDF_WT_bias.tab', header=FALSE)
dsxD.all <- unique(rbind(dsxD,dsxDwt))
names(dsxD) <- 'dsxD Bias'
names(dsxDwt) <- 'dsxD WT Female Bias'
names(dsxD.all) <- 'dsxD'

dsxNullF <- read.csv('/home/jfear/mclab/arbeitman/arbeitman_dsx/deg_output/dsxNullF_Null_bias.tab', header=FALSE)
dsxNullFwt <- read.csv('/home/jfear/mclab/arbeitman/arbeitman_dsx/deg_output/dsxNullF_WT_bias.tab', header=FALSE)
dsxNullF.all <- unique(rbind(dsxNullF,dsxNullFwt))
names(dsxNullF) <- 'dsxNullF Bias'
names(dsxNullFwt) <- 'dsxNull WT Bias'
names(dsxNullF.all) <- 'dsxNullF'

dsxNullM <- read.csv('/home/jfear/mclab/arbeitman/arbeitman_dsx/deg_output/dsxNullM_Null_bias.tab', header=FALSE)
dsxNullMwt <- read.csv('/home/jfear/mclab/arbeitman/arbeitman_dsx/deg_output/dsxNullM_WT_bias.tab', header=FALSE)
dsxNullM.all <- unique(rbind(dsxNullM,dsxNullMwt))
names(dsxNullM) <- 'dsxNullM Bias'
names(dsxNullMwt) <- 'dsxNull WT Bias'
names(dsxNullM.all) <- 'dsxNullM'


# Plot Venn Diagrams
odir <- paste0(mclab, "/cegs_sem_sd_paper/analysis_output/dsx/")

png(paste0(odir,'dsxd_summary_venn.png'))
plot(Venn(c(ctrl.F, ctrl.M, dsxD, dsxDwt)),type='ellipses')
dev.off()

png(paste0(odir,'dsxnullF_summary_venn.png'))
plot(Venn(c(ctrl.F, ctrl.M, dsxNullF, dsxNullFwt)),type='ellipses')
dev.off()

png(paste0(odir,'dsxnullM_summary_venn.png'))
plot(Venn(c(ctrl.F, ctrl.M, dsxNullM, dsxNullMwt)),type='ellipses')
dev.off()


png(paste0(odir,'dsxnullM_dsxnullF_summary_venn.png'))
plot(Venn(c(dsxNullF, dsxNullFwt, dsxNullM, dsxNullMwt)),type='ellipses')
dev.off()

png(paste0(odir,'full_null_summary_venn.png'))
plot(Venn(c(ctrl.F, ctrl.M,dsxNullF, dsxNullFwt, dsxNullM, dsxNullMwt)),type='battle')
dev.off()

png(paste0(odir,'dsxd_nullF_summary_venn.png'))
plot(Venn(c(dsxNullF.all, dsxD.all)))
dev.off()

