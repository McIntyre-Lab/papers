library(Vennerable)

MCLAB <- Sys.getenv('MCLAB')

setwd(paste(MCLAB,"/arbeitman_fru_network/reports_external/venn_diagram",sep=""))

male   <- read.csv(paste(MCLAB,"/arbeitman_fru_network/reports_external/flag_x_induced_repressed_male.csv",sep=""),header=TRUE)
null   <- read.csv(paste(MCLAB,"/arbeitman_fru_network/reports_external/flag_x_induced_repressed_null_male_v2.csv",sep=""),header=TRUE)
female <- read.csv(paste(MCLAB,"/arbeitman_fru_network/reports_external/flag_x_induced_repressed_female.csv",sep=""),header=TRUE)

male.ind.lst <- list(FruMA = male[male$flag_a_ind == 1,]$primary_fbgn, FruMB = male[male$flag_b_ind == 1,]$primary_fbgn, FruMC = male[male$flag_c_ind == 1,]$primary_fbgn)
vm.male.ind <- Venn(male.ind.lst)
svg("male_ind.svg")
plot(vm.male.ind)
dev.off()

male.rep.lst <- list(FruMA = male[male$flag_a_rep == 1,]$primary_fbgn, FruMB = male[male$flag_b_rep == 1,]$primary_fbgn, FruMC = male[male$flag_c_rep == 1,]$primary_fbgn)
vm.male.rep <- Venn(male.rep.lst)
svg("male_rep.svg")
plot(vm.male.rep)
dev.off()

female.ind.lst <- list(FruMA = female[female$flag_a_ind == 1,]$primary_fbgn, FruMB = female[female$flag_b_ind == 1,]$primary_fbgn, FruMC = female[female$flag_c_ind == 1,]$primary_fbgn)
vm.female.ind <- Venn(female.ind.lst)
svg("female_ind.svg")
plot(vm.female.ind)
dev.off()

female.rep.lst <- list(FruMA = female[female$flag_a_rep == 1,]$primary_fbgn, FruMB = female[female$flag_b_rep == 1,]$primary_fbgn, FruMC = female[female$flag_c_rep == 1,]$primary_fbgn)
vm.female.rep <- Venn(female.rep.lst)
svg("female_rep.svg")
plot(vm.female.rep)
dev.off()

null.lst <- list(induced = null[null$flag_induced == 1,]$primary_fbgn, repressed = null[null$flag_repressed == 1,]$primary_fbgn)
vm.null <- Venn(null.lst)
svg("null_v3.svg")
plot(vm.null)
dev.off()


ind.lst <- list( male.FruMA = male[male$flag_a_ind == 1,]$primary_fbgn, male.FruMB = male[male$flag_b_ind == 1,]$primary_fbgn, male.FruMC = male[male$flag_c_ind == 1,]$primary_fbgn, female.FruMA = female[female$flag_a_ind == 1,]$primary_fbgn, female.FruMB = female[female$flag_b_ind == 1,]$primary_fbgn, female.FruMC = female[female$flag_c_ind == 1,]$primary_fbgn, null = null[null$flag_induced == 1,]$primary_fbgn)
vm.ind <- Venn(ind.lst)
svg("male_female_null_ind_v3.svg")
plot(vm.ind,type="battle")
dev.off()

svg("male_female_ind.svg")
plot(vm.ind[,c("male.FruMA", "male.FruMB", "male.FruMC", "female.FruMA", "female.FruMB", "female.FruMC")],type="battle")
dev.off()

svg("male_female_FruMA_ind.svg")
plot(vm.ind[,c("male.FruMA","female.FruMA")],type="circles")
dev.off()

svg("male_female_FruMB_ind.svg")
plot(vm.ind[,c("male.FruMB","female.FruMB")],type="circles")
dev.off()

svg("male_female_FruMC_ind.svg")
plot(vm.ind[,c("male.FruMC","female.FruMC")],type="circles")
dev.off()

svg("male_null_ind_v3.svg")
plot(vm.ind[,c("male.FruMA", "male.FruMB", "male.FruMC", "null")],type="ellipses")
dev.off()

svg("female_null_ind_v3.svg")
plot(vm.ind[,c("female.FruMA", "female.FruMB", "female.FruMC", "null")],type="ellipses")
dev.off()

rep.lst <- list( male.FruMA = male[male$flag_a_rep == 1,]$primary_fbgn, male.FruMB = male[male$flag_b_rep == 1,]$primary_fbgn, male.FruMC = male[male$flag_c_rep == 1,]$primary_fbgn, female.FruMA = female[female$flag_a_rep == 1,]$primary_fbgn, female.FruMB = female[female$flag_b_rep == 1,]$primary_fbgn, female.FruMC = female[female$flag_c_rep == 1,]$primary_fbgn, null = null[null$flag_repressed == 1,]$primary_fbgn)
vm.rep <- Venn(rep.lst)
svg("male_female_null_rep_v3.svg")
plot(vm.rep,type="battle")
dev.off()

svg("male_female_rep.svg")
plot(vm.rep[,c("male.FruMA", "male.FruMB", "male.FruMC", "female.FruMA", "female.FruMB", "female.FruMC")],type="battle")
dev.off()

svg("male_female_FruMA_rep.svg")
plot(vm.rep[,c("male.FruMA","female.FruMA")],type="circles")
dev.off()

svg("male_female_FruMB_rep.svg")
plot(vm.rep[,c("male.FruMB","female.FruMB")],type="circles")
dev.off()

svg("male_female_FruMC_rep.svg")
plot(vm.rep[,c("male.FruMC","female.FruMC")],type="circles")
dev.off()

svg("male_null_rep_v3.svg")
plot(vm.rep[,c("male.FruMA", "male.FruMB", "male.FruMC", "null")],type="ellipses")
dev.off()

svg("female_null_rep_v3.svg")
plot(vm.rep[,c("female.FruMA", "female.FruMB", "female.FruMC", "null")],type="ellipses")
dev.off()



ind2.lst <- list(male.ind = male[male$flag_induced == 1,]$primary_fbgn, female.ind = female[female$flag_induced == 1,]$primary_fbgn, null.ind=null[null$flag_induced == 1,]$primary_fbgn)
vm.ind2 <- Venn(ind2.lst)

svg("male_female_ind.svg")
plot(vm.ind2[,c("male.ind","female.ind")])
dev.off()

svg("male_null_ind_v3.svg")
plot(vm.ind2[,c("male.ind","null.ind")])
dev.off()

svg("female_null_ind_v3.svg")
plot(vm.ind2[,c("female.ind","null.ind")])
dev.off()

rep2.lst <- list(male.rep = male[male$flag_repressed == 1,]$primary_fbgn, female.rep = male[female$flag_repressed == 1,]$primary_fbgn, null.rep=null[null$flag_repressed == 1,]$primary_fbgn)
vm.rep2 <- Venn(rep2.lst)
svg("male_female_rep.svg")
plot(vm.rep2[,c("male.rep","female.rep")])
dev.off()

svg("male_null_rep_v3.svg")
plot(vm.rep2[,c("male.rep","null.rep")])
dev.off()

svg("female_null_rep_v3.svg")
plot(vm.rep2[,c("female.rep","null.rep")])
dev.off()

