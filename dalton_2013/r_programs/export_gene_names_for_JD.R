
setwd("/home/jfear/mclab/arbeitman_fru_network/reports/venn_diagram")

male <- read.csv("/home/jfear/mclab/arbeitman_fru_network/reports/flag_x_induced_repressed_male.csv",header=TRUE)
attach(male)

fruMA.ind <- male[flag_a_ind == 1,]
fruMA.rep <- male[flag_a_rep == 1,]

fruMB.ind <- male[flag_b_ind == 1,]
fruMB.rep <- male[flag_b_rep == 1,]

fruMC.ind <- male[flag_c_ind == 1,]
fruMC.rep <- male[flag_c_rep == 1,]

induced <- male[flag_a_ind == 1 | flag_b_ind == 1 | flag_c_ind == 1,]
repressed <- male[flag_a_rep == 1 | flag_b_rep == 1 | flag_c_rep == 1,]

null.induced <- male[flag_a_ind == 0 & flag_b_ind == 0 & flag_c_ind == 0,]
null.repressed <- male[flag_a_rep == 0 & flag_b_rep == 0 & flag_c_rep == 0,]

null.ind.rep <- male[flag_a_rep == 0 & flag_b_rep == 0 & flag_c_rep == 0 & flag_a_ind == 0 & flag_b_ind == 0 & flag_c_ind == 0,]

write.csv(fruMA.ind, file="/home/jfear/tmp/male_fruMA_induced.csv", row.names=FALSE, quote=FALSE)
write.csv(fruMB.ind, file="/home/jfear/tmp/male_fruMB_induced.csv", row.names=FALSE, quote=FALSE)
write.csv(fruMC.ind, file="/home/jfear/tmp/male_fruMC_induced.csv", row.names=FALSE, quote=FALSE)

write.csv(fruMA.rep, file="/home/jfear/tmp/male_fruMA_repressed.csv", row.names=FALSE, quote=FALSE)
write.csv(fruMB.rep, file="/home/jfear/tmp/male_fruMB_repressed.csv", row.names=FALSE, quote=FALSE)
write.csv(fruMC.rep, file="/home/jfear/tmp/male_fruMC_repressed.csv", row.names=FALSE, quote=FALSE)

write.csv(induced, file="/home/jfear/tmp/male_total_induced.csv", row.names=FALSE, quote=FALSE)
write.csv(repressed, file="/home/jfear/tmp/male_total_repressed.csv", row.names=FALSE, quote=FALSE)

write.csv(null.induced, file="/home/jfear/tmp/male_total_not_induced.csv", row.names=FALSE, quote=FALSE)
write.csv(null.repressed, file="/home/jfear/tmp/male_total_not_repressed.csv", row.names=FALSE, quote=FALSE)

write.csv(null.ind.rep, file="/home/jfear/tmp/male_total_not_induced_or_repressed.csv", row.names=FALSE, quote=FALSE)

detach(male)

female <- read.csv("/home/jfear/mclab/arbeitman_fru_network/reports/flag_x_induced_repressed_female.csv",header=TRUE)

attach(female)

fruMA.ind <- female[flag_a_ind == 1,]
fruMA.rep <- female[flag_a_rep == 1,]

fruMB.ind <- female[flag_b_ind == 1,]
fruMB.rep <- female[flag_b_rep == 1,]

fruMC.ind <- female[flag_c_ind == 1,]
fruMC.rep <- female[flag_c_rep == 1,]

induced <- female[flag_a_ind == 1 | flag_b_ind == 1 | flag_c_ind == 1,]
repressed <- female[flag_a_rep == 1 | flag_b_rep == 1 | flag_c_rep == 1,]

null.induced <- female[flag_a_ind == 0 & flag_b_ind == 0 & flag_c_ind == 0,]
null.repressed <- female[flag_a_rep == 0 & flag_b_rep == 0 & flag_c_rep == 0,]

write.csv(fruMA.ind, file="/home/jfear/tmp/female_fruMA_induced.csv", row.names=FALSE, quote=FALSE)
write.csv(fruMB.ind, file="/home/jfear/tmp/female_fruMB_induced.csv", row.names=FALSE, quote=FALSE)
write.csv(fruMC.ind, file="/home/jfear/tmp/female_fruMC_induced.csv", row.names=FALSE, quote=FALSE)

write.csv(fruMA.rep, file="/home/jfear/tmp/female_fruMA_repressed.csv", row.names=FALSE, quote=FALSE)
write.csv(fruMB.rep, file="/home/jfear/tmp/female_fruMB_repressed.csv", row.names=FALSE, quote=FALSE)
write.csv(fruMC.rep, file="/home/jfear/tmp/female_fruMC_repressed.csv", row.names=FALSE, quote=FALSE)

write.csv(induced, file="/home/jfear/tmp/female_total_induced.csv", row.names=FALSE, quote=FALSE)
write.csv(repressed, file="/home/jfear/tmp/female_total_repressed.csv", row.names=FALSE, quote=FALSE)

write.csv(null.induced, file="/home/jfear/tmp/female_total_not_induced.csv", row.names=FALSE, quote=FALSE)
write.csv(null.repressed, file="/home/jfear/tmp/female_total_not_repressed.csv", row.names=FALSE, quote=FALSE)

write.csv(null.ind.rep, file="/home/jfear/tmp/female_total_not_induced_or_repressed.csv", row.names=FALSE, quote=FALSE)

detach(female)

