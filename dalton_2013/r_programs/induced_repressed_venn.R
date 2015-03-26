library(VennDiagram)

setwd("/home/jfear/mclab/arbeitman_fru_network/reports/venn_diagram")

male <- read.csv("/home/jfear/mclab/arbeitman_fru_network/reports/flag_x_induced_repressed_male.csv",header=TRUE)
female <- read.csv("/home/jfear/mclab/arbeitman_fru_network/reports/flag_x_induced_repressed_female.csv",header=TRUE)

make_venn_two_way <- function(dataset,letter,main,sub,filename=NULL){

    varname <- c(paste("flag_",letter,"_ind",sep=""), paste("flag_",letter,"_rep",sep="")) 
    induced <- dataset[dataset[[varname[1]]] == 1,]$primary_fbgn
    repressed <- dataset[dataset[[varname[2]]] == 1,]$primary_fbgn

    my.venn <- venn.diagram(list("Induced"=induced, "Repressed"=repressed),
                     filename=filename,
                     scale=FALSE,
                     main = main,
                     sub = sub,
                     cat.dist=0.05,
                     cat.pos=0,
                     fill=c(2,4),
                     main.cex=2,
                     sub.cex=1,
                     sub.fontface="bold",
                     width=1200,
                     height=1200,
                     resolution=300
                     )
    return(my.venn)
}


fru.a.male <- make_venn_two_way(male,"a","Fru A","male",filename="fru_a_male_venn.tiff")
fru.b.male <- make_venn_two_way(male,"b","Fru B","male",filename="fru_b_male_venn.tiff")
fru.c.male <- make_venn_two_way(male,"c","Fru C","male",filename="fru_c_male_venn.tiff")

fru.a.female <- make_venn_two_way(female,"a","Fru A","female",filename="fru_a_female_venn.tiff")
fru.b.female <- make_venn_two_way(female,"b","Fru B","female",filename="fru_b_female_venn.tiff")
fru.c.female <- make_venn_two_way(female,"c","Fru C","female",filename="fru_c_female_venn.tiff")

make_venn_three_way <- function(dataset,letter1,letter2,letter3,type,main,sub=NULL,filename=NULL){

    varname <- c(paste("flag_",letter1,"_",type,sep=""), paste("flag_",letter2,"_",type,sep=""), paste("flag_",letter3,"_",type,sep="")) 
    Fru_A <- dataset[dataset[[varname[1]]] == 1,]$primary_fbgn
    Fru_B <- dataset[dataset[[varname[2]]] == 1,]$primary_fbgn
    Fru_C <- dataset[dataset[[varname[3]]] == 1,]$primary_fbgn

    my.venn <- venn.diagram(list("Fru_A"=Fru_A, "Fru_B"=Fru_B, "Fru_C"=Fru_C),
                     filename=filename,
                     scale=FALSE,
                     main = main,
                     fill=c(2,4,5),
                     cat.dist=0.05,
                     #cat.pos=0,
                     cat.cex=1.1,
                     cat.fontface="bold",
                     main.cex=2,
                     sub.cex=1,
                     sub.fontface="bold",
                     rotation=1,
                     reverse=FALSE,
                     width=1200,
                     height=1200,
                     resolution=300,
                     sep.dist=0.1)

    return(my.venn)
}

fru.ind.male <- make_venn_three_way(male,"a","b","c","ind","Induced",filename="fru_ind_male.tiff")
fru.rep.male <- make_venn_three_way(male,"a","b","c","rep","Repressed",filename="fru_rep_male.tiff")

fru.ind.female <- make_venn_three_way(female,"a","b","c","ind","Induced",filename="fru_ind_female.tiff")
fru.rep.female <- make_venn_three_way(female,"a","b","c","rep","Repressed",filename="fru_rep_female.tiff")
