/* Look at the correlation among the Sxl genes */
proc corr data=SEM.dsrp_sex_det_sbs_combine_sym ;
    var Sxl_:;
    run;

data tmp2;
    retain sample;
    set SEM.dsrp_sex_det_sbs_combine_sym;
    sample = strip(patRil) || '_' || strip(matRil);
    keep sample Sxl_:;
    run;

proc export data=tmp2 outfile='/tmp/mydata_sxl.csv' dbms=csv replace;
    putnames=yes;
    run;

/* Paste the following into R

    library(ggplot2)
    library(reshape)
    mydata <- read.csv('/tmp/mydata_sxl.csv',header=TRUE)
    melted <- melt(mydata, id=c("sample"))
    png("/home/jfear/mclab/cegs_sem_sd_paper/reports/sxl_gene_expression.png",width=1200,height=800)
    qplot(factor(variable),value, data=melted, geom="boxplot", xlab="Isoform" ,ylab="Gene Expression") + geom_point(aes(factor(variable), value),position=position_jitter(width=0.01)) + opts(title="Sxl")
    dev.off()

*/
