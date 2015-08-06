/* The goal here is to identify the differences in gene expression between the different fusions */

data SEMLOCAL.dsx;
    set SEM.sex_stacked;
    where symbol = 'dsx';
    run;

proc means data=SEM.sex_stacked;
    by fusion_id;
    output out=sem_stacked_means mean(apn) = mean_apn mean(rpkm) = mean_rpkm;
    run;


data dsx_fusions;
    set sem_stacked_means;
    if fusion_id = 'S47254_SI' or 
    fusion_id = 'S47253_SI' or 
    fusion_id = 'S47252_SI' or 
    fusion_id = 'S47251_SI' or 
    fusion_id = 'S47250_SI' or 
    fusion_id = 'S47249_SI' or 
    fusion_id = 'F47248_SI' then output;
    run;


proc means data=SEM.sex_stacked noprint;
    by fusion_id line;
    output out=sem_line_means mean(apn) = mean_apn mean(rpkm) = mean_rpkm;
    run;


data dsx_fusions_line;
    set sem_line_means;
    if fusion_id = 'S47254_SI' or 
    fusion_id = 'S47253_SI' or 
    fusion_id = 'S47252_SI' or 
    fusion_id = 'S47251_SI' or 
    fusion_id = 'S47250_SI' or 
    fusion_id = 'S47249_SI' or 
    fusion_id = 'F47248_SI' then output;
    run;
