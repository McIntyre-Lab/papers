/********************************************************************************
* Need to check the distribution of the data after normalization.
* 
* libname dros "!MCLAB/Dros_PB_ChIP/sasdata/RNAseq";
********************************************************************************/


/* fusions */
%macro norm_plots (species) ;

/* Female */
data female_&species.;
    retain sample featureID log_uq_apn;
    set dros.norm_stats_fsn_&species._f;
    sample = trim(genotype) || '_' || trim(sex) || trim(rep);
    keep sample featureID log_uq_apn;
    run;

proc export data=female_&species. outfile="/home/ammorse/TB14/etoh_srna/roz/female_&species..csv" dbms=csv replace;
putnames=yes;
run;

data _null_;
    call system("Rscript /home/ammorse/mclab/SHARE/McIntyre_Lab/Dros_PB_ChIP/r_programs/plot_normalization_distributions_v2.R /home/ammorse/TB14/etoh_srna/roz/female_&species..csv &species._fsn_female");
    run;

/* Male */
data male_&species.;
    retain sample featureID log_uq_apn;
    set dros.norm_stats_fsn_&species._m;
    sample = trim(genotype) || '_' || trim(sex) || trim(rep);
    keep sample featureID log_uq_apn;
    run;

proc export data=male_&species. outfile="/home/ammorse/TB14/etoh_srna/roz/male_&species..csv" dbms=csv replace;
putnames=yes;
run;

data _null_;
    call system("Rscript /home/ammorse/mclab/SHARE/McIntyre_Lab/Dros_PB_ChIP/r_programs/plot_normalization_distributions_v2.R /home/ammorse/TB14/etoh_srna/roz/male_&species..csv &species._fsn_male");
    run;
%mend ;

%norm_plots (mel) ;
%norm_plots (sim) ;

/* fragments */

