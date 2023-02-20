/********************************************************************************
* Need to check the distribution of the data after normalization.
* 
* 
********************************************************************************/

libname drosRNA "!MCLAB/Dros_PB_ChIP/sasdata/RNAseq";


/* fusions 

*** uq_median of 70.5 
 'bad' sample with low q3 (sim_12_m_noEtoh_rep2) not included!!
*/

%macro norm_plots (species) ;

/* all samples for each species */
data both_&species.;
    retain sampleID featureID log_uq_apn;
    set drosRNA.norm_stats_fsn_&species.;
    sample = trim(genotype) || '_' || trim(sex) ||'_'|| trim(rep);
    keep sampleID featureID log_uq_apn ;
    if sampleID = "sim_12_m_noEtoh_rep2" then delete ;
    run;

proc sort data = both_&species ;
by sampleID ;
run ;

proc export data=both_&species. outfile="/home/ammorse/TB14/etoh_srna/roz/both_&species..csv" dbms=csv replace;
putnames=yes;
run;
/*
data _null_;
    call system("Rscript /home/ammorse/mclab/SHARE/McIntyre_Lab/Dros_PB_ChIP/r_programs/plot_normalization_distributions_v2.R /home/ammorse/TB14/etoh_srna/roz/both_&species..csv &species._fsn_F_M_etoh_noEtoh");
*/
%mend ;

%norm_plots (mel) ;
%norm_plots (sim) ;


/* using mapped reads FF */
%macro norm_plots (species) ;

/* all samples for each species */
data both_map_&species.;
    retain sampleID featureID log_uq_apn;
    set drosRNA.norm_stats_fsn_&species.;
    sample = trim(genotype) || '_' || trim(sex) ||'_'|| trim(rep);
    keep sampleID featureID log_map_apn ;
    if sampleID = "sim_12_m_noEtoh_rep2" then delete ;
    run;

proc sort data = both_map_&species ;
by sampleID ;
run ;

proc export data=both_map_&species. outfile="/home/ammorse/TB14/etoh_srna/roz/both_map_&species..csv" dbms=csv replace;
putnames=yes;
run;
/*
data _null_;
    call system("Rscript /home/ammorse/mclab/SHARE/McIntyre_Lab/Dros_PB_ChIP/r_programs/plot_normalization_distributions_v2.R /home/ammorse/TB14/etoh_srna/roz/both_&species..csv &species._fsn_F_M_etoh_noEtoh");
*/
%mend ;

%norm_plots (mel) ;
%norm_plots (sim) ;




