data mated;
    retain fusion_id line mating_status rep mapped_reads junc_mapped_reads total_mapped_reads read_length region_length region_depth reads_in_region apn uq_apn log_uq_apn ;
    set CEGLOCAL.norm_basic_stats_m;
    keep fusion_id line mating_status rep mapped_reads junc_mapped_reads total_mapped_reads read_length region_length region_depth reads_in_region apn uq_apn log_uq_apn ;
    run;

data virgin;
    retain fusion_id line mating_status rep mapped_reads junc_mapped_reads total_mapped_reads read_length region_length region_depth reads_in_region apn uq_apn log_uq_apn ;
    set CEGLOCAL.norm_basic_stats_v;
    keep fusion_id line mating_status rep mapped_reads junc_mapped_reads total_mapped_reads read_length region_length region_depth reads_in_region apn uq_apn log_uq_apn ;
    run;

data uqnorm;
    set mated virgin;
    run;

proc export data=uqnorm outfile='!MCLAB/cegs_sergey/pipeline_output/uq_normalized_data_all_lines.csv' dbms=csv replace;
putnames=yes;
run;

proc datasets nolist;
    delete mated;
    delete virgin;
    delete uqnorm;
    run; quit;
