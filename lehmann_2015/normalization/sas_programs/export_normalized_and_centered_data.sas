data centered;
    retain fusion_id line mating_status rep apn uq_ff uq_apn log_uq_apn uq_log_uq_apn uq_log_uq_center;
    set CEGLOCAL.ccfus_norm_centered;
    keep fusion_id line mating_status rep apn uq_ff uq_apn log_uq_apn uq_log_uq_apn uq_log_uq_center;
    run;

proc export data=centered outfile='!MCLAB/cegs_sergey/pipeline_output/uq_normalized_and_centered_data_all_lines.csv' dbms=csv replace;
putnames=yes;
run;

proc datasets nolist;
    delete centered;
    run; quit;
