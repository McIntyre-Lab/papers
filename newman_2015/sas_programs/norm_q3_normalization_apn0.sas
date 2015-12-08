/********************************************************************************
* Normalizations!
* Redoing normalizations on all data with APN>0. Will be dropping fusions that
* aren't expressed later, but this will give us a more accurately normalized data.
********************************************************************************/

/* need to open up the libraries for the diabetes data */

/* for testing purposes libname diabetes '/media/jrbnewman/SAS_WRK/diabetes/sas_data'; */

libname con '/home/jrbnewman/concannon/sas_data';

/* Recreating full dataset. This is apparently large so we will store it locally for now */
/* First merge with on_Calls_gt_apn0 and filter by flag_iso_all_on0 to make a clean dataset */

proc sort data=con.counts_fusion_key;
  by fusion_id;
  run;

proc sort data=con.fusions_on_gt_apn0;
  by fusion_id;
  run;

data all_cvg_cnts_w_key_apn0;
  merge con.counts_fusion_key (in=in1) con.fusions_on_gt_apn0 (in=in2);
  by fusion_id;
  if in1;
  run;

proc sort data=all_cvg_cnts_w_key_apn0;
  by name;
  run;

data all_cvg_cnt_apn0_filtered;
  set all_cvg_cnts_w_key_apn0;
  if APN gt 0;
  run;

/* quartiles */
proc univariate data=all_cvg_cnt_apn0_filtered noprint;
  by Name;
  var apn;
  output out=quartiles_apn0 Q3=Q3;
  run;

proc freq data=all_cvg_cnt_apn0_filtered;
tables Name;
 run;

/* Calculate the Q3 of Q3 for normalizations */
proc means data=quartiles_apn0 Q3;
            var q3;
            run;

/* Q3 of Q3 =   20.0264901*/

/* Now drop five low coverage samples and recalculate Q3 */

data quartiles_apn0_filtered;
   set quartiles_apn0;
if name =  '2009-PC-0221' then delete; *sample 75 cd8;
if name =  '2009-PC-0144' then delete; *sample 48 cd4;
if name =  '2009-PC-0236' then delete; *sample 80;
if name =  '2009-PC-0237' then delete; *sample 80;
if name =  '2009-PC-0235' then delete; *sample 80; 
run;

/* Calculate the median sum, median and uq */
proc means data=quartiles_apn0_filtered q3;
            var q3;
            run;

/* Q3 of Q3 (filtered) =  20.0358566 */

/* Now we need to merge the stats back into the original dataset */

    proc sort data=all_cvg_cnt_apn0_filtered;
        by name;
        run;

    proc sort data=quartiles_apn0;
        by name;
        run;

/* merge and oops test */

    data all_cvg_cnt_apn0_merge oops_apn0;
        merge all_cvg_cnt_apn0_filtered (in=in1) quartiles_apn0 (in=in2);
        by name;
        if in1 and in2 then output all_cvg_cnt_apn0_merge;
        else output oops_apn0; 
        run;


/* Merging worked. Now I need to calculate the adjustments per sample */
/* Divide the across-sample UQ by the sample UQ and then output */


/* Here I'm goint to calculate several adjustment factors so we can choose the best normalization strategy */
/* Then we want to export this as plot the distributions */

    data con.apn0_q3_norm;
        set all_cvg_cnt_apn0_merge;
        log_apn=log2(apn + 1);

        * Q3 normalization - all samples;
        q3_q3_apn=(apn/q3) * 20.0264901;
        q3_q3_ff=20.0264901/q3;
        log_q3_q3_apn=log2(q3_q3_apn + 1);

        * Q3 normalization - removed low coverage samples;
        q3_q3_apn_filter=(apn/q3) * 20.0358566;
        q3_q3_ff_filter=20.0358566/q3;
        log_q3_q3_apn_filter=log2(q3_q3_apn_filter + 1);
        
        * set low coverage samples to missing for filtered normalization;
        if name = '2009-PC-0221' then do;
             flag_low_coverage=1;
             q3_q3_apn_filter=.;
             q3_q3_ff_filter=.;
             log_q3_q3_apn_filter=.;
             end;
        else if name =  '2009-PC-0144' then do; *sample 48 cd4;
             flag_low_coverage=1;
             q3_q3_apn_filter=.;
             q3_q3_ff_filter=.;
             log_q3_q3_apn_filter=.;
             end;
        else if name =  '2009-PC-0236' then do; *sample 80;
             flag_low_coverage=1;
             q3_q3_apn_filter=.;
             q3_q3_ff_filter=.;
             log_q3_q3_apn_filter=.;
             end;
        else if name =  '2009-PC-0237' then do; *sample 80;
             flag_low_coverage=1;
             q3_q3_apn_filter=.;
             q3_q3_ff_filter=.;
             log_q3_q3_apn_filter=.;
             end;
        else if name =  '2009-PC-0235' then do; *sample 80; 
             flag_low_coverage=1;
             q3_q3_apn_filter=.;
             q3_q3_ff_filter=.;
             log_q3_q3_apn_filter=.;
             end;
        else flag_low_coverage=0;
run;

/* Going to export data to plot for later reference */


    proc sort data=con.apn0_q3_norm;
        by name;
        run;

data q3_norm_plots_apn0;
    set con.apn0_q3_norm;
    keep fusion_id Name apn log_apn cell_type log_q3_q3_apn q3_q3_ff q3_q3_ff_filter log_q3_q3_apn_filter;
    run;

/* Export CSV for plotting */

proc export data=q3_norm_plots_apn0 outfile='/home/jrbnewman/concannon/normalization/fusions_q3_norm_apn0.csv' dbms=csv replace;
    putnames=yes;
    run;


