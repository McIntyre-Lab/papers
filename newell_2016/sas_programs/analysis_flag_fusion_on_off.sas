/* Flag fusions on or off

 *   (1) fusion_on0 = 1 if apn > 0
 *   (2) fusion_on5 =1 if apn ge 5
 *   (3) for each line, flag_fusion_on_apn0 = 1 if fusion is expressed in > 50% of reps
 *                      flag_fusion_on_apn5 = 1 if fusion is expressed in > 50% of reps
 */

libname ribo "!MCLAB/arbeitman/arbeitman_ribotag/sas_data";
*libname dmel "!MCLAB/useful_dmel_data/flybase530/sasdata";


proc contents data = ribo.all_coverage_counts_with_key2;
run;

proc sort data = ribo.all_coverage_counts_with_key2;
by sample_id ;
run;

proc sort data = ribo.design_file;
by sample_id;
run;

data ribo ;
set ribo.all_coverage_counts_with_key2;
if APN > 0 then fusion_on0 = 1;
	else fusion_on0 = 0;
if APN ge 5 then fusion_on5 = 1;
	else fusion_on5 = 0;
run;

proc sort data = ribo ;
by sample_id ;
run;

/* look at counts */
proc means data = ribo noprint ;
by type sex ;
output out = sums sum(fusion_on0)=cnt_apn_gt_0 sum (fusion_on5)=cnt_apn_ge_5;
run;

/* IP FEMALE (315905) : apn0: 135541 (42.9%), apn5: 11256 (3.6%)
   IP MALE (379086) : apn0: 193012 (50.9%), apn5: 24003 (6.3%)
   INPUT FEMALE (315905) : apn0: 240283 (76.1%), apn5: 112124 (35.5%)
   INPUT MALE (379086) : apn0: 306809 (80.9%), apn5: 152541 (40.2%)
*/

/* Going to use apn0 for flagging on/off */

proc sort data= ribo;
by trt fusion_id ;
run;

/* fusion is on if expressed in (apn>0) in 50% of the (reps? There are not reps here, just sex and treatment (ip vs input)) */


proc means data = ribo noprint ;
    by trt fusion_id ;
    var fusion_on0 ;
    output out = trt_on0 mean=trt_percent_on ;
    run ;

proc sort data = trt_on0 ;
    by fusion_id ;
    run;


proc transpose data = trt_on0 out = trt_on_sbys0 ;
    by fusion_id ;
    id trt ;
    var trt_percent_on ;
    run;




data ribo.on_calls_gt_apn0 ;      * if 50% of reps then fusion is expressed ;
    set trt_on_sbys0 ;           * using treatments to determine ;
    by fusion_id ;
    if IPmale > 0.5 then flag_IPmale_on = 1 ; else flag_IPmale_on = 0 ;
    if IPfemale >0.5 then flag_IPfemale_on = 1; else flag_IPfemale_on = 0;
    if InputMale > 0.5 then flag_InputMale_on = 1 ; else flag_InputMale_on = 0 ;
    if InputFemale > 0.5 then flag_InputFemale_on = 1 ; else flag_InputFemale_on = 0 ;

   
    if flag_IPmale_on=0 and flag_IPfemale_on=0 and flag_InputMale_on=0 and flag_InputFemale_on = 0 then flag_fusion_on0 = 0; else flag_fusion_on0 = 1;
  
    if flag_IPmale_on = 1 and flag_IPfemale_on=1 and flag_InputMale_on=1 and flag_InputFemale_on = 1 then flag_fusion_all_on0 =1; else flag_fusion_all_on0=0;
     
    keep fusion_id flag_IPmale_on flag_IPfemale_on flag_InputMale_on flag_InputFemale_on flag_fusion_on0 flag_fusion_all_on0;
    run;


proc contents data = trt_on_sbys0;  
    run;
proc freq data = ribo.on_calls_gt_apn0  ;
    tables flag_: ;
    run;     * flag_fusion_on0: 82.2%, flag_fusion_all_on0: 37.7%;



/* Going to do the same with apn ge 5 too. */

proc means data = ribo noprint ;
    by trt fusion_id ;
    var fusion_on5 ;
    output out = trt_on5 mean=trt_percent_on ;
    run ;

proc sort data = trt_on5 ;
    by fusion_id ;
    run;

proc transpose data = trt_on5  out = trt_on_sbys5 ;
    by fusion_id ;
    id trt ;
    var trt_percent_on ;
    run;

data ribo.on_calls_ge_apn5 ;      * if 50% of reps then fusion is expressed ;
    set trt_on_sbys5 ;           * using treatments to determine ;
    by fusion_id ;
    if IPmale > 0.5 then flag_IPmale_on = 1 ; else flag_IPmale_on = 0 ;
    if IPfemale >0.5 then flag_IPfemale_on = 1; else flag_IPfemale_on = 0;
    if InputMale > 0.5 then flag_InputMale_on = 1 ; else flag_InputMale_on = 0 ;
    if InputFemale > 0.5 then flag_InputFemale_on = 1 ; else flag_InputFemale_on = 0 ;

   
    if flag_IPmale_on=0 and flag_IPfemale_on=0 and flag_InputMale_on=0 and flag_InputFemale_on = 0 then flag_fusion_on5 = 0; else flag_fusion_on5 = 1;
    if flag_IPmale_on = 1 and flag_IPfemale_on=1 and flag_InputMale_on=1 and flag_InputFemale_on = 1 then flag_fusion_all_on5 =1; else flag_fusion_all_on5=0;
     
    keep fusion_id flag_IPmale_on flag_IPfemale_on flag_InputMale_on flag_InputFemale_on flag_fusion_on5 flag_fusion_all_on5;
    run;


proc contents data = trt_on_sbys0;  
    run;
proc freq data = ribo.on_calls_gt_apn0  ;
    tables flag_: ;
    run;     




