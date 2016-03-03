/* Flag genes on or off

(1) fusion_on0 = 1 if apn>0
(2) fusion_on5 = 1 if apn ge 5
(3) for each line, flag_fusion_on_apn* = 1 if isoform is expressed in > 50% of reps
*/

/* for testing purposes libname diabetes '/media/jrbnewman/SAS_WRK/diabetes/sas_data'; */

libname con '/home/jrbnewman/concannon/sas_data';

/* Sort and merge with design file */
proc sort data=con.counts_by_event;
  by Name;
  run;

proc sort data=con.design_by_subject_new;
  by Name;
  run;

data counts_fusion_key;
  merge con.counts_by_event (in=in1) con.design_by_subject_new (in=in2);
  by Name;
  if in1;
  run;

*make permenant for now;

data con.counts_fusion_key;
   set counts_fusion_key;
run;

data fusions_w_onflag;
  set counts_fusion_key;
  if APN > 0 then fusion_on0 = 1;
  	else fusion_on0 = 0;
  if APN ge 5 then fusion_on5 = 1;
	else fusion_on5 = 0;
  run;


/* Look at counts */
proc sort data=fusions_w_onflag;
  by cell_type name;
  run;

proc means data=fusions_w_onflag noprint;
  by cell_type name;
  output out = sums sum(fusion_on0)=cnt_apn_gt_0 sum(fusion_on5)=cnt_apn_ge_5;
  run;

/* Gene is on if expressed in (apn>0) 50% of the samples..*/
proc sort data=fusions_w_onflag;
  by cell_type fusion_id;
  run;

proc means data = fusions_w_onflag noprint;
  by cell_type fusion_id;
  var fusion_on0;
  output out=trt_on0 mean=trt_percent_on;
  run;

proc sort data = trt_on0 ;
    by fusion_id ;
    run;

proc transpose data = trt_on0 out = trt_on_sbys0 ;
    by fusion_id ;
    id cell_type ;
    var trt_percent_on ;
    run;


/* CD19 B-cells	 	CD8 T-cells 		CD4 T-cells*/

data con.fusions_on_gt_apn0 ;     
    set trt_on_sbys0 ;           
    by fusion_id ;
    if CD19 > 0.5 then flag_CD19_on = 1 ; else flag_CD19_on = 0 ;
    if CD8 >0.5 then flag_CD8_on = 1; else flag_CD8_on = 0;
    if CD4 > 0.5 then flag_CD4_on = 1 ; else flag_CD4_on = 0 ;
 
    if flag_CD19_on=0 and flag_CD8_on=0 and flag_CD4_on=0 then flag_fusion_on0 = 0; 
	else flag_fusion_on0 = 1;
  
    if flag_CD19_on = 1 and flag_CD8_on=1 and flag_CD4_on=1 then flag_fusion_all_on0 =1; 
	else flag_fusion_all_on0=0;
     
    keep fusion_id flag_CD19_on flag_CD8_on flag_CD4_on flag_fusion_on0 flag_fusion_all_on0;
    run;


proc contents data = trt_on_sbys0;  
    run;

ods listing;
ods html close;
proc freq data = con.fusions_on_gt_apn0  ;
    tables flag_: ;
    run;    

* flag_fusion_on0 = 1, 58.07% ;
* flag_fusion_all_on0 = 1, 49.68%, 163713/329557 fusions ;




/* Now for apn ge 5 */

proc means data = fusions_w_onflag noprint;
  by cell_type fusion_id;
  var fusion_on5;
  output out=trt_on5 mean=trt_percent_on;
  run;

proc sort data = trt_on5 ;
    by fusion_id ;
    run;


proc transpose data = trt_on5 out = trt_on_sbys5 ;
    by fusion_id ;
    id cell_type ;
    var trt_percent_on ;
    run;

/* CD19 B-cells	 	CD8 T-cells 		CD4 T-cells*/

data con.fusion_on_ge_apn5 ;      * if 50% of reps then fusion is expressed ;
    set trt_on_sbys5 ;           * using treatments to determine ;
    by fusion_id ;
    if CD19 > 0.5 then flag_CD19_on = 1 ; else flag_CD19_on = 0 ;
    if CD8 >0.5 then flag_CD8_on = 1; else flag_CD8_on = 0;
    if CD4 > 0.5 then flag_CD4_on = 1 ; else flag_CD4_on = 0 ;
 
    if flag_CD19_on=0 and flag_CD8_on=0 and flag_CD4_on=0 then flag_fusion_on5 = 0; 
	else flag_fusion_on5 = 1;
  
    if flag_CD19_on = 1 and flag_CD8_on=1 and flag_CD4_on=1 then flag_fusion_all_on5 =1; 
	else flag_fusion_all_on5=0;
     
    keep fusion_id flag_CD19_on flag_CD8_on flag_CD4_on flag_fusion_on5 flag_fusion_all_on5;
    run;


proc contents data = trt_on_sbys5;  
    run;
proc freq data = con.fusion_on_ge_apn5  ;
    tables flag_: ;
    run;    

* flag_fusion_on5 = 1, 28.00% ;
* flag_fusion_all_on5 =1, 22.29%, 73460/329557 ;
