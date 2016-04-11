libname dmel '!MCLAB/useful_dmel_data/flybase551/sasdata' ;
libname ribo '!MCLAB/arbeitman/arbeitman_ribotag/sas_data' ;

data dmel ;
set dmel.fb551_si_fusions_unique_flagged;
if exons_per_fusion = 1 and genes_per_fusion = 1 then flag_single_exon_gene = 1;
     else flag_single_exon_gene = 0;

if fbtrs_per_fusion >1 and Constitutive = 1 then flag_single_iso_gene = 1;
     else flag_single_iso_gene = 0;

keep fusion_id genes_per_fusion fbtrs_per_fusion flag_single_exon_gene flag_single_iso_gene constitutive ;
run ;


data isos ;
set dmel.fb551_si_fusions_unique_flagged ;
where constitutive = 1 and fbtrs_per_fusion > 1 and fbgns_per_fusion = 1 ;
run ;



/* On/Off fusions for apn>0 */
proc sort data = ribo.on_calls_gt_apn0;
by fusion_id;

proc sort data = dmel ;
by fusion_id ;
run;


data to_cnt_0 oops;
merge ribo.on_calls_gt_apn0 (in=in1) dmel (in=in2) ;
by fusion_id;
if in1 and in2 then output to_cnt_0 ;
else output oops ; *525 in oops, but every obs from on_calls_gt_apn0 is present in the to_cnt_0, because I used non-redundant fusions to align and get counts;
run;
*there are 63706 fusions in the dmel.fb551_si_fusions_unique_flagged  included in this file are redundant fusions that have been eliminated. There are only 63181 non-redundant fusions;

data flag_on_all_trts_0 ;
set to_cnt_0;
if flag_IP_on = 1 and flag_Input_on = 1 then flag_on_all_trts = 1;
else flag_on_all_trts =0;
run;

proc freq data = flag_on_all_trts_0  ;
tables flag_on_all_trts  ;
run;

proc freq data = flag_on_all_trts_0  ;
tables flag_single_exon_gene*flag_on_all_trts;
run;

proc freq data = flag_on_all_trts_0  ;
tables flag_single_iso_gene*flag_on_all_trts;
run;


proc sort data =  ribo.on_calls_gt_apn5;
by fusion_id ;
proc sort data = dmel ;
by fusion_id ;

data to_cnt_5 oops ;
merge  ribo.on_calls_gt_apn5 (in=in1) dmel (in=in2) ;
by fusion_id ;
if in1 and in2 then output to_cnt_5 ;
else output oops ;   *0 in oops ;
run;

data flag_on_all_trts_5;
set to_cnt_5 ;
if flag_IP_on = 1  and flag_Input_on = 1 then flag_on_all_trts = 1 ;
else flag_on_all_trts = 0;
run ;
