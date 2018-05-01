ods listing; ods html close;

libname event '!MCLAB/event_analysis/sas_data';
libname eventloc '/mnt/store/event_sandbox/sas_data';

/* Import case-only fragment counts */

    data WORK.COUNTS_BY_FRAGMENT    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile '/mnt/store/event_sandbox/hg19_aceview_fragment_counts.csv' delimiter
= ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
       informat sample_id $31. ;
       informat fragment_id $15. ;
       informat mapped_reads best32. ;
       informat read_length best32. ;
       informat region_length best32. ;
       informat region_depth best32. ;
       informat reads_in_region best32. ;
       informat apn best32. ;
       informat rpkm best32. ;
       informat mean best32. ;
       informat std best32. ;
       informat cv best32. ;
       format sample_id $31. ;
       format fragment_id $15. ;
       format mapped_reads best12. ;
       format read_length best12. ;
       format region_length best12. ;
       format region_depth best12. ;
       format reads_in_region best12. ;
       format apn best12. ;
       format rpkm best12. ;
       format mean best12. ;
       format std best12. ;
       format cv best12. ;
    input
                sample_id $
                fragment_id $
                mapped_reads
                read_length
                region_length
                region_depth
                reads_in_region
                apn
                rpkm
                mean
                std
                cv
    ;
    if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
    run;

/* Sum tech reps */

data design_file;
  set con.design_file;
  keep sample_id cell_type name;
run;

data counts_by_fragment2;
  set counts_by_fragment;
  keep sample_id fragment_id apn;
run;

proc sort data=design_file;
   by sample_id;
proc sort data=counts_by_fragment;
   by sample_id;
run;

data counts_w_key;
    merge design_file (in=in1) counts_by_fragment (in=in2);
   by sample_id;
   if in1 and in2;
run;

proc sort data=counts_w_key;
   by name cell_type fragment_id;
proc means data=counts_w_key noprint;
   by name cell_type fragment_id;
   var apn;
   output out=summed_counts_by_frag sum=;
run;

/* Make permenant */

data eventloc.hg19_counts_by_fragment;
   set summed_counts_by_frag;
   drop _TYPE_ _FREQ_;
run;


/* Flag fragment as "on" if APN>=5 */

data flag_counts;
   set eventloc.hg19_counts_by_fragment;
   if apn ge 5 then flag_frag_on=1;
   else flag_frag_on=0;
run;

/* Flag fragments as "on" if detected at APN >= 5 in at least 50% of samples */

proc sort data=flag_counts;
    by fragment_id cell_type;
proc means data=flag_counts noprint;
   by fragment_id cell_type;
   var flag_frag_on;
   out=perc_on_by_trt mean=;
run;

* split by cell type;

data cd4 cd8 cd19;
   set perc_on_by_trt;
   if cell_type="CD4" then output cd4;
   else if cell_type="CD8" then output cd8;
   else if cell_type="CD19" then output cd19;
run;

data flag_cd4;
  set cd4;
  if flag_frag_on >= 0.5 then flag_cd4_on=1;
  else if flag_frag_on=0 then flag_cd4_on=0;
  else flag_cd4_on=.;
  keep fragment_id flag_cd4_on;
run;

data flag_cd8;
  set cd8;
  if flag_frag_on >= 0.5 then flag_cd8_on=1;
  else if flag_frag_on=0 then flag_cd8_on=0;
  else flag_cd8_on=.;
  keep fragment_id flag_cd8_on;
run;

data flag_cd19;
  set cd19;
  if flag_frag_on >= 0.5 then flag_cd19_on=1;
  else if flag_frag_on=0 then flag_cd19_on=0;
  else flag_cd19_on=.;
  keep fragment_id flag_cd19_on;
run;

proc sort data=flag_cd4;
   by fragment_id;
proc sort data=flag_cd8;
   by fragment_id;
proc sort data=flag_cd19;
   by fragment_id;
run;

data flag_fragment_on_apn5;
   merge flag_cd4 flag_cd8 flag_cd19;
   by fragment_id;
   if flag_cd4_on=1 and flag_cd8_on=1 and flag_cd19_on=1 then flag_fragment_all_on_ge5=1;
   else flag_fragment_all_on_ge5=0;
run;

/* Make permenant */

data event.hg19_flag_fragment_on_apn5;
   set flag_fragment_on_apn5;
run;

