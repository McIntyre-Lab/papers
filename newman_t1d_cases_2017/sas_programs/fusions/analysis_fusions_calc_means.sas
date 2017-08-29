/*** CALCULATING MEANS BY FUSION ***/

/* import libraries */
libname con '/home/jrbnewman/concannon/sas_data';


/*** TRANSPOSE DATA ***/

/* Transpose depth by sample - want "zero" counts back in!! */

data fusion_data_for_tpose;
   set con.apn0_q3_norm;
   keep Name fusion_id log_q3_q3_apn_filter;
run;

proc sort data=fusion_data_for_tpose;
    by fusion_id Name;
run;

proc transpose data=fusion_data_for_tpose out=fusion_apn_data_tpose;
   by fusion_id;
   var log_q3_q3_apn_filter;
   id Name;
run;

data con.fusion_apn_data_tpose;
   set fusion_apn_data_tpose;
   array change _numeric_;
            do over change;
            if change=. then change=0;
            end;
   drop _NAME_;
   run ;



/**** CALCULATE GROUP MEANS ****/

/* UNTRANSPOSE AND CALC GROUP MEAN */

proc transpose name=Name data=con.fusion_apn_data_tpose out=fusion_q3apn_untpose(rename=(col1=log_q3_q3_apn));
     by fusion_id;
     var _2009_PC_:;
run;

data fusion_q3apn_untpose_2;
   set fusion_q3apn_untpose;
   Name=tranwrd(Name, "_2009_PC_", "2009-PC-");
run;


/* make permenant -- need for ANOVAs! */

data con.fusion_q3_norm_data_all;
   set fusion_q3apn_untpose_2;
run;

/* get groups */

data cell_types;
  set con.design_file;
  keep Name cell_type;
run;

proc sort data=fusion_q3apn_untpose_2;
   by Name;
run;

proc sort data=cell_types;
by Name;
run;

data fusion_exp_w_cell oops1 oops2;
   merge cell_types (in=in1) fusion_q3apn_untpose_2 (in=in2);
   by Name;
   if in1 and in2 then output fusion_exp_w_cell;
   else if in1 then output oops1;
   else output oops2;
run;


proc sort data=fusion_exp_w_cell;
  by cell_type fusion_id;
run;

proc means data=fusion_exp_w_cell noprint;
   var log_q3_q3_apn;
   by cell_type fusion_id;
   output out=fusion_mean_counts_by_fusion mean=mean;
run;

/* Split means on cell_type */

data means_cd19 means_cd8 means_cd4 ;
   set fusion_mean_counts_by_fusion;
   if cell_type='CD19' then output means_cd19;
   if cell_type='CD4' then output means_cd4;
   if cell_type='CD8' then output means_cd8;
   drop _TYPE_ _FREQ_;
   run;

data means_cd19_2;
   set means_cd19;
   rename mean=mean_logq3q3apn_cd19;
   drop cell_type;
run;

data means_cd4_2;
   set means_cd4;
   rename mean=mean_logq3q3apn_cd4;
   drop cell_type;
run;

data means_cd8_2;
   set means_cd8;
   rename mean=mean_logq3q3apn_cd8;
   drop cell_type;
run;

data means_merge;
   merge means_cd19_2 means_cd4_2 means_cd8_2;
   by fusion_id;
run;


/* make permenant */

data con.fusion_means_by_celltype;
   set means_merge;
run;

