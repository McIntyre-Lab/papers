ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Set up analysis for differential exon-skipping analysis
   For this I am going to analyze the set of exon-skipping events (detected)

   Then summarize by gene */

/* Get events to analyze */

data exonskip;
   set evspl.splicing_events_annot_refseq;
   where flag_exonskip=1;
   keep event_id gene_id;
run;

/* For each cell type, flag if fragments are detected */

    data WORK.SPLICING_COUNTS    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile '!MCLAB/event_analysis/alignment_output/mm10_refseq_splicing_counts.csv' delimiter
= ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
       informat sample_id $4. ;
       informat event_id $335. ;
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
       format sample_id $4. ;
       format event_id $335. ;
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
               event_id $
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

data set_group;
  length cell_type $3.;
  set splicing_counts;
  if sample_id in ('NSC1','NSC2') then cell_type="NSC";
  if sample_id in ('OLD1','OLD2') then cell_type="OLD";
  keep sample_id cell_type event_id apn;
run;

data flag_event_apn_gt0;
  set set_group;
  if apn > 0 then flag_event_apn_gt0=1;
  else flag_event_apn_gt0=0;
run;


proc sort data=flag_event_apn_gt0;
  by cell_type event_id;
proc means data=flag_event_apn_gt0 noprint;
  by cell_type event_id;
  var flag_event_apn_gt0;
  output out=mean_on mean=mean_gt0;
run;

data nsc old;
  set mean_on;
  if cell_type="NSC" then output nsc;
  else output old;
run;

data flag_on_nsc;
  set nsc;
  if mean_gt0 ge 0.5 then flag_event_nsc_on=1;
  else flag_event_nsc_on=0;
  keep event_id flag_event_nsc_on;
run;

data flag_on_old;
  set old;
  if mean_gt0 ge 0.5 then flag_event_old_on=1;
  else flag_event_old_on=0;
  keep event_id flag_event_old_on;
run;

proc sort data=flag_on_nsc;
   by event_id;
proc sort data=flag_on_old;
   by event_id;
run;

data flag_event_on;
   merge flag_on_nsc (in=in1) flag_on_old (in=in2);
   by event_id;
   if in1 and in2;
run;

data events_for_analysis;
   set flag_event_on;
   where flag_event_nsc_on=1 and flag_event_old_on=1;
   keep event_id;
run; *136131 events analyzable;

proc sort data=exonskip;
  by event_id;
proc sort data=events_for_analysis;
   by event_id;
run;

data events_for_analysis2;
   merge events_for_analysis (in=in1) exonskip (in=in2);
   by event_id;
   if in1 and in2;
run; *13907 exonskip events analyzable;

proc sort data=set_group;
   by event_id;
run;

data counts_for_anova;
   merge events_for_analysis2 (in=in1) set_group (in=in2);
   by event_id;
   if in1 and in2;
   log_apn = log(apn + 1);
run;

proc sort data=counts_for_anova;
   by event_id cell_type;
run;

ods listing close;
proc glimmix data=counts_for_anova;
   by event_id;
   class cell_type;
   model log_apn = cell_type / htype=1;
   output out=resid resid=resid pred=pred student=stu;
ods output tests1=anova;
run;
quit;

/* Flag residuals */
proc univariate data=resid normal noprint;
   by event_id;
   var Resid;
   output out=normtest probn=pnorm;
run;

data flag_resids;
  set normtest;
  if pnorm = . then flag_fail_norm=.;
  else if pnorm le 0.05 then flag_fail_norm=1;
  else flag_fail_norm=0;
run;

proc freq data=flag_Resids noprint;
  tables flag_fail_norm / out=splicing_flag_fail_norm;
run;



ods listing;
proc print data=splicing_flag_fail_norm;
run;

/*
        flag_
        fail_
 Obs     norm    COUNT    PERCENT

  1       .          6      .
  2       0      12500    89.9216
  3       1       1401    10.0784
*/

/* Make permenant */

data event.exonskip_flag_failnorm;
   set splicing_flag_fail_norm;
run;

data event.exonskip_anova;
  set anova;
run;

data event.exonskip_flag_resid;
   set flag_resids;
run;

/* Calculate FDR */

data anova_results;
  set event.exonskip_anova;
  where effect="cell_type";
  keep event_id probf;
run;


proc multtest inpvalues(ProbF)=anova_results fdr
   out=anova_results_w_fdr;
run;


data on_flags;
   set flag_event_on;
   where flag_event_nsc_on=1 or flag_event_old_on=1;
   keep event_id flag_event_nsc_on flag_event_old_on;
run;

proc sort data=on_flags;
  by event_id;
proc sort data=exonskip;
  by event_id;
proc sort data=anova_results_w_fdr;
  by event_id;
run;

data exonskip_on;
   merge exonskip (in=in1) on_flags (in=in2);
   by event_id;
   if in1 and in2;
run;

data exonskip_on_w_fdr;
  merge exonskip_on (in=in1) anova_results_w_fdr (in=in2);
  by event_id;
  if in1;
run;

data flag_fdr;
  set exonskip_on_w_fdr;
  if flag_event_nsc_on=1 and flag_event_old_on=1 then do;
         if fdr_p=. then do;
           flag_exonskip_de=.;
           flag_exonskip_dd=.; end;
         else if fdr_p < 0.05 then do;
           flag_exonskip_de=1;
           flag_exonskip_dd=1; end;
         else do;
           flag_exonskip_de=0;
           flag_exonskip_dd=0; end;
         end;
  else if flag_event_nsc_on=1 and flag_event_old_on=0 then do;
          flag_exonskip_de=.;
          flag_exonskip_dd=1; end;
  else if flag_event_nsc_on=0 and flag_event_old_on=1 then do;
          flag_exonskip_de=.;
          flag_exonskip_dd=1; end;
  else do;
          flag_exonskip_de=.;
          flag_exonskip_dd=0; end;
run;

proc freq data=flag_fdr;
   tables flag_exonskip_de flag_exonskip_dd;
run;


/*
                                                 Cumulative    Cumulative
    flag_exonskip_de    Frequency     Percent     Frequency      Percent
    ---------------------------------------------------------------------
                   0       13899       99.97         13899        99.97
                   1           4        0.03         13903       100.00

                          Frequency Missing = 19555


                                                 Cumulative    Cumulative
    flag_exonskip_dd    Frequency     Percent     Frequency      Percent
    ---------------------------------------------------------------------
                   0       13899       41.55         13899        41.55
                   1       19555       58.45         33454       100.00

                            Frequency Missing = 4
*/

proc sort data=flag_fdr;
   by gene_id;
proc means data=flag_fdr noprint;
  by gene_id;
  var flag_exonskip_de flag_exonskip_dd;
  output out=num_de_dd_exonskip_by_gene sum=;
run;

/* Merge with MISO */

data miso_diff_se;
   set event.miso_diff_se_refseq;
run;

proc sort data=miso_diff_se;
   by gene_id;
proc sort data=num_de_dd_exonskip_by_gene;
   by gene_id;
run;

data miso_to_exonskip;
  merge miso_diff_se (in=in1) num_de_dd_exonskip_by_gene (in=in2);
  by gene_id;
  if not in1 then do;
      num_diff_se_bf10=0;
      num_diff_se_bf5=0; end;
  if not in2 then do;
       flag_exonskip_de=0;
       flag_exonskip_dd=0;
  end;
run;

data flag_genes;
  set miso_to_exonskip;
  if num_diff_se_bf10 > 0 then flag_miso_se_diff_bf10=1; else flag_miso_se_diff_bf10=0;
  if num_diff_se_bf5 > 0 then flag_miso_se_diff_bf5=1; else flag_miso_se_diff_bf5=0;
  if flag_exonskip_de > 1 then flag_event_se_de=1; else flag_event_se_de=0;
  if flag_exonskip_dd > 1 then flag_event_se_dd=1; else flag_event_se_dd=0;
run;

proc freq data=flag_genes;
   tables flag_miso_se_diff_bf10*flag_event_se_de
          flag_miso_se_diff_bf10*flag_event_se_dd
          flag_miso_se_diff_bf5*flag_event_se_de
          flag_miso_se_diff_bf5*flag_event_se_dd;
run;



/*
  flag_miso_se_diff_bf10
            flag_event_se_dd

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |   6831 |   4617 |  11448
           |  59.59 |  40.28 |  99.87
           |  59.67 |  40.33 |
           |  99.84 |  99.91 |
  ---------+--------+--------+
         1 |     11 |      4 |     15
           |   0.10 |   0.03 |   0.13
           |  73.33 |  26.67 |
           |   0.16 |   0.09 |
  ---------+--------+--------+
  Total        6842     4621    11463
              59.69    40.31   100.00

             The SAS System         15:16 M

 flag_miso_se_diff_bf5
           flag_event_se_dd

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |   6828 |   4614 |  11442
          |  59.57 |  40.25 |  99.82
          |  59.67 |  40.33 |
          |  99.80 |  99.85 |
 ---------+--------+--------+
        1 |     14 |      7 |     21
          |   0.12 |   0.06 |   0.18
          |  66.67 |  33.33 |
          |   0.20 |   0.15 |
 ---------+--------+--------+
 Total        6842     4621    11463
             59.69    40.31   100.00

*/

proc freq data=flag_genes;
   where gene_id ne "";
   tables flag_miso_se_diff_bf10*flag_event_se_de
          flag_miso_se_diff_bf10*flag_event_se_dd
          flag_miso_se_diff_bf5*flag_event_se_de
          flag_miso_se_diff_bf5*flag_event_se_dd;
run;


/*
flag_miso_se_diff_bf10
          flag_event_se_dd

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 |   6461 |   4617 |  11078
         |  58.26 |  41.63 |  99.89
         |  58.32 |  41.68 |
         |  99.88 |  99.91 |
---------+--------+--------+
       1 |      8 |      4 |     12
         |   0.07 |   0.04 |   0.11
         |  66.67 |  33.33 |
         |   0.12 |   0.09 |
---------+--------+--------+
Total        6469     4621    11090
            58.33    41.67   100.00

   flag_miso_se_diff_bf5
             flag_event_se_dd

   Frequency|
   Percent  |
   Row Pct  |
   Col Pct  |       0|       1|  Total
   ---------+--------+--------+
          0 |   6458 |   4614 |  11072
            |  58.23 |  41.61 |  99.84
            |  58.33 |  41.67 |
            |  99.83 |  99.85 |
   ---------+--------+--------+
          1 |     11 |      7 |     18
            |   0.10 |   0.06 |   0.16
            |  61.11 |  38.89 |
            |   0.17 |   0.15 |
   ---------+--------+--------+
   Total        6469     4621    11090
               58.33    41.67   100.00


*/


