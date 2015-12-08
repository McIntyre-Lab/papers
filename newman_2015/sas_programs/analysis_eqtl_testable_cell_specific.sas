libname eqtl '/home/jrbnewman/concannon/eqtl_analysis/sas_data';
libname splicing '/mnt/data/splicing';
libname con '/home/jrbnewman/concannon/sas_data';

/* Need the list of testable exons, AS events and IR events */

data eqtl_exon_testable;
   set eqtl.eqtl_results_summary;
   if feature_type='exon';
   flag_testable=1;
   keep feature_id flag_testable;
   run;
   
data eqtl_as_testable;
   set eqtl.eqtl_results_summary;
   if feature_type='AS';
   flag_testable=1;
   keep feature_id  flag_testable;
   run;
      
data eqtl_ir_testable;
   set eqtl.eqtl_results_summary;
   if feature_type='IR';
   flag_testable=1;
   keep feature_id  flag_testable;
   run;
   
proc sort data=eqtl_exon_testable nodup;
   by feature_id;
proc sort data=eqtl_as_testable nodup;
   by feature_id;
proc sort data=eqtl_ir_testable nodup;
   by feature_id;
   run;
   

/* Now get the list of tissue-specific exons */


data cell_specific_exons;
   set con.results_by_fusion_w_flags;
   if flag_cd19_on=1 and flag_cd4_on=0 and flag_cd8_on=0 then flag_cell_specific=1;
   else if flag_cd19_on=0 and flag_cd4_on=1 and flag_cd8_on=0 then flag_cell_specific=1;
   else if flag_cd19_on=0 and flag_cd4_on=0 and flag_cd8_on=1 then flag_cell_specific=1;
   else flag_cell_specific=0;
   if flag_cd19_on=0 and flag_cd4_on=0 and flag_cd8_on=0 then delete; *we don't need off fusions for these counts;
   keep fusion_id flag_cell_specific;
   rename fusion_id=feature_id;
   run;

/* Now get the list of tissue-specific AS events */

data cell_specific_as;
   set splicing.splicing_results_w_annot_fdr;
   if flag_intron_retention=0;
   if flag_cd19_on=1 and flag_cd4_on=0 and flag_cd8_on=0 then flag_cell_specific=1;
   else if flag_cd19_on=0 and flag_cd4_on=1 and flag_cd8_on=0 then flag_cell_specific=1;
   else if flag_cd19_on=0 and flag_cd4_on=0 and flag_cd8_on=1 then flag_cell_specific=1;
   else flag_cell_specific=0;
   if flag_cd19_on=0 and flag_cd4_on=0 and flag_cd8_on=0 then delete; *we don't need off events for these counts;
   keep event_id flag_cell_specific;
   rename event_id=feature_id;
   run;

/* Now get the list of tissue-specific IR events */

data cell_specific_ir;
   set splicing.splicing_results_w_annot_fdr;
   if flag_intron_retention=1;
   if flag_cd19_on=1 and flag_cd4_on=0 and flag_cd8_on=0 then flag_cell_specific=1;
   else if flag_cd19_on=0 and flag_cd4_on=1 and flag_cd8_on=0 then flag_cell_specific=1;
   else if flag_cd19_on=0 and flag_cd4_on=0 and flag_cd8_on=1 then flag_cell_specific=1;
   else flag_cell_specific=0;
   if flag_cd19_on=0 and flag_cd4_on=0 and flag_cd8_on=0 then delete; *we don't need off events for these counts;
   keep event_id flag_cell_specific;
   rename event_id=feature_id;
   run;   



/* Sort and merge */

proc sort data=cell_specific_exons nodup;
   by feature_id;
proc sort data=cell_specific_as nodup;
   by feature_id;
proc sort data=cell_specific_ir nodup;
   by feature_id;
   run;


   
proc sort data=eqtl_exon_testable nodup;
   by feature_id;
proc sort data=eqtl_as_testable nodup;
   by feature_id;
proc sort data=eqtl_ir_testable nodup;
   by feature_id;
   run;
   

data exon_test_specific;
   merge eqtl_exon_testable (in=in1) cell_specific_exons (in=in2);
   by feature_id;
   if in1 and in2 then output;
   else if in1 then output;
   else do;
      flag_testable=0;
      output;
      end;
   run;
   
data as_test_specific;
   merge eqtl_as_testable (in=in1) cell_specific_as (in=in2);
   by feature_id;
   if in1 and in2 then output;
   else if in1 then output;
   else do;
      flag_testable=0;
      output;
      end;
   run;

   
data ir_test_specific;
   merge eqtl_ir_testable (in=in1) cell_specific_ir (in=in2);
     if in1 and in2 then output;
   else if in1 then output;
   else do;
      flag_testable=0;
      output;
      end;
   run;
 
   
   /* Test cell specificity */
   
   proc freq data=exon_test_specific;
   tables flag_testable*flag_cell_specific;
   run;
   

/*
Table of flag_testable by flag_cell_specific

    flag_testable
              flag_cell_specific

    Frequency|
    Percent  |
    Row Pct  |
    Col Pct  |       0|       1|  Total
    ---------+--------+--------+
           0 | 165142 |  13335 | 178477
             |  86.29 |   6.97 |  93.25
             |  92.53 |   7.47 |
             |  93.07 |  95.61 |
    ---------+--------+--------+
           1 |  12299 |    613 |  12912
             |   6.43 |   0.32 |   6.75
             |  95.25 |   4.75 |
             |   6.93 |   4.39 |
    ---------+--------+--------+
    Total      177441    13948   191389
                92.71     7.29   100.00



*/


   proc freq data=as_test_specific;
   tables flag_testable*flag_cell_specific;
   run;
  

/*

  flag_testable
           flag_cell_specific

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 | 141556 |  13040 | 154596
          |  83.67 |   7.71 |  91.38
          |  91.57 |   8.43 |
          |  91.36 |  91.62 |
 ---------+--------+--------+
        1 |  13392 |   1192 |  14584
          |   7.92 |   0.70 |   8.62
          |  91.83 |   8.17 |
          |   8.64 |   8.38 |
 ---------+--------+--------+
 Total      154948    14232   169180
             91.59     8.41   100.00

*/


   
   proc freq data=ir_test_specific;
   tables flag_testable*flag_cell_specific;
   run;


/*
Table of flag_testable by flag_cell_specific

    flag_testable
              flag_cell_specific

    Frequency|
    Percent  |
    Row Pct  |
    Col Pct  |       0|       1|  Total
    ---------+--------+--------+
           0 |  33553 |   6065 |  39618
             |  77.93 |  14.09 |  92.02
             |  84.69 |  15.31 |
             |  92.32 |  90.36 |
    ---------+--------+--------+
           1 |   2790 |    647 |   3437
             |   6.48 |   1.50 |   7.98
             |  81.18 |  18.82 |
             |   7.68 |   9.64 |
    ---------+--------+--------+
    Total       36343     6712    43055
                84.41    15.59   100.00




*/
   


