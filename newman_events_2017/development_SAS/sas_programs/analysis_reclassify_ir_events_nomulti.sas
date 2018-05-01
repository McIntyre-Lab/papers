ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Take reclassified IR events and drop any from genes we aren't analyzing */


data intron_class;
  set event.ir_reclassification_v2;
run;

data genes2keep;
   set event.flag_gene_expressed;
   where flag_gene_expressed=1;
   keep gene_id;
run;

proc sort data=intron_class;
  by gene_id;
proc sort data=genes2keep;
  by gene_id;
run;

data intron_class_nomulti;
   merge intron_class (in=in1) genes2keep (in=in2);
   by gene_id;
   if in1 and in2;
run; *154934 IR events kept;



proc freq data=intron_class_nomulti noprint;
   where flag_fusion_on=1 and flag_splicing_on=1;
   tables flag_low_expressed*flag_possible_novel_donor*flag_possible_ir*flag_possible_unprocessed / out=ir_class;
run;

proc print data=ir_class;
run;

/*
                                   flag_
  flag_low_    flag_possible_    possible_    flag_possible_
  expressed      novel_donor         ir         unprocessed     COUNT    PERCENT

      0               0              0               1          24794    97.7990
      0               0              1               0            200     0.7889
      0               1              0               1            314     1.2386
      0               1              1               0             44     0.1736

(not losing too many here)

24794 possible unprocessed transcript
200 possible IR
314 possible novel donor
44 ambiguous IR

*/

data data_for_plots;
  set intron_class_nomulti;
  keep event_id flag_no_computed_intron mean_apn_intron flag_splicing_on mean_apn_ir
       flag_fusion_on mean_apn_fusion flag_low_expressed flag_possible_novel_donor flag_possible_ir flag_possible_unprocessed ;
run;

proc export data=data_for_plots
    outfile="!MCLAB/event_analysis/analysis_output/event_analysis_reclassified_ir_nomulti.csv"
   dbms=csv replace;
run;


/* Make permenant */

data event.ir_reclassification_nomulti;
   set intron_class_nomulti;
run;
