ods listing; ods html close;
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';

/* Import iReckon results and look at by-bin CV */


   %macro import_data(sample);

proc import datafile="/mnt/store/event_sandbox/ireckon/&sample./result.gtf"
     out=&sample._ireckon dbms=tab replace;
     getnames=no; guessingrows=170000;
run;

/* Parse results GTF */

data &sample._ireckon2;
  set &sample._ireckon;
  where VAR3="transcript";
  length sample_id $4.;
  length gene_id $100.;
  length transcript_id $100.;
  format rpkm best32. ;
  format cov best32. ;
      sample_id="&sample.";
  gene_id=compress(scan(scan(VAR9,1,";"),2,'"'));
  transcript_id=compress(scan(scan(VAR9,2,";"),2,'"'));
  rpkm=compress(scan(scan(VAR9,3,";"),2,'"'))+0;
  cov=compress(scan(scan(VAR9,8,";"),2,'"'))+0;

      log_rpkm=log(rpkm + 1);
      /* Bin transcripts */
      if log_rpkm=0 then rpkm_bin=0;
      else if log_rpkm < 0.1 then rpkm_bin=1;
      else if log_rpkm < 0.6 then rpkm_bin=2;
      else if log_rpkm < 1.6 then rpkm_bin=3;
      else rpkm_bin=4;
      keep sample_id gene_id transcript_id log_rpkm rpkm rpkm_bin VAR1 VAR4 VAR5 VAR7;
  
  drop VAR9;
  rename VAR1=chr VAR4=start VAR5=stop VAR7=strand;
  run;

  %mend;
  %import_data(NSC1);
  %import_data(NSC2);

   data stack_ireckon;
      set NSC1_ireckon2 NSC2_ireckon2;
   run;
   
   data event.ireckon_results_nsc;
     set stack_ireckon;
   run; 

   /* Calc CV by sample and bin */
   proc sort data=stack_ireckon;
      by sample_id rpkm_bin;
   proc means data=stack_ireckon noprint;
      by sample_id rpkm_bin;
      var log_rpkm;
      output out=varstats_by_bin_sample cv=rpkm_cv;
   run;

   /* Now I need to transpose this, as I want bin as rows and samples as columns */

   proc sort data=varstats_by_bin_sample;
      by rpkm_bin sample_id;
   proc transpose data=varstats_by_bin_sample out=cv_sbys;
     by rpkm_bin;
     id sample_id;
     var rpkm_cv;
   run;

   proc transpose data=varstats_by_bin_sample out=count_sbys;
     by rpkm_bin;
     id sample_id;
     var _FREQ_;
   run;

   data cv_sbys2;
      set cv_sbys;
      rename NSC1=NSC1_CV NSC2=NSC2_CV;
      drop _NAME_ ;
   run;

   data count_sbys2;
      set count_sbys;
      rename NSC1=NSC1_count NSC2=NSC2_count;
      drop _NAME_ ;
   run;

   data  event.stats_ireckon;
      length xscript_set $32.;
      merge count_sbys2 cv_sbys2;
      by rpkm_bin;
      xscript_set="iReckon";
   run;

proc print data=event.stats_ireckon;
run;

/*

        xscript_                NSC1_    NSC2_
 Obs      set       rpkm_bin    count    count    NSC1_CV    NSC2_CV

  1     iReckon         1        9200     6146    91.9096    87.1610
  2     iReckon         2        3884     4129    53.8262    52.6500
  3     iReckon         3         874      969    28.9392    29.0458
  4     iReckon         4         626      602    41.2229    40.9589


*/

