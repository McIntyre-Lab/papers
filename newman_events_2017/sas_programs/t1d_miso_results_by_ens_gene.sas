
ods listing; ods html close;
libname con '!PATCON/sas_data';
libname event '!MCLAB/event_analysis/sas_data';

/* Flag gene significant if at least one event is different between 2 or more cell types */


data miso_results;
  length ens_gene_id $16.;
  set event.t1d_all_miso_results_v2;
  if gene_id="NA" then delete;
  do i=1 by 1 while(scan(gene_id,i,",") ^= "");
       ens_gene_id=scan(gene_id,i,",");
       output;
       end;
  drop gene_id i;
run;

proc sort data=miso_results;
  by ens_gene_id;
proc means data=miso_results noprint;
  by ens_gene_id;
  var flag_cd4_on flag_cd8_on flag_cd19_on flag_any_on flag_all_off flag_miso_testable
      flag_cd4cd8_testable flag_cd4cd19_testable flag_cd8cd19_testable flag_1cell_event flag_2cells_event
      flag_sig_cd4cd8_bf5 flag_sig_cd4cd19_bf5 flag_sig_cd8cd19_bf5 flag_sig_bf5
      flag_sig_cd4cd8_bf10 flag_sig_cd4cd19_bf10 flag_sig_cd8cd19_bf10 flag_sig_bf10;
  output out=miso_by_gene max=;
run;

proc freq data=miso_by_gene;
   tables flag_sig_bf5*flag_sig_bf10;
run;

/*
  Table of flag_sig_bf5 by flag_sig_bf10

   flag_sig_bf5     flag_sig_bf10

   Frequency|
   Percent  |
   Row Pct  |
   Col Pct  |       0|       1|  Total
   ---------+--------+--------+
          0 |  15146 |      0 |  15146
            |  99.74 |   0.00 |  99.74
            | 100.00 |   0.00 |
            |  99.84 |   0.00 |
   ---------+--------+--------+
          1 |     25 |     14 |     39
            |   0.16 |   0.09 |   0.26
            |  64.10 |  35.90 |
            |   0.16 | 100.00 |
   ---------+--------+--------+
   Total       15171       14    15185
               99.91     0.09   100.00


*/

/* Make permenant -- next step is to convert to Aceview IDs and match to T1D data */


data event.t1d_flag_miso_results_by_gene2; 
   set miso_by_Gene;
   drop _TYPE_ _FREQ_;
run;

