				CLASHING ID
PB.4412_chr13:1_2	80	PB.4412:1_10|2_10
PB.4412_chr13:1_intron	79	NONE
PB.5261_chr18:1_2	66	PB.5261:1_10|1_11|1_12|1_9|2_10|2_11|2_12|2_9|4_10|4_11|4_12|4_9
PB.5261_chr18:intron_2	65	


data check;
  set splicing_info;
  where (gene_id ? "PB.4412" or gene_id ? "PB.5261") and event_id ? "intron"; 
run;


libname conesa '!MCLAB/conesa_pacbio/sas_data';

/* Import fusion and splicing counts, rename some splicing events, reformat and export

   FORMAT: sample_id fusion_id/event_id apn
   split on condition (NSC, OLD) */

data fusions_nsc;
   set conesa.coverage_fusions;
   if sample_id='NSC1' or sample_id='NSC2';
   keep sample_id fusion_id apn;
run;

data fusions_old;
   set conesa.coverage_fusions;
   if sample_id='OLD1' or sample_id='OLD2';
   keep sample_id fusion_id apn;
run;

data splicing_nsc;
   set conesa.coverage_splicing;
   if sample_id='NSC1' or sample_id='NSC2';
   if event_id="PB.4412:1_2" then event_id="PB.4412_chr13:1_2";
   if event_id="PB.4412:1_intron" then event_id="PB.4412_chr13:1_intron";
   if event_id="PB.5261:1_2" then event_id="PB.5261_chr18:1_2";
   if event_id="PB.5261:intron_2" then event_id="PB.5261_chr18:intron_2";
   keep sample_id fusion_id apn;
   rename fusion_id=event_id;
run;

data splicing_old;
   set conesa.coverage_splicing;
   if sample_id='OLD1' or sample_id='OLD2';
   if event_id="PB.4412:1_2" then event_id="PB.4412_chr13:1_2";
   if event_id="PB.4412:1_intron" then event_id="PB.4412_chr13:1_intron";
   if event_id="PB.5261:1_2" then event_id="PB.5261_chr18:1_2";
   if event_id="PB.5261:intron_2" then event_id="PB.5261_chr18:intron_2";
   keep sample_id fusion_id apn;
   rename fusion_id=event_id;
run;

proc export data=fusions_nsc
     outfile='!MCLAB/conesa_pacbio/analysis_output/fusion_coverage_nsc.csv'
     dbms=csv replace;
run;

proc export data=fusions_old
     outfile='!MCLAB/conesa_pacbio/analysis_output/fusion_coverage_old.csv'
     dbms=csv replace;
run;

proc export data=splicing_nsc
     outfile='!MCLAB/conesa_pacbio/analysis_output/splicing_coverage_nsc.csv'
     dbms=csv replace;
run;

proc export data=splicing_old
     outfile='!MCLAB/conesa_pacbio/analysis_output/splicing_coverage_old.csv'
     dbms=csv replace;
run;

