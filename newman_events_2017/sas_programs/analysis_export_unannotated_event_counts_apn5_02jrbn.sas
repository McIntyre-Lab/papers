ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';


/* Export unannotated event counts */

data unannot_counts;
   set event.unannot_events_by_gene_apn5_v3;
   num_unannotated_events=num_unannotated_junctions+num_ir_events;
   if num_unannotated_events=0 then delete;
run;


proc export data=unannot_counts
            outfile="!MCLAB/event_analysis/analysis_output/event_analysis_number_of_unannotated_events_by_gene_nomulti_apn5.csv"
    dbms=csv replace;
run;



