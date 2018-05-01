ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Export genes, transcripts on/off */

* Genes;
proc export data=event.flag_gene_expressed
     outfile="!MCLAB/event_analysis/analysis_output/event_analysis_flag_gene_expressed_nomulti.csv"
     dbms=csv replace;
run;

* Transcripts;

proc export data=event.flag_xscript_w_gene_on
     outfile="!MCLAB/event_analysis/analysis_output/event_analysis_flag_transcript_w_gene_exp_nomulti.csv"
     dbms=csv replace;
run;

