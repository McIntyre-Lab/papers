/* Libraries */

libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname eventloc '/mnt/store/event_sandbox/sas_data';
ods listing; ods html close;

/* Export BED file for fusion genomic coordinates to BLAST against */

data mm10_fusions;
   set mm10.mm10_refseq_fusion_si_bed_v2;
   drop score strand;
run;

proc export data=mm10_fusions
     outfile="!MCLAB/event_analysis/analysis_output/mm10_fusions_si.bed"
     dbms=tab replace;
     putnames=no;
run;
