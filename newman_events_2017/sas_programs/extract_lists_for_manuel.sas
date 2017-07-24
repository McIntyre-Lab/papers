ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Extract transcript and gene lists:

List 1. The 17,050 transcripts with all unique feature detected, PLUS the 47,805 transcripts with no unique features (64,855 total transcripts)  
List 2. List 1, plus the 985 transcripts with at least 75% of unique features detected (65,840 total transcripts)
List 3. List 2, plus the 8,110 transcripts with at least 50% of unique features detected (73,950 total transcripts)
List 4. List 3, plus the 2,773 transcripts with at least 25% of unique features detected (76,723 total transcripts)
List 5. List 4, plus the 1,668 transcripts with >0% of unique features detected (78,391 total transcripts)
List 6. 30,586 transcripts with at least one unique feature detected.

List 7. List of genes with novel junctions
*/

proc freq data=event.feature_dtct_cnt_by_xscript_exp;
   tables bin_xscript_perc_uniq_dtct;
run;


data list1 list2 list3 list4 list5 list6;
   set event.feature_dtct_cnt_by_xscript_exp;
   if bin_xscript_perc_uniq_dtct = "no unique" then do;
      output list1;
      output list2;
      output list3;
      output list4;
      output list5;
      end;
   else if bin_xscript_perc_uniq_dtct = "0-25%" then do;
      output list5;
      output list6;
      end;
   else if bin_xscript_perc_uniq_dtct = "25-50%" then do;
      output list4;
      output list5;
      output list6;
      end;
   else if bin_xscript_perc_uniq_dtct = "50-75%" then do;
      output list3;
      output list4;
      output list5;
      output list6;
      end;
   else if bin_xscript_perc_uniq_dtct = "75-100%" then do;
      output list2;
      output list3;
      output list4;
      output list5;
      output list6;
      end;
   else if bin_xscript_perc_uniq_dtct = "100%" then do;
      output list1;
      output list2;
      output list3;
      output list4;
      output list5;
      output list6;
      end;
   keep transcript_id;
run;


data list7;
   set event.unannotated_events_by_gene;
   where flag_gene_has_unannotated_junc=1;
   keep gene_id;
run;


   data WORK.GENE2SYM    ;
   %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
   infile '!MCLAB/useful_mouse_data/mm10/gff/mm10_refseq_entrezID_symbol.txt' delimiter='09'x
MISSOVER DSD lrecl=32767 ;
      informat gene_id $12. ;
      informat gene_name $18. ;
      format gene_id $12. ;
      format gene_name $18. ;
   input
               gene_id $
               gene_name $
   ;
   if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
   run;

proc sort data=gene2sym nodup;
   by gene_id gene_name;
proc sort data=list7;
   by gene_id;
run;

data list7_w_sym;
  merge list7 (in=in1) gene2sym (in=in2);
  by gene_id;
  if in1 and in2 ;
  if gene_id="100504703" and gene_name="LOC100504703" then delete;
  if gene_id="102639598" and gene_name="LOC102639598" then delete;
  if gene_id="30963" and gene_name="Ptpla" then delete;
  if gene_id="57874" and gene_name="Ptplad1" then delete;
  if gene_id="59002" and gene_name="Wdr8" then delete;
run;

proc export data=list1 outfile="!MCLAB/event_analysis/analysis_output/list_1_transcripts_nomulti.tsv"
   dbms=tab replace;
run;

proc export data=list2 outfile="!MCLAB/event_analysis/analysis_output/list_2_transcripts_nomulti.tsv"
   dbms=tab replace;
run;

proc export data=list3 outfile="!MCLAB/event_analysis/analysis_output/list_3_transcripts_nomulti.tsv"
   dbms=tab replace;
run;

proc export data=list4 outfile="!MCLAB/event_analysis/analysis_output/list_4_transcripts_nomulti.tsv"
   dbms=tab replace;
run;

proc export data=list5 outfile="!MCLAB/event_analysis/analysis_output/list_5_transcripts_nomulti.tsv"
   dbms=tab replace;
run;

proc export data=list6 outfile="!MCLAB/event_analysis/analysis_output/list_6_transcripts_nomulti.tsv"
   dbms=tab replace;
run;

proc export data=list7_w_sym outfile="!MCLAB/event_analysis/analysis_output/list_7_genes_nomulti.tsv"
   dbms=tab replace;
run;


