ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Check that the set of genes with DE MISO SE are the same set with DE/DD SE in Event analysis */

/* Calculate FDR */

data anova_results;
  set event.exonskip_anova;
  where effect="cell_type";
  keep event_id probf;
run;

proc multtest inpvalues(ProbF)=anova_results fdr noprint
   out=anova_results_w_fdr;
run;


/* Make FDR output permenant */
data event.exonskip_anova_w_fdr;
   set anova_results_w_fdr;
run;


data on_flags;
   set event.flag_splicing_on;
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

data event.num_de_exonskip_by_gene;
   set num_de_dd_exonskip_by_gene;
run;

