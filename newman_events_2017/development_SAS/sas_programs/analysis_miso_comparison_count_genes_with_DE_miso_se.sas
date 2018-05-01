ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Check that the set of genes with DE MISO SE are the same set with DE/DD SE in Event analysis */



/* Get list of genes with MISO events and ES junctions detected/analyzed */

data genes_to_keep;
  set event.miso_refseq_exonskip_cmpr_dtct;
  where flag_has_refseq=1 and flag_has_miso_se_dtct=1 and
     (flag_exonskip_dtct_both=1 or flag_exonskip_dtct_one=1);
  keep gene_id;
run;

/* Merge with MISO */

data miso_diff_se;
   set event.miso_diff_se_refseq;
  where flag_refseq_match=1;
run;

proc sort data=miso_diff_se;
   by gene_id;
proc sort data=event.num_de_exonskip_by_gene;
   by gene_id;
proc sort data=genes_to_keep nodup;
   by gene_id;
run;

data miso_to_exonskip;
  merge miso_diff_se (in=in1) event.num_de_exonskip_by_gene (in=in2);
  by gene_id;
  if not in1 then do;
      num_diff_se_bf10=0;
      num_diff_se_bf5=0; end;
  if not in2 then do;
       flag_exonskip_de=0;
       flag_exonskip_dd=0;
  end;
  if in1 then output;
run;

data miso_to_exonskip2;
  merge miso_to_exonskip (in=in1) genes_to_keep (in=in2);
  by gene_id;
  if in1 and in2;
run; *1591 genes;

data flag_genes;
  set miso_to_exonskip2;
  if num_diff_se_bf10 > 0 then flag_miso_se_diff_bf10=1; else flag_miso_se_diff_bf10=0;
  if num_diff_se_bf5 > 0 then flag_miso_se_diff_bf5=1; else flag_miso_se_diff_bf5=0;
  if flag_exonskip_de > 1 then flag_event_se_de=1; else flag_event_se_de=0;
  if flag_exonskip_dd > 1 then flag_event_se_dd=1; else flag_event_se_dd=0;
run;

data event.miso_refseq_exonskip_cmpr_de;
   set flag_genes;
run;


proc freq data=event.miso_refseq_exonskip_cmpr_de;
   tables flag_miso_se_diff_bf10*flag_event_se_de
          flag_miso_se_diff_bf10*flag_event_se_dd
          flag_miso_se_diff_bf5*flag_event_se_de
          flag_miso_se_diff_bf5*flag_event_se_dd;
run;



/*
flag_miso_se_diff_bf10
          flag_event_se_dd

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 |    962 |    617 |   1579
         |  60.47 |  38.78 |  99.25
         |  60.92 |  39.08 |
         |  99.18 |  99.36 |
---------+--------+--------+
       1 |      8 |      4 |     12
         |   0.50 |   0.25 |   0.75
         |  66.67 |  33.33 |
         |   0.82 |   0.64 |
---------+--------+--------+
Total         970      621     1591
            60.97    39.03   100.00

 flag_miso_se_diff_bf5
           flag_event_se_dd

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |    959 |    614 |   1573
          |  60.28 |  38.59 |  98.87
          |  60.97 |  39.03 |
          |  98.87 |  98.87 |
 ---------+--------+--------+
        1 |     11 |      7 |     18
          |   0.69 |   0.44 |   1.13
          |  61.11 |  38.89 |
          |   1.13 |   1.13 |
 ---------+--------+--------+
 Total         970      621     1591
             60.97    39.03   100.00
*/

