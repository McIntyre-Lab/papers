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
  where flag_has_refseq=1 and flag_has_miso_se_dtct=1 ;
  keep gene_id flag_has_multiple_exons;
run;

/* Merge with MISO */

data miso_diff_se;
   set event.miso_diff_se_refseq;
  where flag_refseq_match=1;
run;

data as_genes;
   set event.flag_gene_alt_spliced;
run;


proc sort data=miso_diff_se nodup;
   by gene_id;
proc sort data=event.num_de_exonskip_by_gene;
   by gene_id;
proc sort data=genes_to_keep nodup;
   by gene_id;
run;

data miso_to_exonskip;
  merge miso_diff_se (in=in1) as_genes (in=in2);
  by gene_id;
  if not in1 then do;
      num_diff_se_bf10=0;
      num_diff_se_bf5=0; end;
  if not in2 then do;
       flag_gene_as=0;
  end;
  if in1 then output;
run;

data miso_to_exonskip2;
  merge miso_to_exonskip (in=in1) genes_to_keep (in=in2);
  by gene_id;
  if in1 and in2;
run; *1129 genes;

data flag_genes;
  set miso_to_exonskip2;
  if num_diff_se_bf10 > 0 then flag_miso_se_diff_bf10=1; else flag_miso_se_diff_bf10=0;
  if num_diff_se_bf5 > 0 then flag_miso_se_diff_bf5=1; else flag_miso_se_diff_bf5=0;
run;

data event.miso_refseq_exonskip_cmpr_as;
   set flag_genes;
run;


proc freq data=event.miso_refseq_exonskip_cmpr_as;
   tables flag_miso_se_diff_bf10*flag_gene_as
          flag_miso_se_diff_bf5*flag_gene_as;
run;



/*
  flag_miso_se_diff_bf10
            flag_gene_as

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |     52 |   1070 |   1122
           |   4.61 |  94.77 |  99.38
           |   4.63 |  95.37 |
           | 100.00 |  99.35 |
  ---------+--------+--------+
         1 |      0 |      7 |      7
           |   0.00 |   0.62 |   0.62
           |   0.00 | 100.00 |
           |   0.00 |   0.65 |
  ---------+--------+--------+
  Total          52     1077     1129
               4.61    95.39   100.00

 flag_miso_se_diff_bf5
           flag_gene_as

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |     52 |   1065 |   1117
          |   4.61 |  94.33 |  98.94
          |   4.66 |  95.34 |
          | 100.00 |  98.89 |
 ---------+--------+--------+
        1 |      0 |     12 |     12
          |   0.00 |   1.06 |   1.06
          |   0.00 | 100.00 |
          |   0.00 |   1.11 |
 ---------+--------+--------+
 Total          52     1077     1129
              4.61    95.39   100.00


*/


proc freq data=event.miso_refseq_exonskip_cmpr_as;
   where flag_has_multiple_exons=1;
   tables flag_miso_se_diff_bf10*flag_gene_as
          flag_miso_se_diff_bf5*flag_gene_as;
run;

/*
 flag_miso_se_diff_bf10
           flag_gene_as

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |     39 |   1069 |   1108
          |   3.50 |  95.87 |  99.37
          |   3.52 |  96.48 |
          | 100.00 |  99.35 |
 ---------+--------+--------+
        1 |      0 |      7 |      7
          |   0.00 |   0.63 |   0.63
          |   0.00 | 100.00 |
          |   0.00 |   0.65 |
 ---------+--------+--------+
 Total          39     1076     1115
              3.50    96.50   100.00


  flag_miso_se_diff_bf5
            flag_gene_as

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |     39 |   1064 |   1103
           |   3.50 |  95.43 |  98.92
           |   3.54 |  96.46 |
           | 100.00 |  98.88 |
  ---------+--------+--------+
         1 |      0 |     12 |     12
           |   0.00 |   1.08 |   1.08
           |   0.00 | 100.00 |
           |   0.00 |   1.12 |
  ---------+--------+--------+
  Total          39     1076     1115
               3.50    96.50   100.00

*/

