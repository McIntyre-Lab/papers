libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/splicing/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';

/* for each gene in the junction catalog, I am counting the number of different junction type combinations
   for selecting to test with ASTA

   I want to limit this to genes with 2-10 transcripts */

data genes2drop;
  set mm10.mm10_fusion_2_gene_id;
  where flag_multigene=1;
  keep gene_id;
run;

proc sort data=genes2drop nodup;
  by gene_id;
run;


data xs2gene;
   set mm10.mm10_fusion_2_gene_id;
    keep gene_id transcript_id;
run;

proc sort data=xs2gene nodup;
  by gene_id transcript_id;
proc freq data=xs2gene noprint;
  tables gene_id / out = xs_per_gene;
run;

data genes2keep;
  set xs_per_gene;
  if count < 2 then delete;
  if count > 10 then delete;
  keep gene_id;
run;


data junc_cat;
   set evspl.splicing_Events_annot_refseq;
   if num_transcripts > 0 or transcript_id ^= '' then flag_junction_annotated=1;
   keep event_id gene_id flag_junction_annotated flag_intron_retention flag_exonskip
        flag_alt_donor flag_alt_acceptor;
run;

proc sort data=junc_cat;
   by gene_id;
proc sort data=genes2keep;
   by gene_id;
proc sort data=genes2drop;
   by gene_id;
run;

data junc_cat_subset;
  merge junc_cat (in=in1) genes2keep (in=in2) genes2drop (in=in3);
  by gene_id;
  if in3 then delete;
  else if in1 and in2 then output;
run;

data junc_types;
   length junc_combo $20.;
   set junc_cat_subset;
   if flag_intron_retention=1 then junc_combo="IR";
   else do;
      if flag_junction_annotated=1 then do;
         if flag_exonskip=1 then do;
            if flag_alt_donor=1 then do;
               if flag_alt_acceptor=1 then junc_combo="JA_ES_AD_AA";
               else junc_combo="JA_ES_AD";
            end;
            else do;
               if flag_alt_acceptor=1 then junc_combo="JA_ES_AA";
               else junc_combo="JA_ES";
            end;
         end;
         else do;
            if flag_alt_donor=1 then do;
               if flag_alt_acceptor=1 then junc_combo="JA_AD_AA";
               else junc_combo="JA_AD";
            end;
            else do;
               if flag_alt_acceptor=1 then junc_combo="JA_AA";
               else junc_combo="JA";
            end;
         end;
      end;
      else do;
         if flag_exonskip=1 then do;
            if flag_alt_donor=1 then do;
               if flag_alt_acceptor=1 then junc_combo="JU_ES_AD_AA";
               else junc_combo="JU_ES_AD";
            end;
            else do;
               if flag_alt_acceptor=1 then junc_combo="JU_ES_AA";
               else junc_combo="JU_ES";
            end;
         end;
         else do;
            if flag_alt_donor=1 then do;
               if flag_alt_acceptor=1 then junc_combo="JU_AD_AA";
               else junc_combo="JU_AD";
            end;
            else do;
               if flag_alt_acceptor=1 then junc_combo="JU_AA";
               else junc_combo="JU";
            end;
         end;
      end;
   end;
run;

proc sort data=junc_types;
  by gene_id junc_combo;
proc freq data=junc_types noprint;
  by gene_id;
  tables junc_combo / out=junc_by_gene;
run;

proc transpose data=junc_by_gene out=junc_type_sbys;
   by gene_id;
   id junc_combo;
   var count;
run;

data junc_type_sbys2;
   set junc_type_sbys;
   array change _numeric_;
   do over change;
     if change=. then change=0;
   end;
run;

data subset;
  set junc_type_sbys2;
  if JA > 0 and JA_AA > 0 and JA_ES > 0 and JA_AD > 0 and JA_AD_AA > 0 and (JA_ES_AD > 0 or JA_ES_AA > 0) and
   JA_ES_AD_AA > 0;
run;

/* 10 genes selected */

proc print data=subset (keep=gene_id);
run;

/*
gene_id
21419
217378
234725
24058
435802
666528
70375
71893
74140
74190
*/

/* okay, need at least one gene with a junction which we call “alternative” yet the donor and acceptor combination is constitutively joined */

data junc_of_interest;
  set evspl.splicing_Events_annot_refseq;
   if num_transcripts > 0 or transcript_id ^= '' then flag_junction_annotated=1;
   if (flag_alt_donor=1 or flag_alt_acceptor=1) and flag_junction_annotated=1;
   keep event_id num_transcripts gene_id flag_alt_donor flag_alt_acceptor;
run;

proc sort data=xs_per_gene;
   by gene_id;
proc sort data=junc_of_interest;
   by gene_id;
run;

data junc_subset2;
  merge xs_per_gene (in=in1) junc_of_interest (in=in2) genes2keep (in=in3) genes2drop (in=in4);
  by gene_id;
  if in4 then delete;
  if in1 and in2 and in3 ;
run;

data alt_junc_constit;
   set junc_subset2;
   if num_transcripts = count;
run;

proc freq data=alt_junc_constit;
   tables flag_alt_donor*flag_alt_acceptor;
run;

/*

flag_alt_donor
          flag_alt_acceptor

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 |      0 |     28 |     28
         |   0.00 |  45.90 |  45.90
         |   0.00 | 100.00 |
         |   0.00 | 100.00 |
---------+--------+--------+
       1 |     33 |      0 |     33
         |  54.10 |   0.00 |  54.10
         | 100.00 |   0.00 |
         | 100.00 |   0.00 |
---------+--------+--------+
Total          33       28       61
            54.10    45.90   100.00

None with both. Take 1 example from each: 

12564
12722

*/

/* gene list for ASTA/MISO testing:
21419
217378
234725
24058
435802
666528
70375
71893
74140
74190
12564
12722
*/


data subset_xs2gene;
  set xs2gene;
  if gene_id in ('21419','217378','234725','24058','435802','666528',
                 '70375','71893','74140','74190','12564','12722');
  keep transcript_id;
run;

proc export data=subset_xs2gene 
     outfile="!MCLAB/event_analysis/design_files/subset_xscripts_for_miso_test.txt" 
     dbms=tab replace;
     putnames=no;
run;



