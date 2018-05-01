libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/splicing/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';

/* for each gene in the junction catalog, I am counting the number of different junction type combinations
   for selecting to test with MISO

   I want to limit this to genes with 2-10 transcripts, and with at least 5 annotated junctions detected */

data junc2keep;
  set eventloc.unique_junction2event_mm10;
  where flag_junction_nsc_on=1;
  keep event_id;
run;

data junc2gene;
   set evspl.splicing_events_annot_refseq;
   if flag_junction_annotated=1 or num_transcripts>0;
   keep gene_id event_id;
run;

proc sort data=junc2keep;
  by event_id;
proc sort data=junc2gene;
  by event_id;
run;

data junc2gene2;
  merge junc2gene (in=in1) junc2keep (in=in2);
  by event_id;
  if in1 and in2;
run;

proc freq data=junc2gene2 noprint;
  tables gene_id / out=genes2keep;
run;

data genes2keep2;
  set genes2keep;
  where count ge 5;
  keep gene_id;
run;


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
proc sort data=genes2keep2;
   by gene_id;
run;

data junc_cat_subset;
  merge junc_cat (in=in1) genes2keep (in=in2) genes2drop (in=in3) genes2keep2 (in=in4);
  by gene_id;
  if in3 then delete;
  else if in1 and in2 and in4 then output;
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
  if JA > 0 and JA_AA > 0 and JA_ES > 0 and JA_AD > 0 and JA_AD_AA > 0 and (JA_ES_AD > 0 or JA_ES_AA > 0 or
   JA_ES_AD_AA) > 0;
run;

/* 10 genes selected */

proc print data=subset (keep=gene_id);
run;

/*
gene_id
 101214
 106581
 108011
 114565
 114863
 11864
 12419
 12830
 142980
 16150
 16468
 170787
 17350
 20509
 217378
 22029
 22276
 225608
 22654
 228960
 231386
 232337
 234366
 234725
 269120
 27416
 28240
 56812
 58175
 66433
 66855
 66940
 67669
 69930
 70380
 70382
 70796
 71893
 72205
 72344
 72999
 73750
 74080
 74140
 74190
 76577
 77613
 77733
 94281
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
  merge xs_per_gene (in=in1) junc_of_interest (in=in2) genes2keep (in=in3) genes2drop (in=in4)  genes2keep2 (in=in5) ;
  by gene_id;
  if in4 then delete;
  if in1 and in2 and in3 and in5;
run;

data alt_junc_constit;
   set junc_subset2;
   if num_transcripts = count;
run;
ods listing;
proc freq data=alt_junc_constit;
   tables flag_alt_donor*flag_alt_acceptor;
run;

/*
 Table of flag_alt_donor by flag_alt_acceptor

     flag_alt_donor
               flag_alt_acceptor

     Frequency|
     Percent  |
     Row Pct  |
     Col Pct  |       0|       1|  Total
     ---------+--------+--------+
            0 |      0 |      9 |      9
              |   0.00 |  64.29 |  64.29
              |   0.00 | 100.00 |
              |   0.00 | 100.00 |
     ---------+--------+--------+
            1 |      5 |      0 |      5
              |  35.71 |   0.00 |  35.71
              | 100.00 |   0.00 |
              | 100.00 |   0.00 |
     ---------+--------+--------+
     Total           5        9       14
                 35.71    64.29   100.00

*/

/* gene list for ASTA/MISO testing:
 101214
 106581
 108011
 114565
 114863
 11864
 12419
 12830
 142980
 16150
 16468
 170787
 17350
 20509
 217378
 22029
 22276
 225608
 22654
 228960
 231386
 232337
 234366
 234725
 269120
 27416
 28240
 56812
 58175
 66433
 66855
 66940
 67669
 69930
 70380
 70382
 70796
 71893
 72205
 72344
 72999
 73750
 74080
 74140
 74190
 76577
 77613
 77733
 94281
*/


data subset_xs2gene;
  set xs2gene;
  if gene_id in ('101214','106581','108011','114565','114863','11864','12419',
                 '12830','142980','16150','16468','170787','17350','20509',
                 '217378','22029','22276','225608','22654','228960','231386',
                 '232337','234366','234725','269120','27416','28240','56812',
                 '58175','66433','66855','66940','67669','69930','70380','70382',
                 '70796','71893','72205','72344','72999','73750','74080','74140',
                 '74190','76577','77613','77733','94281');
  keep transcript_id;
run;

proc export data=subset_xs2gene 
     outfile="!MCLAB/event_analysis/design_files/subset_xscripts_for_miso_test2.txt" 
     dbms=tab replace;
     putnames=no;
run;



