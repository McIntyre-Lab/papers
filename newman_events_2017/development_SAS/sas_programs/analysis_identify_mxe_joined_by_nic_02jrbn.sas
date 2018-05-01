/* Libraries */
libname event '!MCLAB/event_analysis/sas_data';
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
ods listing; ods html close;

/* Need to ID instances were a NIC junction suggests that two mutually-exclusive exons are joined.

   (1) ID constitutive and non-constitutive exons (or better: ID exons with constitutive/non-con splice sites)
   (2) ID detected NIC junctions
   (3) Remove NIC junctions where the donor/acceptor exon is constitutive. Remaining junctions should only be
       between alternative/common transcripts
   (4) Merge in donor and acceptor exons, and their transcripts
   (5) For each donor site, iterate through transcripts and check if any of the transcripts that use that donor
       use the acceptor site
   (6) Repeat (5) but for exons
   (7) Count instances where the NIC junction is detected/supported (APN 0,2,5)
*/

/*** (1) Identify constitutive and non-constitutive splice sites ***/

/* (a) count transcripts per gene */

data xs2gene;
  set mm10.mm10_exons_w_info;
  length transcript_id2 $20.;
  do i=1 by 1 while(scan(transcript_id,i,"|") ^= "");
     transcript_id2=scan(transcript_id,i,"|");
     output;
     end;
  keep gene_id transcript_id2;
run;

proc sort data=xs2gene nodup;
   by gene_id transcript_id2;
proc freq data=xs2gene noprint;
   tables gene_id / out=xs_per_gene;
run;

data xs_per_gene2;
  set xs_per_gene;
  xs_per_gene=COUNT;
  drop COUNT PERCENT;
run;

/* (b) count transcripts per donor and acceptor */

data donor2xs;
   set mm10.mm10_exons_w_info;
  length transcript_id2 $20.;
  do i=1 by 1 while(scan(transcript_id,i,"|") ^= "");
     transcript_id2=scan(transcript_id,i,"|");
     output;
     end;
  keep gene_id stop transcript_id2;
  rename stop=donor_site;
run;

proc sort data=donor2xs nodup;
   by gene_id donor_site  transcript_id2;
proc freq data=donor2xs noprint;
   by gene_id;
   tables donor_site / out=xs_per_donor;
run;

data xs_per_donor2;
  set xs_per_donor;
  xs_per_donor=COUNT;
  drop COUNT PERCENT;
run;

data acceptor2xs;
   set mm10.mm10_exons_w_info;
  length transcript_id2 $20.;
  do i=1 by 1 while(scan(transcript_id,i,"|") ^= "");
     transcript_id2=scan(transcript_id,i,"|");
     output;
     end;
  keep gene_id start transcript_id2;
  rename start=acceptor_site;
run;

proc sort data=acceptor2xs nodup;
   by gene_id acceptor_site transcript_id2;
proc freq data=acceptor2xs noprint;
   by gene_id;
   tables acceptor_site / out=xs_per_acceptor;
run;

data xs_per_acceptor2;
  set xs_per_acceptor;
  xs_per_acceptor=COUNT;
  drop COUNT PERCENT;
run;

/* (c) Flag donors and acceptors as constitutive or not */

proc sort data=xs_per_gene2;
   by gene_id;
proc sort data=xs_per_donor2;
   by gene_id;
proc sort data=xs_per_acceptor2;
   by gene_id;
run;

data flag_constit_donor;
  merge xs_per_donor2 (in=in1) xs_per_gene2 (in=in2);
  by gene_id;
  if in1 and in2;
  if xs_per_donor=xs_per_gene then flag_constit_donor=1;
  else flag_constit_donor=0;
run;

data flag_constit_acceptor;
  merge xs_per_acceptor2 (in=in1) xs_per_gene2 (in=in2);
  by gene_id;
  if in1 and in2;
  if xs_per_acceptor=xs_per_gene then flag_constit_acceptor=1;
  else flag_constit_acceptor=0;
run;

/* (2) ID detected NIC junctions */

data nic_on;
  set event.event2star2pacbio_junc_table;
  where junction_type ? "NIC" and flag_in_catalog=1;
  keep junction_id chr strand donor_stop acceptor_start flag_events_detected
       flag_events_NSC_apn_gt0 flag_events_NSC_apn_ge2 flag_events_NSC_apn_ge5;
run;

data junc2event;
  set eventloc.unique_junction2event_mm10;
  keep event_id junction_id;
run;

proc sort data=nic_on;
  by junction_id;
proc sort data=junc2event;
  by junction_id;
run;

data nic_on2;
  merge nic_on (in=in1) junc2event (in=in2);
  by junction_id;
  if in1 and in2;
  rename donor_stop=donor_site acceptor_start=acceptor_site;
run;

data nic_on3;
   length gene_id $12.;
   set nic_on2;
   gene_id=scan(event_id,1,":");
run;

/* (3) Remove NIC junctions where the donor/acceptor exon is constitutive. Remaining junctions should only be
   between alternative/common transcripts */

data constit_donor;
  set flag_constit_donor;
  where flag_constit_donor=1;
  keep gene_id donor_site;
run;

proc sort data=nic_on3  nodup;
  by gene_id donor_site;
proc sort data=constit_donor nodup;
  by gene_id donor_site;
run;

data nic_alt_donor;
   merge nic_on3 (in=in1) constit_donor (in=in2);
   by gene_id donor_site;
   if in1 and in2 then delete;
   else if in1 then output;
run;

data constit_acceptor;
  set flag_constit_acceptor;
  where flag_constit_acceptor=1;
  keep gene_id acceptor_site;
run;

proc sort data=nic_alt_donor  nodup;
  by gene_id acceptor_site;
proc sort data=constit_acceptor nodup;
  by gene_id acceptor_site;
run;

data nic_alt_acceptor;
   merge nic_alt_donor (in=in1) constit_acceptor (in=in2);
   by gene_id acceptor_site;
   if in1 and in2 then delete;
   else if in1 then output;
run;

/* (4) Merge in donor and acceptor exons, and their transcripts */

data donor_exons;
   set mm10.mm10_exons_w_info;
   length transcript_id2 $20.;
   do i = 1 by 1 while(scan(transcript_id,i,"|") ^= "");
       transcript_id2=scan(transcript_id,i,"|");
       output;
       end;
   keep stop exon_id gene_id transcript_id2;
   rename stop=donor_site transcript_id2=transcript_id;
run;

data acceptor_exons;
   set mm10.mm10_exons_w_info;
   length transcript_id2 $20.;
   do i = 1 by 1 while(scan(transcript_id,i,"|") ^= "");
       transcript_id2=scan(transcript_id,i,"|");
       output;
       end;
   keep start exon_id gene_id transcript_id2;
   rename start=acceptor_site transcript_id2=transcript_id;
run;

/* (a) For each gene*donor/acceptor site, concatenate transcripts */

data donor2xs;
   set donor_exons;
   drop exon_id;
run;

data acceptor2xs;
   set acceptor_exons;
   drop exon_id;
run;

proc sort data=donor2xs nodup;
   by gene_id donor_site transcript_id;
proc freq data=donor2xs noprint;
   by gene_id ;
   tables donor_site / out=xs_per_donor;
proc sort data=xs_per_donor;
   by descending count;
proc print data=xs_per_donor(obs=1);
run; *130 transcripts per donor site;

proc sort data=acceptor2xs nodup;
   by gene_id acceptor_site transcript_id;
proc freq data=acceptor2xs noprint;
   by gene_id ;
   tables acceptor_site / out=xs_per_acceptor;
proc sort data=xs_per_acceptor;
   by descending count;
proc print data=xs_per_acceptor(obs=1);
run; *130 transcripts per acceptor site;

data donor2xs_cat; 
  array xs[130] $ 25.;
  retain xs1-xs130;
  set donor2xs;
  by gene_id donor_site;
  if first.donor_site then do;
     call missing(of xs1-xs130);
     records = 0;
  end;
  records + 1;
  xs[records]=transcript_id;
  if last.donor_site then output;
run;

  *clean up the output file;
data donor2xs_cat2;
  set donor2xs_cat;
  length donor_transcript_id $ 2000.;
  rename records= num_donor_xscript;
          donor_transcript_id= catx("|", OF xs1-xs130);
  drop xs1-xs130 transcript_id;
  run;

data acceptor2xs_cat; 
  array xs[130] $ 25.;
  retain xs1-xs130;
  set acceptor2xs;
  by gene_id acceptor_site;
  if first.acceptor_site then do;
     call missing(of xs1-xs130);
     records = 0;
  end;
  records + 1;
  xs[records]=transcript_id;
  if last.acceptor_site then output;
run;

  *clean up the output file;
data acceptor2xs_cat2;
  set acceptor2xs_cat;
  length acceptor_transcript_id $ 2000.;
  rename records= num_acceptor_xscript;
          acceptor_transcript_id= catx("|", OF xs1-xs130);
  drop xs1-xs130 transcript_id;
  run;

/* (b) Merge in the set of donor and acceptor transcripts */

proc sort data=nic_alt_acceptor;
  by gene_id donor_site;
proc sort data=donor2xs_cat2;
  by gene_id donor_site;
run;

data nic_donor2xs;
  merge nic_alt_acceptor (in=in1) donor2xs_cat2 (in=in2);
  by gene_id donor_site;
  if in1 and in2;
run;

proc sort data=nic_donor2xs;
  by gene_id acceptor_site;
proc sort data=acceptor2xs_cat2;
  by gene_id acceptor_site;
run;

data nic_acceptor2xs;
  merge nic_donor2xs (in=in1) acceptor2xs_cat2 (in=in2);
  by gene_id acceptor_site;
  if in1 and in2;
run;


/* (5) For each donor site, iterate through transcripts and check if any of the transcripts that use that donor
       use the acceptor site */

data flag_shared_xs;
   set nic_acceptor2xs;
   num_match=0;
   do i = 1 by 1 while(scan(donor_transcript_id,i,"|") ^= "");
      do j = 1 by 1 while(scan(acceptor_transcript_id,j,"|") ^= "");
         if scan(donor_transcript_id,i,"|") = scan(acceptor_transcript_id,j,"|") 
           then num_match=num_match +1;
      end;
   end;
  if num_match = 0 then flag_donor_acceptor_mxe=1;
  else flag_donor_acceptor_mxe=0;
run;

/* (6) Count instances where the NIC junction is detected/supported (APN 0,2,5) */

proc freq data=flag_shared_xs;
   tables flag_donor_acceptor_mxe
          flag_donor_acceptor_mxe*flag_events_detected
          flag_donor_acceptor_mxe*flag_events_nsc_apn_gt0
          flag_donor_acceptor_mxe*flag_events_nsc_apn_ge2
          flag_donor_acceptor_mxe*flag_events_nsc_apn_ge5;
run;

/* 

     flag_donor_                             Cumulative    Cumulative
    acceptor_mxe    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0      445145       67.79        445145        67.79
               1      211491       32.21        656636       100.00

211491 mutually-exclusive exon pairs in total. Note that this is ANY mutually exclusive pair, not necessarily those
that are located sequentially 

flag_donor_acceptor_mxe
          flag_events_detected

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 | 443403 |   1742 | 445145
         |  67.53 |   0.27 |  67.79
         |  99.61 |   0.39 |
         |  67.97 |  40.47 |
---------+--------+--------+
       1 | 208929 |   2562 | 211491
         |  31.82 |   0.39 |  32.21
         |  98.79 |   1.21 |
         |  32.03 |  59.53 |
---------+--------+--------+
Total      652332     4304   656636
            99.34     0.66   100.00

2562 NIC junctions joining MXEs detected


 flag_donor_acceptor_mxe
           flag_events_NSC_apn_gt0

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 | 444807 |    338 | 445145
          |  67.74 |   0.05 |  67.79
          |  99.92 |   0.08 |
          |  67.86 |  28.74 |
 ---------+--------+--------+
        1 | 210653 |    838 | 211491
          |  32.08 |   0.13 |  32.21
          |  99.60 |   0.40 |
          |  32.14 |  71.26 |
 ---------+--------+--------+
 Total      655460     1176   656636
             99.82     0.18   100.00

838 NIC junctions joining MXEs with support at APN>0

 flag_donor_acceptor_mxe
           flag_events_NSC_apn_ge2

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 | 445087 |     58 | 445145
          |  67.78 |   0.01 |  67.79
          |  99.99 |   0.01 |
          |  67.80 |  33.92 |
 ---------+--------+--------+
        1 | 211378 |    113 | 211491
          |  32.19 |   0.02 |  32.21
          |  99.95 |   0.05 |
          |  32.20 |  66.08 |
 ---------+--------+--------+
 Total      656465      171   656636
             99.97     0.03   100.00

113 NIC junctions joining MXEs with support at APN>=2

 flag_donor_acceptor_mxe
           flag_events_NSC_apn_ge5

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 | 445125 |     20 | 445145
          |  67.79 |   0.00 |  67.79
          | 100.00 |   0.00 |
          |  67.80 |  29.85 |
 ---------+--------+--------+
        1 | 211444 |     47 | 211491
          |  32.20 |   0.01 |  32.21
          |  99.98 |   0.02 |
          |  32.20 |  70.15 |
 ---------+--------+--------+
 Total      656569       67   656636
             99.99     0.01   100.00

47 NIC junctions joining MXEs with support at APN>=5
*/

/* Make permenant */

data event.nic_mxe_with_support;
  set flag_shared_xs;
run;

