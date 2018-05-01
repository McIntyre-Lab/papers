libname event "!MCLAB/event_analysis/sas_data";
libname mm10 "!MCLAB/useful_mouse_data/mm10/sas_data";
ods listing; ods html close;

/* Import SAM entries for STAR-only junctions to see what they are doing */

proc import datafile="!MCLAB/event_analysis/analysis_output/read_check_output_apn0.csv"
     out=star_sam dbms=csv replace;
     guessingrows=max; getnames=no;
run;

data star_sam2;
  set star_sam;
  rename VAR1=junction_id VAR2=junction_seq VAR3=seq_name VAR4=intron_len
         VAR5=read_id VAR6=sam_flag VAR7=chr VAR8=read_start
         VAR9=cigar_string VAR10=read_seq;
run;

/* Flag if the SAM file has more than 1 gap -- this is indicative that the read
   is likely split over multiple splice sites */

data flag_cigar;
   set star_sam2;
   if count(cigar_string,"N") > 1 then flag_multi_junc=1;
   else flag_multi_junc=0;
run;

proc freq data=flag_cigar;
   tables flag_multi_junc;
run;

/*


                                             Cumulative    Cumulative
 flag_multi_junc    Frequency     Percent     Frequency      Percent
 --------------------------------------------------------------------
               0         435       77.54           435        77.54
               1         126       22.46           561       100.00


A fifth of all extracted reads overlay multiple splice junctions
In terms of an individual junction... */

data multijunc monojunc;
  set flag_cigar;
  if flag_multi_junc=1 then output multijunc;
  else output monojunc;
  keep junction_id;
run;

proc sort data=multijunc nodup;
  by junction_id;
proc sort data=monojunc nodup;
  by junction_id;
run;

data multi_v_mono_junc;
  merge monojunc multijunc (in=in2);
  by junction_id;
  if in2 then flag_multi_junc=1;
  else flag_multi_junc=0;
run;
proc freq data=multi_v_mono_junc;
   tables flag_multi_junc;
run;

/*
                                             Cumulative    Cumulative
 flag_multi_junc    Frequency     Percent     Frequency      Percent
 --------------------------------------------------------------------
               0          61       76.25            61        76.25
               1          19       23.75            80       100.00

19 of 80 STAR-only junctions are due to the read mapping to multiple junction sites
This can only really occur if the "middle" exon is smaller than the read, so this appears
to be a microexon issue */

/* Double check: for multi-junction reads, what is the size distribution of the middle exon? */

data multijunc_reads;
  set flag_cigar;
  where flag_multi_junc=1;
  middle_exon_size=scan(scan(cigar_string,2,"N"),1,"M") + 0;
run;

proc freq data=multijunc_reads;
   tables middle_exon_size;
run;
 
/*


 middle_exon_size    Frequency     Percent
 -------------------------------------------
               11          87       69.05
               17           1        0.79
               18           4        3.17
               20          20       15.87
               21           9        7.14
               22           4        3.17
               23           1        0.79

Microexons. Do they match to known microexons?
*/

data multijunc_exon_coord;
   set multijunc_reads;
   donor1_size=scan(cigar_string,1,"M") + 0;
   intron1_size=scan(scan(cigar_string,1,"N"),2,"M") + 0;
   exon_start=read_start+donor1_size+intron1_size - 1;
   exon_stop=read_start+donor1_size+intron1_size+middle_exon_size - 1;
   keep junction_id read_id chr exon_start exon_stop;
run;

data mm10_exons;
  set mm10.mm10_exons_w_info;
  keep chrom start stop;
  rename chrom=chr start=exon_start stop=exon_stop;
run;

proc sort data=mm10_exons nodup;
   by chr exon_start exon_stop;
proc sort data=multijunc_exon_coord;
   by chr exon_start exon_stop;
run;

data flag_microexon_match;
   merge multijunc_exon_coord (in=in1) mm10_exons (in=in2);
   by chr exon_start exon_stop;
   if in1 and in2 then flag_microexon_match=1;
   else flag_microexon_match=0;
   if in1 then output;
run;

proc freq data=flag_microexon_match;
  tables flag_microexon_match;
run;

/*
  flag_microexon_                             Cumulative    Cumulative
            match    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                1         126      100.00           126       100.00

All match to a known microexon
*/


/* Now for the other junction reads. First check for insertions and deletions */

data monojunc_indel;
  set flag_cigar;
  where flag_multi_junc=0;
  if count(cigar_string,"I") > 0 then flag_insertion=1;
  else flag_insertion=0;
  if count(cigar_string,"D") > 0 then flag_deletion=1;
  else flag_deletion=0;
run;

proc freq data=monojunc_indel;
  tables flag_insertion*flag_deletion;
run;

/*
  flag_insertion     flag_deletion

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |    430 |      4 |    434
           |  98.85 |   0.92 |  99.77
           |  99.08 |   0.92 |
           |  99.77 | 100.00 |
  ---------+--------+--------+
         1 |      1 |      0 |      1
           |   0.23 |   0.00 |   0.23
           | 100.00 |   0.00 |
           |   0.23 |   0.00 |
  ---------+--------+--------+
  Total         431        4      435
              99.08     0.92   100.00

One read with an insertion, 4 with a deletion

Check that these are the only reads for these junctions
*/

data insertion deletion other;
  set monojunc_indel;
  if flaG_insertion=1 then output insertion;
  else if flaG_deletion=1 then output deletion;
  else output other;
  keep junction_id;
run;

proc sort data=insertion nodup;
  by junction_id;
proc sort data=deletion nodup;
  by junction_id;
proc sort data=other nodup;
  by junction_id;
run;

data indel_check;
  merge insertion (in=in1) deletion (in=in2) other (in=in3);
  by junction_id;
  if in1 then flag_insertion=1; else flag_insertion=0;
  if in2 then flag_deletion=1; else flag_deletion=0;
  if in3 then flag_no_indel=1; else flag_no_indel=0;
run;

proc freq data=indel_check noprint;
  tables flag_insertion*flag_deletion*flag_no_indel / out=indel_check2;
proc print data=indel_check2;
run;

/*
   flag_        flag_     flag_no_
 insertion    deletion      indel     COUNT    PERCENT

     0            0           1         71     93.4211
     0            1           1          4      5.2632
     1            0           1          1      1.3158


2 junctions that have reads with indels. We don't allow for this with bowtie, so we won't see these
*/


/* For remaining junction, get coordinates of donor and acceptor positions */

data monojunc_reads;
  set monojunc_indel;
  where flag_insertion=0 and flag_deletion=0;
  donor_size=scan(cigar_string,1,"M");
  intron_size=scan(scan(cigar_string,1,"N"),2,"M");
  donor_stop=read_start+donor_size-1;
  acceptor_start=read_start+donor_size+intron_size-1;
  keep junction_id read_id cigar_string read_start chr donor_stop acceptor_start;
run;

data mm10_donors;
  set mm10.mm10_exons_w_info;
  keep chrom start stop;
  rename chrom=chr start=donor_start stop=donor_stop;
run;

data mm10_acceptors;
  set mm10.mm10_exons_w_info;
  keep chrom start stop;
  rename chrom=chr start=acceptor_start stop=acceptor_stop;
run;

/* Take the maximal coordinates */

proc sort data=mm10_donors;
   by chr donor_stop donor_start;
proc means data=mm10_donors noprint;
   by chr donor_stop;
   var donor_start;
   output out=mm10_donors_max(drop=_TYPE_ _FREQ_) min=;
run;

proc sort data=mm10_acceptors;
   by chr acceptor_start acceptor_stop;
proc means data=mm10_acceptors noprint;
   by chr acceptor_start;
   var acceptor_stop;
   output out=mm10_acceptors_max(drop=_TYPE_ _FREQ_) max=;
run;

proc sort data=monojunc_reads;
   by chr donor_stop;
run;

data monojunc_donor nodonor;
  merge monojunc_reads (in=in1) mm10_donors_max (in=in2);
  by chr donor_stop;
  if in1 and in2 then output monojunc_donor;
  else if in1 then output nodonor; *0 obs;
run;

proc sort data=monojunc_donor;
   by chr acceptor_start;
run;

data monojunc_acceptor noacceptor;
  merge monojunc_donor (in=in1) mm10_acceptors_max (in=in2);
   by chr acceptor_start;
  if in1 and in2 then output monojunc_acceptor;
  else if in1 then output noacceptor; *0 obs;
run;

/* Calc donor and acceptor exon sizes and flag if shorter than the alignment */

data monojunc_exon_size;
  set monojunc_acceptor;
  donor_exon_len=donor_stop-donor_start;
  acceptor_exon_len=acceptor_stop-acceptor_start;
  donor_aln_len=scan(cigar_string,1,"M") + 0;
  acceptor_aln_len=scan(scan(cigar_string,2,"N"),1,"M") + 0;
  if donor_aln_len > donor_exon_len then flag_donor_aln_gt_exon=1;
  else flag_donor_aln_gt_exon=0;
  if acceptor_aln_len > acceptor_exon_len then flag_acceptor_aln_gt_exon=1;
  else flag_acceptor_aln_gt_exon=0;
run;

proc freq data=monojunc_exon_size;
  tables flag_donor_aln_gt_exon * flag_acceptor_aln_gt_exon;
run;

/*
 flag_donor_aln_gt_exon
           flag_acceptor_aln_gt_exon

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |    339 |     55 |    394
          |  78.84 |  12.79 |  91.63
          |  86.04 |  13.96 |
          |  90.40 | 100.00 |
 ---------+--------+--------+
        1 |     36 |      0 |     36
          |   8.37 |   0.00 |   8.37
          | 100.00 |   0.00 |
          |   9.60 |   0.00 |
 ---------+--------+--------+
 Total         375       55      430
             87.21    12.79   100.00

*/



data pb_junc;
  set event.event2star2pacbio_junc_table;
  where flag_in_pacbio=1;
  keep junction_id;
run;

proc sort data=pb_junc nodup;
  by junction_id;
proc sort data=monojunc_exon_size;
  by junction_id;
proc sort data=indel_check;
  by junction_id;
proc sort data=flag_microexon_match;
  by junction_id;
run;

data pb_monoj_ex_size;
  merge pb_junc (in=in1) monojunc_exon_size (in=in2);
  by junction_id;
  if in1 and in2;
run;

data pb_indel;
  merge pb_junc (in=in1) indel_check (in=in2);
  by junction_id;
  if in1 and in2;
run;

data pb_multijunc;
  merge pb_junc (in=in1) flag_microexon_match (in=in2);
  by junction_id;
  if in1 and in2;
run;

proc freq data=pb_monoj_ex_size;
  tables flag_donor_aln_gt_exon * flag_acceptor_aln_gt_exon;
run;

data pb_donor_small pb_acceptor_small;
   set pb_monoj_ex_size;
   if flag_donor_aln_gt_exon=1 then output pb_donor_small;
   if flag_acceptor_aln_gt_exon=1 then output pb_acceptor_small;
  keep junction_id;
run;

data pb_microexon;
  set pb_multijunc;
  keep junction_id;
run;

proc sort data=pb_donor_small nodup;
   by junction_id;
proc sort data=pb_acceptor_small nodup;
   by junction_id;
proc sort data=pb_microexon nodup;
   by junction_id;
run;


data nonpb_monoj_ex_size;
  merge pb_junc (in=in1) monojunc_exon_size (in=in2);
  by junction_id;
  if in1 and in2 then delete;
  else if in2 then output;
run;

data nonpb_indel;
  merge pb_junc (in=in1) indel_check (in=in2);
  by junction_id;
  if in1 and in2 then delete;
  else if in2 then output;
run;

data nonpb_multijunc;
  merge pb_junc (in=in1) flag_microexon_match (in=in2);
  by junction_id;
  if in1 and in2 then delete;
  else if in2 then output;
run;

proc freq data=nonpb_monoj_ex_size;
  tables flag_donor_aln_gt_exon * flag_acceptor_aln_gt_exon;
run;



data nonpb_donor_small nonpb_acceptor_small;
   set nonpb_monoj_ex_size;
   if flag_donor_aln_gt_exon=1 then output nonpb_donor_small;
   if flag_acceptor_aln_gt_exon=1 then output nonpb_acceptor_small;
  keep junction_id;
run;

data nonpb_microexon;
  set nonpb_multijunc;
  keep junction_id;
run;


data nonpb_indel2;
  set nonpb_indel;
  where flag_insertion=1 or flag_deletion=1;
  keep junction_id;
run;


proc sort data=nonpb_donor_small nodup;
   by junction_id;
proc sort data=nonpb_acceptor_small nodup;
   by junction_id;
proc sort data=nonpb_microexon nodup;
   by junction_id;
proc sort data=nonpb_indel2 nodup;
   by junction_id;
run;



data nonpb_juncs_star_apn0;
  set nonpb_donor_small (in=in1) nonpb_acceptor_small (in=in2) nonpb_microexon (in=in3)
      nonpb_indel2 (in=in4);
  if in1 then flag_microdonor=1; else flag_microdonor=0;
  if in2 then flag_microacceptor=1; else flag_microacceptor=0;
  if in3 then flag_microexon=1; else flag_microexon=0;
  if in4 then flag_indel=1; else flag_indel=0;
run;

data nic annot;
 set event.event2star2pacbio_junc_table;
 if junction_type="Annotated" then output annot;
 if junction_type="NIC/Unannotated" then output nic;
 keep junction_id;
run;

proc sort data=nic nodup;
   by junction_id;
proc sort data=annot nodup;
   by junction_id;
proc sort data=nonpb_juncs_star_apn0;
   by junction_id;
run;

data nonpb_junc_star_type_apn0;
  merge nonpb_juncs_star_apn0 (in=in1) nic (in=in2) annot (in=in3);
  by junction_id;
  if in2 then flag_nic=1; else flag_nic=0;
  if in3 then flag_annot=1; else flag_annot=0;
  if in1 then output;
run;

proc freq data=nonpb_junc_star_type_apn0 noprint;
   tables flag_microdonor*flag_microacceptor*flag_microexon*flag_indel*flag_annot*flag_nic/out=junc_check;
proc print data=junc_check;
run;

/*
  flag_          flag_          flag_      flag_    flag_
icrodonor    microacceptor    microexon    indel    annot    flag_nic    COUNT    PERCENT

    0              0              0          1        1          0          3      6.2500
    0              0              1          0        1          0         13     27.0833
    0              1              0          0        0          1          1      2.0833
    0              1              0          0        1          0         16     33.3333
    1              0              0          0        1          0         15     31.2500

*/

