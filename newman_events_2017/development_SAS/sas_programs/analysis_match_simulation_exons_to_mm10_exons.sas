/* Libraries */
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
ods listing; ods html close;

/* For exons from simulated data, match to mm10 exons */

/* Collapse exons with the same coordinates */
data mm10_exons;
  length chrom $12.;
  set mm10.mm10_exons_w_info;
  length coord $30.;
  length transcript_id2 $25.;
  coord=catx(":",chrom,start,stop,strand);
  do i=1 by 1 while(scan(transcript_id,i,"|") ^= "");
      transcript_id2=scan(transcript_id,i,"|");
      output; end;
  keep chrom start stop strand transcript_id2 coord;
  rename transcript_id2=transcript_id;
run;

proc sort data=mm10_exons nodup;
  by chrom start stop strand coord transcript_id;
run;

proc freq data=mm10_exons noprint;
  by chrom start stop strand;
  tables coord / out=coord_count;
proc sort data=coord_count;
  by descending count;
run; *130 max;

data xs_cat;                 
  array xs[130] $ 25. ;    
  retain xs1-xs130;
  set mm10_exons;
  by chrom start stop strand coord;  
  if first.coord then do;
    call missing(of xs1-xs130);
    records = 0;
  end;
  records + 1;
  xs[records]=transcript_id;       
  if last.coord then output;
run;   

data xs_cat2;
  set xs_cat;
  length transcript_cat $ 2000;
  transcript_cat = catx("|",OF xs1-xs130);
  keep chrom start stop strand coord transcript_cat;
  run;

/* merge simulated exons to mm10 exons */

data sim_exons;
  set event.simulated_exons_mm9_to_mm10;
  keep simulation test sim_exon_id chr exon_start exon_stop strand flag_mm10_coord;
  rename chr=chrom exon_start=start exon_stop=stop;
run;

proc sort data=xs_cat2;
   by chrom start stop strand;
proc sort data=sim_exons;
   by chrom start stop strand;
run;

data sim2mm10_exons no_sim no_mm10;
  merge xs_cat2 (in=in1) sim_exons (in=in2);
  by chrom start stop strand;
  if in1 and in2 then output sim2mm10_exons;
  else if in1 then output no_sim;
  else output no_mm10; *should only have exons without mm10 coordinates, and exons on contigs here;
run;

proc freq data=no_mm10;
  tables flag_mm10_coord;
run;
/*

                                               Cumulative    Cumulative
   flag_mm10_coord    Frequency     Percent     Frequency      Percent
   --------------------------------------------------------------------
                 0        2385        0.18          2385         0.18
                 1     1300148       99.82       1302533       100.00

This is a lot. Check what chromosomes they are on
*/

proc freq data=no_mm10;
  where flag_mm10_coord=1;
  tables chrom;
run;

/*
                                   Cumulative    Cumulative
 chrom    Frequency     Percent     Frequency      Percent
 ----------------------------------------------------------
 chr1        77632        5.97         77632         5.97
 chr10       59632        4.59        137264        10.56
 chr11       91890        7.07        229154        17.63
 chr12       52815        4.06        281969        21.69
 chr13       53017        4.08        334986        25.77
 chr14       57542        4.43        392528        30.19
 chr15       47429        3.65        439957        33.84
 chr16       41290        3.18        481247        37.01
 chr17       58350        4.49        539597        41.50
 chr18       35626        2.74        575223        44.24
 chr19       37605        2.89        612828        47.14
 chr1_         619        0.05        613447        47.18
 chr2       107705        8.28        721152        55.47
 chr3        64615        4.97        785767        60.44
 chr4        85426        6.57        871193        67.01
 chr4_          44        0.00        871237        67.01
 chr5        77405        5.95        948642        72.96
 chr5_          99        0.01        948741        72.97
 chr6        68804        5.29       1017545        78.26
 chr7        93091        7.16       1110636        85.42
 chr7_          33        0.00       1110669        85.43
 chr8        63042        4.85       1173711        90.28
 chr9        66061        5.08       1239772        95.36
 chrM          119        0.01       1239891        95.37
 chrUn          10        0.00       1239901        95.37
 chrX        51699        3.98       1291600        99.34
 chrX_         139        0.01       1291739        99.35
 chrY         8409        0.65       1300148       100.00

Could be that these are from non-Refseq exons.. skip these
*/


/* For each simulated transcript, count the number of exons */
data sim_exons2;
  set sim_exons;
  length sim_gene_id $15.;
  sim_gene_id=scan(sim_exon_id,1,":");
run;

proc sort data=sim_exons2;
  by simulation test sim_gene_id;
proc freq data=sim_exons2 noprint;
  by simulation test;
  tables sim_gene_id / out=exons_per_sim_gene;
run;

data sim2mm10_exons2;
  set sim2mm10_exons;
  length sim_gene_id $15.;
  sim_gene_id=scan(sim_exon_id,1,":");
run; 

proc sort data=sim2mm10_exons2;
   by simulation test sim_gene_id;
proc sort data=exons_per_sim_gene;
   by simulation test sim_gene_id;
run;

data sim2mm10_exons3;
  merge sim2mm10_exons2 (in=in1) exons_per_sim_gene (in=in2);
  by simulation test sim_gene_id;
  if in1 and in2;
  drop PERCENT;
  rename count=num_exons_per_sim_gene;
run;

/* For each Refseq transcript, count the number of exons */
proc sort data=mm10_exons;
  by transcript_id;
proc freq data=mm10_exons noprint;
  tables transcript_id / out=exons_per_xs_mm10;
run;

/* Merge and make permenant for now */
data sim2mm10_exons4;
  set sim2mm10_exons3;
  length transcript_id $25.;
  do i=1 by 1 while(scan(transcript_cat,i,"|") ^= "" );
      transcript_id=scan(transcript_cat,i,"|");
      output;
      end;
  drop i transcript_cat;
run;

proc sort data=sim2mm10_exons4;
   by transcript_id;
proc sort data=exons_per_xs_mm10;
   by transcript_id;
run;

data sim2mm10_exons5;
  merge exons_per_xs_mm10 (in=in1) sim2mm10_exons4 (in=in2);
  by transcript_id;
  if in1 and in2;
  drop PERCENT;
  rename count=num_exons_per_xscript;
run;

data eventloc.simulated_exons_to_mm10;
   set sim2mm10_exons5;
run;


