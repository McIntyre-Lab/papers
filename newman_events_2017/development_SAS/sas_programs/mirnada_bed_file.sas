proc import datafile="/mnt/store/isoAnnot/references/mirna_targets/fruitfly_predictions_S_C_aug2010.txt"
   out=mirnada_c dbms=tab replace;
   guessingrows=100000;
run;

proc import datafile="/mnt/store/isoAnnot/references/mirna_targets/fruitfly_predictions_S_0_aug2010.txt"
   out=mirnada_nc dbms=tab replace;
   guessingrows=100000;
run;

proc import datafile="/mnt/store/isoAnnot/references/dmel_fb611_ensembl_BDGP6/fbtr2refseq.csv"
   out=rs2fbtr dbms=csv replace;
   guessingrows=35000; getnames=no;
run;

data mirnada_c2;
  length chr $5.;
  length start_stop $50.;
  set mirnada_c;
  chr=scan(genome_coordinates,2,":");
  start_stop=scan(genome_coordinates,3,":");
  keep mirna_name ext_transcript_id chr start_stop;
run;

data mirnada_nc2;
  length chr $5.;
  length start_stop $50.;
  set mirnada_nc;
  chr=scan(genome_coordinates,2,":");
  start_stop=scan(genome_coordinates,3,":");
  keep mirna_name ext_transcript_id chr start_stop;
run;

data mirnada_c_nc;
   length start_stop2 $50.;
   set mirnada_c2 mirnada_nc2;
   do i=1 by 1 while(scan(start_stop,i,",") ^="");
         start_stop2=scan(start_stop,i,",");
         output;
         end;
run;

data mirnada_c_nc_2;
   set mirnada_c_nc;
   start=scan(start_stop2,1,"-");
   stop=scan(start_stop2,2,"-");
   keep mirna_name ext_transcript_id chr start stop;
run;

proc sort data=mirnada_c_nc_2 nodup;
   by mirna_name ext_transcript_id chr start stop;
run;

data rs2fbtr2;
  length ext_transcript_id $15.;
  set rs2fbtr;
  ext_transcript_id=scan(VAR2,1,".");
run;

proc sort data=rs2fbtr2;
  by ext_transcript_id;
proc sort data=mirnada_c_nc_2;
  by ext_transcript_id;
run;

data mirnada_c_nc_2_fb;
  merge rs2fbtr2 (in=in1) mirnada_c_nc_2 (in=in2);
  by ext_transcript_id;
  if in1 and in2;
run;


data export_data;
  retain chr start stop;
  length mirna_transcript_id $50.;
  set mirnada_c_nc_2_fb;
  mirna_transcript_id=catx('|',VAR1,mirna_name);
  keep chr start stop mirna_transcript_id;
run;

proc export data=export_data outfile="/mnt/store/isoAnnot/references/mirna_targets/mirnada_coord_dm3.bed"
   dbms=tab replace; putnames=no;
run;

