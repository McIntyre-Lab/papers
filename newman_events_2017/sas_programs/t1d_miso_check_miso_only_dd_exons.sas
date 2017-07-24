
ods listing; ods html close;
libname event '!MCLAB/event_analysis/sas_data';

/* Of the genes with with DD MISO events ONLY: why are these only seen in MISO? */


data miso_genes;
   set event.t1d_miso2qds_comparison2;
   where flag_miso_gene_event_dd=1 and flag_cell_by_fus_fdr05=0 and 
         flag_gene_exon_dd=0 and flag_sig_bf5=0;
   keep gene_id;
run;

data ens2av;
  set event.hg19_ens2refseq2av;
  keep ens_gene_id gene_id;
run;

proc sort data=miso_genes;
   by gene_id;
proc sort data=ens2av;
   by gene_id;
run;

data miso_ens2av;
  merge ens2av (in=in1) miso_genes (in=in2);
  by gene_id;
  if in1 and in2;
run;
 
/* import miso gene-to-event annotation */

proc import datafile="/mnt/store/miso_sandbox/hg19/hg19/SE.hg19.gff3_to_ensGene2.txt" out=miso_gene2event
   dbms=tab replace; guessingrows=40000;
run;

data miso_gene2event2;
   set miso_gene2event;
   rename gene_id=ens_gene_id;
run;

data miso_ens2av2;
  set miso_ens2av;
  keep ens_gene_id;
run;

proc sort data=miso_ens2av2 nodup;
   by ens_gene_id;
proc sort data=miso_gene2event2;
  by ens_gene_id;
run;

data miso_events;
  merge miso_ens2av2 (in=in1) miso_gene2event2 (in=in2);
   by ens_gene_id;
  if in1 and in2;
run;

data dd_events;
   set event.t1d_all_miso_results_v2;
   where sum(flag_cd4_on,flag_cd8_on,flag_cd19_on) =1 or sum(flag_cd4_on,flag_cd8_on,flag_cd19_on) =2;
   keep event_name;
   rename event_name=event_id;
run;

proc sort data=dd_events;
  by event_id;
proc sort data=miso_events;
  by event_id;
run;

data miso_dd_events;
  merge miso_events (in=in1) dd_events (in=in2);
  by event_id;
  if in1 and in2;
run;

/* Import MISO SE GFF3 */


    data WORK.MISO_GFF    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile '/mnt/store/miso_sandbox/hg19/SE_2.hg19.gff3' delimiter='09'x MISSOVER DSD
lrecl=32767 firstobs=2 ;
       informat chr $2. ;
       informat event_type $2. ;
       informat feature_type $4. ;
       informat start best32. ;
       informat stop best32. ;
       informat score best32. ;
       informat strand $1. ;
       informat frame best32. ;
       informat attrib $245. ;
       format chr $2. ;
       format event_type $2. ;
       format feature_type $4. ;
       format start best12. ;
       format stop best12. ;
       format score best12. ;
       format strand $1. ;
       format frame best12. ;
       format attrib $245. ;
    input
                chr $
                event_type $
                feature_type $
                start
                stop
                score
                strand $
                frame
                attrib $
    ;
    if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
    run;




data miso_gff2;
  set miso_gff;
  where feature_type="exon";
  length event_id $245.;
  event_id=strip(
           tranwrd(
           tranwrd(
           tranwrd(
           tranwrd(
           tranwrd(
           tranwrd(
           tranwrd(scan(attrib,1,";"),"ID=",""),
           ".A.dn",""),".A.se",""),".A.up",""),".B.dn",""),".B.se",""),".B.up","")
           );
  keep event_id chr start stop;
run;

proc sort data=miso_gff2;
   by event_id;
proc sort data=miso_dd_events;
   by event_id;
run;

data miso_dd_events2 noexon;
  merge miso_dd_events (in=in1) miso_gff2 (in=in2);
  by event_id;
  if in1 and in2 then output miso_dd_events2;
  else if in1 then output noexon;
run;



/* Import aceview exon annotations */

proc import datafile="!MCLAB/junction_annotations/pipeline_output/aceview_hg19/hg19_aceview_exons.bed"
   out=av_exons dbms=tab replace; getnames=no;
run;

    data WORK.AV_EXONS    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile '!MCLAB/junction_annotations/pipeline_output/aceview_hg19/hg19_aceview_exons.bed'
delimiter='09'x MISSOVER DSD lrecl=32767 ;
       informat chr $4. ;
       informat start best32. ;
       informat stop best32. ;
       informat exon_id $39. ;
       informat score best32. ;
       informat strand $1. ;
       format chr $4. ;
       format start best12. ;
       format stop best12. ;
       format exon_id $39. ;
       format score best12. ;
       format strand $1. ;
    input
                chr $
                start
                stop
                exon_id $
                score
                strand $
    ;
    if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
    run;



data av_exons2;
  set av_exons;
  stop2=stop+1;
  keep chr start stop2 exon_id;
  rename stop2=stop;
run;

data av_fus;
  set hg19.hg19_aceview_fusions_si_bed;
  start=fusion_start+1;
  keep fusion_id chr start fusion_stop;
  rename fusion_stop=stop;
run;


/* Match MISO exons to Aceview Exons */



proc sort data=miso_dd_events2;
   by chr start;
proc sort data=av_fus;
    by chr start;
run;

data miso2av_match1 miso_no av_no;
  merge miso_dd_events2 (in=in1) av_fus (in=in2);
   by chr start;
  if in1 and in2 then output miso2av_match1;
   else if in1 then output miso_no;
  else output av_no;
run;

data av_no2;
   set av_no;
   drop start ens_gene_id event_id;
run;


proc sort data=miso_no;
   by chr stop;
proc sort data=av_no2;
    by chr stop;
run;

data miso2av_match2 miso_none;
  merge miso_no (in=in1) av_no2 (in=in2);
   by chr stop;
   if in1 and in2 then output miso2av_match2;
   else if in1 then output miso_none;
run;

data miso_none2;
  set miso_none;
  keep ens_gene_id;
run;

proc sort data=miso_none2 nodup;
  by ens_gene_id;
run;

data miso2av_match3;
   set miso2av_match1 miso2av_match2;
run;


data fus_flags;
  set con.fusion_on_ge_apn5_v2;
  if flag_cd4_on=. or flag_cd8_on=. or flag_cd19_on=. then flag_fusion_exc=1; else flag_fusion_exc=0;
  keep fusion_id flag_fusion_exc flag_fusion_on5 flag_fusion_all_on5;
run;

proc sort data=fus_flags;
  by fusion_id;
proc sort data=miso2av_match3 nodup;
  by fusion_id;
run;

data miso2av_w_flag;
  merge miso2av_match3 (in=in1) fus_flags (in=in2);
  by fusion_id;
  if in1 and in2;
run;

proc freq data=miso2av_w_flag noprint;
   tables flag_fusion_on5*flag_fusion_all_on5*flag_fusion_exc / out=fus_check;
proc print data=fus_check;
run;

data miso_exc miso_off miso_on;
   set miso2av_w_flag;
   if flag_fusion_exc=1 then output miso_exc;
   if flag_fusion_all_on5=1 then output miso_on;
   else if flag_fusion_all_on5=0 then output miso_off;
   keep ens_gene_id;
run;

proc sort data=miso_exc nodup;
  by ens_gene_id;
proc sort data=miso_off nodup;
  by ens_gene_id;
proc sort data=miso_on nodup;
  by ens_gene_id;
run;


   

proc sort data=av_fus_w_flag;
   by gene_id;
proc sort data=miso_genes nodup;
  by gene_id;
run;

data miso_genes_fus_check;
  merge miso_genes (in=in1) av_fus_w_flag (in=in2);
  by gene_id;
  if in1 and in2;
run;

proc export data=miso_genes_fus_check outfile="!MCLAB/event_analysis/hg19_fusion_to_check.csv"
    dbms=csv replace;
run;



