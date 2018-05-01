ods listing; ods html close;

libname con '!PATCON/sas_data';
libname eventloc '!MCLAB/event_analysis/sas_data';
libname event '/mnt/store/event_sandbox/sas_data';
libname hg19 '!PATCON/useful_human_data/aceview_hg19/fusions/sas_data';

/* Formatting fragment coverage/detection data for making heatmaps that show fragment detection, APN, and the proportion of isoforms that the fragment is derivable from

Will have the following columns:
1	gene_id
2	fragment_id
3	num transcripts
4	% transcripts
5	detection CD4
6	detection CD8
7	detection CD19
8	FAM129C_fragments.csvAPN CD4
9	APN CD8
10	APN CD19
11	scaled APN CD4
12	scaled APN CD8
13	scaled APN CD19
14-N	transcripts of gene, flag 0/1 is fragment is assigned to transcript

Test this on FAM129C: do for all transcripts, then for reduced set

I am going to scale the APN based on the highest fragment APN per cell type
*/

%macro fragPlotData(geneID);

data frag2gene;
  set hg19.hg19_aceview_exon_fragment_info;
  length fragment_coord $30.;
  where gene_id ? "&geneID.";
  fragment_coord=catx(":",chr,fragment_start,fragment_end);
  keep fragment_id fragment_start fragment_end fragment_coord gene_id transcript_id num_xscripts_per_fragment;
run;

data iso2frag;
  length transcript_id2 $60.;
  set frag2gene;
  do i=1 by 1 while(scan(transcript_id,i,"|") ^= "");
    transcript_id2=scan(transcript_id,i,"|");
    flag_frag_in_iso=1;
    output;
    end;
  keep fragment_id transcript_id2 flag_frag_in_iso gene_id;
  rename transcript_id2=transcript_id;
run;

proc sort data=iso2frag;
  by fragment_id transcript_id ;
proc transpose data=iso2frag out=iso2frag_sbys;
  by fragment_id;
  id transcript_id;
  var flag_frag_in_iso;
run;

data iso2frag_sbys2;
  set iso2frag_sbys;
  array change _numeric_;
    do over change;
    if change=. then change=0;
    end;
  drop _NAME_;
run;

/* Get fragment detection flags */

data frag_dtct;
  set event.hg19_flag_fragment_on_apn0;
  /* Set "missing"/ambiguous to -1 so I can color it grey on plots later */
  if flag_cd4_on=. then flag_cd4_on=0.5;
  if flag_cd8_on=. then flag_cd8_on=0.5;
  if flag_cd19_on=. then flag_cd19_on=0.5;
  keep fragment_id flag_cd4_on flag_cd8_on flag_cd19_on;
  rename fragment_id=fragment_coord;
run;

/* Get mean APN for each fragment */

data frag_apn;
   set event.t1d_case_only_fragment_counts;
   /* drop low coverage samples */
   if Name="2009-PC-0144" then delete;
   if Name="2009-PC-0221" then delete;
   if Name="2009-PC-0235" then delete;
   if Name="2009-PC-0236" then delete;
   if Name="2009-PC-0237" then delete;
   log_apn=log(apn+1);
   keep name cell_type fragment_id log_apn;
   rename fragment_id=fragment_coord;
run;

proc sort data=frag_apn;
  by fragment_coord cell_type;
proc means data=frag_apn noprint;
  by fragment_coord cell_type;
  var log_apn;
  output out=mean_apn_by_cell mean=;
run;

proc transpose data=mean_apn_by_cell out=mean_apn_sbys;
   by fragment_coord;
   id cell_type;
   var log_apn;
run;

/* Calculate the % of transcripts per fragment */

data xs2gene;
  set iso2frag;
  keep gene_id transcript_id;
run;

proc sort data=xs2gene nodup;
  by gene_id transcript_id;
proc freq data=xs2gene noprint;
  tables gene_id / out=xs_per_gene;
run;

%local xsCnt;
data _null_;
  set xs_per_gene nobs=n;
  call symputx('xsCnt',count);
  stop;
run;

%put &xsCnt.;

data perc_iso;
   set frag2gene;
   frag_perc_iso=num_xscripts_per_fragment/&xsCnt.;
   keep fragment_id frag_perc_iso;
run;

/* Merge all together */

data mean_apn_sbys2;
   set mean_apn_sbys;
   rename cd4=cd4_log_apn cd8=cd8_log_apn cd19=cd19_log_apn;
   drop _NAME_;
run;
   
proc sort data=mean_apn_sbys2;
   by fragment_coord;
proc sort data=frag_dtct;
   by fragment_coord;
proc sort data=frag2gene;
   by fragment_coord;
run;

data frag_dtct_apn;
  merge frag2gene (in=in1) mean_apn_sbys2 frag_dtct ;
  by fragment_coord;
  if in1;
run;

proc sort data=frag_dtct_apn;
   by fragment_id;
proc sort data=iso2frag_sbys2;
   by fragment_id;
proc sort data=perc_iso;
   by fragment_id;
run;

data frag_all_info;
  merge frag_dtct_apn (in=in1) iso2frag_sbys2 perc_iso;
  by fragment_id;
  if in1;
run;

/* Calculate perc of max APN */

proc means data=frag_all_info noprint;
  var cd19_log_apn cd4_log_apn cd8_log_apn;
  output out=max_apn max=;
run;

data _null_;
   set max_apn nobs=n;
   call symputx('maxCD4',cd4_log_apn);
   call symputx('maxCD8',cd8_log_apn);
   call symputx('maxCD19',cd19_log_apn);
   stop;
run;

/* Reorder columns */

data frag_all_info2;
  retain gene_id fragment_id fragment_start fragment_end
         frag_perc_iso flag_cd4_on flag_cd8_on flag_cd19_on
         cd4_log_apn cd8_log_apn cd19_log_apn
         cd4_log_apn_scaled cd8_log_apn_scaled cd19_log_apn_scaled;
  set frag_all_info;
  cd4_log_apn_scaled=cd4_log_apn/&maxCD4.;
  cd8_log_apn_scaled=cd8_log_apn/&maxCD8.;
  cd19_log_apn_scaled=cd19_log_apn/&maxCD19.;
  drop transcript_id fragment_coord;
run;

proc sort data=frag_all_info2;
  by fragment_start fragment_end;
run;

proc export data=frag_all_info2 outfile="!MCLAB/event_analysis/analysis_output/fragment_isoform_heatmaps/&geneID._fragments.csv" dbms=csv replace;
run;
%mend;

%fragPlotData(UBASH3A);

%fragPlotData(FAM129C);
%fragPlotData(IKZF3);
%fragPlotData(CD3DandMPZL2);
%fragPlotData(FCRL2);
%fragPlotData(MAGandCD22);
%fragPlotData(PDLIM1);
%fragPlotData(SPIB);




