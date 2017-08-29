libname con '!PATCON/sas_data';
libname hg19 '!PATCON/useful_human_data/aceview_hg19/fusions/sas_data';


/* Making a summary table for exons from T1D genes from Onengut-Gumuscu et al, 2015 */

data t1d_genes;
 set con.immunogene_flags;
 where flag_diabetes_gene=1;
 keep gene_id;
run;

data gene2fus;
  set hg19.hg19_aceview_fusions_si_info;
  keep gene_id fusion_id;
run;

proc sort data=t1d_genes;
  by gene_id;
proc sort data=gene2fus nodup;
  by gene_id fusion_id;
run;

data t1d_fus;
  merge t1d_genes (in=in1) gene2fus (in=in2);
  by gene_id;
  if in1 and in2;
  flag_fusion=1;
run;

data fus_on;
  set con.fusions_on_gt_apn0;
  drop flag_fusion_on0;
run;

proc sort data=t1d_fus;
  by fusion_id;
proc sort data=fus_on;
  by fusion_id;
run;

data t1d_fus_on;
  merge t1d_fus (in=in1) fus_on (in=in2);
  by fusion_id;
  if in1 and in2;
run;


/* Calc number of Sig Exons per contrast */

data results;
  set con.results_by_fusion_w_flags;
  keep fusion_id flag_cd4cd8_fdr05 flag_cd4cd19_fdr05 flag_cd8cd19_fdr05;
run;

proc sort data=results;
  by fusion_id;
proc sort data=t1d_fus_on;
  by fusion_id;
run;

data t1d_results;
  merge t1d_fus_on (in=in1) results;
  by fusion_id;
  if in1;
run;

data lsmeans;
  set con.lsmean_fusions;
  keep fusion_id cell_type estimate;
run;

proc sort data=lsmeans;
  by fusion_id cell_type;
proc transpose data=lsmeans out=lsmeans_sbys;
  by fusion_id;
  var estimate;
  id cell_type;
run;

proc sort data=lsmeans_sbys;
  by fusion_id;
proc sort data=t1d_results;
  by fusion_id;
run;

data t1d_all_results;
  merge t1d_results (in=in1) lsmeans_sbys;
  by fusion_id;
  if in1;
  drop _NAME_;
run;


/* Calc number of fusions, on fusions, sig fusions */

proc sort data=t1d_all_results;
  by gene_id;
proc means data=t1d_all_results noprint;
  by gene_id;
  var flag_fusion flag_cd19_on flag_cd8_on flag_cd4_on flag_fusion_all_on0
      flag_cd4cd8_Fdr05 flag_cd4cd19_Fdr05 flag_cd8cd19_Fdr05;
  output out=fusion_counts_by_gene
         sum(flag_fusion)=num_exons_total
         sum(flag_cd19_on)=num_exons_CD19_on
         sum(flag_cd8_on)=num_exons_CD8_on
         sum(flag_cd4_on)=num_exons_CD4_on
         sum(flag_fusion_all_on0)=num_exons_all_on
         sum(flag_cd4cd8_Fdr05)=num_exons_cd4cd8_fdr05
         sum(flag_cd4cd19_Fdr05)=num_exons_cd4cd19_fdr05
         sum(flag_cd8cd19_Fdr05)=num_exons_cd8cd19_fdr05;
run;

* calc magnitudes -- use LSmeans for this;

data mag;
  set t1d_all_results;
  if flag_fusion_all_on0=1 then do;
  if CD8 ge CD4 then magnitude_cd4cd8=CD8/CD4;
  else magnitude_cd4cd8=-CD4/CD8;

  if CD19 ge CD4 then magnitude_cd4cd19=CD19/CD4;
  else magnitude_cd4cd19=-CD4/CD19;

  if CD19 ge CD8 then magnitude_cd8cd19=CD19/CD8;
  else magnitude_cd8cd19=-CD8/CD19;

  end;
run;

/* Magnitude ranges by contrast */

proc sort data=mag;
  by gene_id;
proc means data=mag noprint;
  by gene_id;
  where flag_cd4cd8_Fdr05=1;
  var magnitude_cd4cd8;
  output out=cd4cd8_mag
         max(magnitude_cd4cd8)=max_cd4cd8_magnitude
         min(magnitude_cd4cd8)=min_cd4cd8_magnitude;
run;


proc means data=mag noprint;
  by gene_id;
  where flag_cd4cd19_Fdr05=1;
  var magnitude_cd4cd19;
  output out=cd4cd19_mag
         max(magnitude_cd4cd19)=max_cd4cd19_magnitude
         min(magnitude_cd4cd19)=min_cd4cd19_magnitude;
run;

proc means data=mag noprint;
  by gene_id;
  where flag_cd8cd19_Fdr05=1;
  var magnitude_cd8cd19;
  output out=cd8cd19_mag
         max(magnitude_cd8cd19)=max_cd8cd19_magnitude
         min(magnitude_cd8cd19)=min_cd8cd19_magnitude;
run;

proc sort data=cd4cd8_mag;
  by gene_id;
proc sort data=cd4cd19_mag;
  by gene_id;
proc sort data=cd8cd19_mag;
  by gene_id;
run;

data all_mags;
  merge cd4cd8_mag cd4cd19_mag cd8cd19_mag;
  by gene_id;
run;

data mag_ranges;
   set all_mags;
   format min_cd4cd8_str $10.;
   format max_cd4cd8_str $10.;
   format min_cd4cd19_str $10.;
   format max_cd4cd19_str $10.;
   format min_cd8cd19_str $10.;
   format max_cd8cd19_str $10.;
   format range_cd4cd8 $20.;
   format range_cd4cd19 $20.;
   format range_cd8cd19 $20.;

   min_cd4cd8_str=strip(put(min_cd4cd8_magnitude, 10.3));
   max_cd4cd8_str=strip(put(max_cd4cd8_magnitude, 10.3));

   min_cd4cd19_str=strip(put(min_cd4cd19_magnitude, 10.3));
   max_cd4cd19_str=strip(put(max_cd4cd19_magnitude, 10.3));

   min_cd8cd19_str=strip(put(min_cd8cd19_magnitude, 10.3));
   max_cd8cd19_str=strip(put(max_cd8cd19_magnitude, 10.3));

   if min_cd4cd8_magnitude=max_cd4cd8_magnitude then range_cd4cd8=max_cd4cd8_str;
   else range_cd4cd8 = cat(strip(min_cd4cd8_str)," - ",strip(max_cd4cd8_str));

   if min_cd4cd19_magnitude=max_cd4cd19_magnitude then range_cd4cd19=max_cd4cd19_str;
   else range_cd4cd19 = cat(strip(min_cd4cd19_str)," - ",strip(max_cd4cd19_str));

   if min_cd8cd19_magnitude=max_cd8cd19_magnitude then range_cd8cd19=max_cd8cd19_str;
   else range_cd8cd19 = cat(strip(min_cd8cd19_str)," - ",strip(max_cd8cd19_str));

   keep gene_id range_cd4cd8 range_cd4cd19 range_cd8cd19;
run;

proc sort data=mag_ranges;
  by gene_id;
proc sort data=fusion_counts_by_gene;
  by gene_id;
run;

data t1d_exon_summary;
  merge fusion_counts_by_gene (in=in1) mag_ranges;
  by gene_id;
  if in1;
  drop _TYPE_ _FREQ_;
run;

/* Make permenant */

data con.diabetes_gene_exon_summary;
  length chr $10.;
  set t1d_exon_summary;
  if gene_id="PTPN22" then do; chr="1p13.2"; order=1; end;
  if gene_id="IL10" then do; chr="1q32.1"; order=2; end;
  if gene_id="AFF3" then do; chr="2q11.2"; order=3; end;
  if gene_id="IFIH1" then do; chr="2q24.2"; order=4; end;
  if gene_id="CTLA4" then do; chr="2q33.2"; order=5; end;
  if gene_id="CCR5" then do; chr="3p21.31"; order=6; end;
  if gene_id="IL2" then do; chr="4q27"; order=7; end;
  if gene_id="IL21" then do; chr="4q27"; order=8; end;
  if gene_id="IL7R" then do; chr="5p13.2"; order=9; end;
  if gene_id="BACH2" then do; chr="6q15"; order=10; end;
  if gene_id="IKZF1" then do; chr="7p12.2"; order=11; end;
  if gene_id="GLIS3" then do; chr="9p24.2"; order=12; end;
  if gene_id="IL2RA" then do; chr="10p15.1"; order=13; end;
  if gene_id="INS-IGF2" then do; chr="11p15.5"; order=14; end;
  if gene_id="CD69" then do; chr="12p13.31"; order=15; end;
  if gene_id="IKZF4" then do; chr="12q13.2"; order=16; end;
  if gene_id="SH2B3" then do; chr="12q24.12"; order=17; end;
  if gene_id="GPR183" then do; chr="13q32.3"; order=18; end;
  if gene_id="RASGRP1" then do; chr="15q14"; order=19; end;
  if gene_id="CTSH" then do; chr="15q25.1"; order=20; end;
  if gene_id="IL27" then do; chr="16p11.2"; order=21; end;
  if gene_id="DEXI" then do; chr="16p13.13"; order=22; end;
  if gene_id="BCAR1" then do; chr="16q23.1"; order=23; end;
  if gene_id="IKZF3" then do; chr="17q12"; order=24; end;
  if gene_id="ORMDL3" then do; chr="17q12"; order=25; end;
  if gene_id="GSDMB" then do; chr="17q12"; order=26; end;
  if gene_id="CCR7" then do; chr="17q21.2"; order=27; end;
  if gene_id="PTPN2" then do; chr="18p11.21"; order=28; end;
  if gene_id="CD226" then do; chr="18q22.2"; order=29; end;
  if gene_id="TYK2" then do; chr="19p13.2"; order=30; end;
  if gene_id="FUT2" then do; chr="19q13.33"; order=31; end;
  if gene_id="UBASH3A" then do; chr="21q22.3"; order=32; end;
  if gene_id="ICOSLG" then do; chr="21q22.3"; order=33; end;
  if gene_id="C1QTNF6andIL2RB" then do; chr="22q12.3"; order=34; end;
  if gene_id="RAC2" then do; chr="22q12.3"; order=35; end;
  if range_cd4cd8="." then range_cd4cd8="n.d.";
  if range_cd4cd19="." then range_cd4cd19="n.d.";
  if range_cd8cd19="." then range_cd8cd19="n.d.";
run;

proc sort data=con.diabetes_gene_exon_summary;
  by order;
run;

* Export as CSV;

data data_for_export;
  set con.diabetes_gene_exon_summary;
  if gene_id="C1QTNF6andIL2RB" then gene_id="C1QTNF6/IL2RB";
  if num_exons_all_on=0 then delete;
  if num_exons_cd4cd8_fdr05=0 and num_exons_cd8cd19_fdr05=0 and num_exons_cd8cd19_fdr05=0 then delete;
  drop order;
run;

proc export data=data_for_export outfile='!MCLAB/jrbnewman/manuscripts/Newman_T1D_splicing/reviewer_responses/t1d_genes_exon_summary.csv' dbms=csv replace;
run;


