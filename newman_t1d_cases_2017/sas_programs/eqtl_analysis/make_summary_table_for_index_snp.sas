/* Make a summary table of T1D-associated SNPs:
   chrom
   position
   index_snp
   cell type
   number of exon cis-eqtl (magnitude in parenthesis)
   number of splicing cis-eqtl (magnitude in parenthesis)
   gene list
*/

libname eqtl '!PATCON/eqtl_analysis/sas_data';
libname con '!PATCON/sas_data';

/* Get index SNPs */

data results;
  set eqtl.eqtl_results_w_onengut_index2;
  where onengut_snp2eqtl_snp_r2 ge 0.8 and flag_eqtl_sig=1;
run;

data ratios;
  set eqtl.results_summary_table_w_means_v3;
  where flag_eqtl_sig=1;
  keep snp_id gene_id feature_id
  magnitude_htz_hmz1_cd4 magnitude_hmz2_hmz1_cd4
  magnitude_htz_hmz1_cd8 magnitude_hmz2_hmz1_cd8
  magnitude_htz_hmz1_cd19 magnitude_hmz2_hmz1_cd19
  ;
run;

data estimate;
  set ratios;
  if magnitude_hmz2_hmz1_cd4=. then cd4_mag=magnitude_htz_hmz1_cd4;
  else cd4_mag=magnitude_hmz2_hmz1_cd4;

  if magnitude_hmz2_hmz1_cd8=. then cd8_mag=magnitude_htz_hmz1_cd8;
  else cd8_mag=magnitude_hmz2_hmz1_cd8;

  if magnitude_hmz2_hmz1_cd19=. then cd19_mag=magnitude_htz_hmz1_cd19;
  else cd19_mag=magnitude_hmz2_hmz1_cd19;

  keep snp_id gene_id feature_id cd4_mag cd8_mag cd19_mag;
run;

proc sort data=estimate;
   by gene_id snp_id feature_id;
proc sort data=results;
   by gene_id snp_id feature_id;
run;

data results_w_mag;
  merge results (in=in1) estimate (in=in2);
  by gene_id snp_id feature_id;
  if in1 and in2;
run;

data cd4;
  length cell_type $5.;
  length feat_type $8.;
  set results_w_mag;
  if CD4_FDR_P < 0.05;
  cell_type="CD4+";
  flag_sig_eqtl=1;
  if CD4_FDR_P = . then delete;
  if cd4_mag < 1 then magnitude=-(1/cd4_mag);
  else magnitude = cd4_mag;
  if feature_type="exon" then feat_type="exon";
  else feat_type="splicing";
  drop CD8_FDR_P cd8_mag CD19_FDR_P cd19_mag  cd4_mag;
  rename CD4_FDR_P=FDR_P;
run;

data cd8;
  length cell_type $5.;
  length feat_type $8.;
  set results_w_mag;
  if CD8_FDR_P < 0.05;
  cell_type="CD8+";
  flag_sig_eqtl=1;
  if CD8_FDR_P = . then delete;
  if cd8_mag < 1 then magnitude=-(1/cd8_mag);
  else magnitude = cd8_mag;
  if feature_type="exon" then feat_type="exon";
  else feat_type="splicing";
  drop CD4_FDR_P cd4_mag CD19_FDR_P cd19_mag cd8_mag;
  rename CD8_FDR_P=FDR_P;
run;


data cd19;
  length cell_type $5.;
  length feat_type $8.;
  set results_w_mag;
  if CD19_FDR_P < 0.05;
  cell_type="CD19+";
  flag_sig_eqtl=1;
  if CD19_FDR_P = . then delete;
  if cd19_mag < 1 then magnitude=-(1/cd19_mag);
  else magnitude = cd19_mag;
  if feature_type="exon" then feat_type="exon";
  else feat_type="splicing";
  drop CD8_FDR_P cd8_mag CD4_FDR_P cd4_mag cd19_mag; 
  rename CD19_FDR_P=FDR_P;
run;


data stack;
  set cd4 cd8 cd19;
run;

proc sort data=stack;
   by onengut_index_snp cell_type feat_type;
proc freq data=stack noprint;
   by onengut_index_snp cell_type feat_type;
   tables flag_sig_eqtl / out=num_sig_by_cell_feat_snp;
run;

/* For exons, cat genes */

data gene2snp;
  length gene_id2 $35.;
  set stack;
  do i=1 by 1 while (scan(gene_id,i,"|") ^= ' ');
     gene_id2=scan(gene_id,i,"|");
     output;
     end;
  keep onengut_index_snp cell_type feat_type gene_id2;
  rename gene_id2=gene_id;
run;

proc sort data=gene2snp nodup;
  by onengut_index_snp cell_type feat_type gene_id;
proc freq data=gene2snp noprint;
  tables onengut_index_snp*cell_type*feat_type / out=gene_count;
proc sort data=gene_count;
  by descending count;
run; *max 5 genes;

data cat_genes; 
  array gene[5] $ 35.;
  retain gene1-gene5;
  set gene2snp;
  by onengut_index_snp cell_type feat_type;
  if first.feat_type then do;
     call missing(of gene1-gene5);
     records = 0;
  end;
  records + 1;
  gene[records]=gene_id;
  if last.feat_type then output;
run;

  *clean up the output file;
data cat_genes2;
  set cat_genes;
  length gene_list $ 200.;
  rename records= num_genes;
         gene_list= catx(", ", OF gene1-gene5);
  drop gene1-gene5 gene_id;
  run;

/* Get magnitude ranges */

proc sort data=stack;
  by onengut_index_snp cell_type feat_type;
proc means data=stack;
  by onengut_index_snp cell_type feat_type;
  var magnitude;
  output out=mag_min_max min=min max=max;
run;


proc sort data=cat_genes2;
   by onengut_index_snp cell_type feat_type;
proc sort data=num_sig_by_cell_feat_snp;
   by onengut_index_snp cell_type feat_type;
proc sort data=mag_min_max;
   by onengut_index_snp cell_type feat_type;
run;

data summary_by_snp;
  merge num_sig_by_cell_feat_snp (in=in1) cat_genes2 (in=in2) mag_min_max (in=in3);
  by onengut_index_snp cell_type feat_type;
  if in1 and in2 and in3;
run;

data add_mag_range;
  length sig_eqtl $255.;
  set summary_by_snp;
  if count = 1 then sig_eqtl=cat(strip(put(count,10.0))," (",strip(put(max,10.3)), ")");
  else sig_eqtl=cat(strip(put(count,10.0)), " (", strip(put(min,10.3)), " - ", strip(put(max,10.3)), ")");
  drop _TYPE_ _FREQ_ min max num_genes count percent flag_sig_eqtl;
run;

data exons;
  set add_mag_range;
  where feat_type="exon";
  rename sig_eqtl=sig_exon gene_list=gene_exon;
  drop feat_type;
run;

data splicing;
  set add_mag_range;
  where feat_type="splicing";
  rename sig_eqtl=sig_splice gene_list=gene_splice;
  drop feat_type;
run;

proc sort data=exons;
  by onengut_index_snp cell_type;
proc sort data=splicing;
  by onengut_index_snp cell_type;
run;

data summary;
  merge exons splicing;
  by onengut_index_snp cell_type;
run;

data snp_pos;
  informat chr $10.;
  informat position best32. ;
  informat onengut_index_snp $31.;
  input chr $ position onengut_index_snp $ ;
  datalines;
1p13.2 114377568 rs2476601
1q32.1 200814959 rs6691977
1q32.1 206939904 rs3024505
2q11.2 100764087 rs13415583
2q13 111615079 rs4849135
2q24.2 163110536 rs2111485
2q24.2 163124637 rs35667974
2q24.2 163136942 rs72871627
2q33.2 204738919 rs3087243
3p21.31 46457412 rs113010081
4q27 123243596 rs75793288
4q32.3 166574267 rs2611215
5p13.2 35883251 rs11954020
6q15 90976768 rs72928038
6q22.32 126752884 rs1538171
7p12.2 50465830 rs62447205
7p12.1 51028987 rs10277986
9p24.2 4290823 rs6476839
10p15.1 6094697 rs61839660
10p15.1 6108340 rs10795791
10p15.1 6129643 rs41295121
10q23.31 90035654 rs12416116
11p15.5 2182224 rs689
11p15.5 2198665 rs72853903
12p13.31 9905851 rs917911
12q13.2 56435504 rs705705
12q24.12 112007756 rs653178
13q32.3 100081766 rs9585056
14q32.2 98488007 rs1456988
14q32.2 101306447 rs56994090
15q14 38847022 rs72727394
15q25.1 79234957 rs34593439
16p11.2 28505660 rs151234
16p13.13 11194771 rs12927355
16p13.13 11351211 rs193778
16q23.1 75252327 rs8056814
17q12 38053207 rs12453507
17q21.2 38775150 rs757411
17q21.31 44073889 rs1052553
18p11.21 12809340 rs1893217
18p11.21 12830538 rs12971201
18q22.2 67526644 rs1615504
19p13.2 10463118 rs34536443
19p13.2 10469975 rs12720356
19q13.32 47219122 rs402072
19q13.33 49206172 rs516246
20p13 1616206 rs6043409
21q22.3 43825357 rs11203202
21q22.3 45621817 rs6518350
22q12.2 30531091 rs4820830
22q12.3 37587111 rs229533
;
run;

proc sort data=snp_pos;
  by onengut_index_snp;
proc sort data=summary;
  by onengut_index_snp;
run;

data summary_w_pos;
  merge snp_pos (in=in1) summary (in=in2);
  by onengut_index_snp;
  if in1 and in2;
run;

/* Export. I am going to manually finish the table */

proc export data=summary_w_pos outfile="!MCLAB/jrbnewman/manuscripts/Newman_T1D_splicing/reviewer_responses/eqtl_summary_by_index_snp_jrbn2.csv" 
  dbms=csv replace;
run;

