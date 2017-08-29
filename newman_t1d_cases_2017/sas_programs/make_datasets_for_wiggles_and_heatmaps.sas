libname con '!PATCON/sas_data';
libname eqtl '!PATCON/eqtl_analysis/sas_data';
libname av '!PATCON/useful_human_data/aceview_hg19/sas_data';
libname fus '!PATCON/useful_human_data/aceview_hg19/fusions/sas_data';
libname splicing '/mnt/data/splicing';
libname splice '/mnt/data/splice';

/* Data for wiggle plots and coverage heatmaps */

/* I want to make genotype-specific wiggle plots for the following genes * SNP * cell:
   UBASH3A rs1893592 CD4+
   PTPN22 rs1217414 in CD19+ */

/* Import subset mpileups and sum tech reps */

    data WORK.COUNTS_FOR_WIGS    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile '!PATCON/pipeline_output/counts_for_wiggles.csv' delimiter = ',' MISSOVER DSD
lrecl=32767 firstobs=2 ;
       informat sample_id $31. ;
       informat chrom $4. ;
       informat pos best32. ;
       informat count best32. ;
       informat gene $7. ;
       format sample_id $31. ;
       format chrom $4. ;
       format pos best12. ;
       format count best12. ;
       format gene $7. ;
    input
                sample_id $
                chrom $
                pos
                count
                gene $
    ;
    if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
    run;

/* Sum tech reps */

data design;
   set con.design_file;
   if name =  '2009-PC-0221' then delete; *sample 75 cd8;
   if name =  '2009-PC-0144' then delete; *sample 48 cd4;
   if name =  '2009-PC-0236' then delete; *sample 80;
   if name =  '2009-PC-0237' then delete; *sample 80;
   if name =  '2009-PC-0235' then delete; *sample 80;
*  if name = '2009-PC-0101' then delete; *sample 34cd8;
if name = '2009-PC-0104' then delete; *sample35 cd8;
if name =  '2009-PC-0153' then delete; *sample 51 cd4;
if name =  '2009-PC-0200' then delete; *sample 67 cd8;
if name =  '2009-PC-0212' then delete; *sample 72 cd8;
if name = '2009-PC-0215' then delete; *sample 73  cd8;
/*funky heatmap samples*/

if name =  '2009-PC-0083' then delete; 
if name =   '2009-PC-0114' then delete;
if name =   '2009-PC-0224' then delete;
if name =   '2009-PC-0228' then delete;

   keep sample_id subject_id cell_type name;
run;

proc sort data=design nodup;
  by sample_id;
proc sort data=counts_for_wigs;
  by sample_id;
run;

data counts_w_key;
  merge design (in=in1) counts_for_wigs (in=in2);
  by sample_id;
  if in1 and in2;
run;

proc sort data=counts_w_key;
  by subject_id cell_type name gene chrom pos;
proc means data=counts_w_key noprint;
  by subject_id cell_type name gene chrom pos;
  var count;
  output out=counts_summed_w_key sum=;
run;

/* Merge in genotypes */

data ptpn22_snp ubash3a_snp;
  set eqtl.snp_data_w_info;
  if snp_id="rs1893592" then output ubash3a_snp;
  if snp_id="rs1217414" then output ptpn22_snp;
  keep snp_id ref_allele alt_allele subject_id genotype;
run;

data ptpn22_snp2;
  set ptpn22_snp;
  length rs1217414 $3.;
  if genotype=0 then rs1217414=catt(ref_allele,ref_allele);
  else if genotype=1 then rs1217414=catt(ref_allele,alt_allele);
  else if genotype=2 then rs1217414=catt(alt_allele,alt_allele);
  else delete;
  drop genotype  ref_allele alt_allele;
run;

data ubash3a_snp2;
  set ubash3a_snp;
  length rs1893592 $3.;
  if genotype=0 then rs1893592=catt(ref_allele,ref_allele);
  else if genotype=1 then rs1893592=catt(ref_allele,alt_allele);
  else if genotype=2 then rs1893592=catt(alt_allele,alt_allele);
  else delete;
  drop genotype ref_allele alt_allele;
run;

proc sort data=ptpn22_snp2;
  by subject_id;
proc sort data=ubash3a_snp2;
  by subject_id;
proc sort data=counts_summed_w_key;
  by subject_id;
run;

data counts_summed_w_key_gt;
  merge counts_summed_w_key (in=in1) ptpn22_snp2 (in=in2) ubash3a_snp2 (in=in3);
  by subject_id;
  if in1;
  drop _TYPE_ _FREQ_;
run;

/* Format and export counts for wiggles */

/* I want to export these data by genotype individually. I need the following columns for the output:
   chrom, pos, count
*/

%macro exportWigs(celltype,gene,snp,genotype);
data data_for_wigs;
  set counts_summed_w_key_gt;
  where cell_type="&celltype." and gene="&gene." and &snp.="&genotype.";
  count2=0; *some dummy variables, since the R script wants 3 columns of counts I can't figure out how to change this;
  count3=0; *so make a bunch of 0 counts -- these won't show up on the plots, so it won't matter much;
  keep chrom pos count count2 count3;
run;

data data_for_wigs_log;
  set counts_summed_w_key_gt;
  where cell_type="&celltype." and gene="&gene." and &snp.="&genotype.";
  log_count=log(count+1);
  count2=0; *some dummy variables, since the R script wants 3 columns of counts I can't figure out how to change this;
  count3=0; *so make a bunch of 0 counts -- these won't show up on the plots, so it won't matter much;
  keep chrom pos log_count count2 count3;
run;


proc export data=data_for_wigs
     outfile="!PATCON/pipeline_output/wiggle_data/counts_&gene._&celltype._&snp._&genotype._2.csv"
     dbms=csv replace;
run;

proc export data=data_for_wigs_log
     outfile="!PATCON/pipeline_output/wiggle_data/log_counts_&gene._&celltype._&snp._&genotype._2.csv"
     dbms=csv replace;
run;


%mend;

%exportWigs(CD4,UBASH3A,rs1893592,AA);
%exportWigs(CD4,UBASH3A,rs1893592,AC);
%exportWigs(CD4,UBASH3A,rs1893592,CC);

%exportWigs(CD19,PTPN22,rs1217414,GG);
%exportWigs(CD19,PTPN22,rs1217414,GA);
%exportWigs(CD19,PTPN22,rs1217414,AA);

*data for wiggles is done. need to make plots!;

* also now need to figure out coverage heatmaps ;

/* Format and export counts for coverage heatmaps -- overall and by subject */


* Fill in blank positions with 0 (need this for heatmaps so they line up;

data ptpn22_pos;
  length chrom $4.;
  length gene $7.;
  format start_pos best32.;
  gene="PTPN22";
  chrom="1";
  start_pos=114356433;
run;

proc sort data=ptpn22_pos;
  by gene;
run;

data ptpn22_all_pos;
  set ptpn22_pos;
  by gene;
  if first.gene then do;
    pos=start_pos;
    output; end;
  do while(pos < 114414448);
      pos=pos+1;
      output;
  end;
run;

   
data ubash3a_pos;
  length chrom $4.;
  length gene $7.;
  format start_pos best32.;
  gene="UBASH3A";
  chrom="21";
  start_pos=43823557;
run;

proc sort data=ubash3a_pos;
  by gene;
run;

data ubash3a_all_pos;
  set ubash3a_pos;
  by gene;
  if first.gene then do;
    pos=start_pos;
    output; end;
  do while(pos < 43868176);
      pos=pos+1;
      output;
  end;
run;



*split data on genotype so I can make the plots by genotype individually then assemble them later;

%macro exportHeat(celltype,gene,snp,geno1,geno2,geno3);


data geno1 geno2 geno3;
   set counts_summed_w_key_gt;
   if cell_type="&celltype." and gene="&gene." then do;
     if &snp.="&geno1." then output geno1;
     if &snp.="&geno2." then output geno2;
     if &snp.="&geno3." then output geno3;
     end;
   keep subject_id chrom pos count;
run;

proc sort data=geno1;
  by subject_id;
proc sort data=geno2;
  by subject_id;
proc sort data=geno3;
  by subject_id;
run;


%macro scale_counts(genotype);

proc means data=&genotype. noprint;
   by subject_id;
   var count;
   output out=&genotype._max max=max_count;
run;

/* Scale counts by the max per subject */


proc sort data=&genotype.;
   by subject_id;
proc sort data=&genotype._max;
   by subject_id;
run;

data &genotype._w_max;
   merge &genotype. (in=in1) &genotype._max (in=in2);
   by subject_id;
   if in1 and in2;
run;

data &genotype._scaled;
   set &genotype._w_max;
   scaled_count=count/max_count;
   drop count max_count _TYPE_ _FREQ_;
run;

proc sort data=&genotype._scaled;
   by chrom pos;
proc transpose data=&genotype._scaled out=scaled_sbys;
   by chrom pos;
   var scaled_count;
   id subject_id;
run;

proc sort data=scaled_sbys;
  by chrom pos;
proc sort data=&gene._all_pos;
  by chrom pos;
run;

data scaled_sbys_all;
  merge &gene._all_pos (in=in1) scaled_sbys (in=in2) ;
  by chrom pos;
  if in1;
run;

data scaled_sbys_all2;
  set scaled_sbys_all;
  array change _numeric_;
       do over change;
       if change=. then change=0;
       end;
   drop _NAME_;
   run ;

proc transpose name=subject_id data=scaled_sbys_all2 out=counts_w_zeros(rename=(col1=count));
     by chrom pos;
     var M:;
run;

data &genotype._scaled_w_zeros;
   set counts_w_zeros;
run;

%mend;

%scale_counts(geno1);
%scale_counts(geno2);
%scale_counts(geno3);

/* Export */
proc export data=geno1_scaled_w_zeros outfile="!PATCON/pipeline_output/wiggle_data/log_scaled_&celltype._&gene._&snp._&geno1..csv"
   dbms=csv replace;
run;

proc export data=geno2_scaled_w_zeros outfile="!PATCON/pipeline_output/wiggle_data/log_scaled_&celltype._&gene._&snp._&geno2..csv"
   dbms=csv replace;
run;

proc export data=geno3_scaled_w_zeros outfile="!PATCON/pipeline_output/wiggle_data/log_scaled_&celltype._&gene._&snp._&geno3..csv"
   dbms=csv replace;
run;

/* Need to also create a "by-genotype" dataset for making an overall heatmap */

data subset_counts_all;
  set counts_summed_w_key_gt;
   if cell_type="&celltype." and gene="&gene.";
   keep &snp. chrom pos count;
run;

proc sort data=subset_counts_all;
   by chrom pos &snp.;
proc means data=subset_counts_all noprint;
   by chrom pos &snp.;
   var count;
   output out=mean_counts_all mean=;
run;

proc sort data=mean_counts_all;
   by &snp.;
proc means data=mean_counts_all noprint;
   by &snp.;
   var count;
   output out=subset_max max=max_count;
run;

/* Scale counts by the max per subject */


proc sort data=mean_counts_all;
   by &snp.;
proc sort data=subset_max;
   by &snp.;
run;

data subset_w_max;
   merge mean_counts_all (in=in1) subset_max (in=in2);
   by &snp.;
   if in1 and in2;
run;

data subset_scaled;
   set subset_w_max;
   scaled_count=count/max_count;
   drop count max_count _TYPE_ _FREQ_;
run;

proc sort data=subset_scaled;
   by chrom pos;
proc transpose data=subset_scaled out=scaled_sbys;
   by chrom pos;
   var scaled_count;
   id &snp.;
run;

proc sort data=scaled_sbys;
  by chrom pos;
proc sort data=&gene._all_pos;
  by chrom pos;
run;

data scaled_sbys_all;
  merge &gene._all_pos (in=in1) scaled_sbys (in=in2) ;
  by chrom pos;
  if in1;
run;

data scaled_sbys_all2;
  set scaled_sbys_all;
  array change _numeric_;
       do over change;
       if change=. then change=0;
       end;
   drop _NAME_;
   run ;

proc transpose name=&snp. data=scaled_sbys_all2 out=counts_w_zeros(rename=(col1=count));
     by chrom pos;
     var &geno1. &geno2. &geno3. ;
run;

data overall_scaled_w_zeros;
   set counts_w_zeros;
run;



proc export data=overall_scaled_w_zeros outfile="!PATCON/pipeline_output/wiggle_data/log_scaled_&celltype._&gene._&snp._overall.csv"
   dbms=csv replace;
run; 

%mend;

%exportHeat(CD4,UBASH3A,rs1893592,AA,AC,CC);
%exportHeat(CD19,PTPN22,rs1217414,GG,GA,AA);

/* Make and export BED file for splicing events 

I am going to make 4 BED files per gene:
1. All fusions
2. Fusions associated with SNP
3. All detected events in cell type
4. Detected events associated with cell type

*/

%macro exportFeatureBED(gene,celltype,snp);

* fusions;
data fus2coord;
   set fus.hg19_aceview_fusions_si_bed;
   keep chr fusion_start fusion_stop fusion_id;
run;

data fus2keep;
   set fus.hg19_aceview_fusions_si_info;
   where gene_id="&gene.";
   keep fusion_id;
run;

proc sort data=fus2coord;
   by fusion_id;
proc sort data=fus2keep nodup;
   by fusion_id;
run;

data fus2out;
  merge fus2keep (in=in1) fus2coord (in=in2);
  by fusion_id;
  if in1 and in2;
run;

data fus2out_bed;
  retain chr fusion_start fusion_stop fusion_id;
  length score $1.;
  length strand $1.;
  format start2 best32.;
  format stop2 best32.;
  length color $11.;
  format num_blocks best32.;
  format block_length best32.;
  format block_start best32.;
  set fus2out;
  score=".";
  strand=".";
  start2=fusion_start;
  stop2=fusion_stop;
  color="0,0,0";
  num_blocks=1;
  block_length=(fusion_stop-fusion_start)+1;
  block_start=0;
run;
  
*fusions with SNP;
data sig_fus;
   set eqtl.results_summary_table_w_means_v2;
   if gene_id =: "&gene." and snp_id="&snp."
         and &celltype._FDR_P < 0.05 and &celltype._FDR_P ne . and feature_type="exon";
   keep feature_id;
   rename feature_id=fusion_id;
run;

proc sort data=sig_fus nodup;
  by fusion_id;
proc sort data=fus2out_bed;
  by fusion_id;
run;

data sig_fus2out_bed;
  merge fus2out_bed (in=in1) sig_fus (in=in2);
  by fusion_id;
  if in1 and in2;
run;


*events;

data splicing_on;
   set splicing.splicing_results_clean;
   if gene_id="&gene." and flag_&celltype._on=1;
   keep event_id;
run;

proc sort data=splicing_on;
   by event_id;
run;

/* Event annotations */

data splicing_annot;
   set splice.splicing_events_annotations;
   if gene_id="&gene.";
run;

proc sort data=splicing_annot;
   by event_id;
proc sort data=splicing_on;
   by event_id;
run;


data events_w_annot;
   merge splicing_annot (in=in1) splicing_on (in=in2);
   by event_id;
   if in1 and in2;
run;

/* Prepare BED file */

data info_for_bed;
   length score $1.;
   length color $11.;
   set events_w_annot;
   score='.';
   block1_start=0;
   if flag_intron_retention=1 then do;
      num_blocks=1;
      if strand='+' then do;
         block1_length=(feature2_stop +1) - feature1_start;
         totalstart1=feature1_start;
         totalstop1=feature2_stop+1;
         totalstart2=feature1_start;
         totalstop2=feature2_stop+1;
         end;
      if strand='-' then do;
         block1_length=feature2_stop - (feature1_start - 1);
         totalstart1=feature1_start-1;
         totalstop1=feature2_stop;
         totalstart2=feature1_start-1;
         totalstop2=feature2_stop;
         end;
      if flag_eqtl=1 then color='0,255,255';
      else color='0,0,255';
      end;
   else do;
      num_blocks=2;
          totalstart1=feature1_start;
          totalstop1=feature2_stop;
          totalstart2=feature1_start;
          totalstop2=feature2_stop;
      block1_length=(feature1_stop+1)-feature1_start;
      block2_length=(feature2_stop+1)-feature2_start;
      block2_start=feature2_start-feature1_start;
      if flag_eqtl=1 then color='255,0,0';
      else color='255,255,0';
      end;
   keep event_id chr totalstart1 totalstop1 totalstart2 totalstop2 strand score color num_blocks block1_length block2_length block1_start block2_start;
run;

data event_bed;
    retain chr;
    retain totalstart1;
    retain totalstop1;
    retain event_id;
    retain score;
    retain strand;
    retain totalstart2;
    retain totalstop2;
    retain color;
    retain num_blocks;
    length lengths_cat $10.;
    length starts_cat $10.;
    set info_for_bed;
    if num_blocks=2 then do;
        lengths_cat=catx(',', block1_length, block2_length);
        starts_cat=catx(',', block1_start, block2_start);
        end;
    else if num_blocks=1 then do;
        lengths_cat=put(block1_length, 5.);
        starts_cat=put(block1_start, 1.);
        end;
    drop block1_length block2_length block1_start block2_start;
run;

*events with SNP;

data event_eqtl;
   set eqtl.eqtl_results_summary_table;
   where gene_id ? "&gene." and snp_id="&snp." and &celltype._FDR_P lt 0.05;
   keep feature_id;
   rename feature_id=event_id;
run;

proc sort data=event_eqtl nodup;
  by event_id;
proc sort data=event_bed;
  by event_id;
run;

data event_bed_eqtl;
  merge event_bed (in=in1) event_eqtl (in=in2);
  by event_id;
  if in1 and in2;
run;

proc sort data=event_bed;
   by chr totalstart1 totalstop1 strand;
proc sort data=event_bed_eqtl;
   by chr totalstart1 totalstop1 strand;
proc sort data=fus2out_bed;
   by chr fusion_start fusion_stop;
proc sort data=sig_fus2out_bed;
   by chr fusion_start fusion_stop;
run;

proc export data=event_bed
   outfile="!PATCON/pipeline_output/wiggle_data/&gene._&celltype._events_detected.bed"
   dbms=tab replace; putnames=no;
run;

proc export data=event_bed_eqtl
   outfile="!PATCON/pipeline_output/wiggle_data/&gene._&celltype._events_detected_with_eqtl.bed"
   dbms=tab replace; putnames=no;
run;

proc export data=fus2out_bed
   outfile="!PATCON/pipeline_output/wiggle_data/&gene._&celltype._exons.bed"
   dbms=tab replace; putnames=no;
run;

proc export data=sig_fus2out_bed
   outfile="!PATCON/pipeline_output/wiggle_data/&gene._&celltype._exons_with_eqtl.bed"
   dbms=tab replace; putnames=no;
run;

%mend;


%exportFeatureBED(UBASH3A,CD4,rs1893592);

%exportFeatureBED(PTPN22,CD19,rs1217414);
*no significant exons for PTPN22, so the export of the final BED file will fail;


/* Export counts by significant feature */

data sig_features;
   set eqtl.results_summary_table_w_means_v2;
   if gene_id =: "UBASH3A" and snp_id="rs1893592"
         and CD4_FDR_P < 0.05 and CD4_FDR_P ne . then output;
   if gene_id =: "PTPN22" and snp_id="rs1217414"
         and CD19_FDR_P < 0.05 and CD19_FDR_P ne . then output;
   keep gene_id snp_id feature_id feature_type CD4_FDR_P CD19_FDR_P 
        cd4_hmz1_expression cd4_htz_expression cd4_hmz2_expression
        cd19_hmz1_expression cd19_htz_expression cd19_hmz2_expression;
run;

proc export data=sig_features
     outfile="!PATCON/pipeline_output/wiggle_data/sig_eqtl_for_plot_annotation.csv"
     dbms=csv replace;
run;





