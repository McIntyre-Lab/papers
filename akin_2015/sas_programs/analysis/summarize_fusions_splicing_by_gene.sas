/* set libraries */
libname sgrloc '/media/jrbnewman/SAS_WRK1/sugrue/sas_data';
libname splice '/media/jrbnewman/SAS_WRK1/splice';
libname sugrue '/home/jrbnewman/McLab/sugrue/sas_data';


/* Merge DE fusions and DE splicing by gene and export */

/* sort then merge */
proc sort data=sugrue.splicing_by_gene_summary;
   by gene_id;
run;


data genes_de_summary_w_fus;
   set sugrue.results_gene_summary;
   keep gene_id exp_exon_count total_exon_count

proc sort data=sugrue.genes_de_summary_w_fus;
   by gene_id;
run;


data genes_fus_se_summary;
    merge sugrue.genes_de_summary_w_fus (in=in1) sugrue.splicing_by_type_summary (in=in2);
    by gene_id;
    if in1 and in2 then output;
    else if in1 then do;
       detected_se_cnt=0;
       detected_ir_cnt=0;
       analyzed_se_cnt=0;
       analyzed_ir_cnt=0;
       diffexp_se_cnt=0;
       diffexp_ir_cnt=0;
       exons_se_cat=',';
       num_events_de=',';
       events_de_cat=',';
       num_events_detected=',';
       events_detect_cat=',';
       output;
       end;
    else do;
      exon_exp_count=0;
      *This should never be needed if there are splicing events!;
      exon_total_count=0; 
      exon_de_count=0;
      exon_updown_cat=',';
      fusions_cat=',';
      exons_cat=',';
      output;
      end;
run;


/* Add flags */

data genes_fus_se_w_flags;
    set genes_fus_se_summary;
    if exon_de_count > 0 then flag_de_exons=1;
    else flag_de_exons=0;

    if diffexp_se_cnt > 0 then flag_de_splicing=1;
    else flag_de_splicing=0;

    if diffexp_ir_cnt > 0 then flag_de_intron_ret=1;
    else flag_de_intron_ret=0;

    * Flag if there are more DE splicing events than DE exons;
    * Expect DE splicing to be no greater than DE exons;
    * Although because using fusions, possibly multiple SE per exon pair; 
    if diffexp_se_cnt > exon_de_count then flag_splice_exon_diff=1;
    else flag_splice_exon_diff=0;

run;

/* Save permenant and export as CSV */

data sugrue.genes_fus_splice_summary;
   set genes_fus_se_w_flags;
run;


proc export data=sugrue.genes_fus_splice_summary
outfile='/home/jrbnewman/McLab/sugrue/generated_files/genes_fusions_splicing_summary.csv'
        dbms=csv replace;
        run;
