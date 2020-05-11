
libname relapse "/home/ammorse/staph_relapse/sasdata";
libname lmm  "/home/ammorse/staph_relapse/sasdata/data_analysis";

/*
link between sra_sampleID and FQ files  

*** NOTE - not keeping sampleID 3867- this is pair to failed library (sampleID 3938) 
*/

data sra_num ;
set relapse.sampleID_st_cc_ref_kitchen_sink ;
where flag_failed_library  =0;
keep sra_sampleID sampleID sampleNum ;
run;

data fq ;
set relapse.fq_list_w_samplenum_sampleID ;
where flag_failed_library = 0 ;
keep sampleNum fqName sampleID ;
run ;

proc sort data = sra_num ;
by sampleID sampleNum ;
proc sort data = fq ;
by sampleID sampleNum ;
run;

data fq_sra only_fq oops ;
merge fq (in=in1) sra_num (in=in2) ;
by sampleID sampleNum ;
if in1 then output fq_sra ;
else if in1 then output only_fq ;
else output oops ;
run ;  /* 0  in fq_only and 0  in oops */

data relapse.link_fqName_2_sra_sampleID ;
set fq_sra ;
keep fqName sra_sampleID ;
run;

proc export data = relapse.link_fqName_2_sra_sampleID 
outfile = "/home/ammorse/staph_relapse/supplementary_material/link_fqName_2_sra_sampleID.csv"
/* outfile = "!MCLAB/staph/seoung_ho_project/supplementary_material/link_fqName_2_sra_sampleID.csv" */
dbms = tab replace ;
run ;


