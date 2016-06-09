libname cegs "Z:/cegs_ase_paper/sas_data";
filename mymacros "Z:/maize_ozone/2014/sas_analysis/macros";
options SASAUTOS=(sasautos mymacros);

 data WORK.m_only    ;
       %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
       infile 'Z:\cegs_ase_paper\output\gene_categories\M_only_genes.csv' delimiter = ','
 MISSOVER DSD lrecl=32767 firstobs=2 ;
          informat symbol $25. ;
          informat gene_compare_enrich $6. ;
          informat COUNT best32. ;
          informat PERCENT best32. ;
          informat VAR5 $87. ;
          format symbol $25. ;
          format gene_compare_enrich $6. ;
          format COUNT best12. ;
          format PERCENT best12. ;
          format VAR5 $87. ;
       input
                   symbol $
                   gene_compare_enrich $
                   COUNT
                   PERCENT
                   VAR5 $
       ;
       if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
       run;

data WORK.v_only    ;
       %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
       infile 'Z:\cegs_ase_paper\output\gene_categories\V_only_genes.csv' delimiter = ','
 MISSOVER DSD lrecl=32767 firstobs=2 ;
         informat symbol $25. ;
          informat gene_compare_enrich $6. ;
          informat COUNT best32. ;
          informat PERCENT best32. ;
          informat VAR5 $1. ;
          format symbol $25. ;
          format gene_compare_enrich $6. ;
          format COUNT best12. ;
          format PERCENT best12. ;
          format VAR5 $1. ;
       input
                   symbol $
                   gene_compare_enrich $
                   COUNT
                   PERCENT
                   VAR5 $
       ;
       if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
       run;


data WORK.final_results2    ;
       %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
       infile 'Z:\cegs_ase_paper\manuscript\tables\Supplemental_Table3.csv' delimiter = ','
 MISSOVER DSD lrecl=32767 firstobs=2 ;
          informat genotype $4. ;
          informat mating_status $1. ;
          informat exonic_region $9. ;
          informat mean_AI best32. ;
          informat AI_intersection_test best32. ;
          informat AI_qsim_test best32. ;
          informat AI_intersection_qsim_both best32. ;
          informat counts_in_genotype best32. ;
          informat counts_in_tester best32. ;
          informat counts_unassigned best32. ;
          informat counts_total best32. ;
          informat mean_apn best32. ;
          informat analyze_cis_effects best32. ;
          informat cis_effects_genotype best32. ;
          informat trans_effects_genotype best32. ;
          informat direction_cis $1. ;
          informat direction_trans $1. ;
          informat chromosome $2. ;
          informat start_site best32. ;
          informat end_site best32. ;
          informat gene_name $25. ;
          informat Flybase551_FBgn $23. ;
          informat flag_AI_alternate_direction best32. ;
          format genotype $4. ;
          format mating_status $1. ;
          format exonic_region $9. ;
          format mean_AI best12. ;
         format AI_intersection_test best12. ;
         format AI_qsim_test best12. ;
         format AI_intersection_qsim_both best12. ;
         format counts_in_genotype best12. ;
         format counts_in_tester best12. ;
         format counts_unassigned best12. ;
         format counts_total best12. ;
         format mean_apn best12. ;
         format analyze_cis_effects best12. ;
         format cis_effects_genotype best12. ;
         format trans_effects_genotype best12. ;
         format direction_cis $1. ;
         format direction_trans $1. ;
         format chromosome $2. ;
         format start_site best12. ;
         format end_site best12. ;
         format gene_name $25. ;
         format Flybase551_FBgn $23. ;
         format flag_AI_alternate_direction best12. ;
      input
                  genotype $
                  mating_status $
                  exonic_region $
                  mean_AI
                  AI_intersection_test
                  AI_qsim_test
                  AI_intersection_qsim_both
                  counts_in_genotype
                  counts_in_tester
                  counts_unassigned
                  counts_total
                  mean_apn
                  analyze_cis_effects
                  cis_effects_genotype
                  trans_effects_genotype
                  direction_cis $
                  direction_trans $
                  chromosome $
                  start_site
                  end_site
                  gene_name $
                  Flybase551_FBgn $
                  flag_AI_alternate_direction
      ;
      if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
      run;

*remove - . and : from symbol;
data m_only2;
set m_only;
rx=rxparse("$'.:-' to '_'");
mycom=compress(symbol);
call rxchange (rx,999,mycom);
comname = compress(mycom);
drop mycom rx;
run;

data v_only2;
set v_only;
rx=rxparse("$'.:-' to '_'");
mycom=compress(symbol);
call rxchange (rx,999,mycom);
comname = compress(mycom);
drop mycom rx;
run;

data final_results3;
set final_results2;
rx=rxparse("$'.:-' to '_'");
mycom=compress(gene_name);
call rxchange (rx,999,mycom);
new_name = compress(mycom);
drop mycom rx;
run;

/*mated only genes*/
%macro mate (name);
data M_gene_&name;
set final_results3;
where new_name contains "&name" and mating_status = "M";
run;

%mend;
%iterdataset(dataset=m_only2, function=%nrstr(%mate(&comname);));

data stack_m;
set M_gene_: ;
keep chromosome start_site end_site gene_name new_name flybase551_fbgn 
genotype mating_status mean_AI exonic_region;
run;
/* missing: mt_trna_l_cun */

proc datasets noprint;
delete M_gene_: ;
run; quit;

/*virgin_only genes*/
%macro virgin (name);
data V_gene_&name;
set final_results3;
where new_name contains "&name" and mating_status = "V";
run;

%mend;
%iterdataset(dataset=v_only2, function=%nrstr(%virgin(&comname);));

data stack_V;
set V_gene_: ;
keep chromosome start_site end_site gene_name new_name flybase551_fbgn 
genotype mating_status mean_AI exonic_region;
run;
/* all ran! */

proc datasets noprint;
delete V_gene_: ;
run; quit;

proc export data=stack_m
outfile = "Z:/cegs_ase_paper/output/m_only_genes_w_AI_mean.csv"
dbms=csv replace;
putnames=yes;
run;

proc export data=stack_v
outfile= "Z:/cegs_ase_paper/output/v_only_genes_w_AI_mean.csv"
dbms=csv replace;
putnames=yes;
run;




/*find a mean for the gene*/
proc sort data=stack_V;
by new_name;
run;

proc means data=stack_V ;
var mean_AI;
by new_name;
output out=mean_V mean=gene_mean_AI;
run;
