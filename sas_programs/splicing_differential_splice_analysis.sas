/*******************************************************************************
* Filename: misc_test_splice_test.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Run a basic splice model to determine if there is a splicing
* effect for each gene.
*
*******************************************************************************/

/* Libraries
    libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
    libname dmel551 '!MCLAB/useful_dmel_data/flybase551/sas_data';
    filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
    options SASAUTOS=(sasautos mymacros);
*/

/* Create Gene List */
data gene_list;
    set SEM.cegsV_splice_data;
    keep symbol_cat;
    run;

proc sort data=gene_list nodupkey;
    by symbol_cat;
    run;


%macro splice(gene);

goptions reset=all;
ods listing close;
ods pdf file="!MCLAB/cegs_sem_sd_paper/analysis_output/splicing/cegsV_splicing_model/&gene..pdf";
proc glimmix data=sem.cegsV_splice_data plots=studentpanel;
    by symbol_cat;
	where symbol_cat="&gene";
    class fusion_id line rep;
    model uq_log_uq_center = line|fusion_id /htype=1;
    lsmeans line*fusion_id /slice=line slice=fusion_id;
    *output out=resid_by_symbol r=resid p=pred;
    *ods output Tests1=model_ts FitStatistics=model_fs;
    run;

ods pdf close;
ods listing;

%mend;
%iterdataset(dataset=gene_list, function=%nrstr(%splice(&symbol_cat);));

/* Clean up */
proc datasets ;
    delete gene_list;
    run; quit;
