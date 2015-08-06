/*
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

/* Create local copy of main dataset */
    data cegsV;
        set SEM.cegs_virgin_norm_cent;
        run;

/* Pull in sex det corrected symbols */
    data sex;
        set SEM.cegsV_sex_det_stack;
        rename symbol_cat = sex;
        keep fusion_id symbol_cat;
        run;

    proc sort data=sex nodupkey;
        by fusion_id;
        run;

/* Pull in FBgns */
    data gene_name;
        set cegsV;
        keep fusion_id fbgn_cat;
        run;

    proc sort data=gene_name nodupkey;
        by fusion_id;
        run;

    data gene_name2;
        set gene_name;
        primary_fbgn = scan(fbgn_cat,1,'|');
        keep fusion_id primary_fbgn;
        run;

/* For sex Det gene replace fbgn with corrected gene symbol */
    data gene_name3;
        merge sex (in=in1) gene_name2 (in=in2);
        by fusion_id;
        if in1 then sym = sex;
        else sym = primary_fbgn;
        if in2 then output gene_name3;
        drop sex; 
        run;

    * Create a gene to symbol key;
    data cegsV_gene2fbgn;
        retain sym;
        set gene_name3;
        rename sym = gene;
        drop fusion_id;
        run;

    proc sort data=cegsV_gene2fbgn nodupkey;
        by gene;
        run;

    * Two of the FBgn's are wrong, because of multigene fusions, fix them here;
    data SEM.cegsV_gene2fbgn;
        set cegsV_gene2fbgn;
        if gene eq 'tra' then primary_fbgn = 'FBgn0003741';
        if gene eq 'sqd' then primary_fbgn = 'FBgn0263396';
        run;

/* Merge sym to the main dataset */
    proc sort data=gene_name3;
        by fusion_id;
        run;

    proc sort data=cegsV;
        by fusion_id;
        run;

    data cegsV_w_name;
        merge cegsV (in=in1) gene_name3;
        by fusion_id;
        drop primary_fbgn;
        run;

/* Calculate mean expression for each gene */
    proc sort data=cegsV_w_name;
        by line sym;
        run;

    proc means data = cegsV_w_name noprint;
        by line sym;
        output out=means mean(mean_exp)=mean_exp;
        run;

    data SEM.cegsV_by_gene_stack;
        set means;
        drop _type_ _freq_;
        run;

/* Flip and create permanent dataset */
    proc transpose data=means out=flip;
        by line;
        var mean_exp;
        id sym;
        run;

    data SEM.cegsV_by_gene_sbs;
        set flip;
        drop _name_;
        run;

    data SEM.cegsV_by_gene_sex_det_sbs;
        set flip;
        drop _name_ FBgn:;
        run;

/* Make gene list */
    data gene_list;
        set gene_name3;
        keep sym;
        run;

    proc sort data=gene_list nodupkey;
        by sym;
        run;

    data SEM.cegsV_gene_list;
        set gene_list;
        rename sym = symbol;
        run;

    proc export data=SEM.cegsV_gene_list outfile='!MCLAB/cegs_sem_sd_paper/design_file/cegsV_gene_list.csv' dbms=csv replace;
        putnames = no;
        run;

/* Mean Center by Gene */
    proc sort data=SEM.cegsV_by_gene_stack;
        by sym;
        run;

    proc standard data=SEM.cegsV_by_gene_sbs mean=0 out=gene_std ;
        run;

    data SEM.cegsV_by_gene_sbs_mean_center;
        set gene_std;
        run;
