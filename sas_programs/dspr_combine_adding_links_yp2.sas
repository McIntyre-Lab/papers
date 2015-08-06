/*******************************************************************************
* Filename: dspr_combine_adding_links_yp2.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Import BIC score from iteratively adding new paths in the SD
* pathway. Format dataset and sort by BIC.
*******************************************************************************/

/* Libraries
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

libname addgene '!MCLAB/cegs_sem_sd_paper/analysis_output/adding_newlinks/dspr_adding_links_yp2/sas_data';

data addlink;
    set addgene.none;
    drop gene;
    run;

/* Create Path */

data SEM.dspr_al_yp2_model_design_file;
    informat model $14.;
    informat modelnum best12.;
    informat path $30.;
    input model$ modelnum path$;
    datalines;
Model_Baseline   0   Baseline
Model_1    1    Sxl->Yp2
Model_2    2    Sxl->fru
Model_3    3    Yp2->Sxl
Model_4    4    Yp2->fru
Model_5    5    fru->Sxl
Model_6    6    fru->Yp2
Model_7    7    Spf45->Yp2
Model_8    8    Spf45->fru
Model_9    9    Spf45->tra
Model_10   10   fl_2_d->Yp2
Model_11   11   fl_2_d->fru
Model_12   12   her->Sxl
Model_13   13   her->fru
Model_14   14   her->tra
Model_15   15   ix->Sxl
Model_16   16   ix->fru
Model_17   17   ix->tra
Model_18   18   snf->Yp2
Model_19   19   snf->fru
Model_20   20   snf->tra
Model_21   21   tra2->Sxl
Model_22   22   tra2->tra
Model_23   23   vir->Yp2
Model_24   24   vir->fru
Model_25   25   Yp2->Spf45
Model_26   26   fru->Spf45
Model_27   27   tra->Spf45
Model_28   28   Yp2->fl_2_d
Model_29   29   fru->fl_2_d
Model_30   30   Sxl->her
Model_31   31   fru->her
Model_32   32   tra->her
Model_33   33   Sxl->ix
Model_34   34   fru->ix
Model_35   35   tra->ix
Model_36   36   Yp2->snf
Model_37   37   fru->snf
Model_38   38   tra->snf
Model_39   39   Sxl->tra2
Model_40   40   tra->tra2
Model_41   41   Yp2->vir
Model_42   42   fru->vir
Model_43   43   fl_2_d->Spf45
Model_44   44   her->Spf45
Model_45   45   ix->Spf45
Model_46   46   snf->Spf45
Model_47   47   tra2->Spf45
Model_48   48   vir->Spf45
Model_49   49   Spf45->fl_2_d
Model_50   50   her->fl_2_d
Model_51   51   ix->fl_2_d
Model_52   52   snf->fl_2_d
Model_53   53   tra2->fl_2_d
Model_54   54   vir->fl_2_d
Model_55   55   Spf45->her
Model_56   56   fl_2_d->her
Model_57   57   ix->her
Model_58   58   snf->her
Model_59   59   tra2->her
Model_60   60   vir->her
Model_61   61   Spf45->ix
Model_62   62   fl_2_d->ix
Model_63   63   her->ix
Model_64   64   snf->ix
Model_65   65   tra2->ix
Model_66   66   vir->ix
Model_67   67   Spf45->snf
Model_68   68   fl_2_d->snf
Model_69   69   her->snf
Model_70   70   ix->snf
Model_71   71   tra2->snf
Model_72   72   vir->snf
Model_73   73   Spf45->tra2
Model_74   74   fl_2_d->tra2
Model_75   75   her->tra2
Model_76   76   ix->tra2
Model_77   77   snf->tra2
Model_78   78   vir->tra2
Model_79   79   Spf45->vir
Model_80   80   fl_2_d->vir
Model_81   81   her->vir
Model_82   82   ix->vir
Model_83   83   snf->vir
Model_84   84   tra2->vir
;
run;


/* merge on paths */
proc sort data=addlink;
    by model;
    run;

proc sort data=SEM.dspr_al_yp2_model_design_file;
    by model;
    run;

data SEM.dspr_al_yp2_stack_bic;
    retain model modelnum path BIC;
    merge addlink (in=in1) SEM.dspr_al_yp2_model_design_file (in=in2);
    by model;
    run;

/* Sort by BIC*/
proc sort data=SEM.dspr_al_yp2_stack_bic;
    by BIC;
    run;

proc export data=SEM.dspr_al_yp2_stack_bic outfile='!MCLAB/cegs_sem_sd_paper/analysis_output/adding_newlinks/dspr_al_yp2_stack_bic.csv' dbms=csv replace;
    putnames=yes;
    run;

/* Clean Up */
proc datasets nolist;
    delete addlink;
    run; quit;

