/*******************************************************************************
* Filename: cegsV_combine_adding_links_yp2.sas
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

libname addgene '!MCLAB/cegs_sem_sd_paper/analysis_output/adding_newlinks/cegs_adding_links_yp2/sas_data';

data addlink;
    set addgene.none;
    drop gene;
    run;

/* Create Path */

data SEM.cegsV_al_yp2_model_design_file;
    informat model $14.;
    informat modelnum best12.;
    informat path $30.;
    input model$ modelnum path$;
    datalines;
Model_Baseline   0   Baseline
Model_1   1   Sxl->Yp2
Model_2   2   Sxl->dsx
Model_3   3   Sxl->fru
Model_4   4   Yp2->Sxl
Model_5   5   Yp2->fru
Model_6   6   Yp2->tra
Model_7   7   dsx->Sxl
Model_8   8   dsx->fru
Model_9   9   fru->Sxl
Model_10  10  fru->Yp2
Model_11  11  fru->dsx
Model_12  12  tra->Yp2
Model_13  13  Spf45->Yp2
Model_14  14  Spf45->dsx
Model_15  15  Spf45->fru
Model_16  16  Spf45->tra
Model_17  17  fl_2_d->Yp2
Model_18  18  fl_2_d->dsx
Model_19  19  fl_2_d->fru
Model_20  20  her->Sxl
Model_21  21  her->dsx
Model_22  22  her->fru
Model_23  23  her->tra
Model_24  24  ix->Sxl
Model_25  25  ix->dsx
Model_26  26  ix->fru
Model_27  27  ix->tra
Model_28  28  snf->Yp2
Model_29  29  snf->dsx
Model_30  30  snf->fru
Model_31  31  snf->tra
Model_32  32  tra2->Sxl
Model_33  33  tra2->Yp2
Model_34  34  tra2->tra
Model_35  35  vir->Yp2
Model_36  36  vir->dsx
Model_37  37  vir->fru
Model_38  38  Yp2->Spf45
Model_39  39  dsx->Spf45
Model_40  40  fru->Spf45
Model_41  41  tra->Spf45
Model_42  42  Yp2->fl_2_d
Model_43  43  dsx->fl_2_d
Model_44  44  fru->fl_2_d
Model_45  45  Sxl->her
Model_46  46  dsx->her
Model_47  47  fru->her
Model_48  48  tra->her
Model_49  49  Sxl->ix
Model_50  50  dsx->ix
Model_51  51  fru->ix
Model_52  52  tra->ix
Model_53  53  Yp2->snf
Model_54  54  dsx->snf
Model_55  55  fru->snf
Model_56  56  tra->snf
Model_57  57  Sxl->tra2
Model_58  58  Yp2->tra2
Model_59  59  tra->tra2
Model_60  60  Yp2->vir
Model_61  61  dsx->vir
Model_62  62  fru->vir
Model_63  63  fl_2_d->Spf45
Model_64  64  her->Spf45
Model_65  65  ix->Spf45
Model_66  66  snf->Spf45
Model_67  67  tra2->Spf45
Model_68  68  vir->Spf45
Model_69  69  Spf45->fl_2_d
Model_70  70  her->fl_2_d
Model_71  71  ix->fl_2_d
Model_72  72  snf->fl_2_d
Model_73  73  tra2->fl_2_d
Model_74  74  vir->fl_2_d
Model_75  75  Spf45->her
Model_76  76  fl_2_d->her
Model_77  77  ix->her
Model_78  78  snf->her
Model_79  79  tra2->her
Model_80  80  vir->her
Model_81  81  Spf45->ix
Model_82  82  fl_2_d->ix
Model_83  83  her->ix
Model_84  84  snf->ix
Model_85  85  tra2->ix
Model_86  86  vir->ix
Model_87  87  Spf45->snf
Model_88  88  fl_2_d->snf
Model_89  89  her->snf
Model_90  90  ix->snf
Model_91  91  tra2->snf
Model_92  92  vir->snf
Model_93  93  Spf45->tra2
Model_94  94  fl_2_d->tra2
Model_95  95  her->tra2
Model_96  96  ix->tra2
Model_97  97  snf->tra2
Model_98  98  vir->tra2
Model_99  99  Spf45->vir
Model_100 100 fl_2_d->vir
Model_101 101 her->vir
Model_102 102 ix->vir
Model_103 103 snf->vir
Model_104 104 tra2->vir
;
run;


/* merge on paths */
proc sort data=addlink;
    by model;
    run;

proc sort data=SEM.cegsV_al_yp2_model_design_file;
    by model;
    run;

data SEM.cegsV_al_yp2_stack_bic;
    retain model modelnum path BIC;
    merge addlink (in=in1) SEM.cegsV_al_yp2_model_design_file (in=in2);
    by model;
    run;

/* Sort by BIC*/
proc sort data=SEM.cegsV_al_yp2_stack_bic;
    by BIC;
    run;

proc export data=SEM.cegsV_al_yp2_stack_bic outfile='!MCLAB/cegs_sem_sd_paper/analysis_output/adding_newlinks/cegsV_al_yp2_stack_bic.csv' dbms=csv replace;
    putnames=yes;
    run;

/* Clean Up */
proc datasets nolist;
    delete addlink;
    run; quit;

