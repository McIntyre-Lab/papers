/*******************************************************************************
* Filename: cegsV_ag_create_model_design_file.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Models are created using the same iterative process for each
* genes.  So Model 1 will be the same for all genes. I wan't to create a design
* file to relate what location the gene was added to the model number.  CEGS
* will be different from DGPR because I used isoforms in DGPR.
*
*******************************************************************************/

/* Libraries
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

data SEM.cegsV_ag_model_design_file;
    informat model $14.;
    informat modelnum best12.;
    informat path $30.;
    input model$ modelnum path$;
    datalines;
        Model_Baseline   0  .
        Model_1   1    Sxl->gene
        Model_2   2    Yp->gene
        Model_3   3    dsx->gene
        Model_4   4    fru->gene
        Model_5   5    tra->gene
        Model_6   6    Spf45->gene
        Model_7   7    fl_2_d->gene
        Model_8   8    her->gene
        Model_9   9    ix->gene
        Model_10  10   snf->gene
        Model_11  11   tra2->gene
        Model_12  12   vir->gene
        Model_13  13   gene->Sxl
        Model_14  14   gene->Yp
        Model_15  15   gene->dsx
        Model_16  16   gene->fru
        Model_17  17   gene->tra
        Model_18  18   gene->Spf45
        Model_19  19   gene->fl_2_d
        Model_20  20   gene->her
        Model_21  21   gene->ix
        Model_22  22   gene->snf
        Model_23  23   gene->tra2
        Model_24  24   gene->vir
        Model_25  25   dsx->gene->Yp
        Model_26  26   tra->gene->dsx
        Model_27  27   tra->gene->fru
        Model_28  28   Spf45->gene->Sxl
        Model_29  29   fl_2_d->gene->Sxl
        Model_30  30   fl_2_d->gene->tra
        Model_31  31   her->gene->Yp
        Model_32  32   ix->gene->Yp
        Model_33  33   snf->gene->Sxl
        Model_34  34   tra2->gene->dsx
        Model_35  35   tra2->gene->fru
        Model_36  36   vir->gene->Sxl
        Model_37  37   vir->gene->tra
        ;
    run;
