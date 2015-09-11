libname sgrloc '/media/jrbnewman/SAS_WRK1/sugrue/sas_data';
libname splice '/media/jrbnewman/SAS_WRK1/splice';
libname sugrue '/home/jrbnewman/McLab/sugrue/sas_data';

/* Transpose depth by sample */

data fusion_data_for_tpose;
   set sugrue.counts_w_flags;
   keep subject fusion_id apn;
run;

proc sort data=fusion_data_for_tpose;
    by fusion_id subject;
run;

proc transpose data=fusion_data_for_tpose out=sugrue.fusion_data_flip;
   by fusion_id;
   var apn;
   id subject;
run;


