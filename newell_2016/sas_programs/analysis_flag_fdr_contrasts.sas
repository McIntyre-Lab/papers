/* Calculate FDR  and get FDR flags. Need the FDR contrasts to be in a more usable format */

libname ribo '!MCLAB/arbeitman/arbeitman_ribotag/sas_data';


/* First do FDR */

%include '!MCLAB/arbeitman/arbeitman_ribotag/sas_programs/FDR_mult_plus_indicators.sas';

%fdr_mult(ribo.Output_contrasts_apn0, probf, fusion_id, 0.20,0.1,0.05);

data ribo.output_fdr_probf_apn0;
  set fdr_probf;
  run;


/*Now transpose the data; first rename the contrasts*/
proc freq data = ribo.output_fdr_probf_apn0;
  table label ;
  run;

data contrast;
  retain fusion_id label;
  set ribo.output_fdr_probf_apn0;
  label raw_p = 'p_contrast';
  label FValue = 'f_contrast';
    rename raw_p  = p_contrast
           FValue = f_contrast;
     if Label = 'IPmale-IPfemale'                 	then contrast =  1 ;
     if Label = 'IPmale-InputMale'           		then contrast =  2 ;
     if Label = 'IPfemale-InputFemale'             	then contrast =  3 ;
     if Label = 'InputMale-InputFemale'               	then contrast =  4 ;
  run;

proc sort data=contrast;
  by fusion_id contrast;
  run;


/*Transpose the contrasts each column will have the value p_all contrast# but will be labeled with the contrast label.*/

proc transpose data=contrast out= p_trans prefix=p_contrast_ ;
  by fusion_id;
  id contrast;
  var p_contrast;
  run;

data p_trans;
    set p_trans;
    label p_contrast_1  =  'p_contrast_IPmale-IPfemale'  	;
    label p_contrast_2  =  'p_contrast_IPmale-InputMale'	;
    label p_contrast_3  =  'p_contrast_IPfemale-InputFemale'    ;
    label p_contrast_4  =  'p_contrast_InputMale-InputFemale'   ;
  run;

%macro flag_pvals (pval, label);
    data flag_&pval;
        set p_trans;
        if &pval = . then flag_&pval._05 = .;
            else if &pval < 0.05 then flag_&pval._05 = 1;
            else flag_&pval._05 = 0;
        label flag_&pval._05 = "flag_p_contrast_05_&label.";
        keep fusion_id &pval flag_&pval._05;
        run;
    
    proc sort data = flag_&pval;
        by fusion_id;
        run;

%mend;

proc contents data=p_trans varnum;
run;

*flag p values;
%flag_pvals(p_contrast_1,IPmale-IPfemale)		;
%flag_pvals(p_contrast_2,IPmale-InputMale)		;
%flag_pvals(p_contrast_3,IPfemale-InputFemale)		;
%flag_pvals(p_contrast_4,InputMale-InputFemale)		;

/*Calculate FDR for the contrasts, then create a flag if the fdr value is significant */

proc transpose data=contrast out=fdr_trans prefix=fdr_p_contrast_ ;
    by fusion_id;
    id contrast;
    var fdr_probF;
    run;

data fdr_trans;
    set fdr_trans;
    label fdr_p_contrast_1  =  'fdr_p_contrast_IPmale-IPfemale'               	;
    label fdr_p_contrast_2  =  'fdr_p_contrast_IPmale-InputMale'         	;
    label fdr_p_contrast_3  =  'fdr_p_contrast_IPfemale-InputFemale'           	;
    label fdr_p_contrast_4  =  'fdr_p_contrast_InputMale-InputFemale'          	;
  run;

%macro flag_fdrs (fdr, label);

    data flag_&fdr;
        set fdr_trans;

        if &fdr = . then flag_&fdr._20= .;
        else if &fdr < 0.20 then flag_&fdr._20 = 1;
        else flag_&fdr._20 = 0;
        label flag_&fdr._20 = "flag_fdr_20_&label.";

        if &fdr = . then flag_&fdr._10= .;
        else if &fdr < 0.10 then flag_&fdr._10 = 1;
        else flag_&fdr._10 = 0;
        label flag_&fdr._10 = "flag_fdr_10_&label.";

        if &fdr = . then flag_&fdr._05 = .;
        else if &fdr < 0.05 then flag_&fdr._05 = 1;
        else flag_&fdr._05 = 0;
        label flag_&fdr._05 = "flag_fdr_05_&label.";

        keep fusion_id &fdr flag_&fdr._20 flag_&fdr._10 flag_&fdr._05;
        run;

    proc sort data = flag_&fdr;
        by fusion_id;
        run;

%mend ;


proc contents data=fdr_trans varnum;
run;

*flag fdrs;
%flag_fdrs(fdr_p_contrast_1,IPmale-IPfemale)		;
%flag_fdrs(fdr_p_contrast_2,IPmale-InputMale)		;
%flag_fdrs(fdr_p_contrast_3,IPfemale-InputFemale)	;
%flag_fdrs(fdr_p_contrast_4,InputMale-InputFemale)	;


/*Merge everything together*/

data ribo.flag_fdr_contrast_by_fusion;
  merge flag_p_contrast_1 flag_p_contrast_2 flag_p_contrast_3 flag_p_contrast_4 
        flag_fdr_p_contrast_1 flag_fdr_p_contrast_2 flag_fdr_p_contrast_3 flag_fdr_p_contrast_4;
  by fusion_id;
run;

