libname seq "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType/sasdata" ;



/*
ALL genes!

NOTE - check where fragment_id = "F327_SI:12" or fragment_id = "F851_SI:34"

merge on/off flags to data 
drop fragment if off in all CD4 groups (group = cellType by sex by case/control (eg. CD4_f_case)) and off in any 3 of the 4 CD4 groups
same for CD8 and CD19


/********************************************************************************
* Normalization using the log2 of the upper quartile adjusted APN.
* into normalization - fragments
********************************************************************************/

/* Create clean dataset 

UQmult of 57 = q3 across all samples of the median of the reads in region

Also calculating fudge factor based on mapped reads

        mapped_FF = (sum_mapped / (median*10**5))
        map_apn = (apn*mapped_FF); 
        log_map_apn = log(map_apn);
        uq_apn = (apn/q3)*57;       
        uq_ff = 57/q3;
        log_uq_apn = log(UQ_apn);

*/



data design ;
set seq.illumina_RO1_meta_design ;
ID = tranwrd(mclab_sampleID, '-', '_') ;
if flag_female = 1 then sex = "f";
else sex = "m" ;
if flag_case = 0 then t1d = "control";
else t1d = "case" ;

num = scan(sampleID, 1, '-') ;
sex = scan(sampleID, 2, '_') ;
type = scan(sampleID, 3, '_') ;
newID = compress(population||'_'||sex||'_'||type||'_'||num) ;

keep mclab_sampleID sampleID population flag_female flag_case sex t1d newID ;
run;

proc sort data = design ;
by sampleID ;
run;

%macro prepping (cell) ;

proc sort data=seq.cvrg_frag_&cell._stack ;
by sampleID ;
run ;

data clean1_&cell.;
merge seq.cvrg_frag_&cell._stack (in=in1) design (in=in2);
by sampleID ;
if in1 ;
run;

proc sort data = clean1_&cell. ;
by featureID ;
proc sort data=seq.onCall_frags_&cell._apn5 ;
by featureID;
proc sort data = seq.flag_shrt_frags ;
by featureID ;
run;

data clean2_&cell.;
merge clean1_&cell. (in=in1) seq.onCall_frags_&cell._apn5  (in=in2) seq.flag_shrt_frags (in=in3) ;
by featureID ;
if in1 ;
run;

data list_&cell. ;
set clean2_&cell. ;
keep featureID sampleID flag_: ;
run;

%mend ;

%prepping (CD4) ;
%prepping (CD8) ;
%prepping (CD19) ;



%macro dumping (cell) ;
/* drop if feature less than 10 bp */

data clean_&cell. ;
set clean2_&cell. ;
if flag_frag_less_10bp = 1 then delete ;
run ;

/* drop if off in all groups or only on in 1 group */
data close_&cell ;
set clean_&cell ;
&cell._sum = flag_&cell._f_case_on5 + flag_&cell._f_control_on5 + flag_&cell._m_case_on5 + flag_&cell._m_control_on5;
if &cell._sum = 0 or &cell._sum = 1 then delete ;
run;

%mend ;

%dumping (CD4) ;   /* 11,007,669 obs */
%dumping (CD8) ;   /* 9,533,538 obs */
%dumping (CD19) ;  /* 6,282,496 obs */



%macro norming (cell, apn) ;

/* Calculate With-in sample basic statistics */
    proc means data=close_&cell. noprint;
        class sex t1d ;
        var reads_in_region;
        output out=mapped_reads_&cell. sum=sum_mapped q3=q3 median=median;
        run;

    data totals_&cell.;
        set mapped_reads_&cell.;
        where _type_ = 3;
        keep sex t1d sum_mapped q3 median;
        run;

    /* Summarize statistics to experiment level */
        proc sort data=totals_&cell.;
            by sex t1d;
            run;
        title "&cell.";
        proc print data = totals_&cell. ; run;

            /*              CD4         CD8         CD19
                f case      71.11      57.69       56.48
                f control   67.74      69.55       70.13
                m case      63.37      69.67       60.91
                m control   68.18      68.74       69.27  */     

        proc means data=totals_&cell. median ;
            var sum_mapped q3 median;
            run;
            /*  q3 = 67.96   CD4 apn5
                q3 = 69.14   CD8 apn5
                q3 = 65.09   CD19 apn5

            FF = (100 / sum_mapped)  */
    title "";
 /* Merge statistics onto dataset */
    proc sort data=close_&cell;
        by  sex t1d ;
        run;

    proc sort data=totals_&cell.;
        by sex t1d ;
        run;

    data stack_&cell. oops_&cell.;
        merge close_&cell. (in=in1) totals_&cell. (in=in2);
        by sex t1d ;
        if in1 then output stack_&cell.;
        else output oops_&cell.; * 0 obs yay!;
        run;
%mend ;
%norming (CD4, 5) ;
%norming (CD8, 5) ;
%norming (CD19, 5) ;


    /* NOTE MAKE SURE TO CHANGE THE UQ MULTIPLIER BY THE NEW Q3 
        *** 68 for CD4 
        *** 69 for CD8
        *** 65 for CD19

        mapped_FF = (sum_mapped / (median*10**5))
        map_apn = (apn*mapped_FF); 
        log_map_apn = log(map_apn);
        uq_apn = (apn/q3)*100;       
        uq_ff = 100/q3;
        log_uq_apn = log(UQ_apn);

*/

%macro stackem (cell, UQ) ;

    data stack_uq_&cell.;
        retain featureID;
        set stack_&cell.;
        mapped_FF = (sum_mapped / (median*10**5)) ;
        map_apn = (apn*mapped_FF); 
        log_map_apn = log(map_apn); 
        uq_apn = (apn/q3)*&UQ;
        uq_ff = &UQ/q3;
        log_uq_apn = log2(UQ_apn);
        run;

    proc freq data=stack_uq_&cell. noprint;
        table featureID/ out=freqs_&cell.;
        run;

%mend ;
%stackem (CD4, 67) ;
%stackem (CD8, 69) ;
%stackem (CD19, 65) ;


%macro siding (cell) ;

proc sort data = stack_uq_&cell. ;
by featureID sampleID ;
proc sort data = list_&cell. ;
by featureID sampleID ;
run;

data stack2_uq_&cell ;
merge list_&cell. (in=in1) stack_uq_&cell. (in=in2) ;
by featureID sampleID ;
if in1 ;
run ;
%mend ;

%siding (CD4) ;  
%siding (CD8) ;
%siding (CD19) ;  

/* create final datasets */
proc contents data = stack2_uq_CD4 ; run ;

%macro finals (cell) ;

data norm_data_frag_&cell. ;
set stack2_uq_&cell ;
keep featureID apn cv log_: map: q3 read: region: rpkm sampleID std sum: uq: ;
run;

data seq.norm_data_frag_&cell._stack ;
set norm_data_frag_&cell. ;
run ;

proc export data =seq.norm_data_frag_&cell._stack 
outfile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType/quantify_t1d_pacbio_transcripts/norm_data_frag_&cell._stack.csv"
dbms =csv replace ;
run ;

%mend ;

%finals (CD4) ;  
%finals (CD8) ;
%finals (CD19) ;





