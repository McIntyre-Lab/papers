libname seq "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType/sasdata" ;


/*
on - offs for fragment coverage counts for all genes
separately for each cell type

flag a fragment as 'on' if on in 50% of the individuals in the group

group = cellType by sex by case/control (eg. CD4_f_case)


flag on if 50% of individuals are on at apn 5 --> flag_&cell._&sex._&t1d_on5 = 1


*/

data design ;
set seq.illumina_RO1_meta_design ;
ID = tranwrd(mclab_sampleID, '-', '_') ;
if flag_female = 1 then sex = "f";
else sex = "m" ;
if flag_case = 0 then t1d = "control";
else t1d = "case" ;
keep mclab_sampleID sampleID population flag_female flag_case sex t1d enrollment_age ;
run;

proc sort data = design nodups;
by _all_ ;
run;

%macro outer (cell) ;

proc sort data = seq.cvrg_frag_&cell._stack ;
by sampleID ;
proc sort data  = design;
by sampleID ;
run;

data cnts_df_&cell. oops_&cell.  ;
merge seq.cvrg_frag_&cell._stack  (in=in1) design (in=in2) ;
by sampleID ;
if in1 then output cnts_df_&cell. ;
else output oops_&cell. ;
run;  /* 0 in oops */

data cnt_flag_&cell. ;
set cnts_df_&cell.;
if apn > 0 then apn_on0 = 1; else apn_on0 = 0;
if apn > 5 then apn_on5 = 1; else apn_on5 = 0;
run;
%mend ;

%outer (CD4) ;
%outer (CD8) ;
%outer (CD19) ;

%macro onCalls (apn, cell) ;

proc sort data=cnt_flag_&cell.;
  by featureID population sex t1d  ;
  run;

proc means data = cnt_flag_&cell. noprint ;
    by featureID population sex t1d ;
    var apn_on&apn ;
    output out = sample_on&apn._&cell. mean=sample_percent_on&apn ;
    run ;

proc sort data = sample_on&apn._&cell. ;
    by featureID ;
    run;

data sample2_on&apn._&cell. ;
retain featureID sample ;
set  sample_on&apn._&cell. ;
/* frag is on if expressed in in 50% of the reps */
if sample_percent_on&apn > 0.5 then flag_frag_on&apn = 1 ;
    else flag_frag_on&apn = 0 ;
sample = compress(population||'_'||sex||'_'||t1d);
run ;

proc transpose data = sample2_on&apn._&cell. out = sample_flag_on&apn._&cell. prefix = flag_ suffix= _on&apn;
    by featureID ;
    id sample ;
    var flag_frag_on&apn ;
    run;

data onCall_frag_&cell._APN&apn. ;
set sample_flag_on&apn._&cell. ;
drop _name_ ;
run ;

data seq.onCall_frags_&cell._APN&apn. ;
set onCall_frag_&cell._APN&apn. ;
run;

%mend  ;
%onCalls (0, CD4) ; 
%onCalls (0, CD8) ; 
%onCalls (0, CD19) ;

%onCalls (5, CD4) ;
%onCalls (5, CD8) ;
%onCalls (5, CD19) ;



proc freq data = seq.onCall_frags_CD4_apn5 ;
tables flag_CD4_f_case_on5 * flag_CD4_m_case_on5 *flag_CD4_f_control_on5 *flag_CD4_m_control_on5  / out = cnt_cd4;
run; 

proc print data = cnt_cd4 ; run; 

proc freq data = seq.onCall_frags_CD8_apn5 ;
tables flag_CD8_f_case_on5 * flag_CD8_m_case_on5 *flag_CD8_f_control_on5 *flag_CD8_m_control_on5  / out = cnt_cd8;
run; 

proc print data = cnt_cd8 ; run; 

proc freq data = seq.onCall_frags_CD19_apn5 ;
tables flag_CD19_f_case_on5 * flag_CD19_m_case_on5 *flag_CD19_f_control_on5 *flag_CD19_m_control_on5  / out = cnt_cd19;
run; 

proc print data = cnt_cd19 ; run; 

/*

                                     flag_CD4_     flag_CD4_
         flag_CD4_     flag_CD4_    f_control_    m_control_
 Obs    f_case_on5    m_case_on5        on5           on5        COUNT

   1         0             0             0             0         43186
   2         0             0             0             1           369
   3         0             0             1             0           299
   4         0             0             1             1           335
   5         0             1             0             0           321
   6         0             1             0             1           367
   7         0             1             1             0           122
   8         0             1             1             1           420
   9         1             0             0             0           989
  10         1             0             0             1           481
  11         1             0             1             0           377
  12         1             0             1             1          1645
  13         1             1             0             0           124
  14         1             1             0             1           218
  15         1             1             1             0           172
  16         1             1             1             1        102294



                                    flag_CD8_     flag_CD8_
        flag_CD8_     flag_CD8_    f_control_    m_control_
Obs    f_case_on5    m_case_on5        on5           on5       COUNT

  1         0             0             0             0        42973
  2         0             0             0             1          374
  3         0             0             1             0          974
  4         0             0             1             1          879
  5         0             1             0             0          971
  6         0             1             0             1          572
  7         0             1             1             0          696
  8         0             1             1             1         3538
  9         1             0             0             0          171
 10         1             0             0             1           75
 11         1             0             1             0          133
 12         1             0             1             1          454
 13         1             1             0             0           76
 14         1             1             0             1           76
 15         1             1             1             0          185
 16         1             1             1             1        99572


                                   flag_CD19_    flag_CD19_
       flag_CD19_    flag_CD19_    f_control_    m_control_
Obs    f_case_on5    m_case_on5        on5           on5       COUNT

  1         0             0             0             0        41547
  2         0             0             0             1          886
  3         0             0             1             0          616
  4         0             0             1             1         3016
  5         0             1             0             0          696
  6         0             1             0             1          474
  7         0             1             1             0          184
  8         0             1             1             1         2842
  9         1             0             0             0          705
 10         1             0             0             1          251
 11         1             0             1             0          239
 12         1             0             1             1         2701
 13         1             1             0             0          197
 14         1             1             0             1          200
 15         1             1             1             0          139
 16         1             1             1             1        97026




*/


