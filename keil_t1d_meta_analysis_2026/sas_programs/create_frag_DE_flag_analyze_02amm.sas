libname seq "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType/sasdata" ;


/* 
create flag_analyze for DE


flagAnalyze_frag_DE = 1 if 

    (1) drop fragment if off in all CD4 groups (group = cellType by sex by case/control (eg. CD4_f_case)) and off in 1 of the 4 CD4
        &cell._sum = flag_&cell._f_case_on5 + flag_&cell._f_control_on5 + flag_&cell._m_case_on5 + flag_&cell._m_control_on5;
        if &cell._sum = 0 or &cell._sum = 1 then flag_drop_frag = 1  ;
        else flag_drop_frag = 0 
        
            flag_frag_group_off = 1 --> drop

    (2) drop fragments less than 10bp 
            flag_frag_less_10bp = 1 --> drop

    (3) for frags that are on:  # cases ge 4 AND # controls ge 4 
            flag_frag_num_on = 1 --> then both case and control have > 4


input:  
    seq.onCall_frags_&cell._apn5
    seq.flag_shrt_frags
    seq.cvrg_frag_&cell._stack

output:
    seq.flag_analyze_frag_DE_&cell.

        flag_analyze_frag_DE = 1 

*/

%macro callOffs (cell) ;

/* 1 */
data oncalls_&cell. ;
set seq.onCall_frags_&cell._apn5  ;
&cell._sum = flag_&cell._f_case_on5 + flag_&cell._f_control_on5 + flag_&cell._m_case_on5 + flag_&cell._m_control_on5;
if &cell._sum = 0 or &cell._sum = 1 then flag_frag_group_off = 1  ;
else flag_frag_group_off = 0 ;
keep featureID flag_frag_group_off ;
run ;

/* 2 */
data shorts ;
set seq.flag_shrt_frags;
run ;

proc sort data = shorts ;
by featureID ;
run ;

/* 3 */
proc sort data = seq.cvrg_frag_&cell._stack ;
by featureID ;
run;

data addCnts_&cell. ;
set seq.cvrg_frag_&cell._stack ;
var1 = scan(sampleID, 1, '_') ;
var2 = scan(var1, 2, '-') ;
var3 = scan(sampleID, 3, '_') ;
ID = compress(var2||'_'||var3) ;
keep sampleID featureID ID   ;
run;

proc freq data = addCnts_&cell. noprint ;
by featureID ;
tables ID / out = byFeature_&cell.;
run ;

proc transpose data = byFeature_&cell. out = cnts_byFeat_&cell prefix=num_;
by featureID ;
var count ;
id ID ;
run ;

data cnts2_byFeat_&cell ; ;
set cnts_byFeat_&cell ;
if num_&cell._case > 4 and num_&cell._control > 4 then flag_frag_num_on = 1 ;
else flag_frag_num_on = 1 ;
keep featureID flag_frag_num_on ;
run;

/* combine flags */
data flag_&cell. ;
merge oncalls_&cell.  (in=in1) shorts (in=in2)  cnts2_byFeat_&cell (in=in3) ;
by featureID ;
run ;

data seq.flag_analyze_frag_DE_&cell. ;
set flag_&cell. ;
if flag_frag_group_off = 0 and flag_frag_less_10bp = 0 and flag_frag_num_on = 1 then flag_analyze_frag_DE = 1 ;
else  flag_analyze_frag_DE = 0 ;
run;

title "DE flag analyze for &cell." ;
proc freq data = seq.flag_analyze_frag_DE_&cell. ;
tables flag_analyze_frag_DE ;
run;
title "";
%mend ;

%callOffs (CD4) ;   /* 97,413 out of 151,719 (64.2%) are analyzable */
%callOffs (CD8) ;   /* 97,281 out of 151,719 (64.1%) are analyzable */
%callOffs (CD19) ;  /* 98,164 out of 151,719 (64.7) are analyzable */


/* export for checking */
%macro exp (cell) ;
proc export data = seq.flag_analyze_frag_DE_&cell.
outfile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType/quantify_t1d_pacbio_transcripts/flag_analyze_frag_DE_&cell..csv"
dbms = csv replace ;
run;

%mend ;

%exp (CD4) ;
%exp (CD8) ;
%exp (CD19) ;



