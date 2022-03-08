
/*
flag on/off 
tpm



import coverage counts and set all samples together

*/


filename mymacros "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/2015/sas_programs/macros";
options SASAUTOS=(sasautos mymacros);

/* macro from /nfshome/adalena.nanni/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/2015/sas_programs/macros/iterdataset.sas
    added here due to sasautos not working properly for AVN
*/
%macro iterdataset(dataset=,function=);
    %local dsid now total rows cols rc;
    %let dsid = %sysfunc(open(&dataset));
    %let now = 0;
    %let rows = %sysfunc(attrn(&dsid, nobs));
    %let cols = %sysfunc(attrn(&dsid, nvars));

    %do %while(%sysfunc(fetch(&dsid)) = 0); %* outer loop across rows;
        %let now = %eval(&now + 1);

        %do i = 1 %to &cols; %* inner loop across coloumns;
            %local v t;
            %let v=%sysfunc(varname(&dsid,&i));
            %local &v;
            %let t = %sysfunc(vartype(&dsid,&i));
            %let &v = %sysfunc(getvar&t(&dsid,&i));
        %end;

        %unquote(&function);

    %end;
    %let rc = %sysfunc(close(&dsid));
%mend;


data shrt ;
set cvrg_cnts_shrtRead ;
if wt_apn > 0 then wt_apn_on0 = 1; else wt_apn_on0 = 0;
if wt_tpm > 0 then wt_tpm_on0 = 1; else wt_tpm_on0 = 0;

if wt_apn > 5 then wt_apn_on5 = 1; else wt_apn_on5 = 0;
if wt_tpm > 5 then wt_tpm_on5 = 1; else wt_tpm_on5 = 0;
run ;


/* Going to use wt_tpm 5 for flagging on/off */
proc sort data= shrt;
by geno trt gene_id;
run;
 
%macro flag_on_off (var1, var2) ;

proc means data = shrt noprint ;
    by geno trt gene_id;
    var wt_&var1._on&var2.;
    output out = shrt_trt_&var1._on&var2. mean=trt_prcnt_&var1._&var2.;
    run ;

proc sort data = shrt_trt_&var1._on&var2. ;
    by gene_id;
    run;

data shrt2_trt_&var1._on&var2. ;
set shrt_trt_&var1._on&var2.;
var = compress(geno||'_'||trt) ;
run ;

proc transpose data = shrt2_trt_&var1._on&var2.  out = shrt_trt_&var1._on&var2._sbys ;
    by gene_id ;
    id var ;
    var trt_prcnt_&var1._&var2. ;
    run;


/* flag definitions:

*/

data onCalls_shrt_gene_&var1.&var2. ;              * if 50% of reps then is expressed ;
    set shrt_trt_&var1._on&var2._sbys ;           * using geno trt to determine ;

    if b73_amb > 0.5 then flag_b73_amb_on = 1; else flag_b73_amb_on = 0;
    if b73_ele > 0.5 then flag_b73_ele_on = 1; else flag_b73_ele_on = 0;

    if mo17_amb > 0.5 then flag_mo17_amb_on = 1; else flag_mo17_amb_on = 0;
    if mo17_ele > 0.5 then flag_mo17_ele_on = 1; else flag_mo17_ele_on = 0;

    if c123_amb > 0.5 then flag_c123_amb_on = 1; else flag_c123_amb_on = 0;
    if c123_ele > 0.5 then flag_c123_ele_on = 1; else flag_c123_ele_on = 0;

    if hp301_amb > 0.5 then flag_hp301_amb_on = 1; else flag_hp301_amb_on = 0;
    if hp301_ele > 0.5 then flag_hp301_ele_on = 1; else flag_hp301_ele_on = 0;

    if nc338_amb > 0.5 then flag_nc338_amb_on = 1; else flag_nc338_amb_on = 0;
    if nc338_ele > 0.5 then flag_nc338_ele_on = 1; else flag_nc338_ele_on = 0;

    /* flag_analyze_b73 = 1 if flag_b73_amb_on = 1 or flag_b73_ele_on = 1 */
    if flag_b73_amb_on = 1 or flag_b73_ele_on = 1 then flag_analyze_b73 = 1; else flag_analyze_b73 = 0 ;
    if flag_mo17_amb_on = 1 or flag_mo17_ele_on = 1 then flag_analyze_mo17 = 1; else flag_analyze_mo17 = 0 ;
    if flag_c123_amb_on = 1 or flag_c123_ele_on = 1 then flag_analyze_c123 = 1; else flag_analyze_c123 = 0 ;
    if flag_hp301_amb_on = 1 or flag_hp301_ele_on = 1 then flag_analyze_hp301 = 1; else flag_analyze_hp301 = 0 ;
    if flag_nc338_amb_on = 1 or flag_nc338_ele_on = 1 then flag_analyze_nc338 = 1; else flag_analyze_nc338 = 0 ;
    
keep gene_id flag_: ;
run ;


proc export data = onCalls_shrt_gene_&var1.&var2.
outfile = "/nfshome/adalena.nanni/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/2018/PacBio/cvr_cnts_gene_mo17_cau/onCalls_shrt_gene_&var1.&var2..csv"
dbms = csv replace ;
run;

%mend ;

/* %flag_on_off (tpm, 0) ; */
%flag_on_off (tpm, 5) ;
/* %flag_on_off (apn, 0) ; */
/* %flag_on_off (apn, 5) ; */









