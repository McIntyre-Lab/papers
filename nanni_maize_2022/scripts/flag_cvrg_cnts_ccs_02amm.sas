
libname make "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/make_combination_flag_file/sasdata";


/*
import coverage counts for ccs reads



*/



filename mymacros "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/2015/sas_programs/macros";
options SASAUTOS=(sasautos mymacros);

proc datasets library = WORK kill nolist ;
run;



data pb_list ;
set make.list_ccs_samples ;
run ;

proc sort data = pb_list nodups ;
by sample;
run;


%macro rename (oldSample, sample, geno, trt, rep) ; /* oldSample = 21_2_mo17_oz, sample = mo17_oz_2 */

proc import datafile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/make_combination_flag_file/cvr_cnts_ccs_b73_gene/cvrg_cnts_gene_&oldSample..csv"
out = cc_&sample
dbms = csv replace ;
guessingrows = MAX ;
run;

data rename_&sample.  ;
set cc_&sample;
if reads_in_region > 0 then flag_detect_ccs = 1 ;
    else flag_detect_ccs = 0 ;
run ;

proc sql noprint;
   select cats(name,'=',name,"_&geno._&trt._&rep.")
          into :list
          separated by ' '
          from dictionary.columns
          where libname = 'WORK' and memname = upcase("rename_&sample.")  ;
quit;

proc datasets library = work nolist ;
modify rename_&sample. ;
rename &list ;
quit ;

data rename2_&sample ;
set rename_&sample ;
rename primary_FBgn_&geno._&trt._&rep = primary_FBgn ;
keep primary_FBgn_&geno._&trt._&rep flag_detect_: ;
run ;

proc sort data = rename2_&sample. ;
by primary_FBgn ;
run;

%mend ;

/*%rename (mo17_p4_c6_amb, mo17_amb_1, mo17, amb, 1); */   /* oldSample, sample, geno, trt, rep */


%iterdataset(dataset=pb_list, function=%nrstr(%rename(&oldSample, &sample, &geno, &trt, &rep);)); 


%macro by_geno (geno) ;  /* for all geno but mo17 */

/* for each geno_trt, merge by gene */
data cvrg_cnts_ccs_&geno._sbys  ;
merge rename2_&geno._: ;
by primary_FBgn ;
run ;

/* sum across flags - if sum > 0 then at least 1 rep in oz or amb has a read */
data flag_detect_ccs_&geno. ;
retain primary_FBgn flag_detect_ccs_&geno._gt0 ; 
set cvrg_cnts_ccs_&geno._sbys ;

if flag_detect_ccs_&geno._amb_1 > 0 or flag_detect_ccs_&geno._oz_1 > 0 then flag_detect_ccs_&geno._gt0  = 1 ;
    else flag_detect_ccs_&geno._gt0 =  0 ;

rename flag_detect_ccs_&geno._amb_1 = flag_detect_ccs_&geno._amb_gt0;
rename flag_detect_ccs_&geno._oz_1 =  flag_detect_ccs_&geno._ele_gt0 ;

keep primary_FBgn flag_detect_ccs_&geno._gt0 flag_detect_ccs_&geno._amb_1 flag_detect_ccs_&geno._oz_1;
run ;

%mend ;

%by_geno (b73) ;
%by_geno (c123) ;
%by_geno (hp301) ;
%by_geno (nc338) ;



/* for mo17_trt, merge by gene */
data flag_detect_mo17  ;
merge rename2_mo17_: ;
by primary_FBgn ;
run ;

/* sum across flags - if sum > 0 then at least 1 rep in oz or amb has a read */ 
data flag_detect_ccs_mo17 ;
retain primary_FBgn flag_detect_ccs_mo17_gt0 ; 
set flag_detect_mo17 ;

flag_detect_sum_ele = sum(of flag_detect_ccs_mo17_oz_1-flag_detect_ccs_mo17_oz_2) ;

if flag_detect_ccs_mo17_amb_1 > 0 then flag_detect_ccs_mo17_amb_gt0 = 1 ;     else flag_detect_ccs_mo17_amb_gt0 = 0 ;

if flag_detect_sum_ele > 0 then flag_detect_ccs_mo17_ele_gt0 = 1 ;     else flag_detect_ccs_mo17_ele_gt0 = 0 ;

if flag_detect_ccs_mo17_amb_gt0 = 1 or flag_detect_ccs_mo17_ele_gt0 = 1 then flag_detect_ccs_mo17_gt0 = 1 ;
    else flag_detect_ccs_mo17_gt0 = 0;

keep primary_FBgn flag_detect_ccs_mo17_amb_gt0 flag_detect_ccs_mo17_ele_gt0 flag_detect_ccs_mo17_gt0 ;
run ;


data make.flag_detect_cvrg_cnts_ccs ;
merge flag_detect_ccs_:  ;
by primary_FBgn ;
run;

proc export data = make.flag_detect_cvrg_cnts_ccs 
outfile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/make_combination_flag_file/flag_detect_cvrg_cnts_ccs.csv"
dbms = csv replace ;
run;

proc contents data = make.flag_detect_cvrg_cnts_ccs ; run;




