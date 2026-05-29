libname seq "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType/sasdata" ;


filename mymacros "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/2018/PacBio/sas_programs/macros";
options SASAUTOS=(sasautos mymacros);



/*


** netanya giving me coverage counts summed across techreps! 

import cvrg cnts summed across TR

*/





/* list of samples to loop over */

data dsgn_lp ;
retain mclab_sampleID sampleID ID;
set seq.illumina_RO1_meta_design ;
ID = tranwrd(mclab_sampleID, '-', '_') ;
*where population ne "CD19";
rename mclab_sampleID = mclabID ;

num = scan(sampleID, 1, '-') ;
sex = scan(sampleID, 2, '_') ;
type = scan(sampleID, 3, '_') ;
newID = compress(population||'_'||sex||'_'||type||'_'||num) ;
run;

proc freq data = dsgn_lp ;
tables analytic_id / out = cnt_individs ;
run;  /* 113 individuals */

proc freq data = dsgn_lp ;
tables flag_case *flag_female * population / out=cnts ;
run;

proc print data = cnts ; run;

/*
flag_     flag_
 case    female    population    COUNT

  0         0         CD19         31
  0         0         CD4          46
  0         0         CD8          38
  0         1         CD19         21
  0         1         CD4          39
  0         1         CD8          35
  1         0         CD19          7
  1         0         CD4          13
  1         0         CD8          13
  1         1         CD19          5
  1         1         CD4          15
  1         1         CD8          12


*/

data dsgn_loop ;
set dsgn_lp ;
keep mclabID sampleID ID population newID ;
run ;

data test ;
set dsgn_loop;
where mclabID = "40010403-CD4";
run ;


/* import frag coverage counts - note techreps are summed */
%macro importing (mclabID, sampleID, ID, population, newID) ;

     data WORK.IN_&ID. ;
     %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
     infile "/TB14/TB14/concannon/coverage_cnts_isoseq_fragments_all_genes_combine_reps/sum_cvg_counts_&mclabID..csv"
     delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
        informat fusion_id $28. ;
        informat mapped_reads best32. ;
        informat read_length best32. ;
        informat region_length best32. ;
        informat region_depth best32. ;
        informat reads_in_region best32. ;
        informat apn best32. ;
        informat rpkm best32. ;
        informat mean best32. ;
        informat std best32. ;
        informat cv best32. ;
        format fusion_id $28. ;
        format mapped_reads best12. ;
        format read_length best12. ;
        format region_length best12. ;
        format region_depth best12. ;
        format reads_in_region best12. ;
        format apn best12. ;
        format rpkm best12. ;
        format mean best12. ;
        format std best12. ;
        format cv best12. ;
     input
                 fusion_id  $
                 mapped_reads
                 read_length
                 region_length
                 region_depth
                 reads_in_region
                 apn
                 rpkm
                 mean
                 std
                 cv   ;
     if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
     run;

data in2_&newID. ;
length sampleID $32. ;
set in_&ID. ;
sampleID = "&sampleID." ;
cellType = "&population." ;
run;

%mend ;

%iterdataset(dataset=dsgn_loop, function=%nrstr(%importing(&mclabID, &sampleID, &ID, &population, &newID);)); 


%macro pops (population) ;

data cvrg_frag_&population. ;
length sampleID $32. ;
length cellType $4. ;
set in2_&population._: ;
run ; 

%mend ;
%pops (CD4) ;
%pops (CD8) ;
%pops (CD19) ;  



/* run some count checks */
proc freq data = cvrg_frag_cd4 ;
tables fusion_id / out = cnts ;
run ;
proc freq data = cnts ;
tables fusion_id / out = cntsing;
run;
data ck2 ;
set cntsing ;
where count ne 1 ;
run;  /* all good */


%macro split (cell) ;

data seq.cvrg_frag_&cell._stack ;
set cvrg_frag_&cell.;
rename fusion_id = featureID ;
drop cellType ;
run;

proc contents data = seq.cvrg_frag_&cell._stack ; run;

data dsgn_&cell. ;
retain mclab_sampleID sampleID ID;
set seq.illumina_RO1_meta_design ;
ID = tranwrd(mclab_sampleID, '-', '_') ;
where population = "&cell.";
rename mclab_sampleID = mclabID ;
run;

proc sort data = dsgn_&cell. ;
by sampleID ;
proc sort data = seq.cvrg_frag_&cell._stack ;
by sampleID ;
run;

data cnts_&cell ;
merge seq.cvrg_frag_&cell._stack (in=in1) dsgn_&cell. (in=in2) ;
by sampleID ;
run ;

title "cnts &cell.";
proc freq data = cnts_&cell ;
tables flag_case flag_female  ;
run;

%mend ;

%split (CD4) ;   /*  obs */
%split (CD8) ;   /*  obs */
%split (CD19) ;  /*  obs */

/*
cellType    num_obs     flag_case=1     flag_female=1
CD19       10,188,864    1,910,412       4,139,226
CD8        15,601,698    3,980,025       7,482,447 
CD4        17,989,713    4,457,628       8,596,854                

*/

%macro cnts (cell) ;   

title "featureID counts for &cell. ";
proc freq data =  cnts_&cell. noprint ;
tables featureID / out = ck_&cell._numID ;
run;

proc sql ;
select count(*) as N from ck_&cell._numID ;
quit ;

title "sample counts for &cell. ";
proc freq data = cnts_&cell. noprint ;
tables sampleID / out = ck_&cell._numSample ;
run ;

proc sql ;
select count(*) as N from ck_&cell._numSample ;
quit ;
title "";
%mend ;

%cnts (cd4) ;  /* 113 samples, 159,201 features */
%cnts (cd8) ;  /*  98 samples, 159,201 features */
%cnts (cd19) ;  /* 64 samples, 159,201 features */




